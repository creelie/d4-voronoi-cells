"""c27_hull_sub.py -- c28_hull_case.py with some bins made types of their own (sub-types)
for the thresholds: the best convex combination of existing combined certificates for the
count vectors of a case file, with one threshold per kind of pair or triple of sub-types.

env SUB=<bins> (0-indexed, comma-separated): each listed bin becomes a sub-type whose range
is the bin, while the kernel and the three-point polynomials stay those of its type.  A
pair or triple of centres is then bounded with the constant of the sub-types it occupies,
as REFINE_PAIRS, REFINE_TRIPLES and SKIP_COVERED=1 do in multi_cap/combo_case_check.py, so
the combination should be checked with every pair and triple of bins that contains a
listed bin refined.  The npz written has c2 and c3 per kind of types (largest over the
sub-kinds), for comparison only.

c28_hull_case.py -- the best convex combination of existing combined certificates for the
count vectors of a case file at 28 centres.

A convex combination of certificates is a certificate: the positivity conditions are
convex and every inequality is linear in the kernel, with thresholds recomputed for the
combination.  With only the weights as unknowns the programme is a small LP, so it can
take far more samples than the full programme, and the bulges of single certificates
between their samples can cancel.  Float throughout; proves nothing.

The case file (env CASE, a path relative to this directory; default case28_lo16.json)
may narrow the types to the bins one count vector occupies, e.g. a close type [2, 2.0161]
and a middle type equal to one bin, leaving gaps between the type ranges.  Such a case
is for this programme only: the combination is then checked on the case with the usual
bins (remap_bins.py recomputes its bin constants), with REFINE_TRIPLES, REFINE_PAIRS and
SKIP_COVERED=1 in multi_cap/combo_case_check.py, which bound the inner products by a at
the upper ends of the occupied bins, as here.

usage: python3 c28_hull_case.py TARGET ITERS cert.npz [cert.npz ...]
TARGET: all (every count vector of the case), or one of SETS for case28_lo16.json
env: CASE; NP, NT (random samples per kind); OUT (prefix of the npz of the best combination);
MG_AA, MG2, MG3, MG_FFF (margins added to the thresholds in the reported bound); FACE=n
(an n-by-n grid on each face of a triple domain where one pair touches); CLIMB_N (climbing
starts per kind and iteration, default 60); PEAKS=<c28_peak2.py log>[,...] (the points
where its ascents ended, added to the triple samples).
"""
import sys, os, time, itertools, json
import numpy as np
from math import comb
from scipy.optimize import linprog
import scipy.sparse as sp
CM = os.path.dirname(os.path.abspath(__file__)) + '/'
TARGETS, iters, certs = sys.argv[1], int(sys.argv[2]), sys.argv[3:]
CASEF = os.environ.get('CASE', 'case28_lo16.json')
sys.argv = ['x', '8', '1', CM + CASEF]
sys.path.insert(0, CM)
import combo_gen2 as G
import typed3pt as T3
from cert_embed import embed, blocks, to_x3
t0 = time.time()
TYPES0 = G.TYPES; ti = {s: i for i, s in enumerate(TYPES0)}
BT = [b[0] for b in G.BINS]; nb = len(G.BINS)
SUB = [int(a) for a in os.environ.get('SUB', '').split(',') if a]
SBT = [('%s%d' % (BT[b], b + 1)) if b in SUB else BT[b] for b in range(nb)]   # sub-type of each bin
PAR = {}
for b in range(nb):
    PAR[SBT[b]] = BT[b]
    if b in SUB:
        G.TRANGE[SBT[b]] = (G.BINS[b][1], G.BINS[b][2])
TYPES = list(dict.fromkeys(SBT))   # in bin order, so the parent types stay in the order of the case
for s in TYPES0:
    PAR.setdefault(s, s)
DEG = {2052: 6, 4098: 8, 7188: 10, 11538: 12}
Zs = [np.load(c if c.startswith('/') else CM + c) for c in certs]
D = max(DEG[len(Z['x3'])] for Z in Zs)
K = len(Zs)


def clip(M):
    """the nearest positive semidefinite matrix: negative eigenvalues set to 0, as the exact
    check does before it rounds (solutions that stopped on a numerical error can have
    eigenvalues near -1e-6, which move the triple polynomials by up to 1e-5)"""
    M = (np.array(M, float) + np.array(M, float).T) / 2
    w, V = np.linalg.eigh(M)
    return (V * np.maximum(w, 0)) @ V.T


def projected(Z):
    d = DEG[len(Z['x3'])]
    B = {name: clip(M) for name, M in blocks(np.array(Z['x3'], float), 3, d).items()}
    A = np.array(Z['A'], float)
    nA = A.shape[1]
    Bd = np.zeros((nA + 1, nA + 1))
    Bd[:nA, :nA] = A[0]; Bd[:nA, nA] = Bd[nA, :nA] = np.array(Z['z'], float); Bd[nA, nA] = float(Z['t'])
    Bd = clip(Bd)
    A = np.stack([Bd[:nA, :nA]] + [clip(a) for a in A[1:]])
    return embed(to_x3(B, 3, d), 3, d, D), A, Bd[:nA, nA], Bd[nA, nA]


PJ = [projected(Z) for Z in Zs]
X3 = np.stack([p[0] for p in PJ], 1)                                                         # n3 x K
AV = np.stack([p[1] for p in PJ])                                                            # K x 13 x 5 x 5
ZV = np.stack([p[2] for p in PJ])                                                            # K x 5
TV = np.array([p[3] for p in PJ])
L = T3.Layout(3, D)
Bp = T3.Builder(L)
for s in TYPES0:
    T3.add_point(Bp, np.array([ti[s]]), ti[s])
PT = Bp.matrix(3) @ X3                                                                       # 3 x K
rng = np.random.default_rng(int(os.environ.get('SEED', '3')))
NP, NT = int(os.environ.get('NP', '60000')), int(os.environ.get('NT', '60000'))
CLIMB_N = int(os.environ.get('CLIMB_N', '60'))   # climbing starts per kind and iteration
PEAKPTS = {}   # PEAKS=<c28_peak2 log>[,...]: the points where its ascents ended, added to the triple samples
for _f in [f for f in os.environ.get('PEAKS', '').split(',') if f]:
    import re as _re
    for _l in open(_f):
        _m = _re.match(r'(\w+): sampled max .* at \[\s*([-\d.e]+)\s+([-\d.e]+)\s+([-\d.e]+)\s*\]', _l)
        if _m:
            PEAKPTS.setdefault(_m.group(1), []).append([float(_m.group(i)) for i in (2, 3, 4)])
PEAKPTS = {k: np.array(v) for k, v in PEAKPTS.items()}
FACE = int(os.environ.get('FACE', '0'))   # n > 0: an n-by-n grid on each face where one pair touches
print('%d certificates, three-point degree %d, n3 %d [%.0f s]' % (K, D, L.n, time.time() - t0), flush=True)


def pair_vals(st, X):
    out = []
    for i0 in range(0, len(X), 10000):
        x = X[i0:i0 + 10000]; p, q, u = x[:, 0], x[:, 1], x[:, 2]
        U, Pp, Pq = G.ubasis(u, G.D2), G.pbasis(p, G.R2), G.pbasis(q, G.R2)
        Kv = np.einsum('nk,na,ikab,nb->ni', U, Pp, AV, Pq, optimize=True)
        B = T3.Builder(L); T3.add_pair(B, np.arange(len(u)), ti[PAR[st[0]]], ti[PAR[st[1]]], u)
        out.append((Kv + B.matrix(len(u)) @ X3 - G.pair(p / 2, q / 2, u)[:, None]).astype(np.float32))
    return np.concatenate(out)


def trip_vals(cb, X):
    out = []
    for i0 in range(0, len(X), 10000):
        x = X[i0:i0 + 10000]
        B = T3.Builder(L); T3.add_triple(B, np.arange(len(x)), tuple(ti[PAR[c]] for c in cb), x[:, 0], x[:, 1], x[:, 2])
        out.append((B.matrix(len(x)) @ X3).astype(np.float32))
    return np.concatenate(out)


def bin_vals(d):
    Pd = G.pbasis(d, G.R2); U1 = G.ubasis(np.ones(1), G.D2)[0]
    Kv = np.einsum('k,na,ikab,nb->ni', U1, Pd, AV, Pd)
    return G.S(d)[:, None] + 0.5 * Kv - Pd @ ZV.T


def rand_gram(T, n):
    outs = []
    for dim in (4, 3, 2):
        V = rng.normal(size=(8 * n, 3, dim)); V /= np.linalg.norm(V, axis=2, keepdims=True)
        g = np.stack([(V[:, 0] * V[:, 1]).sum(1), (V[:, 0] * V[:, 2]).sum(1), (V[:, 1] * V[:, 2]).sum(1)], 1)
        outs.append(g[(g[:, 0] <= T[0]) & (g[:, 1] <= T[1]) & (g[:, 2] <= T[2])][:n])
    return np.concatenate(outs)


def tc(v):
    return {s: sum(n for n, b in zip(v, SBT) if b == s) for s in TYPES}


def npair(c, st):
    return comb(c[st[0]], 2) if st[0] == st[1] else c[st[0]] * c[st[1]]


def ntrip(c, cb):
    out = 1
    for s in set(cb):
        out *= comb(c[s], cb.count(s))
    return out


ALLV = [tuple(v) for v in G._counts()]
SETS = {'all': lambda v: True}
for _n in range(10, 25):
    SETS['k%d' % _n] = (lambda n: lambda v: v[0] == n)(_n)
def _sel(t):
    if t.startswith('v'):   # one count vector, e.g. v14-3-7-0-0-1-2
        w = tuple(int(a) for a in t[1:].split('-')); return lambda v: tuple(v) == w
    return SETS[t]
VS_ALL = [v for v in ALLV if any(_sel(t)(v) for tg in TARGETS.split(',') for t in tg.split('+'))]
pkeys = [st for st in itertools.combinations_with_replacement(TYPES, 2) if any(npair(tc(v), st) for v in VS_ALL)]
tkeys = [cb for cb in itertools.combinations_with_replacement(TYPES, 3) if any(ntrip(tc(v), cb) for v in VS_ALL)]
print('sub-types %s; %d pair kinds, %d triple kinds over %d count vectors' % (TYPES, len(pkeys), len(tkeys), len(VS_ALL)), flush=True)
PS, TS = {}, {}
for st in pkeys:
    P, Q, W = G.pair_samples(*st, n_d=17, n_u=60)
    p, q, u = G.rand_pairs(st[0], st[1], NP)
    X = np.r_[np.stack([P, Q, W], 1), np.stack([p, q, u], 1)]
    PS[st] = [X, pair_vals(st, X)]
for cb in tkeys:
    T = (G.tbound(cb[0], cb[1]), G.tbound(cb[0], cb[2]), G.tbound(cb[1], cb[2]))
    ax = [np.linspace(-1, t, int(os.environ.get('GRID', '48'))) for t in T]
    Gd = np.stack(np.meshgrid(*ax, indexing='ij'), -1).reshape(-1, 3)
    Gd = Gd[G.gram_ok(Gd)]
    X = np.r_[rand_gram(T, NT), Gd, G.coplanar_grid(*T, 64), G.symmetric_grid(cb, *T, 64)]
    if FACE:
        for i in range(3):
            j, k = [a for a in range(3) if a != i]
            A_, B_ = np.meshgrid(np.linspace(-1, T[j], FACE), np.linspace(-1, T[k], FACE), indexing='ij')
            Fg = np.zeros((FACE * FACE, 3)); Fg[:, i] = T[i]; Fg[:, j] = A_.ravel(); Fg[:, k] = B_.ravel()
            X = np.r_[X, Fg[G.gram_ok(Fg)]]
    X = np.r_[X, PEAKPTS.get(''.join(cb), np.zeros((0, 3)))]
    TS[cb] = [X, trip_vals(cb, X)]
BD = [np.linspace(lo, hi, 600) for _, lo, hi in G.BINS]
BV = [bin_vals(d) for d in BD]
print('samples: pairs %d, triples %d [%.0f s]' % (sum(len(v[0]) for v in PS.values()), sum(len(v[0]) for v in TS.values()),
                                                   time.time() - t0), flush=True)


# variable layout: lam (K), c2 (pair kinds), c3 (triple kinds), m (nb), eta
n2, n3 = len(pkeys), len(tkeys)
iL = np.arange(K); i2 = K + np.arange(n2); i3 = K + n2 + np.arange(n3); iM = K + n2 + n3 + np.arange(nb); iE = K + n2 + n3 + nb
NV = iE + 1
def lp(active):
    rows, b = [vrows], [np.zeros(len(vrows))]
    for j, st in enumerate(pkeys):
        V = PS[st][1][active['p'][st]]
        R = np.zeros((len(V), NV)); R[:, iL] = V.astype(float); R[:, i2[j]] = -1; rows.append(R); b.append(np.zeros(len(V)))
    for j, cb in enumerate(tkeys):
        V = TS[cb][1][active['t'][cb]]
        R = np.zeros((len(V), NV)); R[:, iL] = V.astype(float); R[:, i3[j]] = -1; rows.append(R); b.append(np.zeros(len(V)))
    for bb in range(nb):
        V = BV[bb]
        R = np.zeros((len(V), NV)); R[:, iL] = V.astype(float); R[:, iM[bb]] = -1; rows.append(R); b.append(np.zeros(len(V)))
    A = sp.vstack([sp.csr_matrix(R) for R in rows], format='csr'); bv = np.concatenate(b); rows.clear()
    Aeq = np.zeros((1, NV)); Aeq[0, iL] = 1
    cobj = np.zeros(NV); cobj[iE] = 1
    bounds = [(0, None)] * K + [(None, None)] * (NV - K)
    r = linprog(cobj, A_ub=A, b_ub=bv, A_eq=Aeq, b_eq=[1], bounds=bounds, method='highs')
    return r


def top_active(lam=None, keep=4000):
    act = {'p': {}, 't': {}}
    for st in pkeys:
        V = PS[st][1]; s = V.max(1) if lam is None else V @ lam.astype(np.float32)
        act['p'][st] = np.argsort(s)[-keep:]
    for cb in tkeys:
        V = TS[cb][1]; s = V.max(1) if lam is None else V @ lam.astype(np.float32)
        act['t'][cb] = np.argsort(s)[-keep:]
    return act


def pair_values(st, X, x3v, Av):
    p, q, u = X[:, 0], X[:, 1], X[:, 2]
    Kv = np.einsum('nk,na,kab,nb->n', G.ubasis(u, G.D2), G.pbasis(p, G.R2), Av, G.pbasis(q, G.R2))
    B = T3.Builder(L); T3.add_pair(B, np.arange(len(u)), ti[PAR[st[0]]], ti[PAR[st[1]]], u)
    return Kv - G.pair(p / 2, q / 2, u) + B.matrix(len(u)) @ x3v


def triple_values(cb, X, x3v):
    B = T3.Builder(L); T3.add_triple(B, np.arange(len(X)), tuple(ti[PAR[x]] for x in cb), X[:, 0], X[:, 1], X[:, 2])
    return B.matrix(len(X)) @ x3v


def climb_max(lam, x3, Av):
    """local maxima of the combination from its best samples: compass search"""
    newp, newt, cp, ct = {}, {}, {}, {}
    for st in pkeys:
        X, V = PS[st]; s = (V @ lam.astype(np.float32)).astype(float); st_ = X[np.argsort(s)[-CLIMB_N:]]
        f = lambda Y: pair_values(st, Y, x3, Av)
        Y, val = G.climb(f, st_, G.pair_proj(st), lambda Y: np.ones(len(Y), bool))
        newp[st] = Y; cp[st] = max(val.max(), s.max())
    for cb in tkeys:
        X, V = TS[cb]; s = (V @ lam.astype(np.float32)).astype(float); st_ = X[np.argsort(s)[-CLIMB_N:]]
        f = lambda Y: triple_values(cb, Y, x3)
        Y, val = G.climb(f, st_, G.triple_proj(cb), G.gram_ok)
        newt[cb] = Y; ct[cb] = max(val.max(), s.max())
    return newp, newt, cp, ct


def bound_with(lam, c2, c3, m, vs):
    best = (-1, None)
    for v in vs:
        c = tc(v)
        val = lam @ (TV / 2) + sum(v[b] * (m[b] + PT[ti[BT[b]]] @ lam) for b in range(nb)) \
            + sum(npair(c, st) * c2[st] for st in pkeys) + sum(ntrip(c, cb) * c3[cb] for cb in tkeys)
        best = max(best, (val, v))
    return best


MG2 = {st: float(os.environ.get('MG2', '1e-6')) for st in pkeys}
if ('A', 'A') in MG2: MG2[('A', 'A')] = float(os.environ.get('MG_AA', '5e-6'))
MG3 = {cb: float(os.environ.get('MG3', '5e-8')) for cb in tkeys}
if ('F', 'F', 'F') in MG3: MG3[('F', 'F', 'F')] += float(os.environ.get('MG_FFF', '1e-5'))
for target in TARGETS.split(','):
    VS = [v for v in ALLV if any(_sel(t)(v) for t in target.split('+'))]
    vrows = []
    for v in VS:
        c = tc(v); r = np.zeros(NV)
        r[iL] = TV / 2 + sum(v[b] * PT[ti[BT[b]]] for b in range(nb))
        r[i2] = [npair(c, st) for st in pkeys]; r[i3] = [ntrip(c, cb) for cb in tkeys]
        r[iM] = v; r[iE] = -1
        vrows.append(r)
    vrows = np.unique(np.array(vrows), axis=0)
    print('%d count vectors in %s (%d distinct rows)' % (len(VS), target, len(vrows)), flush=True)


    act = top_active()
    for it in range(iters):
        r = lp(act)
        if r.x is None:
            print('LP failed: status %d, %s' % (r.status, r.message), flush=True)
            bad = {k: int((~np.isfinite(PS[k][1])).sum()) for k in pkeys}
            bad.update({k: int((~np.isfinite(TS[k][1])).sum()) for k in tkeys})
            print('non-finite sample values:', {k: v for k, v in bad.items() if v}, flush=True)
            mx = {str(k): float(np.abs(PS[k][1]).max()) for k in pkeys}
            print('largest |pair value| per kind:', mx, flush=True)
            print('largest |bin value|:', [float(np.abs(BV[b]).max()) for b in range(nb)], flush=True)
            break
        lam = np.clip(r.x[iL], 0, None); lam /= lam.sum()
        x3 = X3 @ lam; Av = np.einsum('i,ikab->kab', lam, AV)
        newp, newt, cp, ct = climb_max(lam, x3, Av)
        mtrue = [max((BV[b] @ lam).max(), 0) if False else (BV[b] @ lam).max() for b in range(nb)]
        est = bound_with(lam, {st: cp[st] + MG2[st] for st in pkeys}, {cb: ct[cb] + MG3[cb] for cb in tkeys},
                         [mm + 1e-7 for mm in mtrue], VS)
        gaps = {''.join(k): ct[k] - r.x[i3[j]] for j, k in enumerate(tkeys)}
        gaps.update({''.join(k): cp[k] - r.x[i2[j]] for j, k in enumerate(pkeys)})
        print('iter %d: LP %.6f; with climbed maxima and the margins of the check %.6f at %s; weights %s; largest gaps %s [%.0f s]'
              % (it, r.x[iE], est[0], est[1], {os.path.basename(certs[i]): round(float(lam[i]), 4) for i in range(K) if lam[i] > 1e-4},
                 sorted(((round(v, 9), k) for k, v in gaps.items()), reverse=True)[:4], time.time() - t0), flush=True)
        # add the climbed points and every sample now above the LP thresholds
        for st in pkeys:
            PS[st][0] = np.r_[PS[st][0], newp[st]]; PS[st][1] = np.r_[PS[st][1], pair_vals(st, newp[st])]
        for cb in tkeys:
            TS[cb][0] = np.r_[TS[cb][0], newt[cb]]; TS[cb][1] = np.r_[TS[cb][1], trip_vals(cb, newt[cb])]
        a2 = top_active(lam)
        for st in pkeys:
            a2['p'][st] = np.union1d(np.union1d(act['p'][st], a2['p'][st]), np.arange(len(PS[st][0]) - len(newp[st]), len(PS[st][0])))
        for cb in tkeys:
            a2['t'][cb] = np.union1d(np.union1d(act['t'][cb], a2['t'][cb]), np.arange(len(TS[cb][0]) - len(newt[cb]), len(TS[cb][0])))
        act = a2
        if os.environ.get('OUT'):
            mfine = [(bin_vals(np.linspace(lo, hi, 20001)) @ lam).max() + 2e-7 for _, lo, hi in G.BINS]
            np.savez(os.environ['OUT'] + '_' + target.replace('+', '_') + '.npz', lam=lam, certs=np.array(certs), x3=x3, A=Av, z=lam @ ZV, t=lam @ TV,
                     m=np.array(mfine), d3=D, bound=est[0],
                     c2=np.array([max([cp[st] for st in pkeys if tuple(sorted(PAR[s] for s in st)) == tuple(sorted(k)) ] or [np.nan])
                                  for k in itertools.combinations_with_replacement(TYPES0, 2)]),
                     c3=np.array([max([ct[cb] for cb in tkeys if tuple(sorted(PAR[s] for s in cb)) == tuple(sorted(k))] or [np.nan])
                                  for k in itertools.combinations_with_replacement(TYPES0, 3)]),
                     sub=np.array(SUB), c2sub=np.array([cp[st] for st in pkeys]), c3sub=np.array([ct[cb] for cb in tkeys]),
                     pkeys=np.array([''.join(st) for st in pkeys]), tkeys=np.array(['-'.join(cb) for cb in tkeys]))
