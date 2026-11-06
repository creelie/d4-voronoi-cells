"""c27_hull_mt.py -- c27_hull_sub.py with every bin a type of its own and each certificate
used under one or more typings of the bins (mixed typing).

The three-point part of a certificate is positive for every assignment of its types to the
points, so a certificate may type the centres of a bin as any of its types A, B, F; only the
domains on which its pair and triple polynomials are bounded change.  A convex combination
of (certificate, typing) columns is therefore a certificate: per bin the point term is the
weighted sum of the point terms of the type each column gives the bin, per pair or triple of
bins the polynomial is the weighted sum of the polynomials of the types each column gives
those bins, and the kernel of the distances (A, z, t) is shared.  Thresholds are one per
pair or triple of bins.  Check with multi_cap/combo_mixtype_check.py on a case file with one
type per bin (case_bins.py).

usage: python3 c27_hull_mt.py TARGET ITERS cert.npz[@T1[+T2...]] ...
TYPINGS=T1,T2,... (env): the typings each certificate is used under unless its argument
lists its own; a typing is a string of one type letter (A, B, F) per bin, e.g. AABBBFF
types the bin (2.2, 2.35] as F.  Default: the typing of the case file.
Other env as in c27_hull_sub.py (CASE, NP, NT, GRID, KEEP, CLIMB_N, PEAKS, FACE, OUT, MG_*).
The npz written holds lam, the columns, x3g (one combined three-point vector per typing),
gtyp (the type index of each bin in each typing), A, z, t, m, and c2, c3 per pair and
triple of bins.  Float throughout; proves nothing.
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
SBT = ['%s%d' % (BT[b], b + 1) for b in range(nb)]   # every bin a type of its own
BINOF = {s: b for b, s in enumerate(SBT)}
for b in range(nb):
    G.TRANGE[SBT[b]] = (G.BINS[b][1], G.BINS[b][2])
TYPES = list(SBT)
TYPINGS = [t for t in os.environ.get('TYPINGS', ''.join(BT)).split(',') if t]
SPEC = []
for c in certs:
    f_, _, ty_ = c.partition('@')
    SPEC.append((f_, ty_.split('+') if ty_ else TYPINGS))
TYP = list(dict.fromkeys(t for _, tys in SPEC for t in tys))
for t in TYP:
    assert len(t) == nb and set(t) <= set(TYPES0), t
COLS = [(k, TYP.index(t)) for k, (_, tys) in enumerate(SPEC) for t in tys]
NC = len(COLS)
CK = np.array([k for k, _ in COLS]); CT = np.array([g for _, g in COLS])
TC = [(np.where(CT == g)[0], CK[CT == g]) for g in range(len(TYP))]   # per typing: its columns and their certificates
DEG = {2052: 6, 4098: 8, 7188: 10, 11538: 12}
Zs = [np.load(f_ if f_.startswith('/') else CM + f_) for f_, _ in SPEC]
D = max(DEG[len(Z['x3'])] for Z in Zs)
K = len(Zs)
CNAME = ['%s@%s' % (os.path.basename(SPEC[k][0]), TYP[g]) for k, g in COLS]


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
PTC = np.array([[PT[ti[TYP[g][b]], k] for k, g in COLS] for b in range(nb)])                # nb x NC
TVC, AVC, ZVC = TV[CK], AV[CK], ZV[CK]
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
        Kv = np.einsum('nk,na,ikab,nb->ni', U, Pp, AV, Pq, optimize=True) - G.pair(p / 2, q / 2, u)[:, None]
        o = np.empty((len(u), NC), np.float32)
        for g, (cols, ks) in enumerate(TC):
            B = T3.Builder(L); T3.add_pair(B, np.arange(len(u)), ti[TYP[g][BINOF[st[0]]]], ti[TYP[g][BINOF[st[1]]]], u)
            o[:, cols] = Kv[:, ks] + B.matrix(len(u)) @ X3[:, ks]
        out.append(o)
    return np.concatenate(out)


def trip_vals(cb, X):
    out = []
    for i0 in range(0, len(X), 10000):
        x = X[i0:i0 + 10000]
        o = np.empty((len(x), NC), np.float32)
        for g, (cols, ks) in enumerate(TC):
            B = T3.Builder(L)
            T3.add_triple(B, np.arange(len(x)), tuple(ti[TYP[g][BINOF[c]]] for c in cb), x[:, 0], x[:, 1], x[:, 2])
            o[:, cols] = B.matrix(len(x)) @ X3[:, ks]
        out.append(o)
    return np.concatenate(out)


def bin_vals(d):
    Pd = G.pbasis(d, G.R2); U1 = G.ubasis(np.ones(1), G.D2)[0]
    Kv = np.einsum('k,na,ikab,nb->ni', U1, Pd, AV, Pd)
    return (G.S(d)[:, None] + 0.5 * Kv - Pd @ ZV.T)[:, CK]


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
iL = np.arange(NC); i2 = NC + np.arange(n2); i3 = NC + n2 + np.arange(n3); iM = NC + n2 + n3 + np.arange(nb); iE = NC + n2 + n3 + nb
NV = iE + 1
def lp(active):
    rows, b = [vrows], [np.zeros(len(vrows))]
    def block(V, col):
        n = len(V)
        R = sp.hstack([sp.csr_matrix(V.astype(float)), sp.csr_matrix((n, NV - NC))], format='lil')
        R[:, col] = -1
        return R.tocsr()
    rows = [sp.csr_matrix(vrows)]
    for j, st in enumerate(pkeys):
        V = PS[st][1][active['p'][st]]; rows.append(block(V, i2[j])); b.append(np.zeros(len(V)))
    for j, cb in enumerate(tkeys):
        V = TS[cb][1][active['t'][cb]]; rows.append(block(V, i3[j])); b.append(np.zeros(len(V)))
    for bb in range(nb):
        V = BV[bb]; rows.append(block(V, iM[bb])); b.append(np.zeros(len(V)))
    A = sp.vstack(rows, format='csr'); bv = np.concatenate(b); rows.clear()
    Aeq = np.zeros((1, NV)); Aeq[0, iL] = 1
    cobj = np.zeros(NV); cobj[iE] = 1
    bounds = [(0, None)] * NC + [(None, None)] * (NV - NC)
    r = linprog(cobj, A_ub=A, b_ub=bv, A_eq=Aeq, b_eq=[1], bounds=bounds, method='highs')
    return r


def top_active(lam=None, keep=int(os.environ.get('KEEP', '4000'))):
    act = {'p': {}, 't': {}}
    for st in pkeys:
        V = PS[st][1]; s = V.max(1) if lam is None else V @ lam.astype(np.float32)
        act['p'][st] = np.argsort(s)[-keep:]
    for cb in tkeys:
        V = TS[cb][1]; s = V.max(1) if lam is None else V @ lam.astype(np.float32)
        act['t'][cb] = np.argsort(s)[-keep:]
    return act


def pair_values(st, X, x3g, Av):
    """x3g: the combined three-point vector of each typing"""
    p, q, u = X[:, 0], X[:, 1], X[:, 2]
    out = np.einsum('nk,na,kab,nb->n', G.ubasis(u, G.D2), G.pbasis(p, G.R2), Av, G.pbasis(q, G.R2)) - G.pair(p / 2, q / 2, u)
    for g, x3v in enumerate(x3g):
        B = T3.Builder(L); T3.add_pair(B, np.arange(len(u)), ti[TYP[g][BINOF[st[0]]]], ti[TYP[g][BINOF[st[1]]]], u)
        out = out + B.matrix(len(u)) @ x3v
    return out


def triple_values(cb, X, x3g):
    out = np.zeros(len(X))
    for g, x3v in enumerate(x3g):
        B = T3.Builder(L)
        T3.add_triple(B, np.arange(len(X)), tuple(ti[TYP[g][BINOF[x]]] for x in cb), X[:, 0], X[:, 1], X[:, 2])
        out = out + B.matrix(len(X)) @ x3v
    return out


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
        val = lam @ (TVC / 2) + sum(v[b] * (m[b] + PTC[b] @ lam) for b in range(nb)) \
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
        r[iL] = TVC / 2 + sum(v[b] * PTC[b] for b in range(nb))
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
        x3 = [X3[:, ks] @ lam[cols] for cols, ks in TC]; Av = np.einsum('i,ikab->kab', lam, AVC)
        newp, newt, cp, ct = climb_max(lam, x3, Av)
        mtrue = [max((BV[b] @ lam).max(), 0) if False else (BV[b] @ lam).max() for b in range(nb)]
        est = bound_with(lam, {st: cp[st] + MG2[st] for st in pkeys}, {cb: ct[cb] + MG3[cb] for cb in tkeys},
                         [mm + 1e-7 for mm in mtrue], VS)
        gaps = {''.join(k): ct[k] - r.x[i3[j]] for j, k in enumerate(tkeys)}
        gaps.update({''.join(k): cp[k] - r.x[i2[j]] for j, k in enumerate(pkeys)})
        print('iter %d: LP %.6f; with climbed maxima and the margins of the check %.6f at %s; weights %s; largest gaps %s [%.0f s]'
              % (it, r.x[iE], est[0], est[1], {CNAME[i]: round(float(lam[i]), 4) for i in range(NC) if lam[i] > 1e-4},
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
            np.savez(os.environ['OUT'] + '_' + target.replace('+', '_') + '.npz', lam=lam, cols=np.array(CNAME),
                     certs=np.array([f_ for f_, _ in SPEC]), typings=np.array(TYP), x3=np.zeros(1), x3g=np.array(x3),
                     gtyp=np.array([[ti[c] for c in t] for t in TYP]), A=Av, z=lam @ ZVC, t=lam @ TVC,
                     m=np.array(mfine), d3=D, bound=est[0],
                     c2=np.array([cp.get(k, np.nan) for k in itertools.combinations_with_replacement(TYPES, 2)]),
                     c3=np.array([ct.get(k, np.nan) for k in itertools.combinations_with_replacement(TYPES, 3)]),
                     pkeys=np.array([''.join(st) for st in pkeys]), tkeys=np.array(['-'.join(cb) for cb in tkeys]))
