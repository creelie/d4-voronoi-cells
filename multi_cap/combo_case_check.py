#!/usr/bin/env python3
"""
combo_case_check.py -- combo30_check.py for a case given in a JSON file, as written for
gap_closure/CM/combo_gen2.py: statement (C) at M centres, with the two-point kernel
labelled by distance and the three-point kernel on directions typed by distance range.

The case file gives M, the types and their distance ranges, the bins (each inside the
range of one type) and constraints [bins, lo, hi] on the number of centres in groups
of bins, each of which must be a proved count (prop:C-radial, thm:kissing-stable).
The check is that of combo30_check.py: exact positivity, the pair inequalities by
tensor Bernstein bounds, the bin brackets, the triple inequalities by Taylor forms,
and then, for every integer vector of bin counts that the constraints allow, the bound
    sum_b n_b (m_b + p_type(b)) + t/2 + sum N_st c2_st + sum N_str c3_str
in exact arithmetic, the largest of which must lie below 9 pi^2/8 - 8.

Usage: python3 combo_case_check.py case.json cert.npz d3 [margin2 margin3 marginm]

With PRECHECK=1 in the environment it stops after the thresholds: it samples and
refines every inequality as the check does, skips the branch and bounds, and prints
the bound over the count vectors that those thresholds would give.  That is floating
point and proves nothing; it tells whether a certificate is worth the full check.

With BORDER=1 the matrix [[A_0, z], [z^T, t]] is rounded to an exact positive definite
matrix as a whole, eigenvalues clipped at 0 and shifted as for the other blocks, instead
of rounding A_0 and z and setting t = z^T A_0^{-1} z rounded up; when A_0 is nearly
singular this keeps t near the solver's value.  Either way the bordered matrix is
checked positive semidefinite by exact LDL^T.

With LEVEL=<decimal> the final comparison is with that number instead of 9 pi^2/8 - 8:
the check then proves the density level U(Y) < LEVEL for every packing set of the case,
whose constraints hold whenever U(Y) >= 9 pi^2/8 - 8, so that U(Y) < max(9 pi^2/8 - 8,
LEVEL) for every packing set of M centres.

With TRIPLE_EXTRA=<kind>:<x>[,<kind>:<x>...], e.g. TRIPLE_EXTRA=FFF:1e-5, the threshold
c3 of each named triple kind is raised by x beyond the sampled and refined maximum plus
margin3.  Sampling and SLSQP can stop below a narrow peak of a triple polynomial, and the
branch and bound then fails on that kind alone; the extra room lets it close while every
other threshold stays as it was.  The final bound uses the raised thresholds.

With ONLY=<kind>[,<kind>...], e.g. ONLY=FFF, the thresholds are computed exactly as in
the full check (the sampling sequence is the same), but only the branch and bound of the
named triple kinds runs: the pair inequalities, the bins and the other triple kinds are
skipped.  It ends with a line PARTIAL and proves the named inequalities only.

With EXTRA_CASES=<case.json>[,<case.json>...] the bound with the same thresholds is
also evaluated, in exact arithmetic, over the count vectors of each listed case file
(same M, types and bins, other constraints) and compared with the target.  The pair
and triple inequalities, the bins and the positivity do not depend on the count
vectors, so a passing full check of the main case together with these comparisons
proves the listed cases as well.  The thresholds depend on the main case file only
through the triple kinds it allows (the sampling sequence is fixed), so with
PRECHECK=1 and the main case and arguments of a passing full check the run reproduces
that check's thresholds, printed for comparison, and settles the listed cases at the
cost of the sampling alone.

With REFINE_DIVERSE=<n> the local ascent that refines each triple threshold starts from up
to n samples spread over the domain (taken by decreasing value, each at max-norm distance
at least 0.03 from those before) instead of the 40 largest, which can all lie in the
basin of one peak while a higher one elsewhere is missed (the split check of the
28-centre vector (23, 0, 0, 0, 0, 1, 4) failed on such an interior peak of TRIPLE3_AAF).
Only the floating-point thresholds change; the branch and bound proves them either way.
That peak is a narrow ridge: it falls by 3e-6 within 0.01 of its top, so spread-out
samples miss its basin as well.  With REFINE_GRID=<n> the ascent also starts from every
local maximum of an n^3 grid on the domain of each triple inequality; n = 121 (spacing
about 0.013) finds that peak, n = 41 does not, at about two minutes per inequality.

With REFINE_TRIPLES=<b>-<b'>-<b''>[,...] (bins numbered from 1, e.g. 1-3-7) each listed
combination of bins gets a threshold of its own: a triple of centres in the bins b, b',
b'' has inner products at most a(hi_b, hi_b'), a(hi_b, hi_b''), a(hi_b', hi_b''), so
TRIPLE3 of its kind is sampled, refined and checked by branch and bound on that smaller
set, and in the bound the triples of those bins are counted with the smaller threshold
instead of that of the kind.  These thresholds are sampled with a generator of their
own, so the other thresholds are unchanged.  TRIPLE_EXTRA takes such a combination as
well (e.g. 1-3-7:1e-6).  REFINE_PAIRS=<b>-<b'>[,...] does the same for the pair
inequality on the pairs with distances in the bins b and b', on that rectangle of
distances (sampled with a generator of their own as well).  With ONLY=NONE the branch and bounds of the kinds are skipped
and only those of the listed combinations run; the run then ends PARTIAL, and it proves
the case together with a passing full check of the same case, certificate, arguments
and options, whose thresholds it reproduces.

With PAIR_MARGIN=<kind>:<x>[,...], e.g. PAIR_MARGIN=AA:5e-6, the pair inequalities of the
named kinds (and the pairs of REFINE_PAIRS of those kinds) take the margin x in place of
margin2.  A kind whose function is nearly flat on a whole face of its domain needs many
more boxes at a small margin; this keeps its margin while the others shrink.

With SKIP_COVERED=1, a pair kind or a triple kind all of whose pairs or triples, in every
count vector of the case, lie in combinations of bins listed in REFINE_PAIRS or
REFINE_TRIPLES, is counted in the bound with the thresholds of those combinations alone:
the threshold of the kind is raised to the largest of them where it lies below it, so
that the bound subtracts the full difference, it is not used, and its branch and bound is
skipped; the listed combinations are checked as usual.  This is for certificates built for narrower distance ranges than the
types of the case, whose thresholds on the whole types are large and of no use.  The
count vectors of EXTRA_CASES must then be covered in the same way, which is checked.

With REGION_TRIPLES=<b>-<b>-<p>:<s1>:<s2>[,...] (a combination of REFINE_TRIPLES with two
centres in one bin b and the third, the pivot, in another bin p, e.g. 1-1-7:-0.7944:0.58) the
triples of that combination are split by the inner products u, u' of the directions of its
two centres of bin b with that of its pivot.  The polynomial is checked to be symmetric in u
and u' (exact coefficients), and two more thresholds are sampled, refined and proved by
branch and bound, on sub-boxes of the domain: the middle one with u in [s1, s2] (by the
symmetry, it bounds every triple with u or u' in [s1, s2]) and the mixed one with u <= s1 and
u' >= s2.  Two counts hold for every pivot y, with c0 = a(hi_b, hi_b) the largest inner
product of two directions of bin b and T = a(hi_b, hi_p) that of one of bin b with one of bin
p, when s1 < 0 < s2 < T, s1^2 > c0 and s2^2 < c0:
  at most Ka = floor((1 - c0)/(s1^2 - c0)) directions x of bin b have <x, y> <= s1, since K
  of them have a sum S with |S|^2 <= K + K(K-1) c0 and <S, -y> >= K |s1|;
  at most Kt have <x, y> >= s2, with Kt the largest K such that (1 + g)/2 >= (1 - 2/K)^2,
  g = (c0 - s2^2)/(1 - T^2): their projections to the unit sphere of the hyperplane
  orthogonal to y have pairwise inner products at most g, so the caps of angular radius
  arccos(g)/2 about them are disjoint and K (1 - sqrt((1 + g)/2)) <= 2.
Of the C(n_b, 2) pairs of centres of bin b, for a pivot with ka centres at <= s1 and kt at
>= s2, C(ka, 2) + C(kt, 2) take the threshold of the bins, ka kt the mixed one and the rest
the middle one; the bound takes the largest such sum over ka <= Ka, kt <= Kt.  TRIPLE_EXTRA
takes <b>-<b>-<p>/mid and <b>-<b>-<p>/mix for the two new thresholds.

With SKIP_EMPTY=1 a kind of pairs that occurs in no count vector of the case gets no
threshold and no branch and bound, as the kinds of triples never do; its count is 0 in
every bound.  This is for case files with one type per bin (combo_mixtype_check.py),
where the bins no count vector occupies are types of their own.

With CHECKPOINT=<file> the result of each branch and bound on pairs or triples is
appended to the file, and printed with a line "done:", as soon as it is known, keyed by
the SHA-256 of the pickled function and arguments (polynomial, threshold, domain,
symmetry).  A later run with the same file takes the result of every job whose key is
there and passed from the file instead of running it again, and says so in its "[ok]"
line, with the time of the run that proved it; the logs of the runs together are then
the record of the check.  A job whose arguments differ in any way gets another key and
runs again.  This lets a check that a restart of the machine interrupted resume.  The
float sampling that sets each triple threshold is stored the same way in <file>.setup,
keyed by its polynomial, domain and samples (drawn before the lookup, so that the random
generators advance as in a run without the store), and a resumed run takes it from there.
"""
import itertools
import json
import os
import sys
import time
from fractions import Fraction as Fr
from math import comb

import numpy as np
from flint import arb

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import combo30_check as C  # noqa: E402
from radial_count_check import A_, ldl_psd, solve_exact  # noqa: E402
from radial_case_check import bracket_check  # noqa: E402
from radial_count_sdp import C1, C2  # noqa: E402
from certificate_check import ldl_positive  # noqa: E402
from truncated_search import pair as pair_float  # noqa: E402

CHECKPOINT = os.environ.get('CHECKPOINT')


def _job_key(fn, args):
    import hashlib
    import pickle
    return hashlib.sha256(pickle.dumps((fn.__module__, fn.__name__, args), protocol=4)).hexdigest()


def _run_one(item):
    i, fn, args = item
    return i, fn(*args)


_setup_store = None
_setup_hits = [0]


def setup_cached(tag, parts, g, compute):
    """compute(), or with CHECKPOINT=<file> the value stored for the same parts and samples g in
    <file>.setup by an earlier run (the sampling of a threshold is deterministic given its samples,
    which are drawn before this is called, so the generators advance as without the store)"""
    global _setup_store
    if not CHECKPOINT:
        return compute()
    import base64
    import hashlib
    import pickle
    h = hashlib.sha256(pickle.dumps((tag, parts), protocol=4))
    h.update(np.ascontiguousarray(g).tobytes())
    key, path = h.hexdigest(), CHECKPOINT + '.setup'
    if _setup_store is None:
        _setup_store = {}
        if os.path.exists(path):
            for line in open(path):
                try:
                    e = json.loads(line)
                except ValueError:
                    continue
                _setup_store[e['key']] = e['value']
    if key in _setup_store:
        _setup_hits[0] += 1
        return pickle.loads(base64.b64decode(_setup_store[key]))
    r = compute()
    v = base64.b64encode(pickle.dumps(r, protocol=4)).decode()
    with open(path, 'a') as f:
        f.write(json.dumps({'key': key, 'value': v}) + '\n')
        f.flush()
        os.fsync(f.fileno())
    _setup_store[key] = v
    return r


def run_jobs(fn, arglist, labels, procs):
    """fn(*a) for each a of arglist, in order, on a pool of procs workers; with CHECKPOINT,
    each result is appended to that file as soon as it is known, and a result that a
    previous run stored under the same key, and that passed, is taken from the file"""
    from multiprocessing import Pool
    keys = [_job_key(fn, a) for a in arglist] if CHECKPOINT else [None] * len(arglist)
    stored = {}
    if CHECKPOINT and os.path.exists(CHECKPOINT):
        for line in open(CHECKPOINT):
            try:
                e = json.loads(line)
            except ValueError:   # a line cut short by the interruption
                continue
            if e.get('result') and e['result'][0] is True:
                stored[e['key']] = e
    res, todo = [None] * len(arglist), []
    for i in range(len(arglist)):
        e = stored.get(keys[i])
        if e is not None:
            r = list(e['result'])
            r[1] = '%s [from the checkpoint %s, proved by the run of %s]' % (r[1], os.path.basename(CHECKPOINT), e['when'])
            res[i] = tuple(r)
        else:
            todo.append(i)
    if todo:
        with Pool(min(procs, len(todo))) as pool:
            for i, r in pool.imap_unordered(_run_one, [(i, fn, arglist[i]) for i in todo]):
                res[i] = r
                print('  done: %s: %s %s' % (labels[i], 'ok' if r[0] else 'FAILED', r[1]), flush=True)
                if CHECKPOINT:
                    with open(CHECKPOINT, 'a') as f:
                        f.write(json.dumps({'key': keys[i], 'label': labels[i], 'result': list(r),
                                            'when': time.strftime('%Y-%m-%d %H:%M:%S UTC', time.gmtime())}) + '\n')
                        f.flush()
                        os.fsync(f.fileno())
    return res

CASE = json.load(open(sys.argv[1]))
PRE = os.environ.get('PRECHECK') == '1'
BORDER = os.environ.get('BORDER') == '1'
EXTRA = {k: float(v) for k, v in (e.split(':') for e in os.environ.get('TRIPLE_EXTRA', '').split(',') if e)}
ONLY = [e for e in os.environ.get('ONLY', '').split(',') if e]
REFINE = [tuple(sorted(int(x) - 1 for x in e.split('-'))) for e in os.environ.get('REFINE_TRIPLES', '').split(',') if e]
REFINE_PAIRS = [tuple(sorted(int(x) - 1 for x in e.split('-'))) for e in os.environ.get('REFINE_PAIRS', '').split(',') if e]
PAIR_MARGIN = {k: float(v) for k, v in (e.split(':') for e in os.environ.get('PAIR_MARGIN', '').split(',') if e)}
SKIP_COVERED = os.environ.get('SKIP_COVERED') == '1'
SKIP_EMPTY = os.environ.get('SKIP_EMPTY') == '1'
DIVERSE = int(os.environ.get('REFINE_DIVERSE', '0'))
GRIDN = int(os.environ.get('REFINE_GRID', '0'))
REGION = [(tuple(sorted(int(x) - 1 for x in e.split(':')[0].split('-'))), Fr(e.split(':')[1]), Fr(e.split(':')[2]))
          for e in os.environ.get('REGION_TRIPLES', '').split(',') if e]


def num(x):
    return C.DMAX if x == 'sqrt6' else Fr(str(x))


def counts():
    M, nb = CASE['M'], len(CASE['bins'])
    out = []

    def rec(prefix, left):
        if len(prefix) == nb - 1:
            v = prefix + [left]
            if all(lo <= sum(v[i] for i in bs) <= hi for bs, lo, hi in CASE['constraints']):
                out.append(tuple(v))
            return
        for k in range(left + 1):
            rec(prefix + [k], left - k)
    rec([], M)
    return out


C.TYPES = list(CASE['types'])
C.TRANGE = {k: (num(a), num(b)) for k, (a, b) in CASE['trange'].items()}
C.EDGES = [num(CASE['bins'][0][1])] + [num(b[2]) for b in CASE['bins']]
C.BINTYPE = [b[0] for b in CASE['bins']]
C.COUNTS = counts()
for i, (ty, lo, hi) in enumerate(CASE['bins']):
    assert num(lo) == C.EDGES[i] and C.TRANGE[ty][0] <= num(lo) and num(hi) <= C.TRANGE[ty][1], CASE['bins'][i]


def tcounts(nb):
    out = {s: 0 for s in C.TYPES}
    for ty, n in zip(C.BINTYPE, nb):
        out[ty] += n
    return out


C.tcounts = tcounts


def Nbins(nb, bs):
    """number of pairs or triples of centres with one centre in each listed bin"""
    out = 1
    for b in set(bs):
        out *= comb(nb[b], bs.count(b))
    return out


def covered2(st, c2r, vecs):
    rs = [bs for bs in c2r if (C.BINTYPE[bs[0]], C.BINTYPE[bs[1]]) == st]
    return bool(rs) and all(C.Npair(tcounts(nb), *st) == sum(Nbins(nb, bs) for bs in rs) for nb in vecs), rs


def covered3(combo, c3r, vecs):
    rs = [bs for bs in c3r if tuple(C.BINTYPE[b] for b in bs) == combo]
    return bool(rs) and all(C.Ntriple(tcounts(nb), list(combo)) == sum(Nbins(nb, bs) for bs in rs) for nb in vecs), rs


def main():
    t0 = time.time()
    path, d3 = sys.argv[2], int(sys.argv[3])
    mg2 = float(sys.argv[4]) if len(sys.argv) > 4 else 2e-5
    mg3 = float(sys.argv[5]) if len(sys.argv) > 5 else 2e-6
    mgm = float(sys.argv[6]) if len(sys.argv) > 6 else 2e-6
    Z = np.load(path)
    rng = np.random.default_rng(7)
    if os.environ.get('LEVEL'):
        target, tname = A_(Fr(os.environ['LEVEL'])), os.environ['LEVEL']
    else:
        target, tname = 9 * arb.pi() ** 2 / 8 - 8, '9 pi^2/8 - 8'
    D2, R2 = C.D2, C.R2
    print('statement (C) at %d centres, case %s: combined certificate %s, two-point degree %d/%d, three-point degree %d; '
          '%d count vectors' % (CASE['M'], os.path.basename(sys.argv[1]), os.path.basename(path), D2, R2, d3, len(C.COUNTS)),
          flush=True)
    C.check('dmax exceeds sqrt 6', C.DMAX ** 2 > 6)
    # (a) positivity
    Af = np.array(Z['A'], float)
    if BORDER:
        # round the bordered matrix [[A_0, z], [z^T, t]] as a whole: when A_0 is nearly
        # singular, z^T A_0^{-1} z after rounding A_0 alone can exceed the solver's t
        nA = len(Af[0])
        Bd = np.zeros((nA + 1, nA + 1))
        Bd[:nA, :nA] = Af[0]; Bd[:nA, nA] = Bd[nA, :nA] = np.array(Z['z'], float); Bd[nA, nA] = float(Z['t'])
        Zm = C.psd_exact(Bd, Fr(1, 2 ** 30))
        A0 = [row[:nA] for row in Zm[:nA]]
        z = [Zm[i][nA] for i in range(nA)]
        t = Zm[nA][nA]
    else:
        A0 = C.psd_exact(Af[0], Fr(1, 2 ** 30))
        z = [C.dyad(v) for v in np.array(Z['z'], float)]
    A = [A0] + [C.psd_exact(Af[k], Fr(1, 2 ** 30)) for k in range(1, D2 + 1)]
    good = all(ldl_psd(A[k])[0] for k in range(1, D2 + 1))
    g0, piv = ldl_psd(A0)
    if not BORDER:
        wv = solve_exact(A0, z)
        tq = sum(a * b for a, b in zip(z, wv))
        t = Fr(-(-tq.numerator * 2 ** 48 // tq.denominator), 2 ** 48)
        Zm = [row[:] + [z[i]] for i, row in enumerate(A0)] + [z + [t]]
    gz, _ = ldl_psd(Zm)
    C.check('A_1..A_D and [[A_0, z], [z^T, t]] positive semidefinite (exact LDL^T)',
            good and g0 and all(p > 0 for p in piv) and gz, 't = %.9f, float t %.9f' % (float(t), float(Z['t'])))
    Bf = C.blocks_from_x3(np.array(Z['x3'], float), d3)
    Bk = {name: C.psd_exact(M) for name, M in Bf.items()}
    C.check('three-point blocks positive definite (exact LDL^T)', all(ldl_positive(M) for M in Bk.values()),
            '%d blocks' % len(Bk))
    pt = {s: C.point_term(Bk, d3, i) for i, s in enumerate(C.TYPES)}
    print('  point terms: %s' % {s: '%.6f' % float(v) for s, v in pt.items()}, flush=True)
    # (b) pairs
    Af2 = [np.array([[float(v) for v in row] for row in a]) for a in A]

    def pair_setup(i, j, box, rng_):
        """sample and refine K + PAIR3 - Pi on the pairs with distances in box"""
        Pu = C.pair3_poly(Bk, d3, i, j)
        (a0, a1), (b0, b1) = box
        n = 400000
        p = float(a0) + float(a1 - a0) * rng_.random(n); q = float(b0) + float(b1 - b0) * rng_.random(n)
        k6 = n // 6
        p[:k6] = float(a0); q[k6:2 * k6] = float(b0); p[2 * k6:3 * k6] = float(a1); q[3 * k6:4 * k6] = float(b1)
        top = (p * p + q * q - 4) / (2 * p * q)
        u = -1 + (top + 1) * rng_.random(n) ** 0.5
        u[4 * k6:5 * k6] = top[4 * k6:5 * k6] - 3e-3 * rng_.random(k6)
        gd, ge = np.meshgrid(np.linspace(float(a0), float(a1), 25), np.linspace(float(b0), float(b1), 25), indexing='ij')
        gd, ge = gd.ravel(), ge.ravel()
        gt = (gd * gd + ge * ge - 4) / (2 * gd * ge)
        sv = np.r_[0.0, np.linspace(0, 1, 80) ** 2, 1.0]
        gp = np.repeat(gd, len(sv)); gq = np.repeat(ge, len(sv))
        gu = -1 + (np.repeat(gt, len(sv)) + 1) * np.tile(1 - sv[::-1], len(gd))
        p, q, u = np.r_[p, gp], np.r_[q, gq], np.r_[u, gu]
        Kv = np.zeros(len(u))
        x = (2 * p - float(C1)) / float(C2); y = (2 * q - float(C1)) / float(C2)
        Tx = [np.ones_like(x), x]; Ty = [np.ones_like(y), y]
        for _ in range(2, R2 + 1):
            Tx.append(2 * x * Tx[-1] - Tx[-2]); Ty.append(2 * y * Ty[-1] - Ty[-2])
        Uu = [np.ones_like(u), 2 * u]
        for _ in range(2, D2 + 1):
            Uu.append(2 * u * Uu[-1] - Uu[-2])
        for k in range(D2 + 1):
            Kv += Uu[k] / (k + 1) * np.einsum('na,ab,nb->n', np.stack(Tx, 1), Af2[k], np.stack(Ty, 1))
        P3 = np.polyval([float(c) for c in Pu[::-1]], u)
        v = Kv + P3 - pair_float(p / 2, q / 2, u)
        w = np.argsort(v)[-40:]
        best = C.refine_pair(A, Pu, np.stack([p[w], q[w], u[w]], 1), box)
        return C.above(max(float(v.max()), best) + PAIR_MARGIN.get(C.TYPES[i] + C.TYPES[j], mg2)), float(v.max()), Pu

    c2 = {}
    for kidx, (s, tt) in enumerate(itertools.combinations_with_replacement(C.TYPES, 2)):
        if SKIP_EMPTY and all(C.Npair(tcounts(nb), s, tt) == 0 for nb in C.COUNTS):
            continue
        i, j = C.TYPES.index(s), C.TYPES.index(tt)
        c2[(s, tt)], vmax, _ = pair_setup(i, j, (C.TRANGE[s], C.TRANGE[tt]), rng)
        print('  pair %s%s: float largest K + PAIR3 - Pi %.6e (file c2 %.6e); threshold %.6e'
              % (s, tt, vmax, float(Z['c2'][kidx]), float(c2[(s, tt)])), flush=True)
    jobs = [(A, C.pair3_poly(Bk, d3, C.TYPES.index(s), C.TYPES.index(tt)), c2[(s, tt)], (C.TRANGE[s], C.TRANGE[tt]), s == tt)
            for (s, tt) in c2]
    labels = ['K + PAIR3_%s%s <= Pi + c2 on the admissible pairs' % st for st in c2]
    if PRE or ONLY:
        jobs, labels = [], []
    # pairs refined by bins, with a generator of their own
    c2r = {}
    prng = np.random.default_rng(13)
    for bs in REFINE_PAIRS:
        s, tt = C.BINTYPE[bs[0]], C.BINTYPE[bs[1]]
        box = ((C.EDGES[bs[0]], C.EDGES[bs[0] + 1]), (C.EDGES[bs[1]], C.EDGES[bs[1] + 1]))
        c2r[bs], vmax, Pu = pair_setup(C.TYPES.index(s), C.TYPES.index(tt), box, prng)
        name = 'bins ' + '-'.join(str(b + 1) for b in bs)
        print('  pair %s%s on %s: float largest K + PAIR3 - Pi %.6e; threshold %.6e against %.6e for the kind'
              % (s, tt, name, vmax, float(c2r[bs]), float(c2[(s, tt)])), flush=True)
        if not PRE:
            jobs.append((A, Pu, c2r[bs], box, bs[0] == bs[1]))
            labels.append('K + PAIR3_%s%s <= Pi + its threshold on the pairs of %s' % (s, tt, name))
    skip2 = set()
    if SKIP_COVERED:
        for st in c2:
            cov, rs = covered2(st, c2r, C.COUNTS)
            if cov:
                skip2.add(st)
                c2[st] = max([c2[st]] + [c2r[bs] for bs in rs])
                print('  pair %s%s: in every count vector each of its pairs lies in a listed pair of bins, so the bound counts '
                      'them with the thresholds of those bins alone and the branch and bound of the kind is skipped'
                      % (st[0], st[1]), flush=True)
        if skip2 and not (PRE or ONLY):
            sts = list(c2)
            keep = [k for k in range(len(jobs)) if k >= len(sts) or sts[k] not in skip2]
            jobs, labels = [jobs[k] for k in keep], [labels[k] for k in keep]
    if jobs:
        res = run_jobs(C.pair_box_check, jobs, labels, int(os.environ.get('PAIR_PROCS', '3')))
        for label, (ok, msg) in zip(labels, res):
            C.check(label, ok, msg)
    # (c) bins
    mf = np.array(Z['m'], float)
    assert len(mf) == len(C.BINTYPE)
    m = [C.above(v + mgm) for v in mf]
    if not PRE and not ONLY:
        okb, msg = bracket_check(A, z, D2, R2, C1, C2, C.DMAX, C.EDGES, m)
        C.check('f <= m_b on every bin', okb, msg)
    # (d) triples
    def refine_starts(g, vals, P_GRID=None):
        """the starts of the local ascent: the 40 largest samples, or with REFINE_DIVERSE=<n> up to n
        samples taken in decreasing order of value and at max-norm distance at least 0.03 from those
        already taken, so that a peak whose basin holds none of the 40 largest samples is still climbed"""
        order = np.argsort(vals)[::-1]
        if not DIVERSE:
            return g[order[:40]] if P_GRID is None else np.r_[g[order[:40]], grid_local_maxima(*P_GRID)]
        taken = []
        for i in order:
            if all(np.max(np.abs(g[i] - g[j])) >= 0.03 for j in taken):
                taken.append(i)
                if len(taken) >= DIVERSE:
                    break
        out = g[np.array(taken)]
        return out if P_GRID is None else np.r_[out, grid_local_maxima(*P_GRID)]

    def grid_local_maxima(P, T, n, B=(-1, -1, -1)):
        """the points of an n^3 grid on [B_1, T_1] x [B_2, T_2] x [B_3, T_3] (B_i = -1 unless given),
        restricted to the admissible triples, at which P is at least its value at every admissible
        neighbour (26-neighbourhood)"""
        ax = [np.linspace(float(b), float(t), n) for b, t in zip(B, T)]
        X = np.stack(np.meshgrid(*ax, indexing='ij'), -1).reshape(-1, 3)
        ok = 1 + 2 * X[:, 0] * X[:, 1] * X[:, 2] - (X ** 2).sum(1) >= 0
        V = np.full(len(X), -np.inf)
        idx = np.where(ok)[0]
        for s0 in range(0, len(idx), 20000):
            j = idx[s0:s0 + 20000]
            V[j] = C.peval(P, X[j, 0], X[j, 1], X[j, 2])
        V3 = V.reshape(n, n, n)
        pad = np.pad(V3, 1, constant_values=-np.inf)
        loc = np.isfinite(V3)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    if dx or dy or dz:
                        loc &= V3 >= pad[1 + dx:n + 1 + dx, 1 + dy:n + 1 + dy, 1 + dz:n + 1 + dz]
        return X[loc.reshape(-1)]

    def triple_setup(P, T, labels, rng_, extra):
        """sample and refine TRIPLE3 = P on {u_12 <= T[0], u_13 <= T[1], u_23 <= T[2]}, and lay it
        out for the branch and bound with the symmetry of equal labels"""
        T12, T13, T23 = T
        g = C.random_gram(float(T12), float(T13), float(T23), 200000, rng_)
        ax = [np.r_[-1.0, np.linspace(-1, float(T), 40), float(T)] for T in (T12, T13, T23)]
        G = np.stack(np.meshgrid(*ax, indexing='ij'), -1).reshape(-1, 3)
        G = G[1 + 2 * G[:, 0] * G[:, 1] * G[:, 2] - (G ** 2).sum(1) >= 0]
        g = np.r_[g, G]

        def compute():
            vals = C.peval(P, g[:, 0], g[:, 1], g[:, 2])
            P_GRID = (P, (T12, T13, T23), GRIDN) if GRIDN else None
            best = C.refine_triple(P, refine_starts(g, vals, P_GRID), (float(T12), float(T13), float(T23)))
            thr = C.above(max(float(vals.max()), best) + mg3 + extra)
            perm, sym = C.triple_layout(labels)
            tops = [None] * 3
            for src, val in enumerate((T12, T13, T23)):
                tops[perm.index(src)] = float(C.above(val))
            Pr = {}
            for e, co in P.items():
                ee = [0, 0, 0]
                for src in range(3):
                    ee[perm.index(src)] += e[src]
                Pr[tuple(ee)] = Pr.get(tuple(ee), 0) + co
            return thr, Pr, tops, sym, float(vals.max()), len(g), best
        return setup_cached('triple', (sorted(P.items()), [str(x) for x in T], tuple(labels), extra, mg3, GRIDN,
                                       DIVERSE), g, compute)

    c3 = {}
    tjobs = []
    for combo in itertools.combinations_with_replacement(C.TYPES, 3):
        if all(C.Ntriple(tcounts(nb), list(combo)) == 0 for nb in C.COUNTS):
            continue
        P = C.triple3_poly(Bk, d3, combo)
        T = (C.tmax(combo[0], combo[1]), C.tmax(combo[0], combo[2]), C.tmax(combo[1], combo[2]))
        c3[combo], Pr, tops, sym, vmax, ns, best = triple_setup(P, T, combo, rng, EXTRA.get(''.join(combo), 0.0))
        print('  triple %s: float largest %.6e on %d samples; threshold %.6e%s; %d monomials, symmetry %d'
              % (''.join(combo), vmax, ns, float(c3[combo]),
                 ' (refined %.6e, extra %g)' % (best, EXTRA[''.join(combo)]) if ''.join(combo) in EXTRA else '',
                 len(P), sym), flush=True)
        if ONLY and ''.join(combo) not in ONLY:
            continue
        tjobs.append((''.join(combo), (Pr, tops, sym, float(c3[combo]), ''.join(combo))))
    # triples refined by bins, sampled with a generator of their own so that the
    # thresholds above stay those of the check without refinement
    c3r = {}
    rrng = np.random.default_rng(11)
    for bs in REFINE:
        combo = tuple(C.BINTYPE[b] for b in bs)
        assert combo in c3, bs
        hi = [C.EDGES[b + 1] for b in bs]
        T = (C.amaxq(hi[0], hi[1]), C.amaxq(hi[0], hi[2]), C.amaxq(hi[1], hi[2]))
        name = 'bins ' + '-'.join(str(b + 1) for b in bs)
        c3r[bs], Pr, tops, sym, vmax, ns, best = triple_setup(C.triple3_poly(Bk, d3, combo), T, bs, rrng,
                                                               EXTRA.get('-'.join(str(b + 1) for b in bs), 0.0))
        print('  triple %s on %s (inner products at most %.6f, %.6f, %.6f): float largest %.6e on %d samples; '
              'threshold %.6e against %.6e for the kind; symmetry %d'
              % (''.join(combo), name, *[float(x) for x in T], vmax, ns, float(c3r[bs]), float(c3[combo]), sym), flush=True)
        tjobs.append((''.join(combo) + ' on ' + name, (Pr, tops, sym, float(c3r[bs]), ''.join(combo) + ' ' + name)))
    # triples of a refined combination of bins (b, b, p), split by the inner products of the two centres
    # of bin b with the centre of bin p (REGION_TRIPLES); sampled with a generator of their own as well
    c3g = {}
    grng = np.random.default_rng(13)
    PAIRIDX = {(0, 1): 0, (0, 2): 1, (1, 2): 2}

    def below(x):
        """the largest double <= x, as a rational"""
        y = float(x)
        if Fr(y) > Fr(x):
            y = float(np.nextafter(y, -np.inf))
        return Fr(y)

    def region_setup(P, T, B, rng_, extra):
        """sample and refine TRIPLE3 = P on the sub-box B <= u <= T of its domain"""
        Tf, Bf_ = [float(x) for x in T], [float(x) for x in B]
        g = C.random_gram(*Tf, 200000, rng_)
        g = g[np.all(g >= np.array(Bf_), axis=1)]
        ax = [np.r_[np.linspace(b, t, 40), t] for b, t in zip(Bf_, Tf)]
        G = np.stack(np.meshgrid(*ax, indexing='ij'), -1).reshape(-1, 3)
        g = np.r_[g, G[1 + 2 * G[:, 0] * G[:, 1] * G[:, 2] - (G ** 2).sum(1) >= 0]]

        def compute():
            vals = C.peval(P, g[:, 0], g[:, 1], g[:, 2])
            best = C.refine_triple(P, refine_starts(g, vals, (P, Tf, GRIDN, Bf_) if GRIDN else None), Tf, Bf_)
            return C.above(max(float(vals.max()), best) + mg3 + extra), float(vals.max()), len(g), best
        return setup_cached('region', (sorted(P.items()), [str(x) for x in T], [str(x) for x in B], extra, mg3, GRIDN,
                                       DIVERSE), g, compute)

    for bs, s1, s2 in REGION:
        assert bs in c3r and bs not in c3g, bs
        i, j, k = (0, 1, 2) if bs[0] == bs[1] != bs[2] else (1, 2, 0)
        assert bs[i] == bs[j] != bs[k], bs
        xi, yi = PAIRIDX[tuple(sorted((i, k)))], PAIRIDX[tuple(sorted((j, k)))]
        b, p = bs[i], bs[k]
        combo = tuple(C.BINTYPE[x] for x in bs)
        name = 'bins ' + '-'.join(str(x + 1) for x in bs)
        P = C.triple3_poly(Bk, d3, combo)

        def swap(e):
            e = list(e)
            e[xi], e[yi] = e[yi], e[xi]
            return tuple(e)
        C.check('TRIPLE3 on %s symmetric under exchange of the two centres of bin %d (exact coefficients)' % (name, b + 1),
                all(P.get(swap(e), 0) == co for e, co in P.items()))
        hi = [C.EDGES[x + 1] for x in bs]
        T = [C.amaxq(hi[0], hi[1]), C.amaxq(hi[0], hi[2]), C.amaxq(hi[1], hi[2])]
        c0, Tx = C.amaxq(C.EDGES[b + 1], C.EDGES[b + 1]), T[xi]
        assert s1 < 0 < s2 < Tx < 1 and s1 * s1 > c0 and s2 * s2 < c0, (s1, s2, c0, Tx)
        # at most Ka directions of bin b have inner product <= s1 with a direction of bin p: K of them
        # have sum S with |S|^2 <= K + K(K-1) c0 and <S, -y> >= K |s1|, so K s1^2 <= 1 + (K-1) c0
        Ka = int((1 - c0) / (s1 * s1 - c0))
        # at most Kt have inner product >= s2: projected to the unit sphere of the hyperplane
        # orthogonal to y, two of them have inner product at most g, so the caps of angular radius
        # arccos(g)/2 about them are disjoint, and Kt (1 - sqrt((1 + g)/2)) <= 2
        g = (c0 - s2 * s2) / (1 - Tx * Tx)
        Kt = 2
        while (1 + g) / 2 >= (1 - Fr(2, Kt + 1)) ** 2:
            Kt += 1
        Bm, Tm = [Fr(-1)] * 3, list(T)
        Bm[xi], Tm[xi] = below(s1), C.above(s2)
        Bx, Tx_ = [Fr(-1)] * 3, list(T)
        Tx_[xi], Bx[yi] = C.above(s1), below(s2)
        ex = '-'.join(str(x + 1) for x in bs)
        cm, vm, nm, bm = region_setup(P, Tm, Bm, grng, EXTRA.get(ex + '/mid', 0.0))
        cx, vx, nx, bx = region_setup(P, Tx_, Bx, grng, EXTRA.get(ex + '/mix', 0.0))
        c3g[bs] = (b, p, cm, cx, Ka, Kt)
        print('  triple %s on %s split by the inner products of its two centres of bin %d with the centre of bin %d: '
              'some in [%s, %s], float largest %.6e on %d samples, threshold %.6e; one <= %s and the other >= %s, '
              'float largest %.6e on %d samples, threshold %.6e; against %.6e for the bins; with each centre of bin %d '
              'at most %d centres of bin %d at inner product <= %s and at most %d at >= %s'
              % (''.join(combo), name, b + 1, p + 1, s1, s2, max(vm, bm), nm, float(cm), s1, s2, max(vx, bx), nx,
                 float(cx), float(c3r[bs]), p + 1, Ka, b + 1, s1, Kt, s2), flush=True)
        for lab, Bq, Tq, cq in (('middle', Bm, Tm, cm), ('mixed', Bx, Tx_, cx)):
            tjobs.append((''.join(combo) + ' on ' + name + ' ' + lab,
                          (P, [float(x) for x in Tq], 0, float(cq), ''.join(combo) + ' ' + name + ' ' + lab, 1e-6, 150000,
                           [float(x) for x in Bq])))
    if _setup_hits[0]:
        print('  (the samples and thresholds of %d of these triple setups were taken from %s.setup, stored by an '
              'earlier run with the same samples)' % (_setup_hits[0], os.path.basename(CHECKPOINT)), flush=True)
    skip3 = set()
    if SKIP_COVERED:
        for combo in c3:
            cov, rs = covered3(combo, c3r, C.COUNTS)
            if cov:
                skip3.add(''.join(combo))
                c3[combo] = max([c3[combo]] + [c3r[bs] for bs in rs])
                print('  triple %s: in every count vector each of its triples lies in a listed combination of bins, so the '
                      'bound counts them with the thresholds of those combinations alone and the branch and bound of the '
                      'kind is skipped' % ''.join(combo), flush=True)
        tjobs = [j for j in tjobs if j[0] not in skip3]
    if not PRE and tjobs:
        res = run_jobs(C.triple_box_check, [a for _, a in tjobs], [lab for lab, _ in tjobs],
                       int(os.environ.get('TRIPLE_PROCS', '3')))
        nbox = 0
        for (label, _), (ok, msg, nd) in zip(tjobs, res):
            C.check('TRIPLE3_%s <= its threshold on the admissible triples' % label if ' on ' in label
                    else 'TRIPLE3_%s <= c3 on the admissible triples' % label, ok, msg)
            nbox += nd
        print('  triple boxes in all: %d' % nbox, flush=True)

    def refined_gain(nb):
        return (sum(Nbins(nb, bs) * max(c3[tuple(C.BINTYPE[b] for b in bs)] - c3r[bs], 0) for bs in c3r)
                + sum(Nbins(nb, bs) * max(c2[(C.BINTYPE[bs[0]], C.BINTYPE[bs[1]])] - c2r[bs], 0) for bs in c2r))
    def region_gain(nb):
        """for each split combination (b, b, p): every centre of bin p sees at most C(n_b, 2) pairs of
        centres of bin b, of which, with ka of them at inner product <= s1 and kt at >= s2 with it, at most
        C(ka, 2) + C(kt, 2) have both on one side (threshold of the bins), ka kt one on each side (mixed
        threshold) and the rest some in [s1, s2] (middle threshold); the largest such sum over the allowed
        ka <= Ka, kt <= Kt replaces C(n_b, 2) times the threshold of the bins"""
        out = 0
        for bs, (b, p, cm, cx, Ka, Kt) in c3g.items():
            n, ch = nb[b], c3r[bs]
            if nb[p] == 0 or n < 2:
                continue
            cm, cx, tot, best = min(cm, ch), min(cx, ch), comb(n, 2), None
            for ka in range(min(Ka, n) + 1):
                for kt in range(min(Kt, n - ka) + 1):
                    hot, mix = comb(ka, 2) + comb(kt, 2), ka * kt
                    v = hot * ch + mix * cx + (tot - hot - mix) * cm
                    best = v if best is None or v > best else best
            out += nb[p] * (tot * ch - best)
        return out
    # the bound, over every count vector
    cache = {}
    worst, wvec = None, None
    for nb in C.COUNTS:
        tc = tcounts(nb)
        key = tuple(tc[s] for s in C.TYPES)
        if key not in cache:
            cache[key] = (t / 2 + sum(C.Npair(tc, s, tt) * c2[(s, tt)] for (s, tt) in c2)
                          + sum(C.Ntriple(tc, list(cb)) * c3[cb] for cb in c3))
        val = cache[key] + sum(nb[b] * (m[b] + pt[C.BINTYPE[b]]) for b in range(len(nb)) if nb[b])
        if c3r or c2r:
            val -= refined_gain(nb)
        if c3g:
            val -= region_gain(nb)
        if worst is None or val > worst:
            worst, wvec = val, nb
    print('  %d count vectors, %d type-count vectors; largest bound %.6f at %s' % (len(C.COUNTS), len(cache), float(worst), wvec),
          flush=True)
    for extra in [e for e in os.environ.get('EXTRA_CASES', '').split(',') if e]:
        E = json.load(open(extra))
        assert (E['M'], E['types'], E['trange'], E['bins']) == (CASE['M'], CASE['types'], CASE['trange'], CASE['bins']), extra
        main_constraints, CASE['constraints'] = CASE['constraints'], E['constraints']
        vecs = counts()
        CASE['constraints'] = main_constraints
        if skip2 or skip3:
            C.check('the skipped kinds are covered by the listed bins in the count vectors of %s' % os.path.basename(extra),
                    all(covered2(st, c2r, vecs)[0] for st in skip2)
                    and all(covered3(tuple(k), c3r, vecs)[0] for k in skip3))
        ew, ev = None, None
        for nb in vecs:
            tc = tcounts(nb)
            key = tuple(tc[s] for s in C.TYPES)
            if key not in cache:
                cache[key] = (t / 2 + sum(C.Npair(tc, s, tt) * c2[(s, tt)] for (s, tt) in c2)
                              + sum(C.Ntriple(tc, list(cb)) * c3[cb] for cb in c3))
            val = cache[key] + sum(nb[b] * (m[b] + pt[C.BINTYPE[b]]) for b in range(len(nb)) if nb[b])
            if c3r or c2r:
                val -= refined_gain(nb)
            if c3g:
                val -= region_gain(nb)
            if ew is None or val > ew:
                ew, ev = val, nb
        C.check('with these thresholds, largest bound over the %d count vectors of %s below %s'
                % (len(vecs), os.path.basename(extra), tname), A_(ew) < target,
                '%.6f at %s < %s' % (float(ew), ev, target.str(8)))
    if ONLY and not PRE:
        kinds = ([k for k in ONLY if k in [''.join(c) for c in c3]] + [lab for lab, _ in tjobs if ' on ' in lab]
                 + ['the pairs on bins ' + '-'.join(str(b + 1) for b in bs) for bs in c2r])
        print('PARTIAL: only the inequalities %s were checked by branch and bound; with the thresholds above the bound would be '
              '%.6f against %s = %s [%.0f s]' % (', '.join(kinds) or 'none', float(worst), tname, target.str(8), time.time() - t0))
        return
    if PRE:
        print('PRECHECK (floating point, no branch and bound): the thresholds give %.6f against %s = %s [%.0f s]'
              % (float(worst), tname, target.str(8), time.time() - t0))
        return
    C.check('largest bound over the count vectors below %s' % tname, A_(worst) < target,
            '%.6f < %s' % (float(worst), target.str(8)))
    if os.environ.get('LEVEL'):
        print('PASS: in case %s every packing set of %d centres has U(Y) <= %.6f < %s, so U(Y) < max(9 pi^2/8 - 8, %s) '
              'for every packing set of %d centres [%.0f s]'
              % (os.path.basename(sys.argv[1]), CASE['M'], float(worst), tname, tname, CASE['M'], time.time() - t0))
    else:
        print('PASS: in case %s every packing set of %d centres has U(Y) <= %.6f < 9 pi^2/8 - 8, so T(Y) > 8 [%.0f s]'
              % (os.path.basename(sys.argv[1]), CASE['M'], float(worst), time.time() - t0))


if __name__ == '__main__':
    main()
