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
"""
import itertools
import json
import os
import sys
import time
from fractions import Fraction as Fr

import numpy as np
from flint import arb

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import combo30_check as C  # noqa: E402
from radial_count_check import A_, ldl_psd, solve_exact  # noqa: E402
from radial_case_check import bracket_check  # noqa: E402
from radial_count_sdp import C1, C2  # noqa: E402
from certificate_check import ldl_positive  # noqa: E402
from truncated_search import pair as pair_float  # noqa: E402

CASE = json.load(open(sys.argv[1]))


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


def main():
    t0 = time.time()
    path, d3 = sys.argv[2], int(sys.argv[3])
    mg2 = float(sys.argv[4]) if len(sys.argv) > 4 else 2e-5
    mg3 = float(sys.argv[5]) if len(sys.argv) > 5 else 2e-6
    mgm = float(sys.argv[6]) if len(sys.argv) > 6 else 2e-6
    Z = np.load(path)
    rng = np.random.default_rng(7)
    target = 9 * arb.pi() ** 2 / 8 - 8
    D2, R2 = C.D2, C.R2
    print('statement (C) at %d centres, case %s: combined certificate %s, two-point degree %d/%d, three-point degree %d; '
          '%d count vectors' % (CASE['M'], os.path.basename(sys.argv[1]), os.path.basename(path), D2, R2, d3, len(C.COUNTS)),
          flush=True)
    C.check('dmax exceeds sqrt 6', C.DMAX ** 2 > 6)
    # (a) positivity
    Af = np.array(Z['A'], float)
    A0 = C.psd_exact(Af[0], Fr(1, 2 ** 30))
    A = [A0] + [C.psd_exact(Af[k], Fr(1, 2 ** 30)) for k in range(1, D2 + 1)]
    z = [C.dyad(v) for v in np.array(Z['z'], float)]
    good = all(ldl_psd(A[k])[0] for k in range(1, D2 + 1))
    g0, piv = ldl_psd(A0)
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
    c2 = {}
    for (s, tt) in itertools.combinations_with_replacement(C.TYPES, 2):
        i, j = C.TYPES.index(s), C.TYPES.index(tt)
        Pu = C.pair3_poly(Bk, d3, i, j)
        (a0, a1), (b0, b1) = C.TRANGE[s], C.TRANGE[tt]
        n = 400000
        p = float(a0) + float(a1 - a0) * rng.random(n); q = float(b0) + float(b1 - b0) * rng.random(n)
        k6 = n // 6
        p[:k6] = float(a0); q[k6:2 * k6] = float(b0); p[2 * k6:3 * k6] = float(a1); q[3 * k6:4 * k6] = float(b1)
        top = (p * p + q * q - 4) / (2 * p * q)
        u = -1 + (top + 1) * rng.random(n) ** 0.5
        u[4 * k6:5 * k6] = top[4 * k6:5 * k6] - 3e-3 * rng.random(k6)
        gd, ge = np.meshgrid(np.linspace(float(a0), float(a1), 25), np.linspace(float(b0), float(b1), 25), indexing='ij')
        gd, ge = gd.ravel(), ge.ravel()
        gt = (gd * gd + ge * ge - 4) / (2 * gd * ge)
        sv = np.r_[0.0, np.linspace(0, 1, 80) ** 2, 1.0]
        gp = np.repeat(gd, len(sv)); gq = np.repeat(ge, len(sv))
        gu = -1 + (np.repeat(gt, len(sv)) + 1) * np.tile(1 - sv[::-1], len(gd))
        p, q, u = np.r_[p, gp], np.r_[q, gq], np.r_[u, gu]
        Af2 = [np.array([[float(v) for v in row] for row in a]) for a in A]
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
        best = C.refine_pair(A, Pu, np.stack([p[w], q[w], u[w]], 1), (C.TRANGE[s], C.TRANGE[tt]))
        c2[(s, tt)] = C.above(max(float(v.max()), best) + mg2)
        print('  pair %s%s: float largest K + PAIR3 - Pi %.6e (file c2 %.6e); threshold %.6e'
              % (s, tt, v.max(), float(Z['c2'][len(c2) - 1]), float(c2[(s, tt)])), flush=True)
    jobs = [(A, C.pair3_poly(Bk, d3, C.TYPES.index(s), C.TYPES.index(tt)), c2[(s, tt)], (C.TRANGE[s], C.TRANGE[tt]), s == tt)
            for (s, tt) in c2]
    from multiprocessing import Pool
    with Pool(min(int(os.environ.get('PAIR_PROCS', '3')), len(jobs))) as pool:
        res = pool.starmap(C.pair_box_check, jobs)
    for (s, tt), (ok, msg) in zip(c2, res):
        C.check('K + PAIR3_%s%s <= Pi + c2 on the admissible pairs' % (s, tt), ok, msg)
    # (c) bins
    mf = np.array(Z['m'], float)
    assert len(mf) == len(C.BINTYPE)
    m = [C.above(v + mgm) for v in mf]
    okb, msg = bracket_check(A, z, D2, R2, C1, C2, C.DMAX, C.EDGES, m)
    C.check('f <= m_b on every bin', okb, msg)
    # (d) triples
    c3 = {}
    tjobs = []
    for combo in itertools.combinations_with_replacement(C.TYPES, 3):
        if all(C.Ntriple(tcounts(nb), list(combo)) == 0 for nb in C.COUNTS):
            continue
        P = C.triple3_poly(Bk, d3, combo)
        T12, T13, T23 = C.tmax(combo[0], combo[1]), C.tmax(combo[0], combo[2]), C.tmax(combo[1], combo[2])
        g = C.random_gram(float(T12), float(T13), float(T23), 200000, rng)
        ax = [np.r_[-1.0, np.linspace(-1, float(T), 40), float(T)] for T in (T12, T13, T23)]
        G = np.stack(np.meshgrid(*ax, indexing='ij'), -1).reshape(-1, 3)
        G = G[1 + 2 * G[:, 0] * G[:, 1] * G[:, 2] - (G ** 2).sum(1) >= 0]
        g = np.r_[g, G]
        vals = C.peval(P, g[:, 0], g[:, 1], g[:, 2])
        best = C.refine_triple(P, g[np.argsort(vals)[-40:]], (float(T12), float(T13), float(T23)))
        c3[combo] = C.above(max(float(vals.max()), best) + mg3)
        perm, sym = C.triple_layout(combo)
        tops = [None] * 3
        for src, val in enumerate((T12, T13, T23)):
            tops[perm.index(src)] = float(C.above(val))
        Pr = {}
        for e, co in P.items():
            ee = [0, 0, 0]
            for src in range(3):
                ee[perm.index(src)] += e[src]
            Pr[tuple(ee)] = Pr.get(tuple(ee), 0) + co
        print('  triple %s: float largest %.6e on %d samples; threshold %.6e; %d monomials, symmetry %d'
              % (''.join(combo), vals.max(), len(g), float(c3[combo]), len(P), sym), flush=True)
        tjobs.append((combo, (Pr, tops, sym, float(c3[combo]), ''.join(combo))))
    with Pool(int(os.environ.get('TRIPLE_PROCS', '3'))) as pool:
        res = pool.starmap(C.triple_box_check, [a for _, a in tjobs])
    nbox = 0
    for (combo, _), (ok, msg, nd) in zip(tjobs, res):
        C.check('TRIPLE3_%s <= c3 on the admissible triples' % ''.join(combo), ok, msg)
        nbox += nd
    print('  triple boxes in all: %d' % nbox, flush=True)
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
        if worst is None or val > worst:
            worst, wvec = val, nb
    print('  %d count vectors, %d type-count vectors; largest bound %.6f at %s' % (len(C.COUNTS), len(cache), float(worst), wvec),
          flush=True)
    C.check('largest bound over the count vectors below 9 pi^2/8 - 8', A_(worst) < target,
            '%.6f < %s' % (float(worst), target.str(8)))
    print('PASS: in case %s every packing set of %d centres has U(Y) <= %.6f < 9 pi^2/8 - 8, so T(Y) > 8 [%.0f s]'
          % (os.path.basename(sys.argv[1]), CASE['M'], float(worst), time.time() - t0))


if __name__ == '__main__':
    main()
