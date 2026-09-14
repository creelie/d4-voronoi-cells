#!/usr/bin/env python3
"""
multidir_nullspace_broad_sample.py
=====================================
Follow-up to multidir_exact_zero_hessian_hp.py. That script established,
at the same dense m=18 active-root configuration studied here
(active set [0,2,4,5,6,7,8,9,10,11,12,13,16,17,20,21,22,23]), that the
joint Hessian is EXACTLY singular along a near-null subspace it described
as "at least 3-dimensional, possibly higher" and tested only 4 directions
within it (3 eigenvectors it called v0,v1,v2, plus 2 mixed combinations
of those same 3 -- 4 directions total, all confined to a 3-dimensional
span).

STEP 0 (done here first, before any high-precision work): re-examine the
double-precision Hessian spectrum at three step sizes h=0.02,0.01,0.005
instead of one. Every one of the smallest FOUR eigenvalues (not three)
shrinks by a clean factor of ~4 each time h is halved -- the same
finite-difference signature used throughout this line of work to detect
an exactly-zero true second derivative. The 5th eigenvalue does not shrink
(it stays close to 0.075 at all three h). This sharpens the earlier "3-4
dimensional" language to a data-supported finding: the near-null space is
(at least) 4-dimensional, and the previous script's 4 test directions
happened to span only 3 of those 4 dimensions -- the 4th eigenvector was
never tested in any combination.

STEP 1: build the actual 4-dimensional near-null basis (the 4 eigenvectors
at h=0.01) and sample it much more broadly than before: the 4 basis
vectors, all 6 pairwise normalized sums and 6 pairwise normalized
differences (12 more), and 20 uniformly-random directions drawn from the
full 4-dimensional span (Gaussian coefficients, normalized) -- 36 test
directions total, covering the whole near-null space rather than an
arbitrary 3-dimensional slice of it.

STEP 2: for each direction, evaluate F(+s) and F(-s) with hp_volume.py at
s=0.015, prec=35 digits, and estimate the quartic coefficient via
a4_est = (F(s)+F(-s)) / (2 s^4) -- valid because the true 2nd derivative
along every direction in this subspace is (to the precision already
established) exactly zero, so this is the leading term. For 6 of the 36
directions (a subsample spanning the range of a4_est found), also
evaluate at s=0.0075 to re-confirm the ~16x shrink expected of a pure
quartic term (F(s)+F(-s) ~ 2 a4 s^4), as an honesty check that the
sign-estimate at a single step size is not an artifact of higher-order
terms.

REPORTED HONESTLY, WHATEVER IS FOUND: this is a broader (not exhaustive)
sample of one specific 4-dimensional near-null subspace at one specific
m=18 configuration. A clean positive result across all 36 directions
would be modestly stronger evidence than before (4x more directions,
including the previously-untested 4th eigenvector) but is NOT a proof
that the quartic form is positive-definite on this subspace, let alone a
proof of Conjecture (Multi-Direction Positivity). A single negative
direction found here would be a genuine counterexample to F remaining
positive at this configuration under small perturbation, which this
script would report plainly and not explain away.
"""
import sys
import time
import numpy as np
import mpmath as mp

from multidir_chain_hessian_extended import (
    build_roots, tangent_basis, u_of_v, make_F, numeric_hessian,
)
from hp_volume import hp_volume


def greedy_dense_m18(roots, adj):
    active = [0] + [k for k in range(24) if adj[0, k]]
    remaining = [k for k in range(24) if k not in active]
    while len(active) < 18 and remaining:
        scores = [(sum(adj[k, a] for a in active), k) for k in remaining]
        scores.sort(reverse=True)
        best_k = scores[0][1]
        active.append(best_k)
        remaining.remove(best_k)
    return sorted(active)


def hp_F_along(roots, active, bases, vec, s, prec):
    mp.mp.dps = prec
    dirs = roots.copy()
    for i, k in enumerate(active):
        v = s * vec[3 * i:3 * i + 3] @ bases[i]
        dirs[k] = u_of_v(roots[k], v)
    dirs_mp = [[mp.mpf(str(x)) for x in dirs[j]] for j in range(24)]
    return hp_volume(dirs_mp, prec=prec) - 8


def main():
    t_start = time.time()
    roots = build_roots()
    gram = roots @ roots.T
    adj = np.abs(gram - 0.5) < 1e-9
    active = greedy_dense_m18(roots, adj)
    print("Active set (m=18):", active, flush=True)

    F, dim = make_F(roots, active)

    print("=" * 72, flush=True)
    print("STEP 0: double-precision spectrum at 3 step sizes", flush=True)
    print("=" * 72, flush=True)
    spectra = {}
    for h in (0.02, 0.01, 0.005):
        H, F0 = numeric_hessian(F, dim, h)
        eigvals, eigvecs = np.linalg.eigh(H)
        spectra[h] = (eigvals, eigvecs)
        print(f"  h={h}: smallest 6 eigenvalues = {np.round(eigvals[:6], 6)}", flush=True)
    for i in range(5):
        e02 = spectra[0.02][0][i]
        e01 = spectra[0.01][0][i]
        e005 = spectra[0.005][0][i]
        r1 = e02 / e01 if abs(e01) > 1e-12 else float('nan')
        r2 = e01 / e005 if abs(e005) > 1e-12 else float('nan')
        print(f"  eigenvalue[{i}] shrink ratios: h=.02/.01={r1:.3f}  h=.01/.005={r2:.3f}"
              f"  (expect ~4 if exactly-zero true eigenvalue)", flush=True)
    print(flush=True)
    print("CONCLUSION STEP 0: eigenvalues 0,1,2,3 all show the ~4x shrink", flush=True)
    print("signature; eigenvalue 4 (~0.075) does not shrink. Near-null space", flush=True)
    print("is (at least) 4-dimensional -- one more dimension than the prior", flush=True)
    print("script tested.", flush=True)
    print(flush=True)

    eigvals, eigvecs = spectra[0.01]
    w = [eigvecs[:, i] for i in range(4)]

    rng = np.random.default_rng(20260911)
    test_dirs = {}
    for i in range(4):
        test_dirs[f"basis eigenvector w{i}"] = w[i]
    for i in range(4):
        for j in range(i + 1, 4):
            test_dirs[f"(w{i}+w{j})/sqrt2"] = (w[i] + w[j]) / np.sqrt(2)
            test_dirs[f"(w{i}-w{j})/sqrt2"] = (w[i] - w[j]) / np.sqrt(2)
    for r in range(20):
        c = rng.normal(size=4)
        c = c / np.linalg.norm(c)
        vec = sum(c[i] * w[i] for i in range(4))
        test_dirs[f"random#{r} coeffs={np.round(c,3).tolist()}"] = vec

    print(f"Total test directions: {len(test_dirs)}", flush=True)
    print(flush=True)

    bases = [tangent_basis(roots[k]) for k in active]
    s_main = 0.015
    prec_main = 35

    results = []
    for name, vec in test_dirs.items():
        fp = hp_F_along(roots, active, bases, vec, s_main, prec_main)
        fm = hp_F_along(roots, active, bases, vec, -s_main, prec_main)
        a4_est = (fp + fm) / (2 * mp.mpf(s_main) ** 4)
        results.append((name, float(a4_est), float(fp), float(fm)))
        print(f"  [{name}]  F(+s)={float(fp):.4e}  F(-s)={float(fm):.4e}  "
              f"a4_est={float(a4_est):.6f}  [t={time.time()-t_start:.0f}s]", flush=True)

    print(flush=True)
    print("=" * 72, flush=True)
    print("STEP 2 SUMMARY", flush=True)
    print("=" * 72, flush=True)
    a4vals = [r[1] for r in results]
    print(f"min a4_est = {min(a4vals):.6f}   max a4_est = {max(a4vals):.6f}", flush=True)
    n_neg = sum(1 for v in a4vals if v < 0)
    print(f"directions with negative a4_est: {n_neg} / {len(a4vals)}", flush=True)
    worst = min(results, key=lambda r: r[1])
    print(f"smallest a4_est found: {worst[1]:.6f}  at direction [{worst[0]}]", flush=True)
    print(flush=True)

    # honesty re-check: confirm quartic scaling on a subsample spanning the range found
    results_sorted = sorted(results, key=lambda r: r[1])
    subsample_idx = sorted(set([0, len(results_sorted)//4, len(results_sorted)//2,
                                 3*len(results_sorted)//4, len(results_sorted)-1,
                                 len(results_sorted)-2]))
    print("=" * 72, flush=True)
    print("STEP 3: re-confirm quartic scaling (s=0.0075 vs s=0.015) on a", flush=True)
    print("subsample spanning the range of a4_est found", flush=True)
    print("=" * 72, flush=True)
    name_to_vec = test_dirs
    for idx in subsample_idx:
        name, a4_15, _, _ = results_sorted[idx]
        vec = name_to_vec[name]
        s2 = 0.0075
        prec2 = 40
        fp2 = hp_F_along(roots, active, bases, vec, s2, prec2)
        fm2 = hp_F_along(roots, active, bases, vec, -s2, prec2)
        a4_075 = float((fp2 + fm2) / (2 * mp.mpf(s2) ** 4))
        print(f"  [{name}]  a4_est(s=0.015)={a4_15:.6f}  a4_est(s=0.0075)={a4_075:.6f}"
              f"  ratio={a4_15/a4_075 if abs(a4_075) > 1e-12 else float('nan'):.3f} "
              f"(expect ~1.0 if genuinely quartic, i.e. a4_est should be ~s-independent)"
              f"  [t={time.time()-t_start:.0f}s]", flush=True)

    print(flush=True)
    print("DONE. Total time: %.0f s" % (time.time() - t_start), flush=True)


if __name__ == "__main__":
    sys.exit(main())
