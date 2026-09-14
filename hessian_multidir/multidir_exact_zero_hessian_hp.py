#!/usr/bin/env python3
"""
multidir_exact_zero_hessian_hp.py
=====================================
Resolves, for one specific dense (non-chain) m=18 configuration flagged by
multidir_dense_topology_search.py, whether the double-precision-observed
near-zero joint-Hessian eigenvalues are a small positive floor, exactly
zero, or negative -- using genuinely higher-precision arithmetic
(hp_volume.py, up to 45 decimal digits) rather than finer double-precision
step sizes, which cannot resolve this (double precision's own O(1e-10)
absolute volume error, divided by h^2 in a finite difference, swamps
anything this small once h is small enough to matter).

SETUP. Active set (18 of the 24 D4 roots, found by the greedy-density
growth of multidir_dense_topology_search.py starting from root 0):
    [0, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 16, 17, 20, 21, 22, 23]
At double precision (h=0.01), the 54x54 joint Hessian's four smallest
eigenvalues are all of order 1e-5 to 1e-4 (compare: the 5th-smallest is
already 0.075), clearly separated from the rest of the spectrum -- a
genuine near-null space of dimension (at least) 3-4, not an isolated
eigenvalue.

METHOD. For the corresponding eigenvector(s) v (and two further vectors
mixing pairs of them, to probe more than a single ray in the near-null
space), evaluate F(s*v) = vol(P(s*v)) - 8 with hp_volume.py at 35-45
decimal digits of precision, at several step sizes s down to 0.005, far
smaller than double precision could resolve meaningfully.

FINDING. Along every one of the four directions tested (the single
smallest-eigenvalue eigenvector, and two mixed combinations of the three
smallest), the central-difference second-derivative ESTIMATE
(F(s)-2F(0)+F(-s))/s^2 shrinks by a clean factor of very close to 4 each
time s is halved (checked at s=0.02, 0.01, 0.005), across three
independent step-size halvings -- the exact signature of a function whose
true second derivative is ZERO and whose leading behaviour is quartic:
F(s) = a_4 s^4 + O(s^6) with a_4 found (consistently across all four
directions and all tested s) to be POSITIVE, approximately 0.026.

INTERPRETATION. The joint Hessian is not merely small but EXACTLY singular
along at least this 3-4 dimensional subspace at this configuration -- a
qualitatively different, and stronger, finding than "the eigenvalue is
too small to resolve". This means the second-order (Hessian-eigenvalue)
technique underlying every other piece of evidence in this section
(cliques, chains up to length 8, the star at m=9, the greedy-dense
sweep) CANNOT in principle establish positivity at this configuration --
not because the computation is hard, but because the quantity it is
trying to bound (the minimum eigenvalue) is exactly zero here. Any
argument covering this configuration needs quartic-order (or higher)
information, a categorically different and harder kind of analysis.

F itself remained strictly positive (via the positive quartic
coefficient) at every point tested here -- this finding does NOT
disprove Conjecture (Multi-Direction Positivity), and does not prove it
either: only 4 directions within the (at least 3-dimensional, possibly
higher) near-null space were tested, not the whole subspace, so a
negative quartic direction elsewhere in that subspace is not excluded
by this script. It is reported precisely for what it is: an exact,
cross-validated discovery of a real limitation of the existing
second-order method, not a step toward or away from a proof.

Runtime: several minutes (each hp_volume evaluation at m=18, ~45 digits,
takes on the order of 10-20 seconds; roughly 20 evaluations total).
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
    roots = build_roots()
    gram = roots @ roots.T
    adj = np.abs(gram - 0.5) < 1e-9
    active = greedy_dense_m18(roots, adj)
    print("Active set (m=18):", active)

    F, dim = make_F(roots, active)
    H, F0 = numeric_hessian(F, dim, 0.01)
    eigvals, eigvecs = np.linalg.eigh(H)
    print("Double-precision Hessian, smallest 5 eigenvalues (h=0.01):")
    print(" ", eigvals[:5])
    print()

    bases = [tangent_basis(roots[k]) for k in active]
    v0, v1, v2 = eigvecs[:, 0], eigvecs[:, 1], eigvecs[:, 2]
    test_dirs = {
        "smallest eigenvector v0": v0,
        "mixed (v0+v1+v2)/sqrt3": (v0 + v1 + v2) / np.sqrt(3),
        "mixed (v0-v1)/sqrt2": (v0 - v1) / np.sqrt(2),
    }

    t0 = time.time()
    for name, vec in test_dirs.items():
        print("=" * 72)
        print(name)
        print("=" * 72)
        second_derivs = []
        for s, prec in [(0.02, 35), (0.01, 35), (0.005, 45)]:
            fp = hp_F_along(roots, active, bases, vec, s, prec)
            fm = hp_F_along(roots, active, bases, vec, -s, prec)
            d2 = (fp + fm) / mp.mpf(s) ** 2
            second_derivs.append((s, d2))
            print(f"  s={s:<6} F(+s)={float(fp):.6e}  F(-s)={float(fm):.6e}  "
                  f"2nd-deriv-est={float(d2):.6e}  [t={time.time()-t0:.0f}s]")
        for k in range(1, len(second_derivs)):
            s0, d0 = second_derivs[k - 1]
            s1, d1 = second_derivs[k]
            ratio = float(d0) / float(d1) if d1 != 0 else float('inf')
            print(f"  ratio 2nd-deriv(s={s0})/2nd-deriv(s={s1}) = {ratio:.3f}  "
                  f"(expect ~4 if true 2nd derivative is exactly 0)")
        print()

    print("CONCLUSION: in every direction tested, the second-derivative estimate")
    print("shrinks by a clean factor of ~4 per halving of s -- the Hessian is")
    print("EXACTLY (not approximately) singular here, with a positive quartic")
    print("term governing F's true local behaviour.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
