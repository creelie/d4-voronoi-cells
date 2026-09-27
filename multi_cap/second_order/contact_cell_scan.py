#!/usr/bin/env python3
"""
contact_cell_scan.py -- how far the root system stays a minimum of the
contact-cell volume vol{x : <x, w_i> <= 1} among sets of 24 unit directions
(floating point, qhull volumes; exploration).

(a) Neighbourhood: for t = 0.05, 0.1, 0.2, 0.3, 0.4 (the norm of the tilt
    field, rotations removed, in radians), minimise vol - 8 over tilts of that
    norm from 10 random starts (L-BFGS with central differences, the norm
    held by rescaling), and print the least value and (vol - 8)/t^4.
(b) Globally, without the packing constraint: minimise vol over all sets of
    24 unit vectors from 20 random starts, and print the least volume found.
    The packing constraint (pairwise angles at least 60 degrees) is what
    forces the root system; without it a smaller circumscribed polytope with
    24 facets would show that the local statement cannot be globalised.
Usage: python3 contact_cell_scan.py [a|b]
"""
import sys
import time
import numpy as np
from scipy.optimize import minimize
from scipy.spatial import HalfspaceIntersection, ConvexHull
import tilt_quartic as TQ


def vol_dirs(w):
    w = w / np.linalg.norm(w, axis=1, keepdims=True)
    try:
        hs = np.hstack([w, -np.ones((len(w), 1))])
        return ConvexHull(HalfspaceIntersection(hs, np.zeros(4)).intersections).volume
    except Exception:
        return 1e6            # unbounded or degenerate


def part_a():
    basis = np.vstack([TQ.FLAT, TQ.TRANS])          # 66 x 72, orthonormal, rotations removed
    rng = np.random.default_rng(7)
    print('t      least vol - 8     (vol-8)/t^4   [time]')
    for t in (0.05, 0.1, 0.2, 0.3, 0.4):
        t0 = time.time(); best = np.inf
        for s in range(10):
            def f(z):
                c = z @ basis
                c = t * c / np.sqrt(TQ.ip(c, c))
                return TQ.vol(c) - 8
            def g(z, h=1e-6):
                out = np.zeros_like(z)
                for j in range(len(z)):
                    e = np.zeros_like(z); e[j] = h
                    out[j] = (f(z + e) - f(z - e)) / (2 * h)
                return out
            r = minimize(f, rng.normal(size=66), jac=g, method='L-BFGS-B', options={'maxiter': 150})
            best = min(best, r.fun)
        print('%.2f   %.3e        %.4f        [%.0fs]' % (t, best, best / t ** 4, time.time() - t0), flush=True)


def part_b():
    rng = np.random.default_rng(11)
    best = np.inf
    t0 = time.time()
    for s in range(20):
        x0 = rng.normal(size=(24, 4)).ravel()
        f = lambda z: vol_dirs(z.reshape(24, 4))
        r = minimize(f, x0, method='Powell', options={'maxiter': 20000, 'xtol': 1e-6, 'ftol': 1e-10})
        best = min(best, r.fun)
        print('start %2d: volume %.6f   (best %.6f)  [%.0fs]' % (s, r.fun, best, time.time() - t0), flush=True)
    print('least contact-cell volume found over 24 unit directions, no packing constraint: %.6f' % best)


if __name__ == '__main__':
    part = sys.argv[1] if len(sys.argv) > 1 else 'a'
    part_a() if part == 'a' else part_b()
