#!/usr/bin/env python3
"""
ipoly_test.py -- compares Poly.eval of certificate_check.py (numpy) with its C version
(ipoly.c, used with POLY_EVAL=c) bit for bit, on random polynomials with rational
coefficients of the size of the typed triple inequalities (up to 756 monomials, total
degree up to 20), their derivatives and the Gram determinant, and on random boxes:
proper boxes of widths from 1e-9 to 1, point boxes (as for the value and the gradient at
a centre), and boxes centred at 0.

Usage: python3 ipoly_test.py [trials]
"""
import sys
import time
from fractions import Fraction as Fr

import numpy as np

import certificate_check as CC

rng = np.random.default_rng(1)


def randpoly(nm, deg):
    p = {}
    while len(p) < nm:
        e = tuple(int(x) for x in rng.integers(0, deg + 1, 3))
        if sum(e) <= deg:
            p[e] = Fr(int(rng.integers(-10**12, 10**12)), int(rng.integers(1, 10**9)))
    return p


def main():
    trials = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    det = {(0, 0, 0): Fr(1), (1, 1, 1): Fr(2), (2, 0, 0): Fr(-1), (0, 2, 0): Fr(-1), (0, 0, 2): Fr(-1)}
    tot = bad = 0; tn = tc = 0.0
    for trial in range(trials):
        P = randpoly(756 if trial % 2 == 0 else 200, 20)
        for q in (P, CC.derivative(P, 0), CC.derivative(CC.derivative(P, 1), 2), det):
            pp = CC.Poly(q)
            n = 20000
            c = rng.uniform(-1, 0.6, (3, n)); w = rng.uniform(0, 1, (3, n)) * 10.0 ** rng.integers(-9, 0, (3, n))
            w[:, : n // 4] = 0.0
            c[:, n // 4: n // 3] = 0.0
            L, H = c - w, c + w
            CC.POLY_EVAL_C = False
            t = time.time(); a = pp.eval((L, H)); tn += time.time() - t
            CC.POLY_EVAL_C = True
            t = time.time(); b = pp.eval((L, H)); tc += time.time() - t
            same = (a[0] == b[0]) & (a[1] == b[1])
            tot += n; bad += int((~same).sum())
    print('%d boxes, %d differences; numpy %.2f s, C %.2f s' % (tot, bad, tn, tc))
    sys.exit(1 if bad else 0)


if __name__ == '__main__':
    main()
