#!/usr/bin/env python3
"""cap_paths_data.py -- the data of fig:cap-scheme(b).

The cap of the cross-polytope B = {|x_0|+|x_1|+|x_2|+|x_3| <= 2} cut off by the
half-space <x,u> >= 1 has volume (1/3) g[u_0^2, u_1^2, u_2^2, u_3^2]
(lem:cap-formula), g(a) = a (2 sqrt(a) - 1)^4 for a >= 1/4 and 0 below.
Two paths of unit vectors u start at the vertex direction e_0:
  path 1, towards the midpoint of an edge:  u = (cos t, sin t, 0, 0),         0 <= t <= 45 deg;
  path 2, towards the centre of a facet:    u = (cos t, s, s, s), s = sin t/sqrt3, 0 <= t <= 60 deg.
With nodes at 0 or below 1/4, where g vanishes identically, the divided
differences reduce to
  g[a, b, 0, 0] = (G(a) - G(b))/(a - b),  G(a) = g(a)/a^2,   and   g[a, b, b, b] = g(a)/(a - b)^3  (b <= 1/4).
A Monte Carlo estimate of the cap volume at four directions is printed as a check.
Writes cap_paths.dat.
"""
import math
import numpy as np


def g(a):
    return a * (2 * math.sqrt(a) - 1) ** 4 if a > 0.25 else 0.0


def G(a):
    return g(a) / a ** 2


def dG(a, h=1e-6):
    return (G(a + h) - G(a - h)) / (2 * h)


def path1(t):
    a, b = math.cos(t) ** 2, math.sin(t) ** 2
    if abs(a - b) < 1e-9:
        return dG(a) / 3
    if b == 0:
        return g(a) / a ** 2 / 3 if a == 1 else (G(a) - 0) / (a - b) / 3
    return (G(a) - G(b)) / (a - b) / 3


def path2(t):
    a, b = math.cos(t) ** 2, math.sin(t) ** 2 / 3
    return g(a) / (a - b) ** 3 / 3 if a > b else 0.0


def mc(u, n=4_000_000, seed=1):
    rng = np.random.default_rng(seed)
    e = rng.exponential(size=(n, 5))
    x = 2 * e[:, :4] / e.sum(1, keepdims=True) * rng.choice([-1, 1], size=(n, 4))
    return 32 / 3 * np.mean(x @ np.asarray(u) >= 1)   # vol B = 2^4 * 2^4/4! = 32/3


def main():
    for t in (0.0, 0.3, 0.7):
        print('path 1, t = %.1f: formula %.5f, Monte Carlo %.5f' % (t, path1(t), mc((math.cos(t), math.sin(t), 0, 0))))
    t = math.asin(0.6)
    s = 0.6 / math.sqrt(3)
    print('path 2, t = %.4f: formula %.5f, Monte Carlo %.5f' % (t, path2(t), mc((0.8, s, s, s))))
    with open('cap_paths.dat', 'w') as f:
        f.write('t c1 c2\n')
        for d in np.linspace(0, 60, 121):
            t = math.radians(d)
            c1 = '%.6f' % path1(t) if d <= 45 + 1e-9 else 'nan'
            f.write('%.2f %s %.6f\n' % (d, c1, path2(t)))


if __name__ == '__main__':
    main()
