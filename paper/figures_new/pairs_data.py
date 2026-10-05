#!/usr/bin/env python3
"""pairs_data.py -- data for tikz_pairs.tex: the weight omega(gamma) of
prop:pair-budget, (1/4) int_0^{r*} Lambda(r, gamma) d(sec^4 r) with the lens
measure Lambda of eq:lens in closed form, on [60, 72] degrees, and its values at
the pair angles of the deletion of a root and of the two-point minimiser of
prop:two-point-barrier.  Writes pairs_omega.dat."""
from mpmath import mp, mpf, pi, sin, tan, cos, acos, sqrt, quad, radians, degrees

mp.dps = 30
RS = acos(sqrt(mpf(2) / 3))


def lens(r, g):
    if r <= g / 2:
        return mpf(0)
    return 4 * pi * ((r - g / 2) / 2 - (sin(2 * r) - sin(g)) / 4
                     - tan(g / 2) / 2 * (sin(r) ** 2 - sin(g / 2) ** 2))


def omega(g):
    if g >= 2 * RS:
        return mpf(0)
    # d(sec^4 r) = 4 sec^4 r tan r dr
    return quad(lambda r: lens(r, g) * 4 * tan(r) / cos(r) ** 4, [g / 2, RS]) / 4


if __name__ == '__main__':
    print('omega(60) = %s (paper: 0.001445408...)' % mp.nstr(omega(radians(60)), 12))
    print('2 r* = %s degrees' % mp.nstr(degrees(2 * RS), 10))
    dele = 88 * omega(radians(60))
    mini = sum(n * omega(radians(a)) for n, a in ((103.7, 62.3), (90.0, 103.0), (39.2, 122.6), (20.1, 154.1)))
    print('deletion: 88 omega(60) = %s;  minimiser: %s;  target 8 - A* = 0.0928556' % (mp.nstr(dele, 8), mp.nstr(mini, 8)))
    print('omega(62.3) = %s' % mp.nstr(omega(radians(62.3)), 8))
    with open('pairs_omega.dat', 'w') as f:
        f.write('g w\n')
        n = 240
        for i in range(n + 1):
            g = 60 + 12 * mpf(i) / n
            f.write('%.5f %.8e\n' % (float(g), float(omega(radians(g)) * 1000)))
