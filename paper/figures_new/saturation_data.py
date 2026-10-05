#!/usr/bin/env python3
"""saturation_data.py -- the data of fig:saturation (lem:cells-density, lem:redistribution).

A saturated packing of unit discs in the plane, made by random sequential
addition in the square [-10,10]^2 until 400000 consecutive trials fail and then
filled by a grid search for every point at distance >= 2 from all centres, so that
no further disc fits.  Its Voronoi cells are clipped to the square.  Writes
tikz_saturation_data.tex with the cells of the centres within R = 5.5 of the
origin, every centre within R + 3, and the pairs (c, c') with |c| <= R < |c'| and
|c - c'| < r = 2.6, the terms of the proof of lem:redistribution that do not cancel.
"""
import numpy as np
from scipy.spatial import Voronoi

rng = np.random.default_rng(7)
L, R, r = 10.0, 5.5, 2.6
P = []
fails = 0
while fails < 400000:
    p = rng.uniform(-L, L, 2)
    if all((p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2 >= 4 for q in P[-400:]) and \
       (not P or np.min(np.sum((np.array(P) - p) ** 2, 1)) >= 4):
        P.append(p); fails = 0
    else:
        fails += 1
P = np.array(P)
g = np.linspace(-L, L, 801)
for x in g:
    for y in g:
        if np.min(np.sum((P - (x, y)) ** 2, 1)) >= 4:
            P = np.vstack([P, (x, y)])
assert all(np.min(np.sum((P - (x, y)) ** 2, 1)) < 4 for x in g[::4] for y in g[::4])
d = np.sqrt(((P[:, None] - P[None]) ** 2).sum(-1)) + 10 * np.eye(len(P))
assert d.min() >= 2 - 1e-12
far = np.array([[x, y] for x in (-60, 60) for y in (-60, 60)])
vor = Voronoi(np.vstack([P, far]))
out = []
inner = [i for i in range(len(P)) if np.hypot(*P[i]) <= R]
for i in inner:
    reg = vor.regions[vor.point_region[i]]
    assert -1 not in reg
    V = vor.vertices[reg]
    assert np.max(np.hypot(*(V - P[i]).T)) < 2, 'saturated: cells within 2 of their centres'
    out.append('\\filldraw[cellfill, celldraw] ' + ' -- '.join('(%.3f,%.3f)' % tuple(v) for v in V) + ' -- cycle;')
for i in range(len(P)):
    if np.hypot(*P[i]) <= R + 3:
        out.append('\\fill[%s] (%.3f,%.3f) circle (\\dotr);' % ('cinner' if np.hypot(*P[i]) <= R else 'couter', *P[i]))
for i in inner:
    for j in range(len(P)):
        if np.hypot(*P[j]) > R and np.hypot(*(P[i] - P[j])) < r:
            out.append('\\draw[transfer] (%.3f,%.3f) -- (%.3f,%.3f);' % (*P[i], *P[j]))
open('tikz_saturation_data.tex', 'w').write('\n'.join(out) + '\n')
print(len(P), 'centres;', len(inner), 'within R;', sum('transfer' in o for o in out), 'boundary pairs')
