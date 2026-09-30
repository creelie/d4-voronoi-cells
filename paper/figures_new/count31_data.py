#!/usr/bin/env python3
"""
count31_data.py -- data of the figure for Theorem 21.25 (tikz_count31.tex).

Writes count31_bounds.dat: for each count M, the best two-point bound on the
union of the caps (multi_cap/runs/radial_count_scan.log) and the largest union
found by the search of Figure 36 (9 pi^2/8 minus the least T of
multi_cap/runs/count_survey.log); and count31_kernel.dat: the kernel K of the
certificate of Theorem 21.25 (multi_cap/radial_certificates/radial_31.json) and
the pair term Pi at d = d' = 2, both times 10^3, against the inner product u.
Run from this directory.
"""
import json
import math
import re
import sys
from fractions import Fraction as Fr

import numpy as np

sys.path.insert(0, '../../multi_cap')
from truncated_search import pair                   # noqa: E402
from radial_count_sdp import pbasis, ubasis         # noqa: E402

VB = 9 * math.pi ** 2 / 8
scan = {int(m): float(v) for m, v in re.findall(r'M = +(\d+): best two-point bound on the union ([0-9.]+)',
                                                  open('../../multi_cap/runs/radial_count_scan.log').read())}
found = {int(m): VB - float(v) for m, v in re.findall(r'M = (\d+): least T over .*? = ([0-9.]+)',
                                                        open('../../multi_cap/runs/count_survey.log').read())}
with open('count31_bounds.dat', 'w') as f:
    f.write('M bound found\n')
    for M in range(24, 50):
        b = scan.get(M, float('nan'))
        u = found.get(M, float('nan'))
        if not (math.isnan(b) and math.isnan(u)):
            f.write('%d %s %s\n' % (M, 'nan' if math.isnan(b) else '%.5f' % b, 'nan' if math.isnan(u) else '%.5f' % u))

c = json.load(open('../../multi_cap/radial_certificates/radial_31.json'))
D, r = c['D'], c['r']
A = np.array([[[float(Fr(x)) for x in row] for row in a] for a in c['A']])
u = np.concatenate([np.linspace(-1, 0.25, 250), np.linspace(0.25, 0.5, 251)[1:]])
d = np.full_like(u, 2.0)
K = np.einsum('nk,na,kab,nb->n', ubasis(u, D), pbasis(d, r), A, pbasis(d, r))
P = pair(1.0, 1.0, u)
with open('count31_kernel.dat', 'w') as f:
    f.write('u K Pi\n')
    for row in zip(u, 1e3 * K, 1e3 * P):
        f.write('%.5f %.6f %.6f\n' % row)
print('wrote count31_bounds.dat and count31_kernel.dat; largest K - Pi on the curve: %.3e' % (K - P).max())
