#!/usr/bin/env python3
"""
c28_typed_dirs.py m d rounds -- the typed three-point bound of ../C30/typed3pt.py on
the directions of 28 centres with the radial counts of 28, every centre at the largest
distance its count allows: m within 2.0161, 21 - m at 2.1, three at 2.35 (one each in
(2.1, 2.15], (2.15, 2.2], (2.2, 2.35], merged at the looser 2.35) and four at sqrt6.
A value below 0 would show that no such directions exist; floating point.
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'C30'))
import typed3pt as T3

R6 = 6 ** .5
a = T3.a_
m, d, rounds = map(int, sys.argv[1:4])
ds = [2.0161, 2.1, 2.35, R6]
cnt = [m, 21 - m, 3, 4]
T = [[a(x, y) for y in ds] for x in ds]
print('28: %d within 2.0161, %d at 2.1, 3 at 2.35, 4 at sqrt6; T =\n%s' % (m, 21 - m, np.round(T, 5)), flush=True)
T3.solve(cnt, T, d, rounds)
