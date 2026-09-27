#!/usr/bin/env python3
"""
hull_along_rays.py -- the inversion hull along the two rays of statement (ii) of
"What is left" (Section 2.8 of the paper; floating point, exploration).

For the root system with one centre pushed out to 2 + delta, and with all 24 pushed out
evenly to 2 + S/24, it prints vol V(Y) - 8, the bound F(Y) - 8 = vol(V(Y) cap K(Y)) - 8 of
Proposition 2.13, their difference vol(V(Y) \\ K(Y)) (the part far centres may cut, which
part (c) of the proof of Theorem 2.18 bounds by 13958 Theta^4), and the second-order
model (2/3) S - (1/2) S^2 of Section 2.9.
"""
import numpy as np
from near_contact_probe import R, vol, F

print('%-22s %8s %12s %12s %12s %12s %14s' % ('family', 'S', 'volV-8', 'F-8', 'V minus K', 'model', '13958 Theta^4'))
for name in ('one centre', 'all 24 evenly'):
    for S in (0.004, 0.01, 0.02, 0.05, 0.1, 0.155, 0.197):
        d = np.full(24, 2.0)
        if name == 'one centre':
            d[0] += S; Theta = S
        else:
            d += S / 24; Theta = S / 24
        Y = d[:, None] * R
        v, f = vol(Y) - 8, F(Y) - 8
        print('%-22s %8.3f %12.6f %12.6f %12.3e %12.6f %14.3e' % (name, S, v, f, v - f, 2 / 3 * S - S * S / 2, 13958 * Theta**4))
