"""combo_mixtype_check.py -- exact check of a combination of certificates whose three-point
parts type the bins differently (mixed typing), on a case file with one type per bin.

usage: python3 combo_mixtype_check.py CASE.json MIX.npz D3 [margin2 margin3 marginm]
with the options of combo_case_check.py (SKIP_EMPTY=1 is set here).

The three-point inequality of a certificate holds for every assignment of its types A, B, F
to the points: its blocks are positive semidefinite, and for points x_1, ..., x_N on the
sphere with types s_1, ..., s_N the sum over i, j, k of the typed kernel is a sum of
quadratic forms of those blocks.  So a certificate can type the centres of each bin as any
of its types; only the sets on which its pair and triple polynomials must be bounded, and
the point term of each bin, change.  MIX.npz holds the shared kernel of the distances
(A, z, t) and the bin constants m, as for combo_case_check.py, and, for each typing g of
the bins (a row gtyp[g] giving the type index 0, 1, 2 of A, B, F of every bin), the
three-point vector x3g[g]: the weighted sum of the three-point parts of the certificates
used under that typing.  Each group is rounded to exact positive definite blocks on its
own; then, the case types being the bins, the point term of a bin is the sum over the
groups of the point term of the type the group gives it, and the pair and triple
polynomials of a combination of bins are the sums over the groups of the polynomials of
the types the group gives those bins.  Everything else is combo_case_check.py: the
thresholds are sampled, refined and proved by branch and bound per pair and per triple of
bins, and the bound is taken over every count vector.

The case file (case_bins.py makes one from an ordinary case) must give each bin a type of
its own whose range is the bin.
"""
import os
import sys
from fractions import Fraction as Fr

import numpy as np

os.environ['SKIP_EMPTY'] = '1'
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import combo30_check as C  # noqa: E402
import combo_case_check as K  # noqa: E402  (reads the case file of sys.argv[1])

Z = np.load(sys.argv[2])
CT = ['A', 'B', 'F']                       # the types of every certificate, in their order
GT = np.array(Z['gtyp'], int)              # groups x bins
NG = len(GT)
VT = list(C.TYPES)                         # one type per bin
assert [C.BINTYPE.index(s) for s in VT] == list(range(len(VT))) and len(C.BINTYPE) == len(VT) == GT.shape[1], \
    'the case must give each bin a type of its own, in bin order'
TMAP = [[int(GT[g][C.BINTYPE.index(s)]) for s in VT] for g in range(NG)]   # group -> case type -> certificate type
print('mixed typing: %d groups, typings %s' % (NG, [''.join(CT[i] for i in row) for row in GT]), flush=True)

_blocks, _point, _pair, _triple = C.blocks_from_x3, C.point_term, C.pair3_poly, C.triple3_poly


def blocks_from_x3(x3, d3):
    """the blocks of every group, keyed (g,) + name; the argument (a placeholder) is unused"""
    out, saved = {}, C.TYPES
    C.TYPES = CT
    try:
        for g in range(NG):
            for name, M in _blocks(np.array(Z['x3g'][g], float), d3).items():
                out[(g,) + name] = M
    finally:
        C.TYPES = saved
    return out


def groups(Bk):
    return [{name[1:]: M for name, M in Bk.items() if name[0] == g} for g in range(NG)]


def point_term(Bk, d3, i):
    return sum(_point(Bg, d3, TMAP[g][i]) for g, Bg in enumerate(groups(Bk)))


def pair3_poly(Bk, d3, i, j):
    p = [Fr(0)]
    for g, Bg in enumerate(groups(Bk)):
        p = C.ladd(p, _pair(Bg, d3, TMAP[g][i], TMAP[g][j]))
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def triple3_poly(Bk, d3, types):
    out, saved = {}, C.TYPES
    C.TYPES = CT
    try:
        for g, Bg in enumerate(groups(Bk)):
            out = C.padd(out, _triple(Bg, d3, tuple(CT[TMAP[g][VT.index(x)]] for x in types)))
    finally:
        C.TYPES = saved
    return {e: c for e, c in out.items() if c}


C.blocks_from_x3, C.point_term, C.pair3_poly, C.triple3_poly = blocks_from_x3, point_term, pair3_poly, triple3_poly

if __name__ == '__main__':
    K.main()
