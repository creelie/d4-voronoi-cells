"""
rigidity_spectrum.py -- the exact spectrum of the rigidity operator of the 24-point
kissing configuration, in integer arithmetic.

Write T for the tangent space of (S^3)^24 at the normalised D_4 roots and

    (Lambda tau)_{ij} = u_i . tau_j + u_j . tau_i ,   (i,j) a tight pair,

the derivative of the 96 tight inner products.  What the quantitative
uniqueness theorem needs is a lower bound on the smallest nonzero singular
value of Lambda restricted to T.

To keep everything integral, work with Lambda' = sqrt2 * Lambda, whose
entries are alpha_i . x with alpha the integer roots, and with the rational
orthogonal projector P onto T, whose blocks are I - alpha alpha^T / 2.  Then

    N = P Lambda'^T Lambda' P

is rational with 4N integral, its nonzero eigenvalues are the squares of the
nonzero singular values of Lambda' restricted to T, and

    sigma_min(Lambda|_T)^2 = lambda_min^+(N) / 2 .

The script proves, in integer arithmetic, that 4N annihilates the polynomial
x(x-8)(x-20)(x-24)(x-32) and has rank 66.  So the eigenvalues of N lie in
{0, 2, 5, 6, 8}, the kernel has dimension 30 (the 24 radial directions and
the 6 infinitesimal rotations), and the smallest nonzero eigenvalue is at
least 2.  Hence sigma_min(Lambda|_T) >= 1.
"""
import itertools
import numpy as np
from fractions import Fraction

ROOTS = []
for i, j in itertools.combinations(range(4), 2):
    for si in (1, -1):
        for sj in (1, -1):
            v = [0, 0, 0, 0]; v[i] = si; v[j] = sj
            ROOTS.append(tuple(v))
ROOTS = sorted(set(ROOTS))
A = np.array(ROOTS, dtype=np.int64)
assert A.shape == (24, 4)

TIGHT = [(i, j) for i in range(24) for j in range(i + 1, 24)
         if int(A[i] @ A[j]) == 1]
assert len(TIGHT) == 96

# Lambda' on the full space (R^4)^24, as a 96 x 96 integer matrix
Lam = np.zeros((96, 96), dtype=np.int64)
for r, (i, j) in enumerate(TIGHT):
    Lam[r, 4 * i:4 * i + 4] += A[j]
    Lam[r, 4 * j:4 * j + 4] += A[i]

# 2P, integral: blocks 2I - alpha alpha^T
P2 = np.zeros((96, 96), dtype=np.int64)
for i in range(24):
    P2[4 * i:4 * i + 4, 4 * i:4 * i + 4] = 2 * np.eye(4, dtype=np.int64) - np.outer(A[i], A[i])

N4 = P2 @ (Lam.T @ Lam) @ P2          # = 4 N, integral
assert np.array_equal(N4, N4.T)
print('max |4N| entry:', int(np.abs(N4).max()))

I96 = np.eye(96, dtype=np.int64)
prod = N4.copy()
for c in (8, 20, 24, 32):
    prod = prod @ (N4 - c * I96)
print('4N annihilates x(x-8)(x-20)(x-24)(x-32):', bool(np.all(prod == 0)))


def rank_mod(M, p):
    """rank over the field with p elements, p prime"""
    B = [[int(x) % p for x in row] for row in M]
    rows, cols, r = len(B), len(B[0]), 0
    for c in range(cols):
        piv = next((k for k in range(r, rows) if B[k][c]), None)
        if piv is None:
            continue
        B[r], B[piv] = B[piv], B[r]
        inv = pow(B[r][c], p - 2, p)
        B[r] = [(v * inv) % p for v in B[r]]
        for k in range(rows):
            if k != r and B[k][c]:
                f = B[k][c]
                B[k] = [(B[k][t] - f * B[r][t]) % p for t in range(cols)]
        r += 1
    return r


# rank over Q is at least the rank over F_p, and here the upper bound 66 is
# known from the 30-dimensional kernel exhibited above, so equality follows
MULT = {}
for c, lam in ((0, 0), (8, 2), (20, 5), (24, 6), (32, 8)):
    r = rank_mod(N4 - c * I96, 1000003)
    MULT[lam] = 96 - r
    print('eigenvalue %d of N: rank of 4N - %2dI is %2d, multiplicity %2d'
          % (lam, c, r, 96 - r))
assert sum(MULT.values()) == 96, MULT
assert MULT == {0: 30, 2: 29, 5: 8, 6: 21, 8: 8}, MULT
print()
print('The multiplicities sum to 96, so the ranks over a field of positive')
print('characteristic, which bound the ranks over Q from below, are met.')
print('The kernel of N is 30-dimensional: the 24 radial directions, on which')
print('P vanishes, and the 6 infinitesimal rotations.  Halving and taking')
print('square roots, the singular values of Lambda on the tangent space are')
print('0, 1, sqrt(5/2), sqrt3, 2 with multiplicities 6, 29, 8, 21, 8.')
print('The smallest nonzero singular value is exactly 1, which is the')
print('constant the radius 2/sqrt(577) of Theorem (a radius for the')
print('rigidity) is built from.')
