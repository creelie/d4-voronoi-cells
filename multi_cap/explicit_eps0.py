#!/usr/bin/env python3
"""
explicit_eps0.py -- an explicit value for the constant epsilon_0 of the
near-contact theorem (Section 2.8 of the paper, v1.8.0), in exact and ball
arithmetic.

The theorem says: if no centre of a unit-ball packing of R^4 lies at a distance
from c strictly between 2 + epsilon_0 and sqrt 6, then vol(V_c) >= 8, with
equality only at the 24-cell.  Its constant used to come from a compactness
argument.  Here every step is quantitative.

 A. The robust second-level bound.  The certificate of de Laat, Leijenhorst and
    de Muinck Keizer is an exact identity
        A_2K(Q) + sum_j c_j w_j(Q) sigma_j(Q) = rhs_|Q|,    |Q| = 1, ..., 4,
    with sigma_j = v_j^T T X T^T v_j (X positive definite), c_j >= 0, and each
    w_j a constant, a Gram determinant or a principal minor (nonnegative on
    every real configuration) or an elementary symmetric function e_r of the
    pair weights p(u) = (u + 1)(1/2 - u).  Only the last can be negative, and
    only when an inner product exceeds 1/2.  If every inner product is at most
    1/2 + kappa, then p >= -a with a = kappa (3/2 + kappa), and |p| <= b = 9/16,
    so e_r^- <= C(n, r) a b^(r-1).  The script bounds, in ball arithmetic,
        S_{k,r} = sum over the blocks with weight e_r of c_j sup |sigma_j|,
    using lambda_max(X) <= min(trace X, max row sum |X|) and the l1 norm of
    the tensor Chebyshev coefficients of T^T v_j (|inner products| <= 1).
 B. The two-point polynomial.  p_2 = (u+1)(u+1/2)^2 u^2 (u-1/2) q(u) exactly,
    and -q >= q_min > 0 on [-1, 1/2] (mean value form on 30 000 intervals).
 C. The robust root-lattice step.  If the Gram matrix G of 24 unit vectors is
    within delta of a matrix g with off-diagonal entries in {-1,-1/2,0,1/2},
    then every 5 x 5 minor of the integral matrix 2g differs from the
    corresponding minor of 2G, which is 0, by less than 1, so it is 0; every
    principal minor of order <= 4 exceeds -1, so it is >= 0; so 2g is positive
    semidefinite of rank <= 4 and the vectors are near a copy of the
    normalised roots: d(W) <= 2 D / sqrt(6 - D) + sqrt(24) D / 6, D = 23 delta.
 D. The local volume bound near the root system, with explicit constants.
 E. Assembly: kappa*, epsilon_0 = 2 kappa*.

Usage: python3 explicit_eps0.py /path/to/LasserreSphericalCodes/proofs/4_24
Needs python-flint.  Exit status 0 when every check passes.
"""
import os
import sys
import time
from fractions import Fraction as Fr
from itertools import combinations
from math import comb

from flint import arb, arb_poly, fmpq, fmpq_poly, ctx as fctx

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'zonal'))
sys.set_int_max_str_digits(0)
import verify45 as V                       # the data readers of the zonal verification

fctx.prec = 256
OK = True
NV = [1, 1, 3, 6]


def check(cond, msg):
    global OK
    print(('  PASS  ' if cond else '  FAIL  ') + msg, flush=True)
    OK &= bool(cond)


def aq(q):
    return arb(q.p) / arb(q.q)


def up(x):
    """an upper bound of an arb, as an arb with zero radius"""
    return arb(x.mid() + x.rad())


def amax(x, y):
    return x if (x - y) > 0 else y


# ------------------------------------------------------------------ A
_MC = {}


def mono_cheb(n):
    """x^n = sum_k a_k T_k(x), exactly"""
    if n not in _MC:
        out = {}
        for j in range(n + 1):
            k = abs(n - 2 * j)
            out[k] = out.get(k, fmpq(0)) + fmpq(comb(n, j), 2 ** n)
        _MC[n] = out
    return _MC[n]


def cheb_l1(poly):
    """l1 norm of the tensor Chebyshev coefficients: a bound for sup |poly| on [-1,1]^n"""
    out = {}
    for e, c in poly.items():
        terms = {(): c}
        for n in e:
            new = {}
            for key, v in terms.items():
                for k, a in mono_cheb(n).items():
                    kk = key + (k,)
                    new[kk] = new.get(kk, fmpq(0)) + v * a
            terms = new
        for key, v in terms.items():
            out[key] = out.get(key, fmpq(0)) + v
    s = arb(0)
    for v in out.values():
        if v != 0:
            s += abs(aq(v))
    return s


def lambda_max_bound(path):
    with open(path) as f:
        nr, nc = [int(x) for x in f.readline().split()]
        tr = arb(0)
        rowmax = arb(0)
        for i in range(nr):
            rs = arb(0)
            for j, t in enumerate(f.readline().split()):
                a = aq(V.parse_q(t))
                rs += abs(a)
                if j == i:
                    tr += a
            rowmax = amax(up(rs), rowmax)
    return (up(tr) if (tr - rowmax) < 0 else rowmax), nr


def read_T(path):
    with open(path) as f:
        nr, nc = [int(x) for x in f.readline().split()]
        return [[V.parse_q(t) for t in f.readline().split()] for _ in range(nr)]


def weights(k):
    C = V.ctx(NV[k - 1])
    nv = NV[k - 1]
    tr = V.gram(k, nv)
    xs = [C.gen(i) for i in range(nv)]
    one = C.from_dict({(0,) * nv: fmpq(1)})
    half = C.from_dict({(0,) * nv: fmpq(1, 2)})
    out = []
    if k == 2:
        return [('e1_pair', (xs[0] + one) * (half - xs[0]))]
    out.append(('det', V.det(tr, C)))
    pw = [(y + one) * (half - y) for y in xs]
    for r in range(1, len(pw) + 1):
        g = C.from_dict({})
        for S in combinations(range(len(pw)), r):
            t = one
            for i in S:
                t = t * pw[i]
            g = g + t
        out.append(('e%d_pair' % r, g))
    if k == 4:
        mins = [V.det([[tr[a][b] for b in S] for a in S], C) for S in combinations(range(4), 3)]
        for r in range(1, 5):
            g = C.from_dict({})
            for S in combinations(range(4), r):
                t = one
                for i in S:
                    t = t * mins[i]
                g = g + t
            out.append(('e%d_minor' % r, g))
    return out


def classify(pf, W, nv):
    e = V.poly_from_dict(pf, nv)
    terms = dict(e.terms())
    if not terms:
        return ('zero', fmpq(0))
    if V.total_degree(e) == 0:
        return ('const', list(terms.values())[0])
    for name, w in W:
        st = dict(w.terms())
        if set(st) != set(terms):
            continue
        rat = {terms[m] / st[m] for m in st}
        if len(rat) == 1:
            return (name, rat.pop())
    raise RuntimeError('a prefactor that is not a multiple of a domain weight')


def part_A(D):
    print('A. sizes of the sum-of-squares terms whose weight can change sign', flush=True)
    S = {}
    for k in (3, 4):
        W = weights(k)
        nv = NV[k - 1]
        for fn in sorted(os.listdir(D)):
            if not fn.startswith('poly_%d_' % k):
                continue
            name = fn[len('poly_%d_' % k):-4]
            if not os.path.exists(os.path.join(D, 'dense_' + name + '.txt')):
                continue                                   # no block: the term is absent
            pref, vecs = V.load_lowrank(os.path.join(D, fn), nv)
            cls = [classify(pf, W, nv) for pf in pref]
            for c in cls:
                assert c[1] >= 0, 'negative multiple of a weight'
            if not any(c[0].endswith('_pair') for c in cls):
                continue
            lam, m = lambda_max_bound(os.path.join(D, 'dense_' + name + '.txt'))
            T = read_T(os.path.join(D, 'transform_' + name + '.txt'))
            for (wname, c), v in zip(cls, vecs):
                assert wname.endswith('_pair')
                s2 = arb(0)
                for b in range(m):
                    acc = {}
                    for a in range(len(T)):
                        tab = T[a][b]
                        if tab == 0 or not v[a]:
                            continue
                        for e, q in v[a].items():
                            acc[e] = acc.get(e, fmpq(0)) + tab * q
                    l = cheb_l1(acc)
                    s2 += l * l
                key = (k, int(wname[1]))
                S[key] = S.get(key, arb(0)) + aq(c) * lam * s2
    for key in sorted(S):
        S[key] = up(S[key])
        print('    S_{%d,%d} <= %s' % (key[0], key[1], S[key].str(6, radius=False)))
    return S


# ------------------------------------------------------------------ B
def part_B():
    print('B. the two-point polynomial', flush=True)
    c = [fmpq(0)] * 17
    for line in open(os.path.join(HERE, 'llm24_out', 'llm24_p2.txt')):
        if line.startswith('#') or not line.strip():
            continue
        i, v = line.split(None, 1)
        v = v.strip()
        c[int(i)] = fmpq(*[int(x) for x in v.split('/')]) if '/' in v else fmpq(int(v))
    p = fmpq_poly(c)
    x = fmpq_poly([0, 1])
    h = fmpq(1, 2)
    den = (x + 1) * (x + h) ** 2 * x ** 2 * (x - h)
    q, r = divmod(p, den)
    check(r == 0 and p.degree() == 16, 'p_2 = (u+1)(u+1/2)^2 u^2 (u-1/2) q(u) exactly, deg q = %d' % q.degree())
    A = arb_poly([aq(-q[i]) for i in range(q.degree() + 1)])
    dA = A.derivative()
    N = 30000
    lo = None
    for k in range(N):
        a = Fr(-1) + Fr(3, 2) * Fr(k, N)
        b = Fr(-1) + Fr(3, 2) * Fr(k + 1, N)
        m = arb(fmpq((a + b).numerator, (a + b).denominator)) / 2
        rad = arb(fmpq((b - a).numerator, (b - a).denominator)) / 2
        val = A(m) + dA(arb(m.mid(), rad.mid() + rad.rad())) * arb(0, rad.mid() + rad.rad())
        lb = val.mid() - val.rad()
        lo = lb if lo is None or lb < lo else lo
    qmin = arb(lo)
    check(qmin > 0, '-q >= q_min = %s > 0 on [-1, 1/2]' % qmin.str(6, radius=False))
    P = arb_poly([aq(p[i]) for i in range(17)])
    return P, qmin


# ------------------------------------------------------------------ C
def part_C(delta):
    """returns an upper bound for d(W) when |G - g| <= delta entrywise"""
    print('C. the robust root-lattice step at delta = %s' % delta, flush=True)
    d = arb(fmpq(delta.numerator, delta.denominator))
    row = (arb(8) + 3 * (1 + 2 * d) ** 2).sqrt()             # row norm of a 5-column piece of 2G
    ok = True
    for k in (1, 2, 3, 4, 5):
        err = (row + 2 * d * arb(k).sqrt()) ** k - row ** k
        ok &= bool(err < 1)
    check(ok, 'every k x k minor (k <= 5) of 2g is within less than 1 of that of 2G')
    Dn = 23 * d
    dW = 2 * Dn / (6 - Dn).sqrt() + arb(24).sqrt() * Dn / 6
    return up(dW)


# ------------------------------------------------------------------ D
S11 = Fr(33166248, 10 ** 7)                                  # > sqrt 11
CA_EXACT = (S11 / 2) / (1 - (3 * S11 + Fr(1, 2)) / 48)      # > the constant c_a of estimate (a)
CA = CA_EXACT
ONE_D = Fr(100006, 100000)                                   # 1 + delta, delta <= S <= 6e-5
assert S11 * S11 > 11 and CA < Fr(21199, 10000)


def local_bracket(Th):
    """Step 2 of the proof of the near-contact theorem with explicit constants: with e = max e_i,
    delta = max delta_i and e + delta <= Th, vol(V(Y) cap K(Y)) - 8 >= S * bracket - C2 * Th^4,
    S = sum delta_i.  Every quantity below is increasing in Th, so evaluating at an upper bound is safe."""
    s2 = Fr(14143, 10000)                                   # > sqrt 2
    sigma = (Th / 2 + s2 * Th) / (1 - Th)                   # the lift onto the moved hyperplane
    lam0 = 2 * (s2 * Th + sigma * (1 + Th))                 # inner homothety: (1 - lam0) F_i^0 lifts into F_i
    r1 = s2 * (1 + Th / 2) / (1 - s2 * Th)                  # P(t) lies in B(r1)
    eta1 = Th / 2 + r1 * Th                                 # |<x, u_i> - 1| on F_i
    eta2 = Th / 2 + r1 * Th + eta1                          # outer homothety 1 + 2 eta2
    cosphi = (1 - Th) / (1 + Th)
    Gup = Fr(4, 3) * (1 + 2 * eta2) ** 3
    Alow = Fr(4, 3) * (1 - lam0) ** 3
    ring = Fr(4, 3) * ((1 + 2 * eta2) ** 3 - (1 - lam0) ** 3)
    ymax = s2 * (1 + 2 * eta2) + eta1
    dvol = Fr(4, 3) * max((1 + 2 * eta2) ** 3 - 1, 1 - (1 - lam0) ** 3)
    mu = ((1 / cosphi - 1) * Gup * ymax + ring * (1 + 2 * eta2) + dvol + eta1 * Gup) / (1 - Th * Th / 4)
    sumE = Fr(208, 10)                                      # sum e_i <= sqrt 24 c_a 2 (1 + delta) S < 20.8 S
    bracket = Alow / 2 - mu * sumE                          # (delta_i / 2) A_i - e_i mu, summed
    C2 = 24 * 2 * 4 * (Fr(3, 2) + Fr(142, 100)) ** 4        # 24 vertex caps of volume 2 s^4, s = sqrt2 (theta + theta')
    return bracket, C2, lam0 / Th, mu / Th


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(2)
    t0 = time.time()
    S = part_A(sys.argv[1])
    P, qmin = part_B()
    delta = Fr(3, 10000)
    dW = part_C(delta)
    check(dW < arb(1) / 48, 'then d(W) <= %s < 1/48, the radius of the rigidity estimate' % dW.str(6, radius=False))

    # f(u) = (u+1)(u+1/2)^2 u^2 (1/2-u) is unimodal between consecutive roots, so its least value at distance
    # >= delta from the roots is at one of six points
    f = lambda u: (u + 1) * (u + Fr(1, 2)) ** 2 * u ** 2 * (Fr(1, 2) - u)
    fmin = min(f(Fr(-1) + delta), f(Fr(-1, 2) - delta), f(Fr(-1, 2) + delta), f(-delta), f(delta), f(Fr(1, 2) - delta))
    need = qmin * aq(fmpq(fmin.numerator, fmin.denominator))
    print('E. assembly: a pair at distance >= delta from {-1,-1/2,0,1/2} has p_2 >= %s' % need.str(6, radius=False))

    b = arb(9) / 16
    B3 = sum((arb(comb(3, r)) * b ** (r - 1) * S[(3, r)] for r in range(1, 4) if (3, r) in S), arb(0))
    B4 = sum((arb(comb(6, r)) * b ** (r - 1) * S[(4, r)] for r in range(1, 7) if (4, r) in S), arb(0))
    # p_2 on (1/2, 1/2 + kappa]: p_2(1/2) = 0 and |p_2'| <= L2 there (kappa <= 1e-3)
    L2 = P.derivative()(arb(fmpq(1, 2) + fmpq(1, 2000), fmpq(1, 2000)))
    L2 = up(abs(L2))

    def E(m, kap):
        a = kap * (arb(3) / 2 + kap)
        return comb(m, 2) * kap * L2 + comb(m, 3) * a * B3 + comb(m, 4) * a * B4

    # the largest kappa of the form 10^-k * j with E(24, kappa) < need
    kap = None
    for k in range(20, 40):
        for j in (9, 8, 7, 6, 5, 4, 3, 2, 1):
            cand = arb(j) / arb(10) ** k
            if E(24, cand) < need:
                kap = (j, k)
                break
        if kap:
            break
    kstar = arb(kap[0]) / arb(10) ** kap[1]
    print('    B_3 <= %s, B_4 <= %s, |p_2\'| <= %s near 1/2' % (B3.str(5, radius=False), B4.str(5, radius=False), L2.str(5, radius=False)))
    print('    E(24, kappa) <= %s at kappa* = %d e-%d' % (up(E(24, kstar)).str(6, radius=False), kap[0], kap[1]))
    check(E(24, kstar) < need, 'E(24, kappa*) < %s: every inner product of 24 such directions is within delta of the set' % need.str(4, radius=False))
    check(E(25, kstar) < 1, 'E(25, kappa*) < 1: no 25 directions with inner products <= 1/2 + kappa*')
    check(kstar < arb(fmpq(delta.numerator, delta.denominator)), 'kappa* < delta (the pairs above 1/2 are within delta of 1/2)')
    eps0 = 2 * kstar
    # kappa(eps) = 1/2 - 2/(2+eps)^2 <= eps/2
    check(arb(1) / 2 - 2 / (2 + eps0) ** 2 <= kstar, 'epsilon_0 = 2 kappa* = %d e-%d: 1/2 - 2/(2 + epsilon_0)^2 <= kappa*'
          % (2 * kap[0], kap[1]) if 2 * kap[0] < 10 else 'epsilon_0 = 2 kappa*: 1/2 - 2/(2 + epsilon_0)^2 <= kappa*')

    Smax = 24 * Fr(2 * kap[0], 10 ** kap[1])
    TH = Fr(524, 100)                                       # Theta <= 2 c_a (1 + eps) S + S <= 5.24 S
    check(2 * CA * ONE_D + 1 <= TH and Fr(4899, 1000) * 2 * CA * ONE_D < Fr(208, 10),
          'estimate (a): c_a = %.6f < 2.1199, so Theta <= 5.24 S and sum e_i <= %.3f S < 20.8 S' % (float(CA), float(Fr(4899, 1000) * 2 * CA * ONE_D)))
    Thb = TH * Smax
    bracket, C2, lr, mr = local_bracket(Thb)
    print('D. local bound: vol >= 8 + S (bracket - C2 (5.24)^4 S^3), bracket = %.6f, C2 = %.0f; lam0 <= %.2f Theta, mu <= %.1f Theta'
          % (float(bracket), float(C2), float(lr), float(mr)))
    check(bracket - C2 * TH ** 4 * Smax ** 3 > 0,
          'at every S = sum delta_i in (0, 24 epsilon_0] the bracket is positive: vol > 8; at S = 0 the root system')
    # the same bound on a fixed neighbourhood of the root system, independent of the certificate
    Sl = Fr(6, 100000)
    br, C2l, _, _ = local_bracket(TH * Sl)
    check(TH * Sl <= Fr(1, 1000) and br - C2l * TH ** 4 * Sl ** 3 >= Fr(2, 100),
          'the same holds whenever d(W) <= 1/48 and 0 < S <= 6e-5 (bracket %.4f at S = 6e-5)' % float(br - C2l * TH ** 4 * Sl ** 3))
    print('elapsed %.0f s' % (time.time() - t0))
    if OK:
        print('PASS: epsilon_0 = %s works in the near-contact theorem' % eps0.str(3, radius=False))
        sys.exit(0)
    print('FAILED')
    sys.exit(1)


if __name__ == '__main__':
    main()
