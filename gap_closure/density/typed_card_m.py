#!/usr/bin/env python3
"""typed_card_m.py -- the typed three-point bound of multi_cap/typed_cardinality_sdp.py for
n points of type 1 (pairwise <= t1) plus m points of type 2
(<= t2 with type 1, <= t3 among themselves) on S^3; Z < 0 (corrected) excludes the code.
Exploration only (floating point, sampled constraints refined in rounds).

Usage: python3 typed_card_m.py d n m t1 t2 t3 [rounds]"""
import os, sys, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'multi_cap'))
import numpy as np
import cvxpy as cp
from three_point_sdp import leg_coeffs, Smat
from cardinality_sdp import Rows, grid_1d
from typed_cardinality_sdp import grid_111, grid_112, random_triples, F_values, g_values, climb


def coeffs(n, m):
    return dict(A=n + m, g1=n * (n - 1), g2=2 * n * m, g3=m * (m - 1), F1=n * (n - 1) * (n - 2),
                F2=3 * n * (n - 1) * m, F3=3 * n * m * (m - 1), F4=m * (m - 1) * (m - 2))


def doms(t1, t2, t3):
    # name -> (t_uv, t_w) for triples, (lo, hi) for pairs
    return dict(g1=t1, g2=t2, g3=t3, F1=(t1, t1), F2=(t2, t1), F3=(t2, t3), F4=(t3, t3))


def solve(d, C, S, keys):
    LC = [leg_coeffs(k) for k in range(d + 1)]
    R = Rows(d, LC)
    f = cp.Variable(d)
    Fs = [cp.Variable((d - k + 1, d - k + 1), symmetric=True) for k in range(d + 1)]
    bnd = {k: cp.Variable() for k in keys if k != 'g1'}
    cons = [f >= 0] + [F >> 0 for F in Fs]
    for k in keys:
        if k.startswith('g'):
            Af, AF = R.rows_i(S[k])
            expr = Af @ f + sum(AF[j] @ cp.vec(Fs[j], order='C') for j in range(d + 1))
        else:
            AF = R.rows_ii(*S[k])
            expr = sum(AF[j] @ cp.vec(Fs[j], order='C') for j in range(d + 1))
        cons.append(expr <= (-1 if k == 'g1' else bnd[k]))
    one = np.array([1.0])
    S111 = [Smat(k, d, one, one, one, LC)[0] for k in range(d + 1)]
    A = cp.sum(f) + sum(cp.sum(cp.multiply(S111[k], Fs[k])) for k in range(d + 1))
    cons.append(A <= 30)
    Z = C['A'] * A - C['g1'] + sum(C[k] * bnd[k] for k in bnd)
    prob = cp.Problem(cp.Minimize(Z), cons)
    t0 = time.time()
    opts = dict(tol_gap_abs=1e-7, tol_gap_rel=1e-7, tol_feas=1e-7, max_iter=500,
                static_regularization_constant=1e-6, equilibrate_max_iter=50, tol_ktratio=1e-6)
    prob.solve(solver='CLARABEL', **opts)
    return prob, f.value, [F.value for F in Fs], {k: bnd[k].value for k in bnd}, R, time.time() - t0


def evaluate(R, fv, Fv, C, D, keys, n_climb=600):
    d, LC = R.d, R.LC
    out = {}
    for k in keys:
        if k.startswith('g'):
            u = np.r_[np.linspace(-1, D[k], 80001), D[k]]
            val = g_values(fv, Fv, d, LC, u)
            out[k] = (float(val.max()), u[np.argsort(val)[-800:]])
        else:
            tuv, tw = D[k]
            P = (grid_112(tw, tuv, 90, 40) if tuv != tw else grid_111(tuv, 90, 40),
                 random_triples(250000, tuv, tw))
            vals, pts = [], [[], [], []]
            for (uu, vv, ww) in P:
                for s0 in range(0, len(uu), 50000):
                    a, b, c = uu[s0:s0 + 50000], vv[s0:s0 + 50000], ww[s0:s0 + 50000]
                    val = F_values(Fv, d, LC, a, b, c)
                    idx = np.argsort(val)[-3000:]
                    vals.append(val[idx]); pts[0].append(a[idx]); pts[1].append(b[idx]); pts[2].append(c[idx])
            vals = np.concatenate(vals); pts = [np.concatenate(x) for x in pts]
            order = np.argsort(vals)[::-1]
            top = order[:n_climb]
            cv, cp_ = climb(Fv, d, LC, pts[0][top], pts[1][top], pts[2][top], tuv, tw)
            best = max(float(vals.max()), float(cv.max()))
            keep = order[:3000]
            out[k] = (best, tuple(np.r_[pts[j][keep], cp_[j]] for j in range(3)))
    one = np.array([1.0])
    A = float(np.sum(fv)) + sum(float(np.sum(Smat(k, d, one, one, one, LC)[0] * Fv[k])) for k in range(d + 1))
    Z = C['A'] * A + sum(C[k] * out[k][0] for k in keys)
    return Z, A, out


def run(d, n, m, t1, t2, t3, rounds):
    C = coeffs(n, m); D = doms(t1, t2, t3)
    keys = [k for k in ('g1', 'g2', 'g3', 'F1', 'F2', 'F3', 'F4') if C[k] > 0]
    S = {}
    for k in keys:
        if k.startswith('g'):
            S[k] = grid_1d(D[k], 300)
        else:
            tuv, tw = D[k]
            S[k] = grid_112(tw, tuv, 24, 14) if tuv != tw else grid_111(tuv, 24, 14)
    for r in range(rounds):
        prob, fv, Fv, bv, R, el = solve(d, C, S, keys)
        if fv is None:
            print('  solver status', prob.status, flush=True); return
        Z, A, out = evaluate(R, fv, Fv, C, D, keys)
        print('  d=%d n=%d m=%d t=(%.5f,%.5f,%.5f) round %d: %s, sampled Z %.5f, corrected Z %.5f; A %.4f; %s [%.0f s]'
              % (d, n, m, t1, t2, t3, r + 1, prob.status, prob.value, Z, A,
                 ' '.join('%s %.3e' % (k, out[k][0]) for k in keys), el), flush=True)
        if Z < 0:
            print('  corrected Z below 0: excluded in floating point', flush=True); break
        if prob.value > 1e-4:
            print('  sampled Z above 0: no certificate of degree %d' % d, flush=True); break
        for k in keys:
            if k.startswith('g'):
                S[k] = np.r_[S[k], out[k][1]]
            else:
                S[k] = tuple(np.r_[S[k][j], out[k][1][j]] for j in range(3))


if __name__ == '__main__':
    d, n, m = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    t1, t2, t3 = float(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6])
    rounds = int(sys.argv[7]) if len(sys.argv) > 7 else 3
    run(d, n, m, t1, t2, t3, rounds)
