#!/usr/bin/env python3
"""typed_card_hole.py -- three-point bound for 24 points (pairwise <= t1) plus a point y (<= t2 with each),
when exactly mp of the 24 have inner product above t1 with y (the others <= t1).
Z(mp) < 0 (corrected) excludes that case.  Tammes on the link of y: mp <= 9 for
t1 = 0.508, t2 = 0.6141.  Exploration (floating point, sampled, rounds).

Usage: python3 typed_card_hole.py d t1 t2 rounds mp1 mp2 ..."""
import os, sys, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'multi_cap'))
import numpy as np
import cvxpy as cp
from three_point_sdp import leg_coeffs, Smat
from cardinality_sdp import Rows, grid_1d
from typed_cardinality_sdp import grid_111, F_values, g_values

rng = np.random.default_rng(5)
N = 24


def wrange(u, v):
    r = np.sqrt(np.clip((1 - u * u) * (1 - v * v), 0, None))
    return u * v - r, u * v + r


def grid_box(ur, vr, whi, n_uv, n_w, sym):
    U, V, W = [], [], []
    a = np.linspace(ur[0], ur[1], n_uv); b = np.linspace(vr[0], vr[1], n_uv)
    for i, u in enumerate(a):
        for j, v in enumerate(b):
            if sym and v < u:
                continue
            lo, hi = wrange(u, v)
            hi = min(hi, whi)
            if hi < lo:
                continue
            for w in np.linspace(lo, hi, n_w):
                U.append(u); V.append(v); W.append(w)
    return np.array(U), np.array(V), np.array(W)


def random_box(n, ur, vr, whi):
    u = rng.uniform(*ur, n); v = rng.uniform(*vr, n)
    lo, hi = wrange(u, v); hi = np.minimum(hi, whi)
    ok = hi >= lo
    u, v, lo, hi = u[ok], v[ok], lo[ok], hi[ok]
    w = lo + (hi - lo) * rng.random(len(u))
    # also push some to the boundary of the Gram domain and to w = whi
    k = len(u) // 4
    w[:k] = hi[:k]; w[k:2 * k] = lo[k:2 * k]
    return u, v, w


def climb(Fv, d, LC, u, v, w, ur, vr, whi, steps=(3e-2, 1e-2, 3e-3, 1e-3, 3e-4, 1e-4), tries=16, reps=4):
    def to_w(a, b, ph):
        return a * b + np.sqrt(np.clip((1 - a * a) * (1 - b * b), 0, None)) * np.cos(ph)

    def project(a, b, ph):
        a = np.clip(a, *ur); b = np.clip(b, *vr); ph = np.clip(ph, 0, np.pi)
        r = np.sqrt(np.clip((1 - a * a) * (1 - b * b), 0, None))
        ww = a * b + r * np.cos(ph)
        over = ww > whi
        if over.any():
            c = np.where(r[over] > 1e-12, (whi - a[over] * b[over]) / np.maximum(r[over], 1e-12), -1.0)
            ph = ph.copy(); ph[over] = np.arccos(np.clip(c, -1, 1))
            ww = to_w(a, b, ph)
        return a, b, ph, ww
    r = np.sqrt(np.clip((1 - u * u) * (1 - v * v), 1e-300, None))
    ph = np.arccos(np.clip((w - u * v) / r, -1, 1))
    a, b, ph, ww = project(u.copy(), v.copy(), ph)
    best = F_values(Fv, d, LC, a, b, ww)
    m = len(a)
    for st in steps:
        for _ in range(reps):
            D = rng.normal(size=(tries, m, 3)) * st
            ca, cb, cph, cw = project((a + D[..., 0]).ravel(), (b + D[..., 1]).ravel(), (ph + D[..., 2]).ravel())
            ok = cw <= whi + 1e-15
            val = np.where(ok, F_values(Fv, d, LC, ca, cb, cw), -np.inf).reshape(tries, m)
            j = val.argmax(0); vb = val[j, np.arange(m)]
            up = vb > best; idx = j * m + np.arange(m)
            a = np.where(up, ca[idx], a); b = np.where(up, cb[idx], b)
            ph = np.where(up, cph[idx], ph); best = np.where(up, vb, best)
    return best, (a, b, to_w(a, b, ph))


def run(d, t1, t2, mp, rounds, Acap=30.0):
    LC = [leg_coeffs(k) for k in range(d + 1)]
    R = Rows(d, LC)
    # triple kinds: name -> (ur, vr, whi, sym, coefficient)
    kinds = {'b1': ((-1, t1), (-1, t1), t1, True, N * (N - 1) * (N - 2) + 6 * (N - mp) * (N - mp - 1) // 2)}
    if mp >= 1:
        kinds['b10'] = ((t1, t2), (-1, t1), t1, False, 6 * mp * (N - mp))
    if mp >= 2:
        kinds['b11'] = ((t1, t2), (t1, t2), t1, True, 6 * mp * (mp - 1) // 2)
    S1 = grid_1d(t1, 300)
    S2 = np.linspace(t1, t2, 200) if mp >= 1 else None
    ST = {}
    for k, (ur, vr, whi, sym, c) in kinds.items():
        ST[k] = grid_111(t1, 24, 14) if k == 'b1' else grid_box(ur, vr, whi, 20, 12, sym)
    for rnd in range(rounds):
        f = cp.Variable(d)
        Fs = [cp.Variable((d - k + 1, d - k + 1), symmetric=True) for k in range(d + 1)]
        bv = {k: cp.Variable() for k in kinds}
        a2 = cp.Variable()
        cons = [f >= 0] + [F >> 0 for F in Fs]
        Af, AF = R.rows_i(S1)
        cons.append(Af @ f + sum(AF[j] @ cp.vec(Fs[j], order='C') for j in range(d + 1)) <= -1)
        if mp >= 1:
            Af, AF = R.rows_i(S2)
            cons.append(Af @ f + sum(AF[j] @ cp.vec(Fs[j], order='C') for j in range(d + 1)) <= a2)
        for k in kinds:
            AF = R.rows_ii(*ST[k])
            cons.append(sum(AF[j] @ cp.vec(Fs[j], order='C') for j in range(d + 1)) <= bv[k])
        one = np.array([1.0])
        S111 = [Smat(k, d, one, one, one, LC)[0] for k in range(d + 1)]
        A = cp.sum(f) + sum(cp.sum(cp.multiply(S111[k], Fs[k])) for k in range(d + 1))
        cons.append(A <= Acap)
        Z = (N + 1) * A - N * (N - 1) - 2 * (N - mp) + 2 * mp * a2 + sum(kinds[k][4] * bv[k] for k in kinds)
        prob = cp.Problem(cp.Minimize(Z), cons)
        t0 = time.time()
        prob.solve(solver='CLARABEL', tol_gap_abs=1e-7, tol_gap_rel=1e-7, tol_feas=1e-7, max_iter=500,
                   static_regularization_constant=1e-6, equilibrate_max_iter=50, tol_ktratio=1e-6)
        el = time.time() - t0
        if f.value is None:
            print('  status', prob.status, flush=True); return
        fv = f.value; Fv = [F.value for F in Fs]
        # corrected values
        u = np.r_[np.linspace(-1, t1, 80001)]
        g1 = g_values(fv, Fv, d, LC, u); g1m = float(g1.max()); add1 = u[np.argsort(g1)[-600:]]
        corr = {}
        add = {}
        if mp >= 1:
            u2 = np.linspace(t1, t2, 40001)
            g2 = g_values(fv, Fv, d, LC, u2); corr['a2'] = float(g2.max()); add2 = u2[np.argsort(g2)[-400:]]
        for k, (ur, vr, whi, sym, c) in kinds.items():
            P = [grid_box(ur, vr, whi, 70, 30, sym), random_box(300000, ur, vr, whi)]
            vals, pts = [], [[], [], []]
            for (uu, vv, ww) in P:
                for s0 in range(0, len(uu), 20000):
                    aa, bb, cc = uu[s0:s0 + 20000], vv[s0:s0 + 20000], ww[s0:s0 + 20000]
                    val = F_values(Fv, d, LC, aa, bb, cc)
                    idx = np.argsort(val)[-1500:]
                    vals.append(val[idx]); pts[0].append(aa[idx]); pts[1].append(bb[idx]); pts[2].append(cc[idx])
            vals = np.concatenate(vals); pts = [np.concatenate(x) for x in pts]
            order = np.argsort(vals)[::-1]; top = order[:500]
            cv, cpt = climb(Fv, d, LC, pts[0][top], pts[1][top], pts[2][top], ur, vr, whi)
            corr[k] = max(float(vals.max()), float(cv.max()))
            keep = order[:2500]
            add[k] = tuple(np.r_[pts[j][keep], cpt[j]] for j in range(3))
        Av = float(np.sum(fv)) + sum(float(np.sum(Smat(k, d, one, one, one, LC)[0] * Fv[k])) for k in range(d + 1))
        Zc = (N + 1) * Av + N * (N - 1) * g1m + 2 * (N - mp) * g1m + (2 * mp * corr['a2'] if mp >= 1 else 0) \
            + sum(kinds[k][4] * corr[k] for k in kinds)
        print('  d=%d mp=%d round %d: %s, sampled Z %.5f, corrected Z %.5f; A %.4f, g1 max %.3e, %s [%.0f s]'
              % (d, mp, rnd + 1, prob.status, prob.value, Zc, Av, g1m,
                 ' '.join('%s %.4e (var %.4e)' % (k, corr[k], (a2.value if k == 'a2' else bv[k].value)) for k in corr), el), flush=True)
        if Zc < 0:
            print('  corrected Z below 0: case mp=%d excluded in floating point' % mp, flush=True)
            np.savez('hole_d%d_m%d.npz' % (d, mp), d=d, t1=t1, t2=t2, mp=mp, f=fv, **{'F%d' % k: Fv[k] for k in range(d + 1)})
            return Zc
        if prob.value > 1e-4:
            print('  sampled Z above 0: no certificate of degree %d for mp=%d' % (d, mp), flush=True)
            return prob.value
        S1 = np.r_[S1, add1]
        if mp >= 1:
            S2 = np.r_[S2, add2]
        for k in kinds:
            ST[k] = tuple(np.r_[ST[k][j], add[k][j]] for j in range(3))
    return None


if __name__ == '__main__':
    d = int(sys.argv[1]); t1, t2 = float(sys.argv[2]), float(sys.argv[3])
    rounds = int(sys.argv[4])
    for mp in map(int, sys.argv[5:]):
        run(d, t1, t2, mp, rounds)
