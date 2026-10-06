#!/usr/bin/env python3
"""
c28_maxclose.py starts seed -- look for directions of 28 centres under the radial counts of
28 (21 within 2.1, one each at 2.15, 2.2, 2.35, four at sqrt6, each at its bin maximum)
with as many of the 21 as possible within 2.0161.  The 21 distances are variables
d_j = 2.1 - s_j (2.1 - 2.0161), s_j in [0, 1]; minimise -sum s_j + mu * (squared
violations) for increasing mu, then fix the centres with s_j > 0.98 at 2.0161 and the
rest at 2.1 and re-minimise the violations alone.  Floating point.
"""
import sys
import numpy as np
from scipy.optimize import minimize
starts, seed = map(int, sys.argv[1:3])
R6 = 6 ** .5
n = 28; iu = np.triu_indices(n, 1)
lo, hi = 2.0161, 2.1
fixed = np.array([2.15, 2.2, 2.35] + [R6] * 4)
def dist(s):
    return np.r_[hi - s * (hi - lo), fixed]
def amat(ds):
    return (ds[:, None] ** 2 + ds[None, :] ** 2 - 4) / (2 * ds[:, None] * ds[None, :])
def dadd(ds):  # derivative of a(d_i, d_j) in d_i
    return (ds[:, None] ** 2 - ds[None, :] ** 2 + 4) / (2 * ds[:, None] ** 2 * ds[None, :])
def make(mu):
    def f(z):
        X = z[:4 * n].reshape(n, 4); s = z[4 * n:]
        ds = dist(s)
        nr = np.linalg.norm(X, axis=1); W = X / nr[:, None]; G = W @ W.T
        A = amat(ds); V = G - A; np.fill_diagonal(V, -1)
        Vp = np.maximum(V, 0)
        val = -s.sum() + mu * (Vp[iu] ** 2).sum()
        Gr = 2 * mu * Vp  # symmetric, d val / d G_ij for i<j counted once each: use full symmetric /2
        gW = Gr @ W
        gX = (gW - (gW * W).sum(1)[:, None] * W) / nr[:, None]
        # d val / d d_i = sum_j 2 mu Vp_ij * (-dA_ij/dd_i) over j != i
        gd = -(2 * mu * Vp * dadd(ds)).sum(1)
        gs = -1 + gd[:21] * (-(hi - lo))
        return val, np.r_[gX.ravel(), gs]
    return f
def viol_only(ds):
    tv = amat(ds)[iu]
    def f(x):
        X = x.reshape(n, 4); nr = np.linalg.norm(X, axis=1); W = X / nr[:, None]; G = W @ W.T
        v = G[iu] - tv; vp = np.maximum(v, 0); val = (vp ** 2).sum()
        Gr = np.zeros((n, n)); Gr[iu] = 2 * vp; Gr = Gr + Gr.T; gW = Gr @ W
        gX = (gW - (gW * W).sum(1)[:, None] * W) / nr[:, None]
        return val, gX.ravel()
    return f, tv
opt = {'maxiter': 20000, 'gtol': 1e-13, 'ftol': 1e-15}
rng = np.random.default_rng(seed)
bounds = [(None, None)] * (4 * n) + [(0, 1)] * 21
best = {}
for st in range(starts):
    z = np.r_[rng.standard_normal(4 * n), rng.uniform(0, 1, 21)]
    for mu in (1e1, 1e2, 1e3, 1e4, 1e5):
        z = minimize(make(mu), z, jac=True, method='L-BFGS-B', bounds=bounds, options=opt).x
    s = z[4 * n:]
    for thr in (0.98, 0.9, 0.7):
        k = int((s > thr).sum())
        order = np.argsort(-s)
        ds = np.r_[np.where(np.isin(np.arange(21), order[:k]), lo, hi), fixed]
        f, tv = viol_only(ds)
        r = minimize(f, z[:4 * n], jac=True, method='L-BFGS-B', options=opt)
        X = r.x.reshape(n, 4); W = X / np.linalg.norm(X, axis=1)[:, None]
        mv = ((W @ W.T)[iu] - tv).max()
        if mv < 1e-8 and k > max(best, default=-1):
            np.save('maxclose_%d.npy' % k, r.x)
        if mv < 1e-8:
            best[k] = best.get(k, 0) + 1
    if (st + 1) % 50 == 0:
        print('after %d starts: feasible close counts %s' % (st + 1, dict(sorted(best.items()))), flush=True)
print('largest feasible close count: %s' % max(best, default=None), flush=True)
