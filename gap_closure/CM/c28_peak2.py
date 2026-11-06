"""c28_peak2.py -- the largest value of each triple polynomial of a (mixed) certificate, by dense
sampling and local ascent; a faster form of c28_peak.py (power tables, analytic gradients).
Floating point; it only tells how far the thresholds of an exact check must sit above the
sampled maxima (TRIPLE_EXTRA of multi_cap/combo_case_check.py).

usage: python3 c28_peak2.py case.json cert.npz KIND[,KIND...] SAMPLES_MILLIONS [STARTS]
For each kind it prints the largest value on SAMPLES_MILLIONS million admissible triples
(sampled as in the check, from dimensions 4, 3 and 2) and after SLSQP ascent from the
STARTS (default 400) best of them, and the threshold the check computes with its own
650 000 samples and 40 starts, for comparison.

With MIXTYPE=1 the certificate is a mixed-typing combination (c27_hull_mt.py) on a case file
with one type per bin (case_bins.py), read through multi_cap/combo_mixtype_check.py, and a
KIND names its bins' types separated by '-', e.g. A1-B3-F7 (printed as A1B3F7, the key the
PEAKS option of the mixing programmes reads).
"""
import sys, os, time
import numpy as np
from scipy.optimize import minimize
case, npz = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2])
kinds, nm = sys.argv[3].split(','), float(sys.argv[4])
starts = int(sys.argv[5]) if len(sys.argv) > 5 else 400
MC = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'multi_cap')
sys.path.insert(0, MC)
os.chdir(MC)
sys.argv = ['x', case, npz, '10']
if os.environ.get('MIXTYPE') == '1':
    import combo_mixtype_check  # noqa: F401  (patches the three-point polynomials of C)
import combo_case_check as K
C = K.C
t0 = time.time()
Z = np.load(npz)
Bk = {name: C.psd_exact(M) for name, M in C.blocks_from_x3(np.array(Z['x3'], float), 10).items()}


def tables(P):
    E = np.array(list(P.keys())); co = np.array([float(c) for c in P.values()])
    return E, co, E.max()


def ev(E, co, D, g):
    X = [np.ones(len(g))]; Y = [np.ones(len(g))]; W = [np.ones(len(g))]
    for _ in range(D):
        X.append(X[-1] * g[:, 0]); Y.append(Y[-1] * g[:, 1]); W.append(W[-1] * g[:, 2])
    out = np.zeros(len(g))
    for (a, b, c), w in zip(E, co):
        out += w * X[a] * Y[b] * W[c]
    return out


def f_and_grad(E, co, x):
    p = x[None, :] ** E
    val = co @ np.prod(p, 1)
    gr = np.zeros(3)
    for i in range(3):
        Ei = E[:, i]
        dp = np.where(Ei > 0, Ei * x[i] ** np.maximum(Ei - 1, 0), 0.0)
        q = p.copy(); q[:, i] = dp
        gr[i] = co @ np.prod(q, 1)
    return val, gr


for kind in kinds:
    combo = tuple(kind.split('-')) if '-' in kind else tuple(kind)
    P = C.triple3_poly(Bk, 10, combo)
    E, co, D = tables(P)
    T = tuple(float(C.tmax(a, b)) for a, b in ((combo[0], combo[1]), (combo[0], combo[2]), (combo[1], combo[2])))
    rng = np.random.default_rng(4321)
    topg = np.zeros((0, 3)); topv = np.zeros(0); smax = -np.inf; n = 0
    while n < nm * 1e6:
        g = C.random_gram(*T, 100000, rng)
        v = ev(E, co, D, g); n += len(g)
        smax = max(smax, v.max())
        topg = np.r_[topg, g]; topv = np.r_[topv, v]
        o = np.argsort(topv)[-4 * starts:]; topg, topv = topg[o], topv[o]
    det = {'type': 'ineq', 'fun': lambda x: 1 + 2 * x[0] * x[1] * x[2] - x[0] ** 2 - x[1] ** 2 - x[2] ** 2,
           'jac': lambda x: np.array([2 * x[1] * x[2] - 2 * x[0], 2 * x[0] * x[2] - 2 * x[1], 2 * x[0] * x[1] - 2 * x[2]])}
    best, bx = -np.inf, None
    for x0 in topg[np.argsort(topv)[-starts:]]:
        r = minimize(lambda x: tuple(-y for y in f_and_grad(E, co, x)), x0, jac=True, method='SLSQP',
                     bounds=[(-1, T[0]), (-1, T[1]), (-1, T[2])], constraints=[det], options={'maxiter': 300, 'ftol': 1e-16})
        x = np.clip(r.x, -1, T)
        if 1 + 2 * x[0] * x[1] * x[2] - x[0] ** 2 - x[1] ** 2 - x[2] ** 2 >= -1e-12:
            val = f_and_grad(E, co, x)[0]
            if val > best:
                best, bx = val, x
    print('%s: sampled max %.9e on %.1fM triples; ascent from %d starts %.9e at %s [%.0f s]'
          % (''.join(combo), smax, n / 1e6, starts, best, np.round(bx, 6), time.time() - t0), flush=True)
