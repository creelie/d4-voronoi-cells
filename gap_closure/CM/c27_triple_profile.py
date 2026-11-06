"""c27_triple_profile.py -- how flat the triple polynomials of a mixed-typing mix (c27_hull_mt.py) are:
for each combination of bins (KINDS=1-1-2,...), samples of admissible triples, the largest value
after local ascent from the 20 best samples, quantiles, and the share of samples within 1e-5, 1e-4
and 1e-3 of the maximum.  Floating point; evidence only.
usage: MIX=mix.npz KINDS=1-1-1,1-1-2 [N=20000] python3 c27_triple_profile.py"""
import sys, os, numpy as np
CM = '/home/claude/d4-voronoi-cells/gap_closure/CM/'
sys.argv = ['x', '8', '1', CM + 'case27_k15f3g1.json']; sys.path.insert(0, CM)
import combo_gen2 as G, typed3pt as T3
Z = np.load(os.environ['MIX'], allow_pickle=True)
x3g, gtyp = Z['x3g'], Z['gtyp']; L = T3.Layout(3, 10)
nb = len(G.BINS); BT = [b[0] for b in G.BINS]
for b in range(nb): G.TRANGE['%s%d' % (BT[b], b + 1)] = (G.BINS[b][1], G.BINS[b][2])
rng = np.random.default_rng(1)
def tv(bins, X):
    out = np.zeros(len(X))
    for g in range(len(x3g)):
        B = T3.Builder(L); T3.add_triple(B, np.arange(len(X)), tuple(int(gtyp[g][b]) for b in bins), X[:, 0], X[:, 1], X[:, 2])
        out += B.matrix(len(X)) @ x3g[g]
    return out
def samp(T, n):
    V = rng.normal(size=(n, 3, 4)); V /= np.linalg.norm(V, axis=2, keepdims=True)
    g = np.stack([(V[:, 0] * V[:, 1]).sum(1), (V[:, 0] * V[:, 2]).sum(1), (V[:, 1] * V[:, 2]).sum(1)], 1)
    return g[(g[:, 0] <= T[0]) & (g[:, 1] <= T[1]) & (g[:, 2] <= T[2])]
for kind in os.environ['KINDS'].split(','):
    bins = [int(c) - 1 for c in kind.split('-')]
    names = ['%s%d' % (BT[b], b + 1) for b in bins]
    T = (G.tbound(names[0], names[1]), G.tbound(names[0], names[2]), G.tbound(names[1], names[2]))
    Xs, Vs = [], []
    while sum(len(x) for x in Xs) < int(os.environ.get('N', '200000')):
        X = samp(T, 40000)[:5000]; Xs.append(X); Vs.append(tv(bins, X))
    X = np.concatenate(Xs); V = np.concatenate(Vs); i = V.argmax(); mx = V[i]
    # local climb from the best 20
    best = []
    for j in np.argsort(V)[-20:]:
        x = X[j].copy(); v = V[j]; st = 0.02
        for it in range(300):
            y = np.clip(x + rng.normal(size=(64, 3)) * st, -1, T)
            y = y[G.gram_ok(y)]
            if len(y) == 0: st /= 2; continue
            w = tv(bins, y); k = w.argmax()
            if w[k] > v: x, v = y[k], w[k]
            else: st *= 0.7
            if st < 1e-6: break
        best.append((v, tuple(np.round(x, 4))))
    best.sort(reverse=True)
    q = np.quantile(V, [0.5, 0.9, 0.99])
    print(kind, names, 'T', np.round(T, 4), 'max %.7e' % best[0][0], 'at', best[0][1],
          'quantiles 50/90/99 %.3e %.3e %.3e' % tuple(q),
          'frac within 1e-5/1e-4/1e-3 of max: %.4f %.4f %.4f' % tuple((V >= best[0][0] - e).mean() for e in (1e-5, 1e-4, 1e-3)), flush=True)
    for v, x in best[1:6]: print('   other local max %.7e at %s' % (v, x))
