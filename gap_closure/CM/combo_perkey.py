# Per type-count bound of a combo_direct.py certificate on a case (floating point, evidence only).
# Usage: python3 combo_perkey.py cert.npz d3 rounds case.json   (d3, rounds and case as for combo_direct.py)
import sys, os, itertools
import numpy as np
from scipy.optimize import linprog
cert = sys.argv.pop(1)
sys.path.insert(0, os.getcwd())
from combo_gen2 import T3, TYPES, BINS, COUNTS, CASE, D2, R2, tcounts, Npair, Ntriple
d = np.load(cert)
x3 = d['x3']; m = d['m']; t = float(d['t'])
c2keys = list(itertools.combinations_with_replacement(TYPES, 2))
c3keys = [cb for cb in itertools.combinations_with_replacement(TYPES, 3)
          if not all(Ntriple(tcounts(nb), list(cb)) == 0 for nb in COUNTS)]
c2 = dict(zip(c2keys, d['c2'])); c3 = dict(zip(c3keys, d['c3']))
L = T3.Layout(len(TYPES), int(sys.argv[1]))
ti = {s: i for i, s in enumerate(TYPES)}
Bp = T3.Builder(L)
for s in TYPES:
    T3.add_point(Bp, np.array([ti[s]]), ti[s])
Pm = Bp.matrix(len(TYPES)).toarray()
p = {s: float(Pm[ti[s]] @ x3) for s in TYPES}
print('point terms', p, 't/2', t / 2)
nb = len(BINS)
keys = sorted({tuple(tcounts(v)[ty] for ty in TYPES) for v in COUNTS})
rows = []
for key in keys:
    tc = dict(zip(TYPES, key))
    # LP max of m.n over bin counts with the case constraints and these type counts
    Aub, bub = [], []
    for bs, lo, hi in CASE['constraints']:
        r = np.zeros(nb); r[bs] = 1
        Aub.append(r); bub.append(hi); Aub.append(-r); bub.append(-lo)
    Aeq = np.array([[1.0 if BINS[b][0] == ty else 0.0 for b in range(nb)] for ty in TYPES])
    res = linprog(-m, A_ub=np.array(Aub), b_ub=bub, A_eq=Aeq, b_eq=[tc[s] for s in TYPES], bounds=[(0, None)] * nb)
    lp = -res.fun
    pts = sum(tc[s] * p[s] for s in TYPES)
    pr = sum(Npair(tc, *k) * c2[k] for k in c2keys)
    tr = sum(Ntriple(tc, list(k)) * c3[k] for k in c3keys)
    tot = lp + pts + t / 2 + pr + tr
    rows.append((tot, key, lp, pts, pr, tr, res.x))
rows.sort(reverse=True)
for tot, key, lp, pts, pr, tr, x in rows:
    print('%s  %.5f   bins %s  (m.n %.3f pts %.3f pairs %.3f triples %.3f)' % (key, tot, np.round(x, 2), lp, pts, pr, tr))
