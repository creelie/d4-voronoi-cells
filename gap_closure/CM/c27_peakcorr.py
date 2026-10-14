"""c27_peakcorr.py -- the bound of a mix over the count vectors of a case file,
with each triple threshold replaced by max(mix threshold, dense-search maximum from a
c28_peak2.py log); margins as the mixer's.  Floating point; proves nothing.
usage: python3 c27_peakcorr.py mix.npz peak2.log case.json [MG_AA MG2 MG3 MG_FFF]
(vectors from the case file constraints, M from the case)"""
import sys, re, itertools, numpy as np
from math import comb
import os
CM = os.path.dirname(os.path.abspath(__file__)) + '/'
npz, plog, case = sys.argv[1], sys.argv[2], sys.argv[3]
import json; CS = json.load(open(case))
mg = [float(x) for x in sys.argv[4:8]] if len(sys.argv) > 4 else [5e-7, 5e-7, 5e-8, 1e-5]
sys.argv = ['x', '10', '1', CM + 'case28_k23.json']; sys.path.insert(0, CM)
import combo_gen2 as G, typed3pt as T3
TYPES = G.TYPES; BT = [b[0] for b in G.BINS]
L = T3.Layout(3, 10); B = T3.Builder(L)
for i, s in enumerate(TYPES): T3.add_point(B, np.array([i]), i)
Z = np.load(npz); pt = B.matrix(3) @ Z['x3']
pk = list(itertools.combinations_with_replacement(TYPES, 2)); tk = list(itertools.combinations_with_replacement(TYPES, 3))
c2 = {st: Z['c2'][j] + (mg[0] if st == ('A', 'A') else mg[1]) for j, st in enumerate(pk)}
dense = {}
for line in open(plog):
    m = re.match(r'(\w+): sampled max .* ascent from \d+ starts ([\d.e-]+)', line)
    if m: dense[m.group(1)] = float(m.group(2))
c3 = {}
for j, cb in enumerate(tk):
    k = ''.join(cb); c3[cb] = max(Z['c3'][j], dense.get(k, -1)) + mg[2] + (mg[3] if k == 'FFF' else 0)
def vectors():
    out = []
    def rec(prefix, left):
        if len(prefix) == 6:
            out.append(tuple(prefix + [left])); return
        for n in range(left + 1):
            rec(prefix + [n], left - n)
    rec([], CS['M'])
    ok = lambda v: all(lo <= sum(v[i] for i in idx) <= hi for idx, lo, hi in CS['constraints'])
    return [v for v in out if ok(v) and v[0] <= 24]
best = (-1, None)
for v in vectors():
    c = {s: sum(n for n, b in zip(v, BT) if b == s) for s in TYPES}
    np_ = lambda st: comb(c[st[0]], 2) if st[0] == st[1] else c[st[0]] * c[st[1]]
    def nt(cb):
        o = 1
        for s in set(cb): o *= comb(c[s], cb.count(s))
        return o
    val = float(Z['t']) / 2 + sum(v[b] * (Z['m'][b] + 1e-7 + pt[TYPES.index(BT[b])]) for b in range(7)) \
        + sum(np_(st) * c2[st] for st in pk) + sum(nt(cb) * c3[cb] for cb in tk)
    best = max(best, (val, v))
print('stored bound %.6f; with the dense triple maxima %.6f at %s' % (float(Z['bound']), best[0], best[1]))
