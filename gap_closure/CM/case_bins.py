"""case_bins.py -- a case file with one type per bin, for multi_cap/combo_mixtype_check.py.

usage: python3 case_bins.py CASE.json OUT.json
Each bin b (numbered from 1) of type s becomes the type s<b> with the bin as its range;
the count vectors (constraints) are unchanged.  These are the types c27_hull_mt.py uses.
"""
import json
import sys

E = json.load(open(sys.argv[1]))
names = ['%s%d' % (ty, b + 1) for b, (ty, _, _) in enumerate(E['bins'])]
E['tag'] = E['tag'] + '_bins'
E['types'] = names
E['trange'] = {n: [lo, hi] for n, (_, lo, hi) in zip(names, E['bins'])}
E['bins'] = [[n, lo, hi] for n, (_, lo, hi) in zip(names, E['bins'])]
json.dump(E, open(sys.argv[2], 'w'))
