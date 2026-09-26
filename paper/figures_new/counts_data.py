#!/usr/bin/env python3
"""counts_data.py -- the table of tikz_counts.tex from multi_cap/runs/count_survey.log:
for each count M the least T found and how many of its centres sit on the sphere of
radius sqrt 6."""
import re
rows = re.findall(r'M = (\d+): least T over \d+ local minima \(of \d+ starts\) = ([\d.]+), with (\d+) centres',
                  open('../../multi_cap/runs/count_survey.log').read())
with open('counts_data.dat', 'w') as f:
    f.write('M T on inside\n')
    for M, T, on in rows:
        f.write('%s %s %s %d\n' % (M, T, on, int(M) - int(on)))
print(len(rows), 'counts')
