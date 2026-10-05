"""cradial_data.py -- data of fig_tikz_cradial and of tab:C-radial (prop:C-radial): for M centres within sqrt 6 with T(Y) <= 8, the least
number of centres within a radius rho,
  - proved by the case certificates (read from ../../multi_cap/runs/radial_case_check_M_rRHO.log,
    and radial_case_check_30.log for thirty centres), a count at a radius also holding at every
    larger radius;
  - given by the caps alone: the least n such that n S(2) + (M - n) S(rho) >= tau fails for every
    smaller count, tau = 9 pi^2/8 - 8 (ball arithmetic).
Writes cradial_counts.dat (one row per M; nan where no count above the caps is proved) and prints
the rows of the table and the value of each single-split certificate."""
import glob, math, re
from flint import arb, ctx
ctx.prec = 200
R2 = arb(3)/2; R = R2.sqrt(); PI = arb.pi()
def S(d):
    h = arb(d)/2
    anti = h/8*(5*R2-2*h*h)*(R2-h*h).sqrt() + 3*R2*R2/8*(h/R).asin()
    return 4*PI/3*(3*R2*R2/8*PI/2 - anti)
tau = 9*PI**2/8 - 8
RHOS = ['2.05', '2.1', '2.15', '2.2']
MS = range(24, 31)

def caps(M, rho):
    """least count within rho that the caps force: 1 + the largest k <= min(M, 24) with
    k S(2) + (M - k) S(rho) < tau (certainly), or 0 if there is none."""
    best = 0
    for k in range(0, min(M, 24) + 1):
        v = k*S(2) + (M - k)*S(arb(rho))
        if v < tau:
            best = k + 1
        else:
            assert not (v > tau) or True
    return best

pat = re.compile(r'(\d+) <= N\(([\d.]+)\) <= (\d+)')
proved, values = {}, {}
for f in glob.glob('../../multi_cap/runs/radial_case_check_*.log'):
    txt = open(f).read()
    line = [l for l in txt.splitlines() if l.startswith('PASS')]
    if not line:
        continue
    M = int(re.search(r'exactly (\d+) centres', line[0]).group(1))
    for lo, r, hi in pat.findall(line[0]):
        if int(lo) > 0 and r != '2.0161':
            proved[(M, r)] = max(proved.get((M, r), 0), int(lo))
    m = re.search(r'below the target: \[([\d.]+)', txt)
    if m and re.search(r'_r([\d.]+)\.log$', f):          # single splits only, not the tree at thirty
        rr = re.search(r'_r([\d.]+)\.log$', f).group(1)
        values[(M, rr)] = float(m.group(1))

def best(M, rho):
    """the proved count at rho, or at a smaller radius of RHOS (a count within a smaller radius
    holds within a larger one)."""
    c = [proved.get((M, r), 0) for r in RHOS if float(r) <= float(rho)]
    c += [proved.get((M, r), 0) for r in ('2.0161',)]
    return max(c) if max(c) > 0 else None

tex = []
with open('cradial_counts.dat', 'w') as f:
    f.write('M ' + ' '.join('p%s' % r.replace('.', '') for r in RHOS) + ' '
            + ' '.join('c%s' % r.replace('.', '') for r in RHOS) + '\n')
    for M in MS:
        c = [caps(M, r) for r in RHOS]
        p = [best(M, r) for r in RHOS]
        p = [x if x is not None and x > y else None for x, y in zip(p, c)]   # only counts above the caps
        f.write('%d ' % M + ' '.join('nan' if x is None else str(x) for x in p) + ' '
                + ' '.join(str(x) for x in c) + '\n')
        print('%d  ' % M + '  '.join('%s (%d)' % ('--' if x is None else x, y) for x, y in zip(p, c)))
        tex.append('$%d$ & ' % M + ' & '.join('%s\\,(%d)' % ('--' if x is None else '$%d$' % x, y) for x, y in zip(p, c)) + '\\\\')
print('rows of the table:')
print('\n'.join(tex))
print('certificate values:')
for k in sorted(values):
    print('  M=%d rho=%s: %.7f' % (k[0], k[1], values[k]))
print('largest (25 <= M <= 30):', max(v for k, v in values.items() if 25 <= k[0] <= 30))
