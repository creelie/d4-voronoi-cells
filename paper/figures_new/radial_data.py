"""radial_data.py -- data of fig_tikz_radial: for 24 centres with T <= 8, the number of centres
that can lie beyond a radius rho by the cap budget alone, k(rho) = floor((24 S(2) - tau) /
(S(2) - S(rho))) with tau = 9 pi^2/8 - 8 (ball arithmetic), against the counts the certificates of
the proposition "Where the twenty-four centres lie" allow."""
from flint import arb, ctx
ctx.prec = 200
R2 = arb(3)/2; R = R2.sqrt(); PI = arb.pi()
def S(d):
    h = arb(d)/2
    anti = h/8*(5*R2-2*h*h)*(R2-h*h).sqrt() + 3*R2*R2/8*(h/R).asin()
    return 4*PI/3*(3*R2*R2/8*PI/2 - anti)
tau = 9*PI**2/8 - 8
num = 24*S(2) - tau
rows = []
n = 900
prev = None
with open('radial_budget.dat', 'w') as f:
    f.write('rho k\n')
    for i in range(1, n+1):
        rho = 2 + arb('0.3')*i/n
        q = num/(S(2) - S(rho))
        k = int(float(q.mid())) if q.upper() - q.lower() < 0.5 else None
        k = min(k, 24)
        f.write('%.6f %d\n' % (float(rho.mid()), k))
for r in ['2.0161', '2.05', '2.1', '2.15']:
    print(r, int(float((num/(S(2)-S(arb(r)))).mid())))
