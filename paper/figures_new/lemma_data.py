"""lemma_data.py -- data of fig_tikz_lemma: the cap volume S(d) of a centre at distance d, and
the count bound of the lemma "Twenty-four centres within 2.444" (M-23) S(2.444) against the room
T(X) - 8 >= s(D) - (8 - A*) of the twenty-three-centre theorem.  Ball arithmetic at 200 bits."""
from flint import arb, ctx
ctx.prec = 200
R2 = arb(3)/2; R = R2.sqrt(); R4 = R2*R2
def S(d):
    h = arb(d)/2
    if h >= R: return arb(0)
    anti = h/8*(5*R2-2*h*h)*(R2-h*h).sqrt() + 3*R4/8*(h/R).asin()
    return 4*arb.pi()/3*(3*R4/8*arb.pi()/2 - anti)
tau = 9*arb.pi()**2/8 - 8
gap = 8 - (9*arb.pi()**2/8 - 23*S(2))
room = S(2) - S(arb('2.1648')) - gap
with open('lemma_S.dat', 'w') as f:
    f.write('d S\n')
    n = 600
    for i in range(n+1):
        d = 2 + (arb(6).sqrt() - 2 - arb('1e-4'))*i/n
        f.write('%.6f %.8e\n' % (float(d.mid()), float(S(d).mid())))
s = S(arb('2.444'))
with open('lemma_bars.dat', 'w') as f:
    f.write('M v\n')
    for M in range(25, 36):
        f.write('%d %.8e\n' % (M, float(((M-23)*s).mid())))
for d in ['2', '2.0161', '2.1648', '2.4', '2.444']:
    print('S(%s) =' % d, S(arb(d)).str(8))
print('room', room.str(8), ' tau', tau.str(8), ' 22S(2)+11S(2.444)', (22*S(2)+11*s).str(8))
