"""margin_data.py -- data of fig_tikz_margin: the weight w of the margin programme at s = 0.008,
and max w against the floor 36 s set by the root system, for s from 1e-5 to 0.05."""
import numpy as np
def w(u, s): return (u+1)*(u+0.5)**2*u**2*(0.5+s-u)
s0 = 0.008
u = np.linspace(-1, 0.5+s0, 4001)
with open('margin_w.dat', 'w') as f:
    f.write('u w\n')
    for x in u:
        v = w(x, s0)
        if v > 1e-5: f.write('%.6f %.8e\n' % (x, v))
        else: f.write('%.6f nan\n' % x)
with open('margin_floor.dat', 'w') as f:
    f.write('s maxw floor\n')
    for s in np.logspace(-5, np.log10(0.05), 121):
        uu = np.linspace(-1, 0.5+s, 200001)
        f.write('%.6e %.6e %.6e\n' % (s, w(uu, s).max(), 36*s))
ss = np.logspace(-6, -2, 200001)
m = np.array([w(np.linspace(-1, 0.5+s, 20001), s).max() for s in ss[::2000]])
k = np.argmin(np.abs(m - 36*ss[::2000]))
print('crossing near s =', ss[::2000][k], 'maxw', m[k])
print('max w at 0.008:', w(u, s0).max(), 'at u=', u[np.argmax(w(u, s0))])
