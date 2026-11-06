"""remap_bins.py case.json in.npz out.npz: copy a mixed certificate and recompute its per-bin radial
constants m_b on the bins of another case (max on a fine grid + 2e-7, as hull_c.py does), so that a
certificate mixed on a narrowed case can be checked on a case with the usual bins."""
import sys, numpy as np
import os
CM = os.path.dirname(os.path.abspath(__file__)) + '/'
case, src, dst = sys.argv[1], sys.argv[2], sys.argv[3]
sys.argv = ['x', '10', '1', case if case.startswith('/') else CM + case]; sys.path.insert(0, CM)
import combo_gen2 as G
Z = dict(np.load(src))
A, z = np.array(Z['A'], float), np.array(Z['z'], float)
U1 = G.ubasis(np.ones(1), G.D2)[0]


def f(d):
    Pd = G.pbasis(d, G.R2)
    return G.S(d) + 0.5 * np.einsum('k,na,kab,nb->n', U1, Pd, A, Pd) - Pd @ z


old = np.array(Z['m'], float)
Z['m'] = np.array([f(np.linspace(lo, hi, 20001)).max() + 2e-7 for _, lo, hi in G.BINS])
print('bins', [(b[0], b[1], b[2]) for b in G.BINS]); print('old m', np.round(old, 7)); print('new m', np.round(Z['m'], 7))
np.savez(dst, **Z)
