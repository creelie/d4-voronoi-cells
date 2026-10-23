"""largest value of the triple polynomials of a mixed certificate, by dense sampling and local ascent
usage: python3 c28_peak.py case.json cert.npz KINDS(comma) BATCHES
(found the AAF peak of mix28_k22d.npz that the first exact check missed, peak_k22d.log)"""
import sys, os, time
import numpy as np
case, npz, kinds, nb = sys.argv[1], sys.argv[2], sys.argv[3].split(','), int(sys.argv[4])
case, npz = os.path.abspath(case), os.path.abspath(npz)
MC = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'multi_cap')
sys.path.insert(0, MC)
os.chdir(MC)
sys.argv = ['x', case, npz, '10']
import combo_case_check as K
C = K.C
t0 = time.time()
Z = np.load(npz)
Bf = C.blocks_from_x3(np.array(Z['x3'], float), 10)
Bk = {name: C.psd_exact(M) for name, M in Bf.items()}
print('blocks [%.0f s]' % (time.time() - t0), flush=True)
for kind in kinds:
    combo = tuple(kind)
    P = C.triple3_poly(Bk, 10, combo)
    T = tuple(float(C.tmax(a, b)) for a, b in ((combo[0], combo[1]), (combo[0], combo[2]), (combo[1], combo[2])))
    rng = np.random.default_rng(1234)
    topg = np.zeros((0, 3)); topv = np.zeros(0); smax = -1
    for b in range(nb):
        g = C.random_gram(*T, 200000, rng)
        v = C.peval(P, g[:, 0], g[:, 1], g[:, 2])
        smax = max(smax, v.max())
        topg = np.r_[topg, g]; topv = np.r_[topv, v]
        o = np.argsort(topv)[-3000:]; topg, topv = topg[o], topv[o]
    starts = topg[np.argsort(topv)[-300:]]
    best = C.refine_triple(P, starts, T)
    print('%s: tops %s; sampled max %.9e on %d batches; refined from 300 starts %.9e [%.0f s]'
          % (kind, T, smax, nb, best, time.time() - t0), flush=True)
