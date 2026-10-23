"""Embed a typed three-point vector x3 of degree d into the layout of degree D >= d
(blocks of the extra degrees zero, Chebyshev indices beyond d - k zero)."""
import numpy as np


def blocks(x3, m, d):
    names = [('f', k) for k in range(1, d + 1)] + [('F', P, k) for P in range(m) for k in range(d + 1)]
    out, off = {}, 0
    for name in names:
        s = m if name[0] == 'f' else m * (d - name[2] + 1)
        n = s * (s + 1) // 2
        v = x3[off:off + n]; off += n
        M = np.zeros((s, s)); M[np.arange(s), np.arange(s)] = v[:s]
        iu = np.triu_indices(s, 1); M[iu] = v[s:]; M[(iu[1], iu[0])] = v[s:]
        out[name] = M
    assert off == len(x3)
    return out


def to_x3(B, m, D):
    names = [('f', k) for k in range(1, D + 1)] + [('F', P, k) for P in range(m) for k in range(D + 1)]
    parts = []
    for name in names:
        s = m if name[0] == 'f' else m * (D - name[2] + 1)
        M = B.get(name, np.zeros((s, s)))
        iu = np.triu_indices(s, 1)
        parts += [np.diag(M), M[iu]]
    return np.concatenate(parts)


def embed(x3, m, d, D):
    if d == D:
        return np.array(x3, float)
    B = blocks(np.array(x3, float), m, d)
    out = {}
    for name, M in B.items():
        if name[0] == 'f':
            out[name] = M
        else:
            k = name[2]; mk, MK = d - k + 1, D - k + 1
            idx = np.array([Q * MK + i for Q in range(m) for i in range(mk)])
            E = np.zeros((m * MK, m * MK)); E[np.ix_(idx, idx)] = M
            out[name] = E
    return to_x3(out, m, D)
