"""
fig_c_certificate.py -- the certificate of thm:certificate (fig:certpair).

(a) The slack of the condition (eq:C), Q / 1000 = omega(u) + omega(v) + omega(t) - P / 1000,
    on the slice t = 1/2 of the admissible domain, u <= v <= 1/2, from
    multi_cap/continuation_out/certificate_d8.npz and the closed form of omega.  The slice
    is bounded by the diagonal u = v and by the two arcs of coplanar triples, which meet
    the diagonal at (-sqrt3/2, -sqrt3/2), where the slack is least over the whole domain
    (multi_cap/runs/certificate_tight.log).  The squares are the triples of the root system
    with one root removed that have a pair at 60 degrees.
(b) The values of the relaxations against 8 - A_* (prop:two-point-barrier and
    multi_cap/runs/three_point_sdp_d6.log, three_point_sdp_d8.log), the certificate, and
    the pair sum of the root system with one root removed.
"""
import itertools
import math
import os
import sys

import numpy as np
from matplotlib.colors import LinearSegmentedColormap

from figstyle import *

HERE = os.path.dirname(os.path.abspath(__file__))
MC = os.path.join(HERE, "..", "..", "multi_cap")
sys.path.insert(0, MC)
from fractions import Fraction as Fr  # noqa: E402
from certificate_check import SCALE, build_P  # noqa: E402

Z = np.load(os.path.join(MC, "continuation_out", "certificate_d8.npz"))
d = len(Z["f"]) - 1
f = [Fr(float(x)) for x in Z["f"]]
f = [f[0]] + [max(x, Fr(0)) for x in f[1:]]
F = [[[Fr(float(x)) for x in row] for row in Z["F%d" % k]] for k in range(d + 1)]
P, _ = build_P(d, f, F)
E = np.array(list(P.keys()))
C = np.array([float(c) for c in P.values()])
TS = 1 / math.sqrt(2)


def omega(u):
    u = np.clip(np.asarray(u, float), -1 + 1e-12, 1)
    tau = np.sqrt((1 - u) / (1 + u))
    w = 4 * math.pi * (9 / 32 * np.arctan((TS - tau) / (1 + TS * tau)) - (4 * tau ** 3 - 24 * tau + 11 * math.sqrt(2)) / 96)
    return np.where(u > 1 / 3, w, 0.0)


def slack(X):
    X = np.atleast_2d(X)
    out = np.empty(len(X))
    for i in range(0, len(X), 20000):
        Y = X[i:i + 20000]
        Pv = (C[None, :] * np.prod(Y[:, None, :] ** E[None, :, :], axis=2)).sum(1)
        out[i:i + 20000] = omega(Y).sum(1) - Pv / SCALE
    return out


n = 361
ax_u = np.linspace(-1, 0.5, n)
U, V = np.meshgrid(ax_u, ax_u, indexing="xy")
inside = (U <= V + 1e-12) & (U * U - U * V + V * V <= 0.75 + 1e-12)
S = np.full(U.shape, np.nan)
pts = np.stack([U[inside], V[inside], np.full(inside.sum(), 0.5)], 1)
S[inside] = slack(pts)
umin = -math.sqrt(3) / 2
smin = float(slack(np.array([umin, umin, 0.5]))[0])
assert np.nanmin(S) >= smin - 1e-9

# the triples of the root system less one root with a pair at 60 degrees, sorted u <= v <= 1/2
R = np.array([np.array(a) / math.sqrt(2) for a in itertools.product([-1, 0, 1], repeat=4) if sum(x * x for x in a) == 2])
W = R[1:]
G = W @ W.T
sq = set()
for i, j, k in itertools.combinations(range(23), 3):
    x = sorted(np.round([G[i, j], G[i, k], G[j, k]], 6))
    if abs(x[2] - 0.5) < 1e-9:
        sq.add((x[0], x[1]))
sq = sorted(sq)

fig = plt.figure(figsize=(7.2, 3.3))
gs = fig.add_gridspec(1, 2, width_ratios=[1.12, 1.0], left=0.085, right=0.985, top=0.9, bottom=0.14, wspace=0.42)

ax = fig.add_subplot(gs[0])
cmap = LinearSegmentedColormap.from_list("slack", ["#eef3fb", BLUE, INK])
im = ax.imshow(1e5 * S, origin="lower", extent=(-1, 0.5, -1, 0.5), cmap=cmap, vmin=0, vmax=40,
               interpolation="nearest", aspect="equal", zorder=1)
# the coplanar triples of the slice, on the ellipse u^2 - u v + v^2 = 3/4
a_ = np.linspace(math.radians(60), math.radians(180), 1201)
arc1 = np.stack([np.cos(a_), np.cos(a_ - math.pi / 3)], 1)           # directions at 0, a and 60 degrees
arc2 = np.stack([np.cos(a_), np.cos(a_ + math.pi / 3)], 1)
for arc in (arc1, arc2):
    arc = arc[(arc[:, 0] <= arc[:, 1] + 1e-12) & (arc[:, 1] <= 0.5 + 1e-12)]
    ax.plot(arc[:, 0], arc[:, 1], color=ORANGE, lw=1.2, zorder=3)
ax.plot([umin, 0.5], [umin, 0.5], color=INK3, lw=0.8, ls=(0, (3, 2)), zorder=3)
ax.plot([-0.5, 0.5], [0.5, 0.5], color=INK3, lw=0.8, zorder=2)
ax.scatter([p[0] for p in sq], [p[1] for p in sq], marker="s", s=22, color=YELLOW, edgecolor=INK, linewidth=0.5, zorder=4)
ax.scatter([umin], [umin], s=26, color=ORANGE, edgecolor=INK, linewidth=0.6, zorder=5)
ax.annotate(r"least, $%.2f\times10^{-5}$" % (1e5 * smin), xy=(umin, umin), xytext=(-0.62, -0.93),
            fontsize=7, color=INK, ha="left", va="center",
            arrowprops=dict(arrowstyle="-", color=INK2, lw=0.6),
            bbox=dict(facecolor=SURFACE, edgecolor="none", pad=0.6), zorder=6)
ax.text(-0.99, 0.43, "coplanar\ntriples", fontsize=6.8, color=ORANGE, ha="left", va="center",
        bbox=dict(facecolor=SURFACE, edgecolor="none", pad=0.6), zorder=6)
ax.text(0.15, -0.12, "$u=v$", fontsize=7, color=INK2, rotation=45, ha="center", va="center", zorder=6)
ax.annotate("root system\nless one root", xy=(-0.5, -0.5), xytext=(0.02, -0.72), fontsize=6.8, color=INK2,
            ha="center", va="center", arrowprops=dict(arrowstyle="-", color=INK2, lw=0.6),
            bbox=dict(facecolor=SURFACE, edgecolor="none", pad=0.6), zorder=6)
ax.set_xlim(-1.02, 0.52); ax.set_ylim(-1.02, 0.52)
ticks = [-1, -0.5, 0, 0.5]
tl = [r"$-1$", r"$-\frac{1}{2}$", "0", r"$\frac{1}{2}$"]
ax.set_xticks(ticks); ax.set_xticklabels(tl); ax.set_yticks(ticks); ax.set_yticklabels(tl)
ax.set_xlabel(r"$u=\langle w_i,w_j\rangle$")
ax.set_ylabel(r"$v=\langle w_i,w_k\rangle$")
ax.set_title(r"slack of (C) where $\langle w_j,w_k\rangle=\frac{1}{2}$", fontsize=7.8, color=INK2, loc="right")
cb = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.03)
cb.set_label(r"$10^5\times$ slack", fontsize=7.2)
cb.ax.tick_params(labelsize=7)
panel_label(ax, "(a)", x=-0.2, y=1.03)

ax2 = fig.add_subplot(gs[1]); clean_axes(ax2)
TARGET = 0.0928555703
DELETION = float(omega((W @ W.T)[np.triu_indices(23, 1)]).sum())
bars = [("pairs only", 0.0738), ("triples\n$d=6$", 0.09011), ("triples\n$d=8$", 0.09523)]
x = np.arange(3)
ax2.bar(x, [b for _, b in bars], width=0.56, color=BLUE, alpha=0.85, zorder=3)
for i, (_, b) in enumerate(bars):
    ax2.text(i, 0.0615, "%.4f" % b, ha="center", va="bottom", fontsize=7, color=SURFACE, zorder=4)
ax2.axhline(TARGET, color=ORANGE, lw=1.0, ls=(0, (4, 2)), zorder=2)
ax2.text(-0.42, TARGET + 0.0012, r"target $8-A_*=%.6f$" % TARGET, fontsize=7, color=ORANGE, ha="left", va="bottom",
         bbox=dict(facecolor=SURFACE, edgecolor="none", pad=0.4), zorder=4)
ax2.axhline(DELETION, color=YELLOW, lw=1.0, ls=(0, (1, 1.5)), zorder=2)
ax2.text(-0.42, DELETION - 0.0012, "root system less one root, %.6f" % DELETION, fontsize=7, color=INK2, ha="left", va="top",
         bbox=dict(facecolor=SURFACE, edgecolor="none", pad=0.4), zorder=4)
ax2.scatter([2], [0.0929], s=28, color=AQUA, edgecolor=INK, linewidth=0.6, zorder=5)
ax2.text(2.33, 0.0905, "certificate\n0.0929", fontsize=7, color=INK2, ha="left", va="top",
         bbox=dict(facecolor=SURFACE, edgecolor="none", pad=0.4), zorder=5)
ax2.set_xticks(x); ax2.set_xticklabels([n_ for n_, _ in bars], fontsize=7.2)
ax2.set_xlim(-0.5, 2.95)
ax2.set_ylim(0.06, 0.135)
ax2.set_ylabel(r"lower bound for $\Sigma_{i<j}\,\omega(\gamma_{ij})$")
panel_label(ax2, "(b)", x=-0.24, y=1.03)
save(fig, "fig_c_certificate", outdir="../figures")
