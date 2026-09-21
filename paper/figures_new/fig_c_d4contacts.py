"""
fig_c_d4contacts.py -- the twenty-four contact directions and their four
inner products.

(a) The 24 normalised roots of D_4 in an oblique projection of R^4 to R^3
    (the first three coordinates, with the fourth folded in along a fixed
    direction so that no two roots coincide).  One root is singled out; the
    other 23 are coloured by their inner product with it: eight at 1/2, six
    at 0, eight at -1/2 and the antipode at -1.  The segments join the chosen
    root to its eight neighbours at 60 degrees, the contact graph of the
    kissing configuration.
(b) The Gram matrix of the 24 directions, rows and columns ordered by the
    orbit structure about the chosen root, taking only the values 1, 1/2, 0,
    -1/2, -1.
(c) The 276 pairs by inner product: 96 at 1/2, 72 at 0, 96 at -1/2, 12 at
    -1, which is what Theorem (twenty-four points) forces of every contact
    configuration of size 24, up to an orthogonal map.
"""
import itertools

import numpy as np
import matplotlib.pyplot as plt
from figstyle import *

roots = []
for i, j in itertools.combinations(range(4), 2):
    for si in (1, -1):
        for sj in (1, -1):
            v = np.zeros(4)
            v[i], v[j] = si, sj
            roots.append(v / np.sqrt(2))
R = np.array(roots)
assert R.shape == (24, 4)
G = R @ R.T
vals, counts = np.unique(np.round(G[np.triu_indices(24, 1)], 6), return_counts=True)
assert dict(zip(vals, counts)) == {-1.0: 12, -0.5: 96, 0.0: 72, 0.5: 96}, dict(zip(vals, counts))

# orthogonal projection onto the hyperplane n^perp for a generic unit vector n,
# so that the image of S^3 lies inside the unit ball of R^3 and no two roots coincide
n = np.array([1.0, 2.0, 3.0, 5.0]); n /= np.linalg.norm(n)
Q, _ = np.linalg.qr(np.column_stack([n, np.eye(4)[:, :3]]))
B3 = Q[:, 1:4]                                    # 4 x 3 orthonormal, orthogonal to n
X = R @ B3
d = np.linalg.norm(X[:, None, :] - X[None, :, :], axis=2) + np.eye(24)
assert d.min() > 0.15, d.min()

k0 = 0                                            # the chosen root
ip = np.round(G[k0], 6)
classes = [(0.5, BLUE, r"$\langle r_0,r\rangle=\frac{1}{2}$ (8)"),
           (0.0, AQUA, r"$\langle r_0,r\rangle=0$ (6)"),
           (-0.5, ORANGE, r"$\langle r_0,r\rangle=-\frac{1}{2}$ (8)"),
           (-1.0, VIOLET, r"$\langle r_0,r\rangle=-1$ (1)")]

fig = plt.figure(figsize=(7.2, 3.3))
gs = fig.add_gridspec(1, 3, width_ratios=[1.45, 1.0, 0.95], left=0.0, right=0.99, top=0.9, bottom=0.17, wspace=0.55)
ax = fig.add_subplot(gs[0], projection="3d")
blank_3d(ax)
ax.computed_zorder = False
ax.view_init(elev=16, azim=-52)
# the unit sphere of R^3, as a faint wire frame, for the scale of the projection
uu, vv = np.meshgrid(np.linspace(0, 2 * np.pi, 37), np.linspace(0, np.pi, 19))
ax.plot_wireframe(np.cos(uu) * np.sin(vv), np.sin(uu) * np.sin(vv), np.cos(vv), color=GRID, lw=0.25, zorder=0)
for a, b in itertools.combinations(range(24), 2):
    if abs(G[a, b] - 0.5) < 1e-9:
        ax.plot(*zip(X[a], X[b]), color=INK3, lw=0.45, alpha=0.55, zorder=1)
for a in range(24):
    if a != k0 and abs(G[k0, a] - 0.5) < 1e-9:
        ax.plot(*zip(X[k0], X[a]), color=BLUE, lw=1.3, zorder=3)
for val, col, lab in classes:
    sel = [a for a in range(24) if a != k0 and abs(ip[a] - val) < 1e-9]
    ax.scatter(X[sel, 0], X[sel, 1], X[sel, 2], s=30, color=col, depthshade=False, zorder=5, label=lab, edgecolors=SURFACE, linewidths=0.4)
ax.scatter([X[k0, 0]], [X[k0, 1]], [X[k0, 2]], s=70, color=INK, depthshade=False, zorder=6, label=r"$r_0$", edgecolors=SURFACE, linewidths=0.5)
lim = 1.02
ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_zlim(-lim, lim)
ax.set_box_aspect((1, 1, 1))
ax.text2D(0.02, 0.97, "(a)", transform=ax.transAxes, fontsize=10, fontweight="bold")
ax.text2D(0.5, 0.0, "orthogonal projection of the 24 roots to a 3-space;\nthin edges join every pair at inner product 1/2", transform=ax.transAxes,
          ha="center", va="top", fontsize=6.3, color=INK2)
ax.legend(loc="upper left", fontsize=5.6, frameon=True, framealpha=1, edgecolor="none",
          handletextpad=0.3, borderpad=0.4, labelspacing=0.25, bbox_to_anchor=(0.62, 1.0))

# (b) Gram matrix ordered by orbit about r_0
order = [k0] + [a for val, _, _ in classes for a in range(24) if a != k0 and abs(ip[a] - val) < 1e-9]
Gp = G[np.ix_(order, order)]
ax2 = fig.add_subplot(gs[1])
from matplotlib.colors import ListedColormap, BoundaryNorm
cmap = ListedColormap([VIOLET, ORANGE, "#e9e8e4", BLUE, INK])
norm = BoundaryNorm([-1.25, -0.75, -0.25, 0.25, 0.75, 1.25], cmap.N)
im = ax2.imshow(Gp, cmap=cmap, norm=norm, interpolation="nearest")
ax2.set_xticks([0, 4.5, 11.5, 18.5, 23]); ax2.set_xticklabels([r"$r_0$", r"$\frac{1}{2}$", "0", r"$-\frac{1}{2}$", r"$-1$"], fontsize=6.5)
ax2.set_yticks([0, 4.5, 11.5, 18.5, 23]); ax2.set_yticklabels([r"$r_0$", r"$\frac{1}{2}$", "0", r"$-\frac{1}{2}$", r"$-1$"], fontsize=6.5)
ax2.tick_params(length=0)
for sp in ax2.spines.values():
    sp.set_visible(False)
for b in (0.5, 8.5, 14.5, 22.5):
    ax2.axhline(b, color=SURFACE, lw=1.2); ax2.axvline(b, color=SURFACE, lw=1.2)
cb = fig.colorbar(im, ax=ax2, fraction=0.046, pad=0.03, ticks=[-1, -0.5, 0, 0.5, 1])
cb.ax.set_yticklabels([r"$-1$", r"$-\frac{1}{2}$", "0", r"$\frac{1}{2}$", "1"], fontsize=6.5)
cb.outline.set_visible(False)
ax2.set_title("Gram matrix, ordered about $r_0$", fontsize=7.5, pad=4, loc="right")
ax2.text(-0.48, 1.06, "(b)", transform=ax2.transAxes, fontsize=10, fontweight="bold")

# (c) pair counts
ax3 = fig.add_subplot(gs[2])
clean_axes(ax3)
labels = [r"$-1$", r"$-\frac{1}{2}$", "0", r"$\frac{1}{2}$"]
cnt = [12, 96, 72, 96]
cols = [VIOLET, ORANGE, AQUA, BLUE]
bars = ax3.bar(range(4), cnt, color=cols, width=0.62, zorder=3)
for i, c in enumerate(cnt):
    ax3.text(i, c + 3, str(c), ha="center", fontsize=7, color=INK2,
             bbox=dict(facecolor=SURFACE, edgecolor="none", pad=0.4))
ax3.set_xticks(range(4)); ax3.set_xticklabels(labels)
ax3.set_ylim(0, 150)
ax3.set_xlabel("inner product of a pair")
ax3.set_ylabel("number of pairs (of 276)")
ax3.text(1.5, 147, "the only values a 24-point contact\nconfiguration can take", ha="center", va="top",
         fontsize=6.4, color=INK2, bbox=dict(facecolor=SURFACE, edgecolor="none", pad=0.8))
panel_label(ax3, "(c)", x=-0.22)
save(fig, "fig_c_d4contacts", outdir="../figures")
