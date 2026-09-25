"""
fig_closure.py -- the further reductions of Section 2.8.

(a) A vertex of the 24-cell and its three orthogonal root pairs, in the linear
    projection of R^4 that keeps the vertex direction e_1 vertical and sends
    e_2, e_3, e_4 to three horizontal unit vectors at 120 degrees.  The six
    normalised roots at the vertex v = sqrt2 e_1 lie on a hexagon at height
    1/sqrt2; the points 2u lie on a hexagon at height sqrt2, and v is the
    midpoint of each of the three chords 2u -- 2u' of an orthogonal pair
    (Corollary 2.14).  One root is missing (hollow): its chord is gone, and v
    is still the midpoint of the other two.  The projection is linear, so
    midpoints are exact.
(b) The three-dimensional analogue of the inversion hull (Proposition 2.13):
    the Voronoi cell of the twelve contacts of the face-centred cubic packing,
    a rhombic dodecahedron (blue), inside the convex hull of the points 2u
    (the cuboctahedron, orange), touching its boundary only at its six
    four-valent vertices, as the 24-cell touches the hull of the 2u at its 24
    vertices.
(c) First order near the root system (Theorem 2.18): the exact volume
    excess of the root system pushed out uniformly, 8(1 + delta/2)^4 - 8, and
    of the deletion with one centre held at 2 + delta on the deleted axis,
    1/3 - (1/3)(1 - delta/2)^4, against sum delta_i, with the tangent
    (2/3) sum delta_i that the proof of Theorem 2.18 finds at the root system.
"""
import itertools
import math

import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
from figstyle import *

S2 = math.sqrt(2)


def faces_of(P):
    """planar faces (vertex cycles) and unit normals of the convex hull of P"""
    H = ConvexHull(P)
    groups = {}
    for simplex, eq in zip(H.simplices, H.equations):
        key = tuple(np.round(eq, 6))
        groups.setdefault(key, set()).update(simplex.tolist())
    out = []
    for key, idx in groups.items():
        n = np.array(key[:3])
        pts = P[list(idx)]
        c = pts.mean(0)
        u = pts[0] - c
        u /= np.linalg.norm(u)
        w = np.cross(n, u)
        ang = np.arctan2((pts - c) @ w, (pts - c) @ u)
        out.append((pts[np.argsort(ang)], n))
    return out


def edges_of(faces):
    E = set()
    for pts, _ in faces:
        for k in range(len(pts)):
            a, b = tuple(np.round(pts[k], 6)), tuple(np.round(pts[(k + 1) % len(pts)], 6))
            E.add(tuple(sorted([a, b])))
    return [np.array(e) for e in E]


LIGHT = np.array([0.35, -0.55, 0.76])
LIGHT /= np.linalg.norm(LIGHT)

fig = plt.figure(figsize=(7.2, 2.75))
gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.0, 1.05], left=0.0, right=0.985,
                      top=0.95, bottom=0.14, wspace=0.08)

# ------------------------------------------------------------------ (a)
ax = fig.add_subplot(gs[0], projection="3d")
blank_3d(ax)
ax.computed_zorder = False
ax.view_init(elev=27, azim=-62)
dirs = [np.array([math.cos(t), math.sin(t), 0.0]) for t in (math.pi / 2, 7 * math.pi / 6, 11 * math.pi / 6)]
up = np.array([0, 0, 1.0])
cols = [BLUE, ORANGE, AQUA]
v = S2 * up
# the hexagon of the points 2u, a slice of the hull, as a pale sheet
hexpts = []
for k in range(6):
    t = math.pi / 2 + k * math.pi / 3
    hexpts.append(S2 * np.array([math.cos(t), math.sin(t), 0]) + S2 * up)
ax.add_collection3d(Poly3DCollection([hexpts], facecolors=[(0.93, 0.91, 0.86, 0.6)],
                                     edgecolors=[INK3], linewidths=0.5, zorder=2))
missing = (1, +1)                        # the orange pair loses one root
for j, (a, c) in enumerate(zip(dirs, cols)):
    for s in (+1, -1):
        u = (s * a + up) / S2
        P2 = 2 * u
        ax.plot(*zip(np.zeros(3), P2), color=INK3, lw=0.45, zorder=1)
        if (j, s) == missing:
            ax.scatter(*P2, s=26, facecolors="white", edgecolors=c, linewidths=1.1, zorder=6, depthshade=False)
        else:
            ax.scatter(*P2, s=22, color=c, zorder=6, depthshade=False)
        ax.scatter(*u, s=7, color=c, zorder=6, depthshade=False)
    A, B = S2 * (dirs[j] + up), S2 * (-dirs[j] + up)
    ls = (0, (2.5, 2)) if j == missing[0] else "-"
    ax.plot(*zip(A, B), color=c, lw=1.5, ls=ls, zorder=5)
ax.scatter(*v, s=34, color=INK, zorder=8, depthshade=False)
ax.scatter(0, 0, 0, s=18, color=INK, zorder=8, depthshade=False)
ax.text(0.12, 0.0, -0.22, "0", fontsize=8, color=INK)
ax.text(0.16, 0.05, S2 + 0.16, "$v$", fontsize=9, color=INK, zorder=9)
ax.set_xlim(-1.25, 1.25); ax.set_ylim(-1.25, 1.25); ax.set_zlim(-0.05, 1.6)
fig.text(0.01, 0.94, "(a)", fontsize=10, fontweight="bold", ha="left", va="top", color=INK)

# ------------------------------------------------------------------ (b)
ax = fig.add_subplot(gs[1], projection="3d")
blank_3d(ax)
ax.computed_zorder = False
ax.view_init(elev=19, azim=-58)
U = []
for i, j in itertools.combinations(range(3), 2):
    for si in (1, -1):
        for sj in (1, -1):
            w = np.zeros(3); w[i] = si; w[j] = sj
            U.append(w / S2)
U = np.array(U)
cubo = 2 * U
rd = [S2 * s * np.eye(3)[i] for i in range(3) for s in (1, -1)]
rd += [np.array(p) / S2 for p in itertools.product((1, -1), repeat=3)]
rd = np.array(rd)
F_rd = faces_of(rd)
F_cu = faces_of(cubo)
view = np.array([math.cos(math.radians(19)) * math.cos(math.radians(-58)),
                 math.cos(math.radians(19)) * math.sin(math.radians(-58)),
                 math.sin(math.radians(19))])
# painter's order by hand (computed_zorder is off): hull back faces, back edges and
# back vertices; then the visible faces of the cell; then the hull's front edges,
# front vertices and the visible contact vertices of the cell
E_cu = edges_of(F_cu)
back_e = [e for e in E_cu if e.mean(0) @ view < 0]
front_e = [e for e in E_cu if e.mean(0) @ view >= 0]
for pts, n in F_cu:
    if n @ view < 0:
        lam = 0.55 + 0.45 * max(0.0, -(n @ LIGHT))
        col = np.clip(0.55 + 0.45 * np.array(to_rgb(ORANGE)) * lam, 0, 1)
        ax.add_collection3d(Poly3DCollection([pts], facecolors=[(*col, 0.35)], edgecolors="none", zorder=1))
ax.add_collection3d(Line3DCollection(back_e, colors=[(*to_rgb(ORANGE), 0.45)], linewidths=0.7, zorder=2))
bk = cubo[cubo @ view < 0]
ax.scatter(*bk.T, s=7, color=ORANGE, alpha=0.5, depthshade=False, zorder=2)
for pts, n in F_rd:
    if n @ view <= 0:
        continue
    lam = max(0.0, n @ LIGHT)
    col = np.clip(np.array(to_rgb(BLUE)) * (0.62 + 0.45 * lam), 0, 1)
    ax.add_collection3d(Poly3DCollection([pts], facecolors=[col], edgecolors=[(1, 1, 1, 0.9)],
                                         linewidths=0.5, zorder=3))
for pts, n in F_cu:
    if n @ view >= 0:
        ax.add_collection3d(Poly3DCollection([pts], facecolors=[(*to_rgb(ORANGE), 0.06)],
                                             edgecolors="none", zorder=4))
ax.add_collection3d(Line3DCollection(front_e, colors=[ORANGE], linewidths=0.95, zorder=5))
fr = cubo[cubo @ view >= 0]
ax.scatter(*fr.T, s=9, color=ORANGE, depthshade=False, zorder=6)
touch = rd[:6]
vis = touch[touch @ view > 0]
ax.scatter(*vis.T, s=15, color=INK, depthshade=False, zorder=7)
ax.set_xlim(-2.0, 2.0); ax.set_ylim(-2.0, 2.0); ax.set_zlim(-2.0, 2.0)
fig.text(0.345, 0.94, "(b)", fontsize=10, fontweight="bold", ha="left", va="top", color=INK)

# ------------------------------------------------------------------ (c)
ax = fig.add_subplot(gs[2]); clean_axes(ax)
x = np.linspace(0, 0.6, 300)
ax.plot(x, 8 * (1 + x / 48) ** 4 - 8, color=BLUE, lw=1.4, label="root system, all 24 out")
ax.plot(x, 1 / 3 - (1 / 3) * (1 - x / 2) ** 4, color=ORANGE, lw=1.4, label="deletion, one centre on the axis")
ax.plot(x, (2 / 3) * x, color=INK3, lw=0.9, ls=(0, (4, 2)), label=r"tangent $\frac{2}{3}\sum\delta_i$")
ax.set_xlim(0, 0.6); ax.set_ylim(0, 0.43)
ax.set_xticks([0, 0.2, 0.4, 0.6]); ax.set_yticks([0, 0.1, 0.2, 0.3, 0.4])
ax.tick_params(labelsize=6.6)
ax.set_xlabel(r"total excess distance $\delta_1+\delta_2+\cdots$", fontsize=7.2)
ax.set_ylabel(r"$\mathrm{vol}(V_c)-8$", fontsize=7.4)
ax.legend(loc="lower right", fontsize=6.0, handlelength=1.6, borderpad=0.2)
panel_label(ax, "(c)", x=-0.13, y=0.97)

save(fig, "fig_closure", outdir="../figures")
