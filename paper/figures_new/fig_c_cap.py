"""
fig_c_cap.py -- the cap of the cross-polytope.

(a) The section of the reference cell V, of the body Q (the twenty-three
    contact half-spaces that stay at roots) and of the cross-polytope B by
    the plane spanned by the base direction u_0 = (1,1,0,0)/sqrt2 and the
    transverse direction v = (0,0,1,0).  In the frame of the four pairwise
    orthogonal roots (1,1,0,0)/sqrt2, (1,-1,0,0)/sqrt2, (0,0,1,1)/sqrt2,
    (0,0,1,-1)/sqrt2 the point a u_0 + b v has coordinates (a, 0, b/sqrt2,
    b/sqrt2), so the section of B = {|x_0|+|x_1|+|x_2|+|x_3| <= 2} is the
    rhombus |a| + sqrt2 |b| <= 2 with vertices +-2u_0 and +-sqrt2 v.  Every
    polygon is computed as a half-plane intersection; nothing is sketched.
    The shaded region is the section of the cap cut from Q by a half-space at
    unit distance whose normal is tilted 30 degrees from u_0 towards v.
(b) The volume of the cap cut from B by a half-space at unit distance,
    vol(B and {<x,u> >= 1}), along three great circles leaving the vertex
    direction e_0: towards the edge centre (e_0+e_1)/sqrt2, towards the facet
    centre (e_0+e_1+e_2)/sqrt3, and towards (1,1,1,1)/2.  The closed form
    (1/3) g[a_0,..,a_3] of the cap theorem is a divided difference and is
    ill conditioned where two of the a_k coincide, so the curves are the direct
    polytope volumes, and the closed form is checked against them wherever
    its nodes are separated.
"""
import itertools

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from scipy.optimize import linprog
from scipy.spatial import ConvexHull, HalfspaceIntersection
from figstyle import *

roots = []
for i, j in itertools.combinations(range(4), 2):
    for si in (1, -1):
        for sj in (1, -1):
            r = np.zeros(4); r[i], r[j] = si, sj
            roots.append(r / np.sqrt(2))
R = np.array(roots)
u0 = np.array([1, 1, 0, 0]) / np.sqrt(2)
v = np.array([0, 0, 1.0, 0])
k0 = int(np.argmin(np.linalg.norm(R - u0, axis=1)))
assert np.allclose(R[k0], u0)


def section(normals, extra=()):
    """the polygon {(a,b): <a u0 + b v, n> <= 1} and extra half-planes (c, d, e): c a + d b <= e"""
    hs = [[n @ u0, n @ v, -1.0] for n in normals] + [[c, d, -e] for c, d, e in extra]
    A = np.array([h[:2] for h in hs]); b = -np.array([h[2] for h in hs])
    res = linprog([0, 0, -1], A_ub=np.hstack([A, np.linalg.norm(A, axis=1)[:, None]]), b_ub=b, bounds=[(None, None)] * 2 + [(0, None)])
    H = HalfspaceIntersection(np.array(hs), res.x[:2])
    P = H.intersections
    c = P.mean(axis=0)
    return P[np.argsort(np.arctan2(P[:, 1] - c[1], P[:, 0] - c[0]))]


V = section(R)
Q = section([R[k] for k in range(24) if k != k0])
frame = np.array([[1, 1, 0, 0], [1, -1, 0, 0], [0, 0, 1, 1], [0, 0, 1, -1]]) / np.sqrt(2)
signs = np.array(list(itertools.product((1, -1), repeat=4)), float)
B = section([s @ frame / 2 for s in signs])            # |x|_1 <= 2 in the frame: <x, s.frame> <= 2
assert abs(np.abs(B[:, 0]).max() - 2) < 1e-9 and abs(np.abs(B[:, 1]).max() - np.sqrt(2)) < 1e-9
phi = np.radians(30)
u = np.cos(phi) * u0 + np.sin(phi) * v
cap = section([R[k] for k in range(24) if k != k0], extra=[(-(u @ u0), -(u @ v), -1.0)])   # <x,u> >= 1

fig = plt.figure(figsize=(7.2, 3.1))
gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.05], left=0.02, right=0.985, top=0.9, bottom=0.17, wspace=0.28)
ax = fig.add_subplot(gs[0])
ax.set_aspect("equal"); ax.set_axis_off()
ax.add_patch(Polygon(B, closed=True, facecolor=to_rgb(BLUE) + (0.10,), edgecolor=BLUE, lw=0.9, zorder=1))
ax.add_patch(Polygon(Q, closed=True, facecolor=to_rgb(YELLOW) + (0.18,), edgecolor=YELLOW, lw=0.9, zorder=2))
ax.add_patch(Polygon(V, closed=True, facecolor=to_rgb(AQUA) + (0.16,), edgecolor=AQUA, lw=1.0, zorder=3))
ax.add_patch(Polygon(cap, closed=True, facecolor=to_rgb(ORANGE) + (0.45,), edgecolor=ORANGE, lw=0.8, zorder=4))
t = np.linspace(0, 2 * np.pi, 400)
ax.plot(np.cos(t), np.sin(t), color=INK2, lw=0.6, ls=(0, (3, 2)), zorder=3)
# the cutting line <x,u> = 1 in the plane: a cos(phi) + b sin(phi) = 1
L = np.array([[np.cos(phi) - 1.2 * np.sin(phi), np.sin(phi) + 1.2 * np.cos(phi)], [np.cos(phi) + 0.9 * np.sin(phi), np.sin(phi) - 0.9 * np.cos(phi)]])
ax.plot(L[:, 0], L[:, 1], color=ORANGE, lw=1.1, zorder=5)
ax.annotate("", xy=(np.cos(phi), np.sin(phi)), xytext=(0, 0), arrowprops=dict(arrowstyle="-|>", color=ORANGE, lw=1.0, mutation_scale=8), zorder=6)
ax.scatter([0], [0], s=8, color=INK, zorder=6)
for x, y, s, col, ha, va in ((2.08, 0, r"$2u_0$", YELLOW, "left", "center"), (-2.08, 0, r"$-2u_0$", BLUE, "right", "center"),
                              (0, np.sqrt(2) + 0.08, r"$\sqrt{2}\,v$", BLUE, "center", "bottom"), (0, -np.sqrt(2) - 0.08, r"$-\sqrt{2}\,v$", BLUE, "center", "top")):
    ax.text(x, y, s, fontsize=7, color=col, ha=ha, va=va)
ax.text(-0.55, 0.12, r"$\mathrm{V}$", fontsize=8, color=AQUA, ha="center", va="center")
ax.text(1.9, -0.62, r"$Q$", fontsize=7.5, color=YELLOW, ha="center", va="center")
ax.plot([1.78, 1.42], [-0.55, -0.3], color=YELLOW, lw=0.5)
ax.text(-1.5, 0.0, r"$B$", fontsize=8, color=BLUE, ha="center", va="center")
ax.text(0.66, 0.16, r"$u$", fontsize=7, color=ORANGE, ha="center", va="center")
ax.text(0.98, 1.38, r"$\langle x,u\rangle=1$", fontsize=6.8, color=ORANGE, ha="left", va="center")
ax.text(-0.78, -1.02, r"$|x|=1$", fontsize=6.5, color=INK2, ha="right", va="center")
ax.plot([-0.76, -0.55], [-1.0, -0.84], color=INK2, lw=0.5)
ax.set_xlim(-2.75, 2.75); ax.set_ylim(-1.8, 1.75)
ax.text(0.02, 0.98, "(a)", transform=ax.transAxes, fontsize=10, fontweight="bold", ha="left", va="top")

# (b) the cap volume along three great circles
SIGNS = signs


def cap_formula(w):
    a = np.abs(np.asarray(w, float)).copy()
    for i in range(4):
        for j in range(i + 1, 4):
            if abs(a[i] - a[j]) < 1e-9:
                a[j] += 1e-9 * (j + 1)
    a /= np.linalg.norm(a)
    s = 0.0
    for k in range(4):
        if 2 * a[k] <= 1:
            continue
        d = 1.0
        for l in range(4):
            if l != k:
                d *= a[k] ** 2 - a[l] ** 2
        s += (2 * a[k] - 1) ** 4 * a[k] ** 2 / d
    return s / 3


def cap_direct(w):
    A = np.vstack([SIGNS, -np.asarray(w, float).reshape(1, 4)])
    b = np.concatenate([2 * np.ones(16), [-1.0]])
    nr = np.linalg.norm(A, axis=1)
    cheb = linprog(np.r_[np.zeros(4), -1.0], A_ub=np.hstack([A, nr[:, None]]), b_ub=b, bounds=[(None, None)] * 4 + [(0, 50)], method="highs")
    if not cheb.success or cheb.x[-1] <= 1e-11:
        return 0.0
    hs = HalfspaceIntersection(np.hstack([A, -b[:, None]]), cheb.x[:4])
    return ConvexHull(hs.intersections, qhull_options="QJ").volume


e0 = np.array([1.0, 0, 0, 0])
targets = [(np.array([1, 1, 0, 0]) / np.sqrt(2), "towards an edge centre", BLUE),
           (np.array([1, 1, 1, 0]) / np.sqrt(3), "towards a facet centre", AQUA),
           (np.array([1, 1, 1, 1]) / 2, r"towards $\frac{1}{2}(1,1,1,1)$", YELLOW)]
ax2 = fig.add_subplot(gs[1]); clean_axes(ax2)
worst = 0.0
for tgt, lab, col in targets:
    amax = np.arccos(e0 @ tgt)
    ang = np.linspace(0, amax, 121)
    e1 = tgt - (tgt @ e0) * e0; e1 /= np.linalg.norm(e1)
    vol = []
    for i, a in enumerate(ang):
        w = np.cos(a) * e0 + np.sin(a) * e1
        vol.append(cap_direct(w))
        sq = np.sort(w ** 2)
        if np.diff(sq).min() > 1e-3:          # the divided difference is well conditioned away from coincident nodes
            worst = max(worst, abs(vol[-1] - cap_formula(w)))
    ax2.plot(np.degrees(ang), vol, color=col, lw=1.3, label=lab)
assert worst < 1e-7, worst
assert abs(cap_direct(e0) - 1 / 3) < 1e-9
ax2.axhline(1 / 3, color=INK2, lw=0.6, ls=(0, (3, 2)))
ax2.scatter([0], [1 / 3], s=18, color=ORANGE, zorder=5)
ax2.set_xlim(-1, 62); ax2.set_ylim(0, 0.36)
ax2.set_xticks([0, 15, 30, 45, 60])
ax2.set_yticks([0, 0.1, 0.2, 0.3, 1 / 3])
ax2.set_yticklabels(["0", "0.1", "0.2", "0.3", r"$\frac{1}{3}$"])
ax2.set_xlabel("angle from the vertex direction (degrees)")
ax2.set_ylabel("cap volume")
ax2.text(2, 0.345, r"the maximum $\frac{1}{3}$, at the vertex direction", fontsize=6.5, color=INK2, ha="left", va="bottom")
ax2.legend(loc="lower left", fontsize=6.6, frameon=True, framealpha=1, edgecolor="none", handlelength=1.6, borderpad=0.3, bbox_to_anchor=(0.0, 0.0))
panel_label(ax2, "(b)", x=-0.14, y=1.0)
save(fig, "fig_c_cap", outdir="../figures")
