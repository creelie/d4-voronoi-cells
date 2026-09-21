"""
fig_free_volume.py -- the free-volume identity, drawn on real geometry.

The identity  vol(V_c) = vol(V) - vol(union of C_j) + vol(union of E_j)  holds
in every dimension for every choice of tilts, so it can be drawn faithfully in
three dimensions, where the pieces can be seen.  The reference cell is the
Voronoi cell of the face-centred cubic contact configuration, the rhombic
dodecahedron  V = {x : <x, w_i> <= 1, i = 1..12}  with the twelve unit
directions w_i = (+-1, +-1, 0)/sqrt 2 and permutations, of volume 4 sqrt 2.

(a) One direction w_j is tilted by phi to w_j'.  The cap C_j = V and
    {<x, w_j'> >= 1} is what the tilted half-space cuts off the reference
    cell; the correction E_j = V_c and {<x, w_j> >= 1} is what the motion frees
    beyond the old facet.  Both are convex polytopes, computed as half-space
    intersections, and drawn to scale inside V.
(b) The three volumes as functions of the tilt, and the identity checked at
    every angle to machine precision.
(c) Two directions tilted at once: the caps overlap, so the correction is the
    volume of a union and not a sum of volumes.
"""
import itertools

import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull, HalfspaceIntersection
from figstyle import *

W = []
for i, j in itertools.combinations(range(3), 2):
    for si in (1, -1):
        for sj in (1, -1):
            v = np.zeros(3); v[i], v[j] = si, sj
            W.append(v / np.sqrt(2))
W = np.array(W)                                  # 12 x 3 unit vectors, fcc contacts


def cell(normals, bound=4.0):
    """the polytope {<x, n> <= 1 for n in normals}, clipped to a box that never binds"""
    hs = [np.append(n, -1.0) for n in normals]
    for k in range(3):
        for s in (1, -1):
            e = np.zeros(3); e[k] = s
            hs.append(np.append(e, -bound))
    H = HalfspaceIntersection(np.array(hs), np.zeros(3))
    return ConvexHull(H.intersections)


def piece(normals, extra_ge):
    """polytope {<x,n> <= 1 for n in normals} and {<x,m> >= 1 for m in extra_ge}; None if empty"""
    hs = [np.append(n, -1.0) for n in normals] + [np.append(-m, 1.0) for m in extra_ge]
    # interior point: solve a small LP by trying candidate points along the extra normals
    from scipy.optimize import linprog
    A = np.array([h[:3] for h in hs]); b = -np.array([h[3] for h in hs])
    # maximise slack t: A x + t <= b
    c = np.zeros(4); c[3] = -1
    Aub = np.hstack([A, np.ones((len(hs), 1))])
    res = linprog(c, A_ub=Aub, b_ub=b, bounds=[(None, None)] * 3 + [(0, None)])
    if not res.success or res.x[3] <= 0:
        return None
    H = HalfspaceIntersection(np.array(hs), res.x[:3])
    return ConvexHull(H.intersections)


def tilt(w, phi, seed=3):
    rng = np.random.default_rng(seed)
    e = rng.normal(size=3); e -= e @ w * w; e /= np.linalg.norm(e)
    return np.cos(phi) * w + np.sin(phi) * e


def draw_hull(ax, hull, face, edge, alpha, lw=0.6, zorder=2):
    pts = hull.points
    faces = []
    # merge coplanar simplices into polygons by their outward normal
    groups = {}
    for simplex, eq in zip(hull.simplices, hull.equations):
        key = tuple(np.round(eq, 6))
        groups.setdefault(key, set()).update(simplex)
    for key, idx in groups.items():
        n = np.array(key[:3])
        P = pts[list(idx)]
        c = P.mean(axis=0)
        u = np.cross(n, [0.3, 0.5, 0.8]); u /= np.linalg.norm(u); v = np.cross(n, u)
        ang = np.arctan2((P - c) @ v, (P - c) @ u)
        faces.append(P[np.argsort(ang)])
    ls = np.array([0.4, -0.6, 0.7]); ls /= np.linalg.norm(ls)
    cols = []
    for key in groups:
        shade = 0.72 + 0.28 * abs(np.array(key[:3]) @ ls)
        cols.append(np.append(np.clip(np.array(to_rgb(face)) * shade, 0, 1), alpha))
    coll = Poly3DCollection(faces, facecolors=cols, edgecolors=[np.append(to_rgb(edge), min(1, alpha + 0.35))] * len(faces), linewidths=lw)
    coll.set_zorder(zorder)
    ax.add_collection3d(coll)


V = cell(W)
volV = V.volume
assert abs(volV - 4 * np.sqrt(2)) < 1e-9, volV

cam = np.array([np.cos(np.radians(20)) * np.cos(np.radians(-47)), np.cos(np.radians(20)) * np.sin(np.radians(-47)), np.sin(np.radians(20))])
j = int(np.argmax(W @ cam))
phi0 = np.radians(22)
wj = tilt(W[j], phi0)
Wc = W.copy(); Wc[j] = wj
Vc = cell(Wc)
C = piece(list(W), [wj])                          # cut off from V by the tilted half-space
E = piece(list(Wc), [W[j]])                       # freed beyond the old facet
assert C is not None and E is not None
ident = Vc.volume - (volV - C.volume + E.volume)
assert abs(ident) < 1e-9, ident

fig = plt.figure(figsize=(7.2, 3.4))
gs = fig.add_gridspec(1, 3, width_ratios=[1.3, 1.0, 1.15], left=0.0, right=0.99, top=0.93, bottom=0.17, wspace=0.28)

ax = fig.add_subplot(gs[0], projection="3d")
blank_3d(ax); ax.computed_zorder = False
ax.view_init(elev=20, azim=-47)
draw_hull(ax, V, "#c9ccd2", INK3, 0.22, lw=0.7, zorder=1)
draw_hull(ax, C, ORANGE, ORANGE, 0.85, lw=0.6, zorder=4)
draw_hull(ax, E, AQUA, AQUA, 0.85, lw=0.6, zorder=4)
# the two planes near the facet, as small quads
for nrm, col in ((W[j], INK2), (wj, ORANGE)):
    u = np.cross(nrm, [0.2, 0.7, 0.4]); u /= np.linalg.norm(u); v = np.cross(nrm, u)
    q = [nrm + 0.95 * (a * u + b * v) for a, b in ((1, 1), (1, -1), (-1, -1), (-1, 1))]
    pc = Poly3DCollection([q], facecolors=[np.append(to_rgb(col), 0.10)], edgecolors=[np.append(to_rgb(col), 0.9)], linewidths=0.8, linestyles=(0, (4, 2)))
    pc.set_zorder(3); ax.add_collection3d(pc)
ax.quiver(0, 0, 0, *W[j], color=INK2, lw=1.0, arrow_length_ratio=0.12, zorder=5)
ax.quiver(0, 0, 0, *wj, color=ORANGE, lw=1.0, arrow_length_ratio=0.12, zorder=5)
lim = 1.42
ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_zlim(-lim, lim); ax.set_box_aspect((1, 1, 1))
ax.text2D(0.02, 0.97, "(a)", transform=ax.transAxes, fontsize=10, fontweight="bold")
ax.text2D(0.5, 0.0, r"$\mathrm{V}$ (grey), $C_j$ (orange), $E_j$ (green); tilt $\varphi=$" + f"{np.degrees(phi0):.0f}" + r"$^\circ$",
          transform=ax.transAxes, ha="center", va="top", fontsize=6.6, color=INK2)
ax.text2D(0.03, 0.88, r"$w_j$", transform=ax.transAxes, fontsize=7.5, color=INK2)
ax.text2D(0.03, 0.82, r"$w_j'=\cos\varphi\,w_j+\sin\varphi\,e$", transform=ax.transAxes, fontsize=7.0, color=ORANGE)

# (b) volumes against the tilt, identity checked at every angle
ax2 = fig.add_subplot(gs[1]); clean_axes(ax2)
phis = np.radians(np.linspace(0.5, 40, 80))
vC, vE, vVc, err = [], [], [], []
for p in phis:
    wp = tilt(W[j], p); Wp = W.copy(); Wp[j] = wp
    Vp = cell(Wp)
    Cp = piece(list(W), [wp]); Ep = piece(list(Wp), [W[j]])     # the two pieces, each as its own polytope
    c_ = Cp.volume if Cp is not None else 0.0
    e_ = Ep.volume if Ep is not None else 0.0
    vC.append(c_); vE.append(e_); vVc.append(Vp.volume); err.append(abs(Vp.volume - (volV - c_ + e_)))
deg = np.degrees(phis)
ax2.plot(deg, vC, color=ORANGE, lw=1.3, label=r"$\mathrm{vol}(C_j)$, cut off")
ax2.plot(deg, vE, color=AQUA, lw=1.3, label=r"$\mathrm{vol}(E_j)$, freed")
ax2.plot(deg, np.array(vVc) - volV, color=BLUE, lw=1.3, label=r"$\mathrm{vol}(V_c)-\mathrm{vol}(\mathrm{V})$")
ax2.axhline(0, color=INK2, lw=0.7)
ax2.set_xlabel(r"tilt $\varphi$ (degrees)")
ax2.set_ylabel("volume")
ax2.legend(fontsize=6.2, frameon=True, framealpha=1, edgecolor="none", loc="upper left")
ax2.text(0.6, 0.182, f"identity at 80 angles:\nlargest residual {max(err):.1e}",
         ha="left", va="top", fontsize=5.8, color=INK2, bbox=dict(facecolor=SURFACE, edgecolor="none", pad=0.5))
panel_label(ax2, "(b)", x=-0.24)

# (c) two tilts: caps overlap
ax3 = fig.add_subplot(gs[2], projection="3d")
blank_3d(ax3); ax3.computed_zorder = False
ax3.view_init(elev=20, azim=-47)
k = [i for i in range(12) if i != j and abs(W[i] @ W[j] - 0.5) < 1e-9][0]     # a neighbour at 60 degrees
def toward(w, target, phi):
    e = target - (target @ w) * w; e /= np.linalg.norm(e)
    return np.cos(phi) * w + np.sin(phi) * e
wj = toward(W[j], W[k], np.radians(24)); wk = toward(W[k], W[j], np.radians(24))
W2 = W.copy(); W2[j] = wj; W2[k] = wk
Cj = piece(list(W), [wj]); Ck = piece(list(W), [wk]); Cjk = piece(list(W), [wj, wk])
draw_hull(ax3, V, "#c9ccd2", INK3, 0.22, lw=0.7, zorder=1)
draw_hull(ax3, Cj, ORANGE, ORANGE, 0.55, lw=0.6, zorder=4)
draw_hull(ax3, Ck, YELLOW, YELLOW, 0.55, lw=0.6, zorder=4)
if Cjk is not None:
    draw_hull(ax3, Cjk, VIOLET, VIOLET, 0.95, lw=0.6, zorder=6)
ax3.set_xlim(-lim, lim); ax3.set_ylim(-lim, lim); ax3.set_zlim(-lim, lim); ax3.set_box_aspect((1, 1, 1))
ax3.text2D(0.02, 0.97, "(c)", transform=ax3.transAxes, fontsize=10, fontweight="bold")
ov = Cjk.volume if Cjk is not None else 0.0
ax3.text2D(0.5, 0.0, f"two caps: overlap (violet) of volume {ov:.3f};\n" + r"$\mathrm{vol}(C_j\cup C_k)=\mathrm{vol}(C_j)+\mathrm{vol}(C_k)-$" + f"{ov:.3f}",
           transform=ax3.transAxes, ha="center", va="top", fontsize=6.4, color=INK2)
save(fig, "fig_free_volume", outdir="../figures")
