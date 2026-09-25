"""
fig_labelled.py -- the labelled certificate of Theorem 2.17 (Section 2.8).

(a) A plane section through c = 0 and two centres y_i, y_j at the packing
    boundary |y_i - y_j| = 2, with d_i = 2.10 and d_j = 2, so that the
    cosine u = amax(d_i, d_j) = 0.525 of their angle exceeds 1/2, which two
    contacts cannot do.  Near the sphere of radius sqrt(3/2): the hyperplanes
    of the two centres at heights h = d/2, the caps beyond them (orange,
    aqua) and their overlap, the pair term (violet).  Moving the hyperplane
    of y_i out to height tau removes, from the pair term, the part of the
    section at height tau beyond the hyperplane of y_j (thick), a share fr of
    that section; the proof bounds it by its value at the other centre's
    height 1, and compares it with 1/11 of the section, which the cap of y_i
    gives back.  To scale.
(b) The bound along two lines of the domain past the face t = 1/2 that the
    contacts end at: Q0 = 1000 sum omega - P alone (dashed) turns negative,
    but the packing makes the centres move out, and Q0 plus the lower bound
    1000 sum Gamma_i of labelled_certificate_check.py (solid) stays
    positive.  Blue: u = v = -sqrt3/2, where (C) comes closest to 0 on the
    face; orange: u = v = t.  The strip is the region II_s of the proof.
"""
import math
import sys
import os

import numpy as np
from fractions import Fraction as Fr
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, FancyArrowPatch, Circle
from figstyle import *

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'multi_cap'))
here = os.getcwd()
os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'multi_cap'))
import certificate_check as CC
import labelled_certificate_check as LCC
from mpmath import iv

R = math.sqrt(1.5)

fig = plt.figure(figsize=(7.2, 2.9))
gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.25], left=0.01, right=0.985, top=0.93, bottom=0.15, wspace=0.16)

# ------------------------------------------------------------------ (a)
ax = fig.add_subplot(gs[0])
ax.set_aspect("equal"); ax.set_axis_off()
di, dj = 2.10, 2.00
hi_, hj_ = di / 2, dj / 2
u = (di * di + dj * dj - 4) / (2 * di * dj)
g = math.acos(u)
thi, thj = g / 2, -g / 2
wi = np.array([math.cos(thi), math.sin(thi)]); wj = np.array([math.cos(thj), math.sin(thj)])


def arc(a0, a1, n=300):
    t = np.linspace(a0, a1, n)
    return np.column_stack([R * np.cos(t), R * np.sin(t)])


def draw_section(ax, lw_scale=1.0, labels=True):
    ai, aj = math.acos(hi_ / R), math.acos(hj_ / R)
    ax.add_patch(Circle((0, 0), R, facecolor=(*to_rgb(BLUE), 0.07), edgecolor="none", zorder=0))
    ax.add_patch(Polygon(arc(thi - ai, thi + ai), closed=True, facecolor=(*to_rgb(ORANGE), 0.25), edgecolor="none", zorder=1))
    ax.add_patch(Polygon(arc(thj - aj, thj + aj), closed=True, facecolor=(*to_rgb(AQUA), 0.25), edgecolor="none", zorder=1))
    p = np.linalg.solve(np.array([wi, wj]), [hi_, hj_])
    lo, hi = max(thi - ai, thj - aj), min(thi + ai, thj + aj)
    ax.add_patch(Polygon([p.tolist()] + arc(lo, hi).tolist(), closed=True, facecolor=(*to_rgb(VIOLET), 0.6),
                         edgecolor=VIOLET, linewidth=0.5 * lw_scale, zorder=3))
    t = np.linspace(-math.pi, math.pi, 600)
    ax.plot(R * np.cos(t), R * np.sin(t), color=BLUE, lw=1.1 * lw_scale, zorder=4)
    for w, h, col, a in ((wi, hi_, ORANGE, ai), (wj, hj_, AQUA, aj)):
        n = np.array([-w[1], w[0]]); half = math.sqrt(R * R - h * h)
        P = h * w[None, :] + np.array([-half, half])[:, None] * n[None, :]
        ax.plot(P[:, 0], P[:, 1], color=col, lw=1.1 * lw_scale, zorder=5)
    return p


p = draw_section(ax)
# the section of y_i at height tau, and its part beyond the hyperplane of y_j
tau = 1.075
n = np.array([-wi[1], wi[0]]); half = math.sqrt(R * R - tau * tau)
c0 = tau * wi
ax.plot(*zip(c0 - half * n, c0 + half * n), color=INK2, lw=0.9, ls=(0, (3, 1.6)), zorder=6)
s_cut = (hj_ - c0 @ wj) / (n @ wj)
ends = sorted([s_cut, -half if n @ wj < 0 else half])
ax.plot(*zip(c0 + ends[0] * n, c0 + ends[1] * n), color=INK, lw=2.6, solid_capstyle="butt", zorder=7)
ax.plot(0, 0, marker="o", ms=3, color=INK)
X0, X1, Y0, Y1 = 0.99, 1.30, -0.215, 0.33
ax.set_xlim(X0, X1); ax.set_ylim(Y0, Y1)
ax.text(1.236, 0.222, r"$\partial B(\sqrt{3/2})$", fontsize=6.8, color=BLUE, ha="left")
ax.text(1.005, 0.305, r"hyperplane of $y_i$, height $h_i=1.05$", fontsize=6.3, color="#c4501f", ha="left", va="top", bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="none", alpha=0.9))
ax.text(1.005, -0.2, r"hyperplane of $y_j$, height $h_j=1$", fontsize=6.3, color="#138a5f", ha="left", va="bottom", bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="none", alpha=0.9))
lab = c0 + (ends[1] + 0.105) * n - 0.017 * wi
ax.text(lab[0], lab[1], r"section at height $\tau$", fontsize=6.1, color=INK2, ha="center", va="center",
        rotation=math.degrees(thi) + 90 - 180, rotation_mode="anchor", bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="none", alpha=0.9))
ax.annotate("pair term", xy=(p[0] + 0.02, p[1] - 0.004), xytext=(1.255, -0.07), fontsize=6.4, color=VIOLET, bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="none", alpha=0.9),
            arrowprops=dict(arrowstyle="-", color=VIOLET, lw=0.5), ha="left", va="center")
ax.annotate(r"share $\mathrm{fr}$", xy=tuple(c0 + (ends[0] + ends[1]) / 2 * n), xytext=(1.255, 0.12), fontsize=6.4,
            color=INK, arrowprops=dict(arrowstyle="-", color=INK, lw=0.5), ha="left", va="center")
# inset: the whole disc, the two directions and the zoom window
axi = ax.inset_axes([0.0, 0.33, 0.33, 0.33])
axi.set_aspect("equal"); axi.set_axis_off()
draw_section(axi, lw_scale=0.6)
for w in (wi, wj):
    axi.add_patch(FancyArrowPatch((0, 0), tuple(1.72 * w), arrowstyle="-|>", mutation_scale=5, color=INK3, lw=0.5, zorder=6))
axi.plot([X0, X1, X1, X0, X0], [Y0, Y0, Y1, Y1, Y0], color=INK, lw=0.5, zorder=8)
axi.plot(0, 0, marker="o", ms=1.6, color=INK)
axi.set_xlim(-1.3, 1.95); axi.set_ylim(-1.3, 1.3)
axi.text(1.72 * wi[0], 1.72 * wi[1] + 0.12, r"$y_i$", fontsize=6.0, color=INK2, ha="center")
axi.text(1.72 * wj[0], 1.72 * wj[1] - 0.3, r"$y_j$", fontsize=6.0, color=INK2, ha="center")
fig.text(0.012, 0.95, "(a)", fontsize=10, fontweight="bold", ha="left", va="top", color=INK)

# ------------------------------------------------------------------ (b)
ax = fig.add_subplot(gs[1]); clean_axes(ax)
d = 8
Z = np.load('continuation_out/certificate_d%d.npz' % d)
f = [Fr(float(x)) for x in Z['f']]
F = [[[Fr(float(x)) for x in row] for row in Z['F%d' % k]] for k in range(d + 1)]
f = [f[0]] + [max(x, Fr(0)) for x in f[1:]]
P, _ = CC.build_P(d, f, F)
mons = np.array(list(P.keys())); coef = np.array([float(c) for c in P.values()])
om0 = CC.closed_forms()[0]


def omega(x):
    if x <= 1 / 3:
        return 0.0
    return float(om0(iv.mpf(math.acos(x))).mid)


def Q0(u, v, t):
    return 1000 * (omega(u) + omega(v) + omega(t)) - float(np.sum(coef * u ** mons[:, 0] * v ** mons[:, 1] * t ** mons[:, 2]))


G = LCC.GammaTables()
aD = float(LCC.amax(LCC.D, LCC.D))
ts = np.linspace(0.5, aD - 1e-9, 90)
lines = {}
for name, fn in (("A", lambda x: (-math.sqrt(3) / 2, -math.sqrt(3) / 2, x)), ("B", lambda x: (x, x, x))):
    q0 = np.array([Q0(*fn(x)) for x in ts])
    L = np.array([fn(x) for x in ts]).T
    lines[name] = (q0, q0 + G.gain(L, L))
ax.axvspan(0.5, 0.51, color=PANEL, zorder=0)
ax.axhline(0, color=INK2, lw=0.7, zorder=1)
ax.plot(ts, lines["A"][0], color=BLUE, lw=1.1, ls=(0, (4, 2)), zorder=3)
ax.plot(ts, lines["A"][1], color=BLUE, lw=1.5, zorder=4)
ax.plot(ts, lines["B"][0], color=ORANGE, lw=1.1, ls=(0, (4, 2)), zorder=3)
ax.plot(ts, lines["B"][1], color=ORANGE, lw=1.5, zorder=4)
ax.set_xlim(0.5, 0.5745); ax.set_ylim(-3.2, 5.2)
ax.set_xticks([0.50, 0.52, 0.54, 0.56]); ax.set_yticks([-3, -2, -1, 0, 1, 2, 3, 4, 5])
ax.tick_params(labelsize=6.6)
ax.set_xlabel(r"largest inner product $t$ of the three directions", fontsize=7.2)
ax.set_ylabel(r"bound on $Q$ (scale $1000$)", fontsize=7.2)
ax.text(0.5055, -2.95, r"II$_s$", fontsize=6.8, color=INK2, ha="center")
iA = 58; iB = 40
ax.text(ts[iA] + 0.001, lines["A"][1][iA] - 0.75, r"$u=v=-\sqrt{3}/2$, with $\Gamma$", fontsize=6.4, color=BLUE, ha="left")
ax.text(ts[iB] - 0.0005, lines["B"][1][iB] + 0.5, r"$u=v=t$, with $\Gamma$", fontsize=6.4, color=ORANGE, ha="right")
ax.text(0.5705, lines["B"][0][-8] + 0.55, r"$Q_0$ alone", fontsize=6.4, color=INK2, ha="center")
panel_label(ax, "(b)", x=-0.1, y=0.98)

os.chdir(here)
save(fig, "fig_labelled", outdir="../figures")
