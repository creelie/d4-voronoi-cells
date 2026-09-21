"""
fig_c_counting.py -- the contact count and the pairwise estimate at 23.

Panels (a) and (b), the covering bound against the count and the map of how
each part of the range is settled, are kept from the earlier composite.
Panels (c) and (d) are computed here from the formulas of the section on the
overlaps: the lens measure Lambda(r, gamma), the per-pair weight

    omega(gamma) = (1/4) int_0^R Lambda(r, gamma) d(sec^4 r),

and the bracket (1/4)[2 pi^2 + int_0^R (2 pi^2 - m C(r)) d(sec^4 r)] at m = 23,
for the two integration limits R = r_23 (where the bracket is the covering
bound (pi m/3) tan^3 r_23) and R = r_* = arccos sqrt(2/3).  The values the
paper quotes, omega(60 deg) = 0.000922... and 0.001445..., the bracket
7.907144... at r_*, the deletion at 88 pairs and the crossings at 90.3 and 64.2
pairs, are asserted before anything is drawn.
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.optimize import brentq
from figstyle import *
from oldpanels import load, split_rows, trim, show

m = 23
r23 = brentq(lambda r: 2 * r - np.sin(2 * r) - 2 * np.pi / m, 0.1, 1.5)
rs = np.arccos(np.sqrt(2 / 3))


def Lam(r, g):
    if r <= g / 2:
        return 0.0
    return 4 * np.pi * quad(lambda t: np.sin(t) ** 2 * (1 - np.tan(g / 2) / np.tan(t)), g / 2, r)[0]


def omega(g, R):
    if R <= g / 2:
        return 0.0
    return 0.25 * quad(lambda r: Lam(r, g) * 4 / np.cos(r) ** 4 * np.tan(r), g / 2, R)[0]


def bracket(R):
    C = lambda r: np.pi * (2 * r - np.sin(2 * r))
    return 0.25 * (2 * np.pi ** 2 + quad(lambda r: (2 * np.pi ** 2 - m * C(r)) * 4 / np.cos(r) ** 4 * np.tan(r), 0, R)[0])


w_s, w_23 = omega(np.pi / 3, rs), omega(np.pi / 3, r23)
A_s, A_23 = bracket(rs), bracket(r23)
assert abs(w_s - 0.0014454) < 1e-6 and abs(w_23 - 0.0009222) < 1e-6
assert abs(A_s - 7.907144) < 1e-6 and abs(A_23 - np.pi * m / 3 * np.tan(r23) ** 3) < 1e-9
n_s, n_23 = (8 - A_s) / w_s, (8 - A_23) / w_23
assert abs(n_s - 64.2) < 0.05 and abs(n_23 - 90.3) < 0.05, (n_s, n_23)
del_s, del_23 = A_s + 88 * w_s, A_23 + 88 * w_23
assert abs(del_s - 8.034340) < 1e-5 and abs(del_23 - 7.997885) < 1e-5, (del_s, del_23)

old = load("old/fig_c_counting.png")
top = trim(split_rows(old, 20)[0])
W_IN = 7.0
top_h = W_IN * top.shape[0] / top.shape[1]

fig = plt.figure(figsize=(W_IN, top_h + 2.6))
gs = fig.add_gridspec(2, 2, height_ratios=[top_h, 2.2], hspace=0.22, wspace=0.32,
                      left=0.09, right=0.985, top=0.99, bottom=0.09)
axt = fig.add_subplot(gs[0, :])
show(axt, top)

# (c) the estimate against the number of pairs at 60 degrees
ax = fig.add_subplot(gs[1, 0]); clean_axes(ax)
n = np.linspace(0, 115, 3)
ax.plot(n, A_s + n * w_s, color=AQUA, lw=1.3, label=r"to $r_*$")
ax.plot(n, A_23 + n * w_23, color=INK3, lw=1.3, label=r"to $r_{23}$")
ax.axhline(8, color=ORANGE, lw=0.8, ls=(0, (4, 2)))
ax.vlines(88, 7.88, 8.06, color=YELLOW, lw=0.7)
ax.vlines(115, 7.88, 7.96, color=INK3, lw=0.6, ls=(0, (1.5, 2)))
ax.scatter([88, 88], [del_s, del_23], s=16, color=YELLOW, zorder=5)
ax.scatter([n_s], [8], s=16, color=AQUA, zorder=5)
ax.scatter([n_23], [8], s=16, color=INK3, zorder=5)
ax.set_xlim(0, 132); ax.set_ylim(7.88, 8.07)
ax.set_xticks([0, 20, 40, 64.2, 88, 115])
ax.set_xticklabels(["0", "20", "40", "64.2", "88", "115"], fontsize=6.3)
ax.set_yticks([7.90, 7.94, 7.98, 8.02, 8.06])
ax.tick_params(labelsize=6.3)
ax.set_xlabel(r"pairs at exactly $60^\circ$", fontsize=7.5)
ax.set_ylabel("second-order estimate", fontsize=7.5)
ax.text(2, 8.003, "target 8", color=ORANGE, fontsize=6.3, ha="left", va="bottom")
ax.text(86, 8.040, "deletion, 8.034340", color=YELLOW, fontsize=6.0, ha="right", va="bottom")
ax.text(92, 7.981, "deletion, 7.997885", color=YELLOW, fontsize=6.0, ha="left", va="top", bbox=dict(facecolor=SURFACE, edgecolor="none", pad=0.3))
ax.text(61, 8.004, "64.2 pairs", color=AQUA, fontsize=6.0, ha="right", va="bottom")
ax.text(93, 7.996, "90.3 pairs", color=INK3, fontsize=6.0, ha="left", va="top")
ax.text(113.5, 7.895, "ceiling", color=INK3, fontsize=6.0, ha="right", va="bottom", rotation=90)
ax.legend(loc="upper left", fontsize=6.3, frameon=True, framealpha=1, edgecolor="none", handlelength=1.5, borderpad=0.25, bbox_to_anchor=(0.0, 1.0))
panel_label(ax, "(c)", x=-0.3, y=1.0)

# (d) the per-pair weight
ax2 = fig.add_subplot(gs[1, 1]); clean_axes(ax2)
gam = np.radians(np.linspace(60, 71, 120))
o_s = np.array([omega(g, rs) for g in gam]) * 1e4
o_23 = np.array([omega(g, r23) for g in gam]) * 1e4
ax2.plot(np.degrees(gam), o_s, color=AQUA, lw=1.3, label=r"to $r_*$")
ax2.plot(np.degrees(gam), o_23, color=INK3, lw=1.3, label=r"to $r_{23}$")
ax2.scatter([60, 60], [w_s * 1e4, w_23 * 1e4], s=14, color=[AQUA, INK3], zorder=5)
ax2.vlines([2 * np.degrees(r23), 2 * np.degrees(rs)], 0, 3, color=INK3, lw=0.5, ls=(0, (1.5, 2)))
ax2.set_xlim(59.5, 72.2); ax2.set_ylim(0, 16)
ax2.set_xticks([60, 62, 64, 66, 68, 2 * np.degrees(r23), 2 * np.degrees(rs)])
ax2.set_xticklabels(["60", "62", "64", "66", "68", "", ""], fontsize=6.3)
ax2.text(2 * np.degrees(r23) - 0.15, 3.4, r"$2r_{23}$", fontsize=6.3, ha="right", va="bottom", color=INK2, bbox=dict(facecolor=SURFACE, edgecolor="none", pad=0.3))
ax2.text(2 * np.degrees(rs) + 0.15, 3.4, r"$2r_*$", fontsize=6.3, ha="left", va="bottom", color=INK2, bbox=dict(facecolor=SURFACE, edgecolor="none", pad=0.3))
ax2.tick_params(labelsize=6.3)
ax2.set_xlabel(r"pair angle $\gamma$ in degrees", fontsize=7.5)
ax2.set_ylabel(r"$\omega(\gamma)$, units of $10^{-4}$", fontsize=7.5)
ax2.text(60.5, 14.7, r"$1.4454\times10^{-3}$", color=AQUA, fontsize=6.3, ha="left", va="bottom")
ax2.text(63.6, 9.0, r"$9.222\times10^{-4}$", color=INK3, fontsize=6.3, ha="left", va="bottom")
ax2.legend(loc="upper right", fontsize=6.5, frameon=True, framealpha=1, edgecolor="none", handlelength=1.6, borderpad=0.3, bbox_to_anchor=(1.0, 1.0))
panel_label(ax2, "(d)", x=-0.3, y=1.0)
save(fig, "fig_c_counting", outdir="../figures")
