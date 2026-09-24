"""
fig_noncontact.py -- neighbours that do not touch the centre (Section 2.7).

(a) Why the radial reduction needs every active neighbour to touch: two
    neighbours at distance 2.4, exactly 2 apart (a packing), at 49.25 degrees;
    pulled in to distance 2 they overlap, so the configuration at the
    all-contact corner is not a packing and not a contact configuration.
    A two-dimensional section, to scale.
(b) The distance criterion (Proposition 2.10): Phi for m contacts and k further
    neighbours at a common distance d, against d; it passes 8 at the
    thresholds that multi_cap/shell_reduction.py certifies in ball arithmetic.
(c) The numerical search of the open case (multi_cap/runs/
    shell_neighbour_search_delta_*.log): with one neighbour held at distance
    2 + delta, the least cell volume found over 60 starts, against the exact
    volume 25/3 - (1/3)(1 - delta/2)^4 of the deletion with that neighbour on
    the deleted root's axis.
"""
import glob
import math
import os
import re

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from mpmath import mp, mpf, pi, sin, cos, acos, tan, quad
from figstyle import *

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.join(HERE, "..", "..", "multi_cap", "runs")

# ---------------------------------------------------------------- (b) data
mp.dps = 20
T = 2 * pi ** 2
C = lambda r: pi * (2 * r - sin(2 * r))


def cap(d, r):
    x = (d / 2) * cos(r)
    return C(acos(x)) if x < 1 else mpf(0)


def Phi(ds):
    a = pi / 4
    pts = sorted(set([mpf(0), a] + [acos(2 / d) for d in ds if d > 2 and acos(2 / d) < a]))
    f = lambda r: max(mpf(0), T - sum(cap(d, r) for d in ds)) * 4 * tan(r) / cos(r) ** 4
    return float((T + quad(f, pts)) / 4)


CASES = [((22, 1), 2.16756, BLUE), ((21, 2), 2.07044, ORANGE),
         ((20, 3), 2.04503, AQUA), ((22, 2), 2.23400, VIOLET)]
dgrid = np.linspace(2.0, 2.3, 61)
curves = {}
for (m, k), t, col in CASES:
    curves[(m, k)] = [Phi([mpf(2)] * m + [mpf(d)] * k) for d in dgrid]

# ---------------------------------------------------------------- (c) data
pts = []
for f in sorted(glob.glob(os.path.join(RUNS, "shell_neighbour_search_delta_*_60_seed7.log"))):
    delta = float(re.search(r"delta_([0-9.]+)_60", f).group(1))
    txt = open(f).read()
    m = re.search(r"least volume with a shell neighbour: ([0-9.]+)", txt)
    if m:
        pts.append((2 + delta, float(m.group(1))))
pts.sort()

# ---------------------------------------------------------------- figure
fig = plt.figure(figsize=(7.2, 2.55))
gs = fig.add_gridspec(1, 3, width_ratios=[0.9, 1.05, 1.05], left=0.02, right=0.985,
                      top=0.87, bottom=0.2, wspace=0.52)

# (a)
ax = fig.add_subplot(gs[0])
ax.set_aspect("equal")
ax.axis("off")
half = math.radians(49.2486) / 2
dirs = [np.array([math.cos(s * half), math.sin(s * half)]) for s in (1, -1)]
ax.add_patch(Circle((0, 0), 1, facecolor=BLUE, alpha=0.25, edgecolor=BLUE, lw=0.8))
for u in dirs:
    ax.add_patch(Circle(tuple(2.4 * u), 1, facecolor="none", edgecolor=INK, lw=0.9))
    ax.add_patch(Circle(tuple(2.0 * u), 1, facecolor="none", edgecolor=ORANGE, lw=0.9,
                        ls=(0, (3, 2))))
# overlap of the pulled-in discs
th = np.linspace(0, 2 * np.pi, 400)
c1, c2 = 2 * dirs[0], 2 * dirs[1]
P = np.stack([np.cos(th), np.sin(th)], 1)
ring = c1 + P
inside = np.linalg.norm(ring - c2, axis=1) <= 1
ax.fill(*(np.vstack([ring[inside], (c2 + P)[np.linalg.norm(c2 + P - c1, axis=1) <= 1][::-1]]).T),
        color=ORANGE, alpha=0.55, lw=0)
ax.plot([0, 2.4 * dirs[0][0]], [0, 2.4 * dirs[0][1]], color=INK3, lw=0.5)
ax.plot([0, 2.4 * dirs[1][0]], [0, 2.4 * dirs[1][1]], color=INK3, lw=0.5)
ax.set_xlim(-1.15, 3.55)
ax.set_ylim(-2.95, 2.3)
ax.text(2.2, 2.05, "distance 2.4,\n2 apart", fontsize=6.2, ha="center", va="bottom", color=INK)
ax.text(1.15, -2.9, "pulled in to 2:\nthey overlap", fontsize=6.2, ha="center", va="bottom", color=ORANGE)
fig.text(0.012, 0.9, "(a)", fontsize=10, fontweight="bold", ha="left", va="bottom", color=INK)

# (b)
ax2 = fig.add_subplot(gs[1]); clean_axes(ax2)
for (m, k), t, col in CASES:
    ax2.plot(dgrid, curves[(m, k)], color=col, lw=1.1, label="%d + %d" % (m, k))
    ax2.plot([t], [8.0], "o", ms=2.8, color=col, zorder=4)
ax2.axhline(8.0, color=INK3, lw=0.7, ls=(0, (4, 2)))
ax2.set_xlim(2.0, 2.3); ax2.set_ylim(7.84, 8.34)
ax2.set_xticks([2.0, 2.1, 2.2, 2.3])
ax2.tick_params(labelsize=6.6)
ax2.set_xlabel("distance $d$ of the non-contacts", fontsize=7.2)
ax2.set_ylabel(r"$\Phi$", fontsize=7.4)
ax2.legend(loc="upper left", fontsize=5.6, frameon=True, framealpha=1, edgecolor="none", ncol=2,
           title="contacts + others", title_fontsize=5.6, handlelength=1.2, borderpad=0.3,
           columnspacing=0.8)
ax2.set_title("the distance criterion", fontsize=7.4, loc="right", color=INK2)
panel_label(ax2, "(b)", x=-0.3, y=1.0)

# (c)
ax3 = fig.add_subplot(gs[2]); clean_axes(ax3)
xs = np.linspace(2.0, 2.85, 200)
ax3.plot(xs, 25 / 3 - (1 / 3) * (1 - (xs - 2) / 2) ** 4, color=INK3, lw=1.0,
         label="on the deleted root's axis")
if pts:
    ax3.plot([p[0] for p in pts], [p[1] for p in pts], "o", ms=3.2, color=BLUE,
             label="least found, 60 starts")
ax3.axhline(8.0, color=INK3, lw=0.7, ls=(0, (4, 2)))
ax3.set_xlim(1.98, 2.87); ax3.set_ylim(7.95, 8.36)
ax3.set_xticks([2.0, 2.2, 2.4, 2.6, 2.8])
ax3.tick_params(labelsize=6.6)
ax3.set_xlabel("distance of the held neighbour", fontsize=7.2)
ax3.set_ylabel("cell volume", fontsize=7.2)
ax3.legend(loc="lower right", fontsize=5.8, frameon=True, framealpha=1, edgecolor="none",
           handlelength=1.4, borderpad=0.3)
ax3.set_title("the open case, searched", fontsize=7.4, loc="right", color=INK2)
panel_label(ax3, "(c)", x=-0.3, y=1.0)

save(fig, "fig_noncontact", outdir="../figures")
