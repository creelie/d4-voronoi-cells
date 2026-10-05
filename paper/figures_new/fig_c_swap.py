"""
fig_c_swap.py -- the volume defect along the swap paths.

The curves are the enumerated cell volumes written by
swap_configs/swap_defect_curves.py (swap_curves.dat, swap_mixed.dat); nothing
here is a schematic.

(a) The three swap paths, m = 2, 3, 6, over the packing-valid range
    (0, theta_*), theta_* = pi/3, on a logarithmic scale; the dashed line is
    the one combinatorial transition at theta = pi/6.
(b) The swap path against the mixed path, which differ in one sign; the
    crossover at theta = 0.6492 and the zero of the swap path at theta_*.
"""
import os

import numpy as np
import matplotlib.pyplot as plt
from figstyle import *

ROOT = "../../swap_configs"
for cand in ("../../swap_configs", ROOT, "../swap_configs", "swap_configs"):
    if os.path.exists(os.path.join(cand, "swap_curves.dat")):
        ROOT = cand
        break
d = np.loadtxt(os.path.join(ROOT, "swap_curves.dat"), skiprows=1)
m = np.loadtxt(os.path.join(ROOT, "swap_mixed.dat"), skiprows=1)
th, s2, s3, s6 = d.T
tm, swap, mixed = m.T
peak2, peak3, peak6 = s2.max(), s3.max(), s6.max()
assert abs(peak2 - 0.2317) < 5e-4 and abs(peak3 - 0.1692) < 5e-4 and abs(peak6 - 6.4e-4) < 1e-5
# the crossover, by linear interpolation of the sign change of swap - mixed
g = swap - mixed
k = np.where((g[:-1] > 0) & (g[1:] <= 0))[0][0]
tc = tm[k] - g[k] * (tm[k + 1] - tm[k]) / (g[k + 1] - g[k])
assert abs(tc - 0.6492) < 1e-3, tc

fig = plt.figure(figsize=(7.2, 3.35))
gs = fig.add_gridspec(1, 2, left=0.085, right=0.985, top=0.93, bottom=0.3, wspace=0.3)

ax = fig.add_subplot(gs[0]); clean_axes(ax)
pos = lambda y: np.where(y > 0, y, np.nan)
ax.semilogy(th, pos(s2), color=BLUE, lw=1.3, label=r"$m=2$")
ax.semilogy(th, pos(s3), color=AQUA, lw=1.3, label=r"$m=3$")
ax.semilogy(th, pos(s6), color=YELLOW, lw=1.3, label=r"$m=6$")
ax.axvline(np.pi / 6, color=INK2, lw=0.7, ls=(0, (3, 2)), zorder=1)
ax.set_xlim(0, np.pi / 3); ax.set_ylim(1e-6, 3)
ax.set_xticks([0, np.pi / 12, np.pi / 6, np.pi / 4, np.pi / 3])
ax.set_xticklabels([r"$0$", r"$\frac{\pi}{12}$", r"$\frac{\pi}{6}$", r"$\frac{\pi}{4}$", r"$\theta_*$"])
ax.set_xlabel(r"swap angle $\theta$")
ax.set_ylabel(r"volume defect $\Delta(\theta)$")
ax.text(np.pi / 6 + 0.03, 0.75, "the one combinatorial\ntransition, at $\\pi/6$", fontsize=6.3, color=INK2, ha="left", va="bottom",
        bbox=dict(facecolor=SURFACE, edgecolor="none", pad=0.4))
ax.text(np.pi / 6 + 0.06, 0.2317 * 1.35, r"$0.2317$", fontsize=6.5, color=BLUE, ha="left", va="bottom",
        bbox=dict(facecolor=SURFACE, edgecolor="none", pad=0.3))
ax.text(np.pi / 6 + 0.06, 6.4e-4 * 1.6, r"$6.4\times10^{-4}$", fontsize=6.5, color=YELLOW, ha="left", va="bottom",
        bbox=dict(facecolor=SURFACE, edgecolor="none", pad=0.3))
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=3, fontsize=7, frameon=False, handlelength=1.8, columnspacing=1.6)
panel_label(ax, "(a)", x=-0.16, y=1.0)

ax2 = fig.add_subplot(gs[1]); clean_axes(ax2)
ax2.semilogy(tm, pos(swap), color=BLUE, lw=1.3, label=r"swap path, $c_1=c_2=+1$")
ax2.semilogy(tm, pos(mixed), color=ORANGE, lw=1.3, label=r"mixed path, $c_1=-1,\ c_2=+1$")
ax2.axvline(tc, color=INK2, lw=0.7, ls=(0, (3, 2)), zorder=1)
ax2.scatter([tc], [np.interp(tc, tm, mixed)], s=14, color=INK, zorder=5)
ax2.set_xlim(0, 1.5); ax2.set_ylim(1e-4, 3)
ax2.set_xticks([0, 0.25, 0.5, 0.75, np.pi / 3, 1.25, 1.5])
ax2.set_xticklabels([r"$0$", r"$0.25$", r"$0.5$", r"$0.75$", r"$\theta_*$", r"$1.25$", r"$1.5$"])
ax2.set_xlabel(r"angle $\theta$")
ax2.set_ylabel(r"volume defect $\Delta(\theta)$")
ax2.text(tc + 0.03, 1.5, f"crossover at ${tc:.4f}$", fontsize=6.5, color=INK2, ha="left", va="bottom",
         bbox=dict(facecolor=SURFACE, edgecolor="none", pad=0.4))
ax2.text(np.pi / 3 + 0.04, 2.2e-4, r"$\Delta(\theta_*)=0$", fontsize=6.5, color=BLUE, ha="left", va="bottom",
         bbox=dict(facecolor=SURFACE, edgecolor="none", pad=0.3))
ax2.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2, fontsize=7, frameon=False, handlelength=1.8, columnspacing=1.6)
panel_label(ax2, "(b)", x=-0.16, y=1.0)
save(fig, "fig_c_swap", outdir="../figures")
