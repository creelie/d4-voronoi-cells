"""
fig_c_verify.py -- the two steps of the certificate verification that the
authors' own code leaves out, as carried out in the directory zonal.

(a) Why step 3 fits in a few hundred megabytes: the counts are the counters
    of the kernel (zonal/runs/triple_counts.log, step3_part*.log), the size of
    its memo table on the largest entry, and the binomial bound C(20,6) on the
    number of coefficients of P(S).
(b) Step 5 for the four-point constraint: the zonal side of the constraint
    filling up over the run, signature by signature and through the three
    splittings of four points into two pairs (zonal/runs/
    step5_constraint_4_zonal.log).  It stops at 53572 coefficients, where the
    sum-of-squares side stops as well; the two agree term by term.
"""
import numpy as np
import matplotlib.pyplot as plt
from figstyle import *

CASCADE = [
    ("monomial triples,\nbefore the parity rule", 2.392e11, r"$2.39\times10^{11}$"),
    ("triples the parity\nrule leaves", 7.854e9, r"$7.85\times10^{9}$"),
    ("distinct canonical\nforms, largest entry", 2505430, r"$2\,505\,430$"),
    ("coefficients of $P(S)$,\nany signature", 38760, r"$38\,760$"),
]
TRACE = [(17, 81), (26, 4104), (302, 20536), (1106, 34552), (1396, 36263),
         (1397, 36267), (1408, 36399), (1705, 38711), (2519, 46378),
         (2810, 47847), (2811, 47851), (2822, 47923), (3122, 48623),
         (3942, 52351), (4236, 53572)]
SPLITS = [(1396, "(01)(23)"), (2810, "(02)(13)"), (4236, "(03)(12)")]
FINAL = 53572

fig = plt.figure(figsize=(7.2, 2.7))
gs = fig.add_gridspec(1, 2, width_ratios=[1.15, 1.0], left=0.2, right=0.985, top=0.9, bottom=0.2, wspace=0.42)

ax = fig.add_subplot(gs[0])
for side in ("top", "right"):
    ax.spines[side].set_visible(False)
ax.grid(True, axis="x", zorder=0); ax.set_axisbelow(True)
names = [n for n, v, l in CASCADE][::-1]
vals = np.array([v for n, v, l in CASCADE][::-1])
labs = [l for n, v, l in CASCADE][::-1]
y = np.arange(4)
ax.barh(y, vals, height=0.58, color=BLUE, zorder=3)
ax.set_yticks(y); ax.set_yticklabels(names, fontsize=6.8)
ax.set_xscale("log"); ax.set_xlim(1e4, 3e13)
ax.set_xticks([1e4, 1e6, 1e8, 1e10, 1e12])
ax.tick_params(labelsize=6.8)
ax.set_xlabel("count", fontsize=7.5)
for i, (v, l) in enumerate(zip(vals, labs)):
    ax.text(v * 2.0, i, l, va="center", ha="left", fontsize=6.8, color=INK2,
            bbox=dict(facecolor=SURFACE, edgecolor="none", pad=0.4))
ax.set_title("what step 3 visits", fontsize=7.8, loc="right", color=INK2)
fig.text(0.012, 0.9, "(a)", fontsize=10, fontweight="bold", ha="left", va="bottom", color=INK)

ax2 = fig.add_subplot(gs[1]); clean_axes(ax2)
t = np.array([a for a, b in TRACE]) / 60.0
n = np.array([b for a, b in TRACE])
ax2.plot(t, n, "-o", color=BLUE, ms=2.6, lw=1.1, label="zonal side, as it fills")
ax2.axhline(FINAL, color=ORANGE, ls=(0, (4, 2)), lw=1.0, label=f"sum-of-squares side, {FINAL} terms")
for k, (sec, name) in enumerate(SPLITS):
    ax2.vlines(sec / 60.0, 22000, 80000, color=INK3, lw=0.6, ls=(0, (1.5, 2)), zorder=1)
    if k < 2:
        ax2.text(sec / 60.0 + 0.8, 64500, name, fontsize=5.8, ha="left", va="bottom", color=INK2, bbox=dict(facecolor=SURFACE, edgecolor="none", pad=0.3))
    else:
        ax2.text(sec / 60.0 - 0.8, 58500, name, fontsize=5.8, ha="right", va="bottom", color=INK2, bbox=dict(facecolor=SURFACE, edgecolor="none", pad=0.3))
ax2.text(23.3 + 0.8, 71500, "the three splittings\ninto two pairs", fontsize=5.8, ha="left", va="bottom", color=INK2)
ax2.set_xlim(0, 76); ax2.set_ylim(0, 80000)
ax2.set_yticks([0, 20000, 40000, 60000, 80000])
ax2.set_yticklabels(["0", "20 000", "40 000", "60 000", "80 000"])
ax2.tick_params(labelsize=6.8)
ax2.set_xlabel("minutes of the run", fontsize=7.5)
ax2.set_ylabel("coefficients carried", fontsize=7.5)
ax2.legend(loc="lower right", fontsize=6.2, frameon=True, framealpha=1, edgecolor="none", handlelength=1.8, borderpad=0.3)
ax2.set_title("step 5, the four-point constraint", fontsize=7.8, loc="right", color=INK2)
panel_label(ax2, "(b)", x=-0.3, y=1.0)
save(fig, "fig_c_verify", outdir="../figures")
