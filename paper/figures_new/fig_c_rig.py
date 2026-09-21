"""
fig_c_rig.py -- the rigidity of the root configuration.

Panels (a) and (b) are the rendered pictures of the neighbour sum in three
and in four dimensions, kept from the earlier composite.  Panels (c) and (d)
are drawn here from the statements they illustrate.

(c) The singular values of the rigidity operator Lambda on the 72-dimensional
    tangent space of (S^3)^24 at the root configuration, with multiplicities
    (the lemma on the rigidity spectrum): 0 with multiplicity 6, the
    infinitesimal rotations; 1 with multiplicity 29; sqrt(5/2) with 8;
    sqrt 3 with 21; 2 with 8.  The multiplicities add to 72.
(d) The inequality |e|^2 <= (577/4)|e|^4 that feasibility forces on the
    displacement e from the orbit of the root system: its two sides cross at
    |e| = 2/sqrt(577), and below that no contact configuration of 24
    directions exists.
"""
import numpy as np
import matplotlib.pyplot as plt
from figstyle import *
from oldpanels import load, split_rows, trim, show

old = load("old/fig_c_rig.png")
top = trim(split_rows(old, 20)[0])
W_IN = 4.55
top_h = W_IN * top.shape[0] / top.shape[1]

fig = plt.figure(figsize=(W_IN, top_h + 2.35))
gs = fig.add_gridspec(2, 2, height_ratios=[top_h, 2.0], hspace=0.28, wspace=0.55,
                      left=0.12, right=0.98, top=0.99, bottom=0.10)
axt = fig.add_subplot(gs[0, :])
show(axt, top)

# (c) the spectrum
ax = fig.add_subplot(gs[1, 0]); clean_axes(ax)
sv = [0.0, 1.0, np.sqrt(2.5), np.sqrt(3.0), 2.0]
mult = [6, 29, 8, 21, 8]
assert sum(mult) == 72
cols = [INK3, ORANGE, BLUE, BLUE, BLUE]
ax.bar(range(5), mult, color=cols, width=0.6, zorder=3)
for i, m in enumerate(mult):
    ax.text(i, m + 1.2, str(m), ha="center", va="bottom", fontsize=7, color=INK2,
            bbox=dict(facecolor=SURFACE, edgecolor="none", pad=0.3))
ax.set_xticks(range(5))
ax.set_xticklabels([r"$0$", r"$1$", r"$\sqrt{5/2}$", r"$\sqrt{3}$", r"$2$"], fontsize=7)
ax.set_ylim(0, 44)
ax.set_yticks([0, 10, 20, 30, 40])
ax.set_xlabel(r"singular value of $\Lambda$", fontsize=7.5)
ax.set_ylabel("multiplicity", fontsize=7.5)
ax.tick_params(labelsize=6.5)
ax.text(1.55, 43, "grey: the kernel,\nthe six rotations;\norange: the least\nnonzero value, exactly 1",
        ha="left", va="top", fontsize=5.8, color=INK2, linespacing=1.25,
        bbox=dict(facecolor=SURFACE, edgecolor="none", pad=0.4))
panel_label(ax, "(c)", x=-0.3, y=1.0)

# (d) the two sides of the inequality
ax2 = fig.add_subplot(gs[1, 1]); clean_axes(ax2)
e = np.linspace(0, 0.135, 400)
e0 = 2 / np.sqrt(577)
ax2.fill_between([0, e0], 0, 0.02, color=PANEL, zorder=0)
ax2.plot(e, e**2, color=BLUE, lw=1.3, label=r"$\|\varepsilon\|^2$")
ax2.plot(e, 577 / 4 * e**4, color=ORANGE, lw=1.3, label=r"$\frac{577}{4}\|\varepsilon\|^4$")
ax2.axvline(e0, color=INK2, lw=0.7, ls=(0, (3, 2)))
ax2.scatter([e0], [e0**2], s=16, color=INK, zorder=5)
ax2.set_xlim(0, 0.135); ax2.set_ylim(0, 0.019)
ax2.set_xticks([0, 0.04, 0.08, 0.12])
ax2.set_yticks([0, 0.005, 0.010, 0.015])
ax2.tick_params(labelsize=6.5)
ax2.set_xlabel(r"displacement $\|\varepsilon\|=d(W)$", fontsize=7.5)
ax2.set_ylabel("the two sides", fontsize=7.5)
ax2.text(0.041, 0.0135, "no contact\nconfiguration\nin this band", ha="center", va="center", fontsize=6, color=INK2)
ax2.text(e0 - 0.003, 0.0178, r"$2/\sqrt{577}$", ha="right", va="center", fontsize=6.5, color=INK2)
ax2.legend(loc="lower right", fontsize=6.5, frameon=True, framealpha=1, edgecolor="none", handlelength=1.6, borderpad=0.3)
panel_label(ax2, "(d)", x=-0.34, y=1.0)
save(fig, "fig_c_rig", outdir="../figures")
