"""
fig_c_boundary.py -- the facet of the reference cell and what the global
routes deliver.

Panel (a), the rendered facet inside the 3-space orthogonal to its contact
direction, is kept from the earlier composite.  Panel (b) is the bar chart of
the lower bounds for vol(V_c) that the boundary form and the other global
routes deliver at the D_4 configuration, where the conjecture is tight:
the true value 8; the contact-cap estimate 7.7990 (the covering bound at
m = 24); the convexity bound 7.7351 for the polar volume integral; the
Mahler-type product bound 16/3 = 5.3333; the inscribed ball, pi^2/2 =
4.9348; and the facet-by-facet bound 8 pi/(3 sqrt 3) = 4.8368.
"""
import numpy as np
import matplotlib.pyplot as plt
from figstyle import *
from oldpanels import load, split_cols, trim, show

old = load("old/fig_c_boundary.png")
left = trim(split_cols(old, 20)[0])

rows = [("true value at $D_4$", 8.0, BLUE),
        ("contact caps", 7.7990, INK3),
        ("convexity of $t^{-4}$", 7.7351, VIOLET),
        ("Mahler-type product", 16 / 3, ORANGE),
        ("inscribed ball $B^4$", np.pi ** 2 / 2, AQUA),
        ("facet by facet", 8 * np.pi / (3 * np.sqrt(3)), YELLOW)]
assert abs(rows[4][1] - 4.9348) < 1e-4 and abs(rows[5][1] - 4.8368) < 1e-4 and abs(rows[3][1] - 5.3333) < 1e-4

W_IN = 7.2
fig = plt.figure(figsize=(W_IN, 3.2))
gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.15], left=0.01, right=0.985, top=0.93, bottom=0.16, wspace=0.55)
axl = fig.add_subplot(gs[0]); show(axl, left)

ax = fig.add_subplot(gs[1])
for side in ("top", "right"):
    ax.spines[side].set_visible(False)
ax.grid(True, axis="x", zorder=0); ax.set_axisbelow(True)
y = np.arange(len(rows))[::-1]
for yy, (name, val, col) in zip(y, rows):
    ax.barh(yy, val, height=0.56, color=col, zorder=3)
    ax.text(val - 0.12, yy, f"{val:.4f}" if val != 8 else "8", va="center", ha="right", fontsize=6.8, color=SURFACE, zorder=5)
ax.axvline(8, color=INK2, lw=0.8, ls=(0, (4, 2)), zorder=4)
ax.set_yticks(y); ax.set_yticklabels([r[0] for r in rows], fontsize=7)
ax.set_xlim(0, 9.6); ax.set_xticks([0, 2, 4, 6, 8])
ax.tick_params(labelsize=6.8)
ax.set_xlabel(r"lower bound delivered for $\mathrm{vol}(V_c)$", fontsize=7.5)
ax.text(8.15, 5.55, "target 8", fontsize=6.5, color=INK2, ha="left", va="center")
panel_label(ax, "(b)", x=-0.38, y=1.0)
save(fig, "fig_c_boundary", outdir="../figures")
