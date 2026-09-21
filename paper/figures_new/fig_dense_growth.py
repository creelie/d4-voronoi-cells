"""
fig_dense_growth.py -- the Richardson-extrapolated minimum eigenvalue of the
joint Hessian H_A under greedy dense growth, m = 9 to 23, the values of the
table in the appendix of numerical tables.  From m = 18 onward the
extrapolated value is below the double-precision noise floor and is plotted
at zero, in the shaded band.
"""
import numpy as np
import matplotlib.pyplot as plt
from figstyle import *

m = np.array([9, 12, 15, 18, 20, 22, 23])
lam = np.array([0.083342, 0.073070, 0.044694, -0.000013, 0.000017, 0.000023, 0.000008])
floor = 1e-4

fig, ax = plt.subplots(figsize=(6.0, 3.4))
fig.subplots_adjust(left=0.13, right=0.98, top=0.95, bottom=0.16)
clean_axes(ax)
ax.axvspan(17.2, 24, color=PANEL, zorder=0)
ax.axhline(0, color=ORANGE, lw=0.8, ls=(0, (4, 2)), zorder=2)
ax.plot(m, lam, "-o", color=BLUE, lw=1.4, ms=4.5, zorder=4)
ax.set_xlim(8, 24); ax.set_ylim(-0.006, 0.095)
ax.set_xticks(m)
ax.set_yticks([0, 0.02, 0.04, 0.06, 0.08])
ax.set_xlabel(r"active-set size $m$ (greedy dense growth)")
ax.set_ylabel(r"$\lambda_{\min}(H_A)$, Richardson extrapolated")
ax.text(9.2, 0.0865, r"$m=9$: the star configuration, $\lambda_{\min}\approx\frac{1}{12}$", fontsize=7, color=BLUE, ha="left", va="bottom")
ax.text(13.6, 0.062, "falling as $m$ grows", fontsize=7, color=BLUE, ha="left", va="bottom")
ax.text(20.5, 0.03, "$m=18$ onward: an exact null\neigenvalue, below the noise\nfloor of double precision;\nat $m=18$ the near-null space\nis four dimensional", fontsize=6.6, color=INK2, ha="center", va="center", linespacing=1.3)
ax.text(23.8, 0.004, r"$\lambda=0$", fontsize=7, color=ORANGE, ha="right", va="bottom")
save(fig, "fig_dense_growth", outdir="../figures")
