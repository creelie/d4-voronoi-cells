"""
fig_c_sigma2.py -- the two-point polynomial of the certificate and the
structure of the certificate kernel.

(a) sigma_2(u) on [-1, 1/2], evaluated exactly from the deposited rational
    coefficients and scaled by its maximum on the interval: it is nonnegative,
    it vanishes at -1 and 1/2 to first order and at -1/2 and 0 to second order,
    and nowhere else.  The insets magnify the two double zeros, where the
    curve touches the axis without crossing.
(b) The same on a logarithmic scale, which shows the four zeros as dips of
    the four orders.
(c) The 127 positive definite blocks of the certificate by size, split into
    the sixty kernel blocks X_lambda, one per signature, and the sixty-seven
    sum-of-squares multipliers of the constraints on two, three and four
    points.
"""
import pickle
import sys
from fractions import Fraction

import numpy as np
import matplotlib.pyplot as plt
from figstyle import *

sys.set_int_max_str_digits(2_000_000)

d = pickle.load(open("p2.pkl", "rb"))
N = max(d)
coef = [Fraction(*d[i]) for i in range(N + 1)]


def p2(x):
    x = Fraction(x)
    acc = Fraction(0)
    for c in reversed(coef):
        acc = acc * x + c
    return acc


# exact evaluation on a rational grid, then one normalisation
xs = [Fraction(-1) + Fraction(3, 2) * Fraction(i, 1200) for i in range(1201)]
vals = [p2(x) for x in xs]
vmax = max(vals)
y = np.array([float(v / vmax) for v in vals])
x = np.array([float(t) for t in xs])

fig = plt.figure(figsize=(7.2, 5.4))
gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.85], hspace=0.42, wspace=0.28,
                      left=0.09, right=0.98, top=0.95, bottom=0.09)
ax = fig.add_subplot(gs[0, :])
ax2 = fig.add_subplot(gs[1, 0])
ax3 = fig.add_subplot(gs[1, 1])
for a in (ax, ax2, ax3):
    clean_axes(a)

ax.plot(x, y, color=BLUE, lw=1.4, zorder=3)
ax.axhline(0, color=INK2, lw=0.7)
zeros = [(-1.0, 1), (-0.5, 2), (0.0, 2), (0.5, 1)]
for z, m in zeros:
    ax.scatter([z], [0], s=26, color=ORANGE, zorder=5)
    ax.text(z, 0.06, f"zero of order {m}", ha="center", va="bottom", fontsize=6.6, color=ORANGE,
            bbox=dict(facecolor=SURFACE, edgecolor="none", pad=0.6))
ax.set_xticks([-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5])
ax.set_xticklabels([r"$-1$", r"$-\frac{3}{4}$", r"$-\frac{1}{2}$", r"$-\frac{1}{4}$", r"$0$", r"$\frac{1}{4}$", r"$\frac{1}{2}$"])
ax.set_xlim(-1.03, 0.53)
ax.set_ylim(-0.06, 1.08)
ax.set_xlabel(r"inner product $u=\langle x,y\rangle$")
ax.set_ylabel(r"$\sigma_2(u)$, scaled to maximum $1$")
ax.text(-0.5, 0.98, r"$\sigma_2 \geq 0$ on $[-1,\frac{1}{2}]$; degree $16$; exact rational coefficients of about $15\,700$ digits",
        fontsize=7, color=INK2, va="top", ha="center", bbox=dict(facecolor=SURFACE, edgecolor="none", pad=1.0))
panel_label(ax, "(a)", x=-0.06)

# insets at the double zeros, placed over the flat part of the curve
for z, (ix0, iy0) in zip((-0.5, 0.0), ((0.12, 0.40), (0.42, 0.40))):
    ins = ax.inset_axes([ix0, iy0, 0.22, 0.34])
    sel = np.abs(x - z) < 0.06
    ins.plot(x[sel], y[sel], color=BLUE, lw=1.2)
    ins.axhline(0, color=INK2, lw=0.6)
    ins.scatter([z], [0], s=14, color=ORANGE, zorder=5)
    ins.set_xticks([z - 0.05, z, z + 0.05])
    ins.set_xticklabels([f"{z-0.05:.2f}", f"{z:g}", f"{z+0.05:.2f}"], fontsize=5.5)
    ins.set_yticks([])
    ins.tick_params(length=2, pad=1)
    for sp in ("top", "right"):
        ins.spines[sp].set_visible(False)
    ins.set_title(f"near $u={z:g}$: touches, does not cross", fontsize=5.8, color=INK2, pad=2)

# (b) logarithmic
pos = y > 0
ax2.semilogy(x[pos], y[pos], color=BLUE, lw=1.2)
for z, m in zeros:
    ax2.axvline(z, color=ORANGE, lw=0.7, ls=(0, (3, 2)))
ax2.set_ylim(1e-12, 3)
ax2.set_xlim(-1.03, 0.53)
ax2.set_xlabel(r"$u$")
ax2.set_ylabel(r"$\sigma_2(u)$ (log scale)")
ax2.text(-0.75, 3e-11, "dips at the four zeros;\nthe grid never lands on one", fontsize=6.4, color=INK2,
         bbox=dict(facecolor=SURFACE, edgecolor="none", pad=0.8))
panel_label(ax2, "(b)", x=-0.2)

# (c) block sizes of the certificate: from the verification log
# the sixty kernel blocks are indexed by signatures |lambda| <= 14; the largest is 350 x 350
# distribution read from llm24_certificate_check.log (blocks and their orders)
import re, os
log = None
for cand in ("llm24_certificate_check.log", "../../multi_cap/runs/llm24_certificate_check.log"):
    if os.path.exists(cand):
        log = open(cand).read()
        break
orders_k, orders_s = [], []
if log:
    for m in re.finditer(r"block \(([^)]*)\)\s+(\d+) x (\d+)", log):
        kind, n = m.group(1).strip("'"), int(m.group(2))
        (orders_s if kind.startswith("sos") else orders_k).append(n)
assert len(orders_k) == 60 and len(orders_s) == 67, (len(orders_k), len(orders_s))
ok_ = sorted(orders_k); os_ = sorted(orders_s)
ax3.plot(range(1, len(ok_) + 1), ok_, color=BLUE, marker="o", ms=2.4, lw=0.9, label=r"kernel blocks $X_\lambda$, one per signature (60)")
ax3.plot(range(1, len(os_) + 1), os_, color=ORANGE, marker="s", ms=2.4, lw=0.9, label="sum-of-squares multipliers (67)")
ax3.set_yscale("log")
ax3.set_xlabel("blocks, sorted by order")
ax3.set_ylabel("order of the block")
ax3.set_ylim(0.7, 700)
ax3.legend(fontsize=6.2, frameon=True, framealpha=1, edgecolor="none", loc="upper left")
ax3.text(2, 130, f"127 blocks of total order {sum(orders_k)+sum(orders_s)},\nthe largest $350\\times350$; every one\npositive definite by a Cholesky\nfactorisation in ball arithmetic",
         fontsize=6.2, color=INK2, va="top", bbox=dict(facecolor=SURFACE, edgecolor="none", pad=0.8))
panel_label(ax3, "(c)", x=-0.2)
save(fig, "fig_c_sigma2", outdir="../figures", vector=True)
