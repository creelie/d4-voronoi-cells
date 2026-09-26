"""
fig_cardinality.py -- the two halves of the twenty-four-point theorem at positive slack
(Section 7.13).

(a) The three-point bound on A(4, 1/2 + s), the number of points of S^3 with inner
    products at most 1/2 + s, from the sampled programmes of degree 8 and 10
    (multi_cap/runs/cardinality_sweep_d*.log), with the certified values of the
    theorem that the kissing number is stable (multi_cap/runs/certify_cardinality_*.log)
    and the slack at which a 25-point code exists (multi_cap/runs/code25_search.log).
(b) The ceiling of the theorem 'Twenty-four points, approximately': with the triple
    and quadruple terms set to zero, the largest slack kappa its proof reaches,
    kappa < q_min f_min(delta) / (276 L_2), against the rounding tolerance delta
    (constants from multi_cap/runs/robust_ceiling.log).
"""
import glob
import os
import re

import numpy as np
import matplotlib.pyplot as plt
from figstyle import *

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.join(HERE, "..", "..", "multi_cap", "runs")


def sweep(d):
    pts = []
    for line in open(os.path.join(RUNS, "cardinality_sweep_d%d.log" % d)):
        m = re.match(r"RESULT d=\d+ t=\S+ s=(\S+)\s+sampled (\S+)", line)
        if m:
            pts.append((float(m.group(1)), float(m.group(2))))
    return sorted(pts)


def certified():
    out = []
    for fn in sorted(glob.glob(os.path.join(RUNS, "certify_cardinality_*.log"))):
        txt = open(fn).read()
        if "PASS" not in txt:
            continue
        s = float(re.search(r"PASS: A\(4, (\S+)\)", txt).group(1)) - 0.5
        b = float(re.search(r"B - 1 = (\S+)", txt).group(1)) + 1
        d = int(re.search(r"certificate degree (\d+)", txt).group(1))
        out.append((s, b, d))
    return sorted(out)


d8, d10, cert = sweep(8), sweep(10), certified()
s25 = float(re.search(r"at most 1/2 \+ (\S+)", open(os.path.join(RUNS, "code25_search.log")).read()).group(1))
cl = open(os.path.join(RUNS, "robust_ceiling.log")).read()
L2 = float(re.search(r"L2 = .*?= (\S+)", cl).group(1))
qmin = float(re.search(r"qmin = .*?= (\S+)", cl).group(1))

fig, (ax, bx) = plt.subplots(1, 2, figsize=(6.5, 2.75), gridspec_kw=dict(wspace=0.42))

# (a)
ax.plot(*zip(*d8), color=BLUE, marker="o", ms=3.2, label="degree 8, sampled")
ax.plot(*zip(*d10), color=ORANGE, marker="o", ms=3.2, label="degree 10, sampled")
ax.plot([c[0] for c in cert], [c[1] for c in cert], ls="none", marker="s", ms=3.0, color=AQUA,
        label="certified", zorder=5)
ax.axhline(25, color=INK3, lw=0.8, ls="--")
ax.axvline(s25, color=MAGENTA, lw=0.9, ls=":")
ax.text(s25 + 0.0025, 26.75, "25 points\nexist from\n$s=%.4f$" % s25, color=MAGENTA, fontsize=6.8,
        va="center", ha="left", bbox=dict(facecolor=SURFACE, edgecolor="none", pad=0.8))
ax.text(0.062, 25.2, "bound $=25$", color=INK3, fontsize=6.8, va="bottom",
        bbox=dict(facecolor=SURFACE, edgecolor="none", pad=0.8))
ax.text(0.0115, 24.3, "certified at $s=%s$" % ", ".join("%g" % c[0] for c in cert), color="#12805a",
        fontsize=6.4, va="center", ha="left", zorder=6, bbox=dict(facecolor=SURFACE, edgecolor="none", pad=0.6))
ax.set_xlim(-0.003, 0.084)
ax.set_ylim(23.6, 34.3)
ax.set_xticks([0, 0.02, 0.04, 0.06, 0.08])
ax.set_xlabel(r"slack $s$ (inner products at most $1/2+s$)", fontsize=7.5)
ax.set_ylabel(r"three-point bound on $A(4,1/2+s)$", fontsize=7.5)
ax.tick_params(labelsize=6.5)
ax.legend(loc="upper left", fontsize=6.6, frameon=True, facecolor=SURFACE, edgecolor="none",
          framealpha=1.0, handlelength=1.6, borderpad=0.3)
clean_axes(ax)
panel_label(ax, "(a)", x=-0.2, y=1.0)

# (b)
f = lambda u: (u + 1) * (u + 0.5) ** 2 * u ** 2 * (0.5 - u)
fmin = lambda d: min(f(-1 + d), f(-0.5 - d), f(-0.5 + d), f(-d), f(d), f(0.5 - d))
kap = lambda d: qmin * fmin(d) / (276 * L2)
ds = np.logspace(-5, np.log10(0.25), 300)
bx.loglog(ds, [kap(d) for d in ds], color=BLUE, label=r"ceiling, $B_3=B_4=0$")
for d in (3e-4, 0.25):
    bx.plot([d], [kap(d)], marker="o", ms=3.8, color=ORANGE, zorder=4)
bx.text(4.2e-4, 1.2e-14, r"$\delta=3\cdot10^{-4}$: $%.1f\cdot10^{-14}$" % (kap(3e-4) * 1e14),
        fontsize=6.6, color=ORANGE, va="top", ha="left")
bx.text(0.17, 1.2e-8, r"$\delta=1/4$: $%.1f\cdot10^{-9}$" % (kap(0.25) * 1e9), fontsize=6.6,
        color=ORANGE, va="bottom", ha="right")
bx.axhline(2e-26, color=INK3, lw=0.8, ls="--")
bx.text(1.4e-5, 2e-25, r"the proof as it stands: $2\cdot10^{-26}$", fontsize=6.6, color=INK2, va="bottom",
        bbox=dict(facecolor=SURFACE, edgecolor="none", pad=0.8))
bx.axhline(1e-2, color=MAGENTA, lw=0.9, ls=":")
bx.text(0.3, 2.5e-2, r"wanted: $10^{-2}$", fontsize=6.6, color=MAGENTA, va="bottom", ha="right",
        bbox=dict(facecolor=SURFACE, edgecolor="none", pad=0.8))
bx.set_xlim(1e-5, 0.4)
bx.set_ylim(1e-28, 1)
bx.set_yticks([1e-28, 1e-24, 1e-20, 1e-16, 1e-12, 1e-8, 1e-4, 1])
bx.set_xlabel(r"rounding tolerance $\delta$", fontsize=7.5)
bx.set_ylabel(r"largest slack $\kappa$ reached", fontsize=7.5)
bx.tick_params(labelsize=6.5)
bx.legend(loc="upper left", bbox_to_anchor=(0.0, 0.875), fontsize=6.6, frameon=True, facecolor=SURFACE,
          edgecolor="none", framealpha=1.0, handlelength=1.6, borderpad=0.3)
clean_axes(bx)
panel_label(bx, "(b)", x=-0.2, y=1.0)

save(fig, "fig_cardinality", outdir=os.path.join(HERE, "..", "figures"))
