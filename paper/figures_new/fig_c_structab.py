"""
fig_c_structab.py -- the root system of D4 and its dual 24-cell seen from a
deep hole, each in its own panel.

Directions on S^3 are drawn in the azimuthal equidistant chart about the hole
h = e_1: a unit vector u at angle a from h goes to a times the unit tangent of
the arc from h to u.  Edges join the pairs at 60 degrees, drawn as segments.
"""
import itertools
import os
import numpy as np
from PIL import Image, ImageDraw, ImageFont
import render3d as r3

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, os.pardir, "figures", "fig_c_structab.png")

h = np.array([1.0, 0, 0, 0])


def chart(u):
    u = u / np.linalg.norm(u)
    c = np.clip(u @ h, -1, 1)
    a = np.degrees(np.arccos(c))
    t = u - c * h
    n = np.linalg.norm(t)
    if n < 1e-12:
        return np.zeros(3)
    return a * t[1:] / n


roots = []
for i, j in itertools.combinations(range(4), 2):
    for si in (1, -1):
        for sj in (1, -1):
            v = np.zeros(4); v[i] = si; v[j] = sj
            roots.append(v / np.sqrt(2))
roots = np.array(roots)

dual = [s * np.eye(4)[i] for i in range(4) for s in (1, -1)]
dual += [np.array(s) / 2 for s in itertools.product((1, -1), repeat=4)]
dual = np.array(dual)


def pairs60(V):
    return [(i, j) for i, j in itertools.combinations(range(len(V)), 2) if abs(V[i] @ V[j] - 0.5) < 1e-9]


R = r3.rotation(np.radians(28), np.radians(-22), np.radians(4))
bounds = ((-146, -146), (146, 146))
size = (900, 900)

# (a) the 24 roots, coloured by angle from the hole
sa = r3.Scene()
col_of = {45: r3.CLAY, 90: r3.TEAL, 135: r3.BLUE}
P = [chart(u) for u in roots]
ang = [int(round(np.degrees(np.arccos(u @ h)))) for u in roots]
for i, j in pairs60(roots):
    same = ang[i] == ang[j]
    sa.tube(P[i], P[j], 1.6 if same else 0.8, col_of[ang[i]] if same else r3.GREY)
for p, a in zip(P, ang):
    sa.sphere(p, 4.2, col_of[a])
sa.sphere(np.zeros(3), 3.0, r3.DARK)
# translucent shells: faces of the two octahedra and of the cuboctahedron
for a in (45, 135):
    for sx, sy, sz in itertools.product((1, -1), repeat=3):
        sa.poly([(a * sx, 0, 0), (0, a * sy, 0), (0, 0, a * sz)], r3.FACE, 0.10)
ima = r3.render(sa, R, size, bounds=bounds)

# (b) the dual 24-cell, coloured by angle from the hole in the same order of
# shells as (a); the antipode of the hole (180 degrees) is left out
sb = r3.Scene()
keep = [k for k in range(len(dual)) if dual[k] @ h > -0.999]
Q = {k: chart(dual[k]) for k in keep}
dang = {k: int(round(np.degrees(np.arccos(np.clip(dual[k] @ h, -1, 1))))) for k in keep}
dcol = {0: r3.DARK, 60: r3.CLAY, 90: r3.TEAL, 120: r3.BLUE}
for i, j in pairs60(dual):
    if i in Q and j in Q:
        same = dang[i] == dang[j]
        sb.tube(Q[i], Q[j], 1.6 if same else 0.8, dcol[dang[i]] if same else r3.GREY)
for k in keep:
    sb.sphere(Q[k], 4.2, dcol[dang[k]])
for sx, sy, sz in itertools.product((1, -1), repeat=3):
    sb.poly([(90 * sx, 0, 0), (0, 90 * sy, 0), (0, 0, 90 * sz)], r3.FACE, 0.10)
imb = r3.render(sb, R, size, bounds=bounds)


def compose(panels, labels, gap=40, font_px=46):
    font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSerif-Bold.ttf", font_px)
    w = sum(p.width for p in panels) + gap * (len(panels) - 1)
    hgt = max(p.height for p in panels) + font_px + 10
    out = Image.new("RGB", (w, hgt), (255, 255, 255))
    d = ImageDraw.Draw(out)
    x = 0
    for p, lab in zip(panels, labels):
        out.paste(p, (x, font_px + 10))
        d.text((x + 4, 0), lab, fill=(11, 11, 11), font=font)
        x += p.width + gap
    return out


if __name__ == "__main__":
    fig = compose([ima, imb], ["(a)", "(b)"])
    fig.save(OUT, dpi=(400, 400))
    print("wrote", os.path.normpath(OUT), fig.size)
