"""
fig_c_cell24.py -- the reference 24-cell drawn from its coordinates.

Panels (a), (b) and (c) are the ray-traced panels of the earlier composite
(the orthogonal projection, the Schlegel diagram and the cap inequality one
dimension down), cut out along their white gutters and placed unchanged.
Panel (d) draws the three 16-cells that triality permutes, each on its own,
in one and the same orthogonal projection of R^4 to R^3: the axis family and
the two halves of the sign family split by the parity of the minus signs.
Each is drawn with its 24 edges and the translucent hull of its image.
"""
import itertools
import os
import numpy as np
from scipy.spatial import ConvexHull
from PIL import Image, ImageDraw, ImageFont
import render3d as r3

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, os.pardir, "figures", "fig_c_cell24.png")
OLD = os.path.join(HERE, "old", "fig_c_cell24.png")

axis = np.array([s * np.eye(4)[i] for i in range(4) for s in (1, -1)])
signs = np.array(list(itertools.product((1, -1), repeat=4))) / 2.0
even = signs[(signs < 0).sum(axis=1) % 2 == 0]
odd = signs[(signs < 0).sum(axis=1) % 2 == 1]
# the sign family scaled to the norm of the axis family, so that all three
# 16-cells are inscribed in one sphere, as they are in the 24-cell
families = [(axis, r3.TEAL), (even * 2 / np.linalg.norm(even[0]) / 2, r3.OCHRE),
            (odd * 2 / np.linalg.norm(odd[0]) / 2, r3.VIOLET)]

# a fixed orthogonal projection R^4 -> R^3, generic for all three families
M = np.array([[0.93, 0.21, -0.27, 0.12],
              [0.17, 0.88, 0.31, -0.33],
              [-0.25, 0.30, 0.79, 0.47]])
Qm, _ = np.linalg.qr(M.T)
P = Qm.T  # rows orthonormal
R = r3.rotation(np.radians(20), np.radians(-14), np.radians(3))


def sixteen_cell(V, col):
    s = r3.Scene()
    X = V @ P.T
    for i, j in itertools.combinations(range(len(V)), 2):
        if abs(V[i] @ V[j] + 1) > 1e-9:  # every pair but the antipodal ones
            s.tube(X[i], X[j], 0.016, col)
    for x in X:
        s.sphere(x, 0.042, col)
    hull = ConvexHull(X)
    for f in hull.simplices:
        s.poly(X[f], col, 0.10)
    return s


def old_quadrants():
    im = Image.open(OLD).convert("RGB")
    a = np.asarray(im).astype(float) / 255
    white = a.mean(axis=2) > 0.985

    def widest_gap(line, lo, hi):
        best = (0, 0); i = lo
        while i < hi:
            if line[i]:
                j = i
                while j < hi and line[j]:
                    j += 1
                if j - i > best[1] - best[0]:
                    best = (i, j)
                i = j
            else:
                i += 1
        return (best[0] + best[1]) // 2

    H, W = white.shape
    rc = widest_gap(white.all(axis=1), H // 3, 2 * H // 3)
    top, bot = im.crop((0, 0, W, rc)), im.crop((0, rc, W, H))
    out = []
    for half in (top, bot):
        h = np.asarray(half).astype(float).mean(axis=2) / 255 > 0.985
        cc = widest_gap(h.all(axis=0), W // 3, 2 * W // 3)
        for box in ((0, 0, cc, half.height), (cc, 0, W, half.height)):
            q = half.crop(box)
            d = ImageDraw.Draw(q)
            d.rectangle((0, 0, 220, 200), fill=(255, 255, 255))  # old panel letter
            out.append(r3.trim(q))
    return out  # (a), (b), (c), (d) of the earlier composite


def row(panels, height, gap):
    panels = [p.resize((round(p.width * height / p.height), height), Image.LANCZOS) for p in panels]
    w = sum(p.width for p in panels) + gap * (len(panels) - 1)
    out = Image.new("RGB", (w, height), (255, 255, 255))
    x = 0
    starts = []
    for p in panels:
        out.paste(p, (x, 0)); starts.append(x)
        x += p.width + gap
    return out, starts


if __name__ == "__main__":
    qa, qb, _, qd = old_quadrants()
    cells = [r3.trim(r3.render(sixteen_cell(V, col), R, (800, 800), bounds=((-1.08, -1.08), (1.08, 1.08))))
             for V, col in families]
    font_px = 46
    font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSerif-Bold.ttf", font_px)
    r1, s1 = row([qa, qb, qd], 620, 50)
    r2, s2 = row(cells, 480, 110)
    W = max(r1.width, r2.width)
    lab = font_px + 12
    fig = Image.new("RGB", (W, lab + r1.height + 40 + lab + r2.height), (255, 255, 255))
    d = ImageDraw.Draw(fig)
    x1 = (W - r1.width) // 2
    fig.paste(r1, (x1, lab))
    for s, t in zip(s1, ["(a)", "(b)", "(c)"]):
        d.text((x1 + s + 2, 0), t, fill=(11, 11, 11), font=font)
    y2 = lab + r1.height + 40
    x2 = (W - r2.width) // 2
    fig.paste(r2, (x2, y2 + lab))
    d.text((x2 + 2, y2), "(d)", fill=(11, 11, 11), font=font)
    fig.save(OUT, dpi=(400, 400))
    print("wrote", os.path.normpath(OUT), fig.size)
