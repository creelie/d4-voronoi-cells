"""
fig_c_structcd.py -- (a) the 600-cell seen from one of its vertices, with an
inscribed 24-cell through that vertex picked out; (b) the cap lemma on a great
2-sphere, a ray-traced panel kept from the earlier composite.

Panel (a) uses the azimuthal equidistant chart about the vertex p = e_1, as in
fig_c_structab: a unit vector at angle a from p goes to a times the unit
tangent of the arc from p to it.  The antipode of p, at 180 degrees, has no
single image and is left out, so the 24-cell shows 23 of its vertices.
"""
import itertools
import os
import numpy as np
from PIL import Image, ImageDraw, ImageFont
import render3d as r3

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, os.pardir, "figures", "fig_c_structcd.png")
OLD = os.path.join(HERE, "old", "fig_c_structcd.png")

phi = (1 + np.sqrt(5)) / 2


def even_perms():
    out = []
    for p in itertools.permutations(range(4)):
        inv = sum(1 for i in range(4) for j in range(i + 1, 4) if p[i] > p[j])
        if inv % 2 == 0:
            out.append(p)
    return out


hurwitz = [s * np.eye(4)[i] for i in range(4) for s in (1, -1)]
hurwitz += [np.array(s) / 2 for s in itertools.product((1, -1), repeat=4)]
icos = []
for p in even_perms():
    for s in itertools.product((1, -1), repeat=3):
        base = np.array([s[0] * phi, s[1] * 1.0, s[2] / phi, 0.0]) / 2
        v = np.zeros(4)
        for k in range(4):
            v[p[k]] = base[k]
        icos.append(v)
V = np.array(hurwitz + icos)
V = np.unique(np.round(V, 12), axis=0)
assert len(V) == 120
is24 = np.array([any(np.allclose(v, w) for w in hurwitz) for v in V])
pole = np.array([1.0, 0, 0, 0])


def chart(x):
    c = np.clip(x @ pole, -1, 1)
    t = x - c * pole
    n = np.linalg.norm(t)
    return np.zeros(3) if n < 1e-12 else np.degrees(np.arccos(c)) * t[1:] / n


s = r3.Scene()
keep = V @ pole > -0.999
for v, h24 in zip(V[keep], is24[keep]):
    s.sphere(chart(v), 4.4 if h24 else 2.6, r3.DARK if h24 else r3.GREY)
H = V[is24 & keep]
for i, j in itertools.combinations(range(len(H)), 2):
    if abs(H[i] @ H[j] - 0.5) < 1e-9:
        s.tube(chart(H[i]), chart(H[j]), 0.9, r3.DARK)
R = r3.rotation(np.radians(28), np.radians(-22), np.radians(4))
ima = r3.trim(r3.render(s, R, (980, 980)))


def old_panel_b():
    a = np.asarray(Image.open(OLD).convert("RGB")).astype(float) / 255
    white = (a.mean(axis=2) > 0.985).all(axis=0)
    # the widest run of white columns in the middle third separates the panels
    n = len(white); best = (0, 0)
    i = n // 3
    while i < 2 * n // 3:
        if white[i]:
            j = i
            while j < n and white[j]:
                j += 1
            if j - i > best[1] - best[0]:
                best = (i, j)
            i = j
        else:
            i += 1
    cut = (best[0] + best[1]) // 2
    im = Image.open(OLD).convert("RGB").crop((cut, 0, n, a.shape[0]))
    # remove the old panel letter in the top-left corner
    d = ImageDraw.Draw(im)
    d.rectangle((0, 0, 140, 110), fill=(255, 255, 255))
    return r3.trim(im)


def compose(panels, labels, gap=60, font_px=46):
    font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSerif-Bold.ttf", font_px)
    hgt = max(p.height for p in panels)
    panels = [p if p.height == hgt else p.resize((round(p.width * hgt / p.height), hgt), Image.LANCZOS) for p in panels]
    w = sum(p.width for p in panels) + gap * (len(panels) - 1)
    out = Image.new("RGB", (w, hgt + font_px + 10), (255, 255, 255))
    d = ImageDraw.Draw(out)
    x = 0
    for p, lab in zip(panels, labels):
        out.paste(p, (x, font_px + 10))
        d.text((x + 4, 0), lab, fill=(11, 11, 11), font=font)
        x += p.width + gap
    return out


if __name__ == "__main__":
    fig = compose([ima, old_panel_b()], ["(a)", "(b)"])
    fig.save(OUT, dpi=(400, 400))
    print("wrote", os.path.normpath(OUT), fig.size)
