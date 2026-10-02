"""
render3d.py -- a small ray-cast renderer for the wireframe panels of the paper.

A scene is a list of opaque spheres and tubes (capsules) and of translucent
planar polygons, seen in orthographic projection.  Opaque primitives are
z-buffered pixel by pixel with Phong shading; the polygons are then laid over
them from back to front wherever they lie in front of the opaque surface.  The
image is computed at SS times the output resolution and averaged down.
"""
import numpy as np
from PIL import Image

TEAL = (0.07, 0.42, 0.42)
BLUE = (0.12, 0.27, 0.52)
CLAY = (0.66, 0.26, 0.18)
OCHRE = (0.76, 0.52, 0.13)
VIOLET = (0.36, 0.26, 0.55)
GREY = (0.62, 0.66, 0.70)
DARK = (0.10, 0.12, 0.16)
FACE = (0.80, 0.84, 0.88)

LIGHT = np.array([-0.45, 0.55, 0.70])
LIGHT = LIGHT / np.linalg.norm(LIGHT)
HALF = LIGHT + np.array([0.0, 0.0, 1.0])
HALF = HALF / np.linalg.norm(HALF)


def rotation(yaw, pitch, roll=0.0):
    """rotation taking world coordinates to view coordinates (x right, y up, z to the viewer)"""
    cy, sy = np.cos(yaw), np.sin(yaw)
    cp, sp = np.cos(pitch), np.sin(pitch)
    cr, sr = np.cos(roll), np.sin(roll)
    Ry = np.array([[cy, 0, sy], [0, 1, 0], [-sy, 0, cy]])
    Rx = np.array([[1, 0, 0], [0, cp, -sp], [0, sp, cp]])
    Rz = np.array([[cr, -sr, 0], [sr, cr, 0], [0, 0, 1]])
    return Rz @ Rx @ Ry


class Scene:
    def __init__(self):
        self.spheres = []   # (centre, radius, colour)
        self.tubes = []     # (a, b, radius, colour)
        self.polys = []     # (vertices, colour, alpha)

    def sphere(self, c, r, col):
        self.spheres.append((np.asarray(c, float), float(r), np.asarray(col, float)))

    def tube(self, a, b, r, col):
        self.tubes.append((np.asarray(a, float), np.asarray(b, float), float(r), np.asarray(col, float)))

    def polyline(self, pts, r, col, closed=False):
        pts = np.asarray(pts, float)
        n = len(pts)
        for i in range(n if closed else n - 1):
            self.tube(pts[i], pts[(i + 1) % n], r, col)
        for p in pts:
            self.sphere(p, r, col)

    def poly(self, verts, col, alpha):
        self.polys.append((np.asarray(verts, float), np.asarray(col, float), float(alpha)))


def _shade(col, nx, ny, nz, amb=0.38, dif=0.62, spec=0.45, shin=40.0):
    lam = np.clip(nx * LIGHT[0] + ny * LIGHT[1] + nz * LIGHT[2], 0, 1)
    h = np.clip(nx * HALF[0] + ny * HALF[1] + nz * HALF[2], 0, 1) ** shin
    out = col[None, :] * (amb + dif * lam)[:, None] + spec * h[:, None]
    return np.clip(out, 0, 1)


def render(scene, R, size, ss=3, margin=0.04, bounds=None, background=(1, 1, 1)):
    """render the scene in the view R; size=(W, H) in output pixels"""
    W, H = size[0] * ss, size[1] * ss
    tf = lambda p: R @ p
    sph = [(tf(c), r, col) for c, r, col in scene.spheres]
    tub = [(tf(a), tf(b), r, col) for a, b, r, col in scene.tubes]
    pol = [(np.array([tf(v) for v in vs]), col, al) for vs, col, al in scene.polys]
    if bounds is None:
        pts = [c for c, r, _ in sph] + [a for a, _, _, _ in tub] + [b for _, b, _, _ in tub]
        pts += [v for vs, _, _ in pol for v in vs]
        pts = np.array(pts)
        rad = max([r for _, r, _ in sph] + [r for *_, r, _ in tub] + [0])
        lo = pts[:, :2].min(0) - rad
        hi = pts[:, :2].max(0) + rad
    else:
        lo, hi = np.array(bounds[0], float), np.array(bounds[1], float)
    span = hi - lo
    scale = min(W * (1 - 2 * margin) / span[0], H * (1 - 2 * margin) / span[1])
    off = np.array([W, H]) / 2 - scale * (lo + hi) / 2

    def px(p):
        return p[:2] * scale + off

    img = np.empty((H, W, 3))
    img[:] = background
    zb = np.full((H, W), -np.inf)

    def window(cx, cy, rr):
        x0 = max(int(np.floor(cx - rr)), 0); x1 = min(int(np.ceil(cx + rr)) + 1, W)
        y0 = max(int(np.floor(cy - rr)), 0); y1 = min(int(np.ceil(cy + rr)) + 1, H)
        return x0, x1, y0, y1

    for c, r, col in sph:
        cx, cy = px(c); rp = r * scale
        x0, x1, y0, y1 = window(cx, cy, rp)
        if x0 >= x1 or y0 >= y1:
            continue
        X, Y = np.meshgrid(np.arange(x0, x1) + 0.5, np.arange(y0, y1) + 0.5)
        dx, dy = (X - cx) / rp, (Y - cy) / rp
        d2 = dx * dx + dy * dy
        m = d2 < 1
        if not m.any():
            continue
        nz = np.sqrt(np.clip(1 - d2, 0, 1))
        z = c[2] + r * nz
        sub = zb[y0:y1, x0:x1]
        w = m & (z > sub)
        if not w.any():
            continue
        sub[w] = z[w]
        img[y0:y1, x0:x1][w] = _shade(col, dx[w], dy[w], nz[w])

    for a, b, r, col in tub:
        pa, pb = px(a), px(b); rp = r * scale
        x0 = max(int(np.floor(min(pa[0], pb[0]) - rp)), 0); x1 = min(int(np.ceil(max(pa[0], pb[0]) + rp)) + 1, W)
        y0 = max(int(np.floor(min(pa[1], pb[1]) - rp)), 0); y1 = min(int(np.ceil(max(pa[1], pb[1]) + rp)) + 1, H)
        if x0 >= x1 or y0 >= y1:
            continue
        X, Y = np.meshgrid(np.arange(x0, x1) + 0.5, np.arange(y0, y1) + 0.5)
        d = pb - pa
        L2 = d @ d
        s = np.zeros_like(X) if L2 < 1e-12 else np.clip(((X - pa[0]) * d[0] + (Y - pa[1]) * d[1]) / L2, 0, 1)
        cx = pa[0] + s * d[0]; cy = pa[1] + s * d[1]
        cz = a[2] + s * (b[2] - a[2])
        dx, dy = (X - cx) / rp, (Y - cy) / rp
        d2 = dx * dx + dy * dy
        m = d2 < 1
        if not m.any():
            continue
        nz = np.sqrt(np.clip(1 - d2, 0, 1))
        z = cz + r * nz
        sub = zb[y0:y1, x0:x1]
        w = m & (z > sub)
        if not w.any():
            continue
        sub[w] = z[w]
        img[y0:y1, x0:x1][w] = _shade(col, dx[w], dy[w], nz[w], spec=0.35)

    order = sorted(range(len(pol)), key=lambda i: pol[i][0][:, 2].mean())
    for i in order:
        vs, col, al = pol[i]
        P = np.array([px(v) for v in vs])
        n3 = np.cross(vs[1] - vs[0], vs[2] - vs[0])
        if abs(n3[2]) < 1e-9:
            continue
        n3 = n3 / np.linalg.norm(n3)
        lam = abs(n3 @ LIGHT)
        fcol = np.clip(col * (0.82 + 0.25 * lam), 0, 1)
        x0 = max(int(np.floor(P[:, 0].min())), 0); x1 = min(int(np.ceil(P[:, 0].max())) + 1, W)
        y0 = max(int(np.floor(P[:, 1].min())), 0); y1 = min(int(np.ceil(P[:, 1].max())) + 1, H)
        if x0 >= x1 or y0 >= y1:
            continue
        X, Y = np.meshgrid(np.arange(x0, x1) + 0.5, np.arange(y0, y1) + 0.5)
        sgn = np.sign((P[1, 0] - P[0, 0]) * (P[2, 1] - P[0, 1]) - (P[1, 1] - P[0, 1]) * (P[2, 0] - P[0, 0]))
        inside = np.ones_like(X, bool)
        k = len(P)
        for j in range(k):
            e = P[(j + 1) % k] - P[j]
            inside &= sgn * (e[0] * (Y - P[j, 1]) - e[1] * (X - P[j, 0])) >= 0
        # depth on the plane, in world units
        wx = (X - off[0]) / scale; wy = (Y - off[1]) / scale
        z = vs[0, 2] - (n3[0] * (wx - vs[0, 0]) + n3[1] * (wy - vs[0, 1])) / n3[2]
        w = inside & (z > zb[y0:y1, x0:x1] - 1e-9)
        if not w.any():
            continue
        blk = img[y0:y1, x0:x1]
        blk[w] = (1 - al) * blk[w] + al * fcol

    img = img[::-1]  # y up
    img = img.reshape(size[1], ss, size[0], ss, 3).mean(axis=(1, 3))
    return Image.fromarray((np.clip(img, 0, 1) * 255 + 0.5).astype(np.uint8))


def trim(im, pad=6, level=250):
    a = np.asarray(im.convert("L"))
    rows = np.where((a < level).any(axis=1))[0]
    cols = np.where((a < level).any(axis=0))[0]
    if len(rows) == 0:
        return im
    return im.crop((max(cols[0] - pad, 0), max(rows[0] - pad, 0),
                    min(cols[-1] + pad + 1, im.width), min(rows[-1] + pad + 1, im.height)))
