#!/usr/bin/env python3
"""
tikz_overlap_check.py -- the collision test of the TikZ figures.

A figure that uses it names every text node and lists the names in \\boxes
(nodes drawn as boxes, which may carry arrows to their borders) and \\labels
(free text).  With \\def\\checkbboxes{} it writes the corners of each named node
and of the whole picture to the log; with \\def\\nolabels{} it draws everything
except the free text.  The test compiles the figure three ways and fails if

  1. two named nodes overlap (boxes and labels alike), or
  2. under any free label the drawing without labels has an edge: a curve, a
     line, a marker or the outline of a shape passing under the text.  A label
     may sit on a fill or a smooth shading, and on nothing else.

In the automatic mode (\\autocheck in the figure, or --auto here) every node is
checked: a node that contains another is a frame and is left out; both tests
look at each node 1.5pt in from its border (inside the default inner sep of
3.33pt), the second a further 0.6pt in, with only the text of the nodes hidden.

Usage: python3 tikz_overlap_check.py [--auto] FIGURE.tex [PAD_PT]
Exit status 0 when both tests pass.
"""
import os
import re
import subprocess
import sys
import tempfile

import numpy as np
from PIL import Image

DPI = 200
EDGE = 12      # grey levels between neighbouring pixels that make an edge
show = None
if '--show' in sys.argv:
    k = sys.argv.index('--show'); show = sys.argv[k + 1]; del sys.argv[k:k + 2]
args = [a for a in sys.argv[1:] if a != '--auto']
force = '--auto' in sys.argv
tex = os.path.abspath(args[0])
pad = float(args[1]) if len(args) > 1 else 1.0
pre = '\\def\\forceauto{}' if force else ''
base = os.path.splitext(os.path.basename(tex))[0]
src = open(tex).read()


def build(prefix, work):
    job = os.path.join(work, base + prefix.replace('\\', '').replace('{', '').replace('}', ''))
    with open(job + '.tex', 'w') as f:
        f.write(prefix + src)
    subprocess.run(['pdflatex', '-interaction=nonstopmode', '-halt-on-error',
                    '-output-directory', work, job + '.tex'],
                   cwd=os.path.dirname(tex), stdout=subprocess.DEVNULL, check=True)
    return job


with tempfile.TemporaryDirectory() as work:
    job = build(pre + '\\def\\checkbboxes{}', work)
    log = open(job + '.log', errors='replace').read().replace('\n', '')
    nodes = {}
    num = r' ([-\d.]+)pt'
    for m in re.findall(r'BBOX (box|label|auto) (\S+)' + num * 8, log):
        kind, name, c = m[0], m[1], [float(v) for v in m[2:]]
        xs, ys = c[0::2], c[1::2]
        # the extent over all four corners: a rotated node has its anchors rotated
        nodes[name] = (kind, min(xs), min(ys), max(xs), max(ys))
    pic = re.search(r'PICBB ([-\d.]+)pt ([-\d.]+)pt ([-\d.]+)pt ([-\d.]+)pt', log)
    border = float(re.search(r'BORDER ([\d.]+)pt', log).group(1))
    px0, py0, px1, py1 = map(float, pic.groups())
    fails = []
    auto = any(k == 'auto' for k, *_ in nodes.values())
    if auto:
        # leave out points (no extent), empty nodes (twice the default inner sep,
        # 6.67pt square, where any node with text is larger) and frames (a node
        # that contains another)
        pts = {n for n, (k, x1, y1, x2, y2) in nodes.items()
               if x2 - x1 < 0.5 or y2 - y1 < 0.5 or (x2 - x1 < 7.5 and y2 - y1 < 7.5)}
        for n in pts:
            del nodes[n]
        frames = set()
        for a, (_, ax1, ay1, ax2, ay2) in nodes.items():
            for b, (_, bx1, by1, bx2, by2) in nodes.items():
                if a != b and ax1 <= bx1 + 0.01 and ay1 <= by1 + 0.01 and ax2 >= bx2 - 0.01 and ay2 >= by2 - 0.01:
                    frames.add(a)
        for n in frames:
            del nodes[n]
        # nodes carry their inner sep (0.33 em unless set); compare them 1.5pt in
        # from their borders, which still checks a margin of about 1.8pt around
        # text with the default
        for n, (k, x1, y1, x2, y2) in list(nodes.items()):
            nodes[n] = (k, x1 + 1.5, y1 + 1.5, x2 - 1.5, y2 - 1.5)
    names = sorted(nodes)
    for a in range(len(names)):
        for b in range(a + 1, len(names)):
            _, ax1, ay1, ax2, ay2 = nodes[names[a]]
            _, bx1, by1, bx2, by2 = nodes[names[b]]
            if min(ax2, bx2) - max(ax1, bx1) > -pad and min(ay2, by2) - max(ay1, by1) > -pad:
                fails.append('nodes %s and %s overlap' % (names[a], names[b]))
    job2 = build(pre + '\\def\\nolabels{}', work)
    subprocess.run(['pdftoppm', '-r', str(DPI), '-png', '-singlefile', job2 + '.pdf', job2], check=True)
    img = np.asarray(Image.open(job2 + '.png').convert('L'), dtype=float)
    H, Wd = img.shape
    s = DPI / 72.27
    bad = []
    for name, (kind, x1, y1, x2, y2) in nodes.items():
        if kind == 'box':
            continue
        if kind == 'auto':
            # a further 0.6pt in: rounded corners of radius 2pt reach that far inside
            x1, y1, x2, y2 = x1 + 0.6, y1 + 0.6, x2 - 0.6, y2 - 0.6
        c1 = int((x1 - px0 + border) * s) + 1; c2 = int((x2 - px0 + border) * s) - 1
        r1 = int((py1 - y2 + border) * s) + 1; r2 = int((py1 - y1 + border) * s) - 1
        patch = img[max(r1, 0):min(r2, H), max(c1, 0):min(c2, Wd)]
        if patch.size == 0:
            fails.append('label %s has an empty pixel box' % name)
            continue
        # background, a fill, or a smooth shading may lie under a label; an edge may
        # not: a line, a curve, a marker or the outline of a shape changes the grey
        # level between neighbouring pixels by far more than any shading does
        jump = np.zeros_like(patch)
        if patch.shape[0] > 1:
            jump[1:, :] = np.maximum(jump[1:, :], np.abs(np.diff(patch, axis=0)))
        if patch.shape[1] > 1:
            jump[:, 1:] = np.maximum(jump[:, 1:], np.abs(np.diff(patch, axis=1)))
        edges = int(np.sum(jump > EDGE))
        if edges:
            fails.append('label %s lies over %d edge pixels' % (name, edges))
            bad.append((name, c1, r1, c2, r2))
    if auto:
        print('%s: %d nodes checked automatically (%d frames and %d points or empty nodes left out)'
              % (base, len(nodes), len(frames), len(pts)))
    else:
        print('%s: %d boxes, %d labels checked' % (base, sum(k == 'box' for k, *_ in nodes.values()),
                                                   sum(k == 'label' for k, *_ in nodes.values())))
    for f in fails:
        print('  FAIL: ' + f)
    print('  PASS' if not fails else '  %d failures' % len(fails))
    if show:
        # the figure as drawn, with the box of every label that failed in red
        from PIL import ImageDraw
        job3 = build(pre, work)
        subprocess.run(['pdftoppm', '-r', str(DPI), '-png', '-singlefile', job3 + '.pdf', job3], check=True)
        im = Image.open(job3 + '.png').convert('RGB'); dr = ImageDraw.Draw(im)
        for name, c1, r1, c2, r2 in bad:
            dr.rectangle([c1, r1, c2, r2], outline=(255, 0, 0)); dr.text((c1, max(r1 - 10, 0)), name, fill=(255, 0, 0))
        im.save(show)
    sys.exit(1 if fails else 0)
