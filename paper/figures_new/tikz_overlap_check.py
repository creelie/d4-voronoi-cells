#!/usr/bin/env python3
"""
tikz_overlap_check.py -- the collision test of the TikZ figures.

A figure that uses it names every text node and lists the names in \\boxes
(nodes drawn as boxes, which may carry arrows to their borders) and \\labels
(free text).  With \\def\\checkbboxes{} it writes the corners of each named node
and of the whole picture to the log; with \\def\\nolabels{} it draws everything
except the free text.  The test compiles the figure three ways and fails if

  1. two named nodes overlap (boxes and labels alike), or
  2. under any free label the drawing without labels is not flat: a curve, a
     line, a marker or the edge of a box passing under the text.  A label may
     sit on a uniform fill, such as a shaded band, and on nothing else.

Usage: python3 tikz_overlap_check.py FIGURE.tex [PAD_PT]
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
tex = os.path.abspath(sys.argv[1])
pad = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0
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
    job = build('\\def\\checkbboxes{}', work)
    log = open(job + '.log', errors='replace').read().replace('\n', '')
    nodes = {}
    num = r' ([-\d.]+)pt'
    for m in re.findall(r'BBOX (box|label) (\S+)' + num * 8, log):
        kind, name, c = m[0], m[1], [float(v) for v in m[2:]]
        xs, ys = c[0::2], c[1::2]
        # the extent over all four corners: a rotated node has its anchors rotated
        nodes[name] = (kind, min(xs), min(ys), max(xs), max(ys))
    pic = re.search(r'PICBB ([-\d.]+)pt ([-\d.]+)pt ([-\d.]+)pt ([-\d.]+)pt', log)
    border = float(re.search(r'BORDER ([\d.]+)pt', log).group(1))
    px0, py0, px1, py1 = map(float, pic.groups())
    fails = []
    names = sorted(nodes)
    for a in range(len(names)):
        for b in range(a + 1, len(names)):
            _, ax1, ay1, ax2, ay2 = nodes[names[a]]
            _, bx1, by1, bx2, by2 = nodes[names[b]]
            if min(ax2, bx2) - max(ax1, bx1) > -pad and min(ay2, by2) - max(ay1, by1) > -pad:
                fails.append('nodes %s and %s overlap' % (names[a], names[b]))
    job2 = build('\\def\\nolabels{}', work)
    subprocess.run(['pdftoppm', '-r', str(DPI), '-png', '-singlefile', job2 + '.pdf', job2], check=True)
    img = np.asarray(Image.open(job2 + '.png').convert('L'), dtype=float)
    H, Wd = img.shape
    s = DPI / 72.27
    for name, (kind, x1, y1, x2, y2) in nodes.items():
        if kind != 'label':
            continue
        c1 = int((x1 - px0 + border) * s) + 1; c2 = int((x2 - px0 + border) * s) - 1
        r1 = int((py1 - y2 + border) * s) + 1; r2 = int((py1 - y1 + border) * s) - 1
        patch = img[max(r1, 0):min(r2, H), max(c1, 0):min(c2, Wd)]
        if patch.size == 0:
            fails.append('label %s has an empty pixel box' % name)
            continue
        # background, or one flat fill (a shaded band): nothing drawn varies under the label
        if patch.max() - patch.min() > 2:
            dark = int(np.sum(np.abs(patch - np.median(patch)) > 2))
            fails.append('label %s lies over %d drawn pixels' % (name, dark))
    print('%s: %d boxes, %d labels checked' % (base, sum(k == 'box' for k, *_ in nodes.values()),
                                               sum(k == 'label' for k, *_ in nodes.values())))
    for f in fails:
        print('  FAIL: ' + f)
    print('  PASS' if not fails else '  %d failures' % len(fails))
    sys.exit(1 if fails else 0)
