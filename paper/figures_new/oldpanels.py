"""
oldpanels.py -- reuse of the rendered panels of a composite figure.

Several composites of the paper combine ray-traced three-dimensional panels,
which are kept, with data panels that have been redrawn here so that no label
touches a curve, a bar or another label.  The rendered panels are cut out of
the composite along its white gutters and placed unchanged on the new grid.
"""
import numpy as np
import matplotlib.pyplot as plt


def load(path):
    a = plt.imread(path)
    if a.ndim == 3 and a.shape[2] == 4:
        al = a[..., 3:4]
        a = a[..., :3] * al + (1.0 - al)
    return a


def white_bands(mask, axis, min_len):
    """maximal runs of all-white rows (axis=0) or columns (axis=1) of at least min_len"""
    line = mask.all(axis=1 - axis)
    runs, start = [], None
    for i, w in enumerate(line):
        if w and start is None:
            start = i
        if not w and start is not None:
            if i - start >= min_len:
                runs.append((start, i))
            start = None
    if start is not None and len(line) - start >= min_len:
        runs.append((start, len(line)))
    return runs


def split_rows(img, min_gap=12):
    """the composite cut into rows along its widest interior white bands"""
    white = img.mean(axis=2) > 0.985
    bands = [b for b in white_bands(white, 0, min_gap) if b[0] > 0 and b[1] < img.shape[0]]
    cuts = [0] + [(a + b) // 2 for a, b in bands] + [img.shape[0]]
    return [img[cuts[i]:cuts[i + 1]] for i in range(len(cuts) - 1)]


def split_cols(img, min_gap=12):
    white = img.mean(axis=2) > 0.985
    bands = [b for b in white_bands(white, 1, min_gap) if b[0] > 0 and b[1] < img.shape[1]]
    cuts = [0] + [(a + b) // 2 for a, b in bands] + [img.shape[1]]
    return [img[:, cuts[i]:cuts[i + 1]] for i in range(len(cuts) - 1)]


def trim(img, pad=4):
    """the image with its uniform white border removed"""
    ink = img.mean(axis=2) < 0.985
    rows, cols = np.where(ink.any(axis=1))[0], np.where(ink.any(axis=0))[0]
    if len(rows) == 0:
        return img
    r0, r1 = max(rows[0] - pad, 0), min(rows[-1] + pad + 1, img.shape[0])
    c0, c1 = max(cols[0] - pad, 0), min(cols[-1] + pad + 1, img.shape[1])
    return img[r0:r1, c0:c1]


def show(ax, img):
    ax.imshow(img, interpolation="lanczos")
    ax.set_axis_off()
