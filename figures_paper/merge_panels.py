"""
merge_panels.py -- the composite figures of the D4 paper.

The manuscript was carrying forty-nine separate single-panel floats, each
asked for at about 0.85 of a column. In a hundred-page two-column document
revtex cannot place that many small floats near the text that cites them, so
it defers them and then dumps them: thirty "float is stuck" notices, and pages
carrying four unrelated figures with a few lines of text squeezed between.

This script groups panels that belong together into one composite each, drawn
at the width it is printed at, and writes them into this directory. Panels are
scaled so that a panel in a composite is at least as large on the page as it
was on its own, so nothing loses legibility in the process.

Writes fig_c_*.png. The originals are left alone; rewrite_figures.py rewires
the manuscript.
"""
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent

COL = 3.4          # revtex column width, inches
TEXT = 7.0         # revtex text width, inches

plt.rcParams.update({'font.family': 'serif', 'font.size': 9,
                     'figure.dpi': 200, 'savefig.dpi': 400,
                     'savefig.bbox': 'tight', 'savefig.pad_inches': 0.02})


def trimmed(name):
    """the panel with its uniform white border removed, and its aspect ratio"""
    a = plt.imread(HERE / (name + '.png'))
    if a.ndim == 3 and a.shape[2] == 4:
        al = a[..., 3:4]
        a = a[..., :3] * al + (1.0 - al)
    g = a.mean(axis=2) if a.ndim == 3 else a
    ink = g < 0.995
    rows, cols = np.any(ink, axis=1), np.any(ink, axis=0)
    if not rows.any():
        return a, a.shape[1] / a.shape[0]
    r0, r1 = np.argmax(rows), len(rows) - np.argmax(rows[::-1])
    c0, c1 = np.argmax(cols), len(cols) - np.argmax(cols[::-1])
    pad = 3
    r0, c0 = max(0, r0 - pad), max(0, c0 - pad)
    r1, c1 = min(a.shape[0], r1 + pad), min(a.shape[1], c1 + pad)
    a = a[r0:r1, c0:c1]
    return a, a.shape[1] / a.shape[0]


def composite(out, names, ncols, width_in, letters=True, hspace=0.10,
              wspace=0.04):
    """Lay the named panels out in a grid `ncols` wide, at `width_in` inches.

    Row heights follow the tallest panel of each row, so a portrait panel
    beside a landscape one does not force the landscape one to grow."""
    imgs = [trimmed(n) for n in names]
    nrows = int(np.ceil(len(imgs) / ncols))
    cell_w = width_in / ncols
    row_h = []
    for r in range(nrows):
        row = imgs[r * ncols:(r + 1) * ncols]
        row_h.append(max(cell_w / ar for _, ar in row))
    fig_h = sum(row_h) * (1 + hspace) + (0.18 if letters else 0.0) * nrows
    fig = plt.figure(figsize=(width_in, fig_h))
    gs = fig.add_gridspec(nrows, ncols, height_ratios=row_h,
                          hspace=hspace, wspace=wspace)
    for i, (img, ar) in enumerate(imgs):
        ax = fig.add_subplot(gs[i // ncols, i % ncols])
        ax.imshow(img, interpolation='lanczos')
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(False)
        if letters:
            ax.set_title('(%s)' % 'abcdefgh'[i], fontsize=9, loc='left',
                         pad=2.5, fontweight='bold')
    fig.savefig(HERE / (out + '.png'))
    plt.close(fig)
    w, h = plt.imread(HERE / (out + '.png')).shape[1::-1]
    print('%-22s %d panels, %d x %d, prints %.1f x %.1f in'
          % (out, len(names), w, h, width_in, width_in * h / w))


# (output, panels, columns, printed width)
GROUPS = [
    ('fig_c_corner', ['fig_corner_a', 'fig_corner_b',
                      'fig_corner_c', 'fig_corner_d'], 2, 0.78 * TEXT),
    ('fig_c_cell24', ['fig_cell24_a', 'fig_cell24_b',
                      'fig_cell24_c', 'fig_cell24_d'], 2, 0.78 * TEXT),
    ('fig_c_boundary', ['fig_boundary_bounds_a', 'fig_boundary_bounds_b'],
     2, TEXT),
    ('fig_c_counting', ['fig_contact_count_a', 'fig_contact_count_b',
                        'fig_pair_budget_a', 'fig_pair_budget_b'], 2,
     0.86 * TEXT),
    ('fig_c_rig', ['fig_rig_a', 'fig_rig_b', 'fig_rig_c', 'fig_rig_d'],
     2, 0.82 * TEXT),
    ('fig_c_rootmeet', ['fig_root_meet_a', 'fig_root_meet_b'], 2, TEXT),
    ('fig_c_structab', ['fig_struct_a', 'fig_struct_b'], 2, 0.86 * TEXT),
    ('fig_c_slack', ['fig_slack_continuation_a', 'fig_slack_continuation_b'],
     2, TEXT),
    ('fig_c_structcd', ['fig_struct_c', 'fig_struct_d'], 2, 0.86 * TEXT),
    ('fig_c_pf1', ['fig_pf_a', 'fig_pf_c'], 2, 0.86 * TEXT),
    ('fig_c_pf2', ['fig_pf_d', 'fig_pf_b'], 2, 0.86 * TEXT),
    ('fig_c_certificate', ['fig_certificate_a', 'fig_certificate_b'], 2, TEXT),
    ('fig_c_domain3d', ['fig_domain3d_a', 'fig_domain3d_b'], 2, TEXT),
    ('fig_c_arcs', ['fig_arc1_breakpoints', 'fig_arc2_breakpoints'], 2, TEXT),
    ('fig_c_certmap', ['fig_certificate_map_a', 'fig_certificate_map_b'],
     2, TEXT),
    ('fig_c_cap', ['fig_cap_geometry_a', 'fig_cap_geometry_b'], 2, TEXT),
    ('fig_c_defect', ['fig_defect_positivity_a', 'fig_defect_positivity_b'],
     2, TEXT),
    ('fig_c_swap', ['fig_swap_curves_a', 'fig_swap_curves_b'], 2, TEXT),
]

if __name__ == '__main__':
    for out, names, ncols, w in GROUPS:
        composite(out, names, ncols, w)
    print('%d composites, replacing %d panels'
          % (len(GROUPS), sum(len(g[1]) for g in GROUPS)))
