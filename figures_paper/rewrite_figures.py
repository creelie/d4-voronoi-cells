"""
rewrite_figures.py -- point the manuscript at the composite figures.

For each group: build one figure* environment carrying the composite, with the
panel captions joined in order and lettered, put it where the first panel of
the group stood, delete the rest, and rewrite every reference to a panel into
a reference to the composite with the panel letter after it.

Run from the directory holding D4.tex:  python3 figures/rewrite_figures.py
"""
import re
import sys
from pathlib import Path

TEX = Path('D4.tex')

#  new label, old labels in panel order, composite file, width as a fraction
#  of \textwidth
GROUPS = [
    ('fig:corner', ['fig:cornera', 'fig:cornerb', 'fig:cornerc', 'fig:cornerd'],
     'fig_c_corner', 0.78),
    ('fig:cell24', ['fig:cell24a', 'fig:cell24b', 'fig:cell24c', 'fig:cell24d'],
     'fig_c_cell24', 0.78),
    ('fig:boundary', ['fig:boundary-bounds', 'fig:boundary-routes'],
     'fig_c_boundary', 1.0),
    ('fig:counting', ['fig:contact-count', 'fig:contact-settled',
                      'fig:pair-budget', 'fig:pair-weight'],
     'fig_c_counting', 0.86),
    ('fig:rig', ['fig:riga', 'fig:rigb', 'fig:rigc', 'fig:rigd'],
     'fig_c_rig', 0.82),
    ('fig:rootmeet', ['fig:root-meet', 'fig:root-meet-size'],
     'fig_c_rootmeet', 1.0),
    ('fig:structab', ['fig:structa', 'fig:structb'], 'fig_c_structab', 0.86),
    ('fig:slack', ['fig:slack-continuation', 'fig:slack-rate'],
     'fig_c_slack', 1.0),
    ('fig:structcd', ['fig:structc', 'fig:structd'], 'fig_c_structcd', 0.86),
    ('fig:pfone', ['fig:pfa', 'fig:pfc'], 'fig_c_pf1', 0.86),
    ('fig:pftwo', ['fig:pfd', 'fig:pfb'], 'fig_c_pf2', 0.86),
    ('fig:certpair', ['fig:certificate', 'fig:certificate-values'],
     'fig_c_certificate', 1.0),
    ('fig:dom3d', ['fig:domain3d', 'fig:domain3d-bb'], 'fig_c_domain3d', 1.0),
    ('fig:arcs', ['fig:arc1-breakpoints', 'fig:arc2-breakpoints'],
     'fig_c_arcs', 1.0),
    ('fig:certmap', ['fig:cert-map', 'fig:cert-regions'],
     'fig_c_certmap', 1.0),
    ('fig:cap', ['fig:cap-geometry', 'fig:cap-volume'], 'fig_c_cap', 1.0),
    ('fig:defect', ['fig:defect-positivity', 'fig:defect-small'],
     'fig_c_defect', 1.0),
    ('fig:swap', ['fig:swap-configs', 'fig:swap-crossover'],
     'fig_c_swap', 1.0),
]

LETTERS = 'abcdefgh'


def find_block(s, label):
    """the whole figure environment carrying this label"""
    for m in re.finditer(r'\\begin\{figure\*?\}.*?\\end\{figure\*?\}\n?', s,
                         re.S):
        if ('\\label{%s}' % label) in m.group(0):
            return m
    raise KeyError(label)


def caption_of(block):
    """the caption text, brace-matched so that nested braces survive"""
    i = block.index('\\caption{') + len('\\caption{')
    depth, j = 1, i
    while depth:
        if block[j] == '{':
            depth += 1
        elif block[j] == '}':
            depth -= 1
        j += 1
    return block[i:j - 1].strip()


def main():
    s = TEX.read_text()
    n_removed = 0

    for new_label, olds, img, frac in GROUPS:
        caps = []
        blocks = []
        for lab in olds:
            m = find_block(s, lab)
            blocks.append(m)
            caps.append(caption_of(m.group(0)))

        merged = ' '.join('(%s) %s' % (LETTERS[i], c)
                          for i, c in enumerate(caps))
        new_env = (
            '\\begin{figure*}[tb]\n'
            '\\centering\n'
            '\\includegraphics[width=%g\\textwidth]{%s}\n'
            '\\caption{%s}\n'
            '\\label{%s}\n'
            '\\end{figure*}\n' % (frac, img, merged, new_label))

        # replace from the last block backwards so earlier spans stay valid
        spans = sorted((m.start(), m.end()) for m in blocks)
        for k, (a, b) in enumerate(reversed(spans)):
            if k == len(spans) - 1:          # the first block in the file
                s = s[:a] + new_env + s[b:]
            else:
                s = s[:a] + s[b:]
                n_removed += 1

    # references: a panel becomes the composite plus its letter
    ref_map = {}
    for new_label, olds, _, _ in GROUPS:
        for i, lab in enumerate(olds):
            ref_map[lab] = (new_label, LETTERS[i])

    def fix(m):
        cmd, body = m.group(1), m.group(2)
        keys = [k.strip() for k in body.split(',')]
        if len(keys) == 1 and keys[0] in ref_map:
            lab, letter = ref_map[keys[0]]
            return '\\%s{%s}(%s)' % (cmd, lab, letter)
        if any(k in ref_map for k in keys):
            print('  multi-key reference left alone: %s' % m.group(0))
        return m.group(0)

    s, n_refs = re.subn(r'\\(cref|Cref|ref|autoref)\{([^}]+)\}',
                        lambda m: fix(m), s), 0
    s = s[0] if isinstance(s, tuple) else s

    TEX.write_text(s)

    remaining = len(re.findall(r'\\begin\{figure\*?\}', s))
    print('removed %d figure environments' % n_removed)
    print('%d float environments remain' % remaining)
    for lab in ref_map:
        if ('\\label{%s}' % lab) in s:
            print('  WARNING: old label still defined: %s' % lab)


if __name__ == '__main__':
    main()
