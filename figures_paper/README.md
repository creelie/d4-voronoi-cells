# Figure composition for the manuscript

Note for v1.3.0. The two scripts here record how the composite figures
`fig_c_*.png` were first assembled from single panels, when the manuscript
was set two-column; the manuscript is now single-column amsart and the
composites are carried in `paper/figures/`. The single-panel sources that
`merge_panels.py` reads are not part of this package. The data panels of
seven composites were redrawn for v1.3.0, and three figures added, by the scripts in
`paper/figures_new/`, which reuse the rendered three-dimensional panels of
the earlier composites (`paper/figures_new/old/`) and pass a label
collision test before writing; those scripts, not these, reproduce the
figures shipped with the paper. What follows is the original description.

The manuscript was carrying forty-nine single-panel floats, each asked for at
about 0.85 of a column. In a hundred-page two-column revtex document that is
more floats than the class can place near the text that cites them: the build
emitted thirty "float is stuck" notices and dumped up to four unrelated
figures onto a page.

These two scripts group the panels that belong together into composites drawn
at the width they are printed at, and rewire the manuscript to use them.

- `merge_panels.py` reads the individual panels in `figures/` and writes
  eighteen composites `fig_c_*.png`, each laid out on a grid with lettered
  panels. Panels are scaled so that a panel inside a composite is at least as
  large on the page as it was on its own, so nothing loses legibility.
- `rewrite_figures.py` replaces each group of figure environments in `D4.tex`
  with one `figure*` carrying the composite and the panel captions joined in
  order, and rewrites every reference to a panel into a reference to the
  composite plus its letter.

Run from the directory holding `D4.tex`:

    python3 figures_paper/merge_panels.py      # from inside figures/
    python3 figures_paper/rewrite_figures.py

Result: forty-nine floats become twenty-three, no float is stuck, no column is
overfull, and twenty-one of the twenty-two figure-bearing pages carry exactly
one composite at the top with text below.
