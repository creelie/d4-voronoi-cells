Compile D4.tex with pdfLaTeX, three passes.

Sections 7.10 and 7.11 prove the uniform annular correction and the
comparison between pair and triple cap intersections. Section 8 gives
the density bound below 0.634625. Appendices I through K give the
constants, interval derivations and the 26-neighbour profile check.
The new diagrams and the revised density figure are drawn in TikZ.

The standard-library checkers can be run from this directory.
  python anc/annular_constants.py
  python anc/cap_comparison.py
  python anc/cap_profiles_26.py
  python anc/exterior_test.py
The corresponding results are in the log files under anc/.

The unrestricted cases (G) and (C) remain open. The title-page imprint includes the verified version DOI. The earlier computational
and formal-verification materials retain their stated scope.
