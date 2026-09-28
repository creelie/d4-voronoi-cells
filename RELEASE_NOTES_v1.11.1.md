# v1.11.1

Supplementary package for *The Sphere Packing Problem in Dimension 4 and the
Twenty-Four-Cell Conjecture*.  This release adds two counts to
Proposition 2.41 and one negative result on the case of statement (C) with
24 centres close to the centre.  Conjecture 1.6 remains open.  (G), and (C)
for 25 to 30 centres, are still the two statements that would close it.
The new code was run, and the certificates found and checked, with the
assistance of Claude, an AI model made by Anthropic.

## What the paper now proves

Everything proved in v1.11.0 still stands.  New, in Proposition 2.41: if
M centres lie within sqrt 6 of the centre and T(Y) <= 8, then

- at M = 26, at least 15 of them lie within 2.05 (the caps alone force 12);
- at M = 29, at least 8 of them lie within 2.0161 (the caps force none).

Each is a single split of the case certificates of Theorem 2.39, at degree
12 in the angle and 4 in the distance, checked like the others: positivity
by exact LDL^T, the inequality K <= Pi by a branch and bound with Pi in
ball arithmetic, the bounds on the bins, and the largest value over the
count vectors by exact dynamic programming.  The largest values are
3.10135 and 3.10157, below 9 pi^2/8 - 8 = 3.10330; the checks take 19
seconds and five minutes on one core.  Table 3 and Figure 17 now show the
count at 26 centres.

Not proved: 17 of 28 centres within 2.05.  Its floating-point certificate
reaches 3.10307, only 2.3e-4 below the level, which leaves too thin a
margin over the 378 pairs, and the exact check did not finish within an
hour.  The table keeps its dash there.

## What the hole programme does not give

Remark 2.37 now records, in floating point, why the labelled three-point
programme cannot settle the case of 24 close centres.  A centre within
sqrt 6 beside 24 centres within 2.0161 needs a direction at inner product
at most 0.6141 with 24 directions of slack 0.008.  The value of a labelled
certificate is a polynomial in the number n of directions, and it excludes
n when it is negative.

- With the further direction, no certificate of degree 7, 10 or 12 excludes
  n = 24, and none of degree 10 or 12 excludes even n = 24.5.
- Without it, none of degree 10 excludes n = 24.9, and none of degree 12
  excludes n = 24.7.

Read as counts, the further direction would have to be worth at least 0.7
of a point, and these certificates value it at least half a point lower.
The values of the certificates at n = 24 do fall with the degree, from
0.127 at degree 10 to 0.0096 at degree 12, which earlier looked like
progress.  They fall as fast when the further direction may lie at inner
product 0.71, where the root system with a direction in one of its deep
holes meets every constraint (0.144 and 0.0116), so the fall says nothing
about exclusion.  The sampling only lowers the values, so a positive value
is no artefact of it.  The statement itself has room: the deepest hole
found in any 24-point code of slack 0.008 is at 0.678.  What is missing is
a certificate that sees single directions.

## New and changed files

- multi_cap/radial_certificates/case_26_r2.05.json, case_29_r2.0161.json
  and their specs, with runs/radial_case_{sdp,check}_{26_r2.05,29_r2.0161}.log.
- multi_cap/hole_labelled_sdp.py: the labelled three-point programme
  (floating point, CVXPY and Clarabel), with runs/hole_labelled_*.log.
- paper/D4.tex, paper/D4.pdf (185 pages), paper/figures/fig_tikz_cradial.png,
  paper/figures_new/cradial_counts.dat: Proposition 2.41, Table 3,
  Figure 17, Remark 2.37 and the appendix on the computations.
- README.md, CITATION.cff and .zenodo.json: v1.11.1.
