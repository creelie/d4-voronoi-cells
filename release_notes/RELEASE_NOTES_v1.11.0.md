# v1.11.0

Supplementary package for *The Sphere Packing Problem in Dimension 4 and the
Twenty-Four-Cell Conjecture*.  This release proves one new result on
statement (C), the bound on the pair terms for 25 or more centres within
sqrt 6.  It shows where the centres of a counterexample must lie: for each
count M from 25 to 30, how many of the M centres lie within 2.05, 2.1, 2.15
and 2.2 of the centre, and, from 27 centres on, a radius below 2.444 within
which 24 of them lie.  Conjecture 1.6 remains open.  (G), and (C) for 25 to
30 centres, are still the two statements that would close it; the new
counts narrow where a counterexample to (C) can sit, but they do not close
any count.  The new code was run, and the certificates found and checked,
with the assistance of Claude, an AI model made by Anthropic.

## What the paper now proves

Everything proved in v1.10.0 still stands.  New:

- **Proposition 2.41 (where the centres lie from twenty-five to
  thirty).**  If M centres lie within sqrt 6 of the centre, 25 <= M <= 30, and
  T(Y) <= 8, then at least the following numbers of them lie within 2.05,
  2.1, 2.15 and 2.2.  The caps alone force the numbers in parentheses, and a
  dash means that no count above them is proved.

        M     2.05      2.1       2.15      2.2
        25    17 (15)   21 (20)   22 (21)   23 (22)
        26    - (12)    20 (18)   22 (21)   23 (22)
        27    14 (9)    20 (17)   21 (20)   22 (21)
        28    - (6)     21 (16)   22 (19)   23 (21)
        29    19 (2)    22 (15)   23 (19)   23 (21)
        30    22 (0)    22 (14)   23 (18)   23 (20)

  Twenty-four of the M centres lie within 2.25 at M = 30, within 2.3 at
  M = 29 and within 2.35 at M = 27 and 28.  Lemma 2.38 gave 2.444.

  Each count is the single split N(rho) <= n - 1 against N(rho) >= n of the
  case certificates of Theorem 2.39, run with M centres in place of thirty,
  at degree 12 in the angle and 4 in the distance.  The first case carries a
  certificate whose largest value over the count vectors of the case lies
  below 9 pi^2/8 - 8 = 3.10330, so T(Y) > 8 there.  The largest of these
  values is 3.10276, at M = 27 and rho = 2.35.  Each certificate is checked like those of Theorem
  2.39.  The check covers positivity by exact LDL^T, the inequality K <= Pi
  by a branch and bound with Pi in ball arithmetic, the bounds on the bins,
  and the largest value over the count vectors by exact dynamic
  programming.  Each check takes at most three minutes on one core.

  Table 3 sets the proved counts against what the caps alone force, which
  is much less where the centres are crowded.  At 29 centres, for example,
  the caps put only two within 2.05 and fifteen within 2.1, against
  nineteen and twenty-two proved.  Figure 17 draws the table.
- **Remark 2.40 (what is left at thirty).**  Theorem 2.39 leaves at most
  24 centres within 2.4, and Proposition 2.41 puts 24 within 2.25.  So the
  case it leaves open is two shells: 24 centres within 2.25, 23 of them
  within 2.15 and 22 within 2.05, and six beyond 2.4, with none in between.
  The introduction and the conclusion now describe the case at thirty this
  way.

## What the counts do not give

- At 2.0161, where the case of 23 close centres is decided, a single split
  gives little.  In floating point its first case already exceeds the level
  when it asks for twelve centres within 2.0161 at 29 centres, for three at
  25 and 28, and for one at 26 and 27.
- Below 27 centres a single split does not bring the 2.444 of Lemma 2.38
  under 2.4.  Its first case stays at 3.1124 and 3.1068 in floating point
  at 25 and 26 centres.
- The counts say nothing about the directions.  The case of 24 close
  centres needs each of its 24 directions within 0.093 of a root system,
  and the case of 23 a certificate that tells a packing from the fictitious
  configurations of Remark 2.37.  Neither follows from counts.

## The margins

A certificate holds its kernel below Pi by a margin at the sampled pairs.
This lets the exact branch and bound close quickly.  With M centres the
margin is given up at each of the M(M - 1)/2 pairs.  At 29 centres a
margin of 5e-5 where the caps are disjoint, which kept the checks of
Proposition 2.44 (the old 2.43) fast, raised the bound by about 0.02 and
put the first three cases tried above the level.  The margins here are therefore chosen
for each case from the room its floating-point optimum leaves, between
2.5e-6 and 5e-5.  Each spec file records the margins tried.

## The paper

- New: Proposition 2.41, Table 3 and Figure 17, placed after Remark 2.40.
  Figure 17 passes the collision test
  (paper/figures_new/tikz_overlap_check.py).
- The introduction and the conclusion mention the new counts and the two
  shells at thirty.  The appendix on the certificates says how the margins
  are chosen.
- Numbering after the new material moves by one.  Theorem 2.41 is now 2.42,
  Corollary 2.42 is 2.43, Proposition 2.43 is 2.44 and Remark 2.44 is 2.45.
  Figures 17 to 43 are now 18 to 44, and Tables 3 to 7 are now 4 to 8.  The
  file map of the README and the header of multi_cap/pushout_check.py
  follow the new numbering.
- 184 pages; no overfull boxes, no undefined references, no missing
  glyphs.

## New and changed files

- multi_cap/radial_certificates/case_M_rRHO.json and case_M_rRHO_spec.json
  for M = 25 to 30: the certificates of Proposition 2.41 and their specs.
- multi_cap/runs/radial_case_sdp_M_rRHO.log and radial_case_check_M_rRHO.log:
  the runs that found them and the exact checks.
- paper/figures_new/tikz_cradial.tex, cradial_data.py and cradial_counts.dat:
  Figure 17.  cradial_data.py reads the check logs, computes the counts
  forced by the caps in ball arithmetic, and prints the rows of Table 3.
- CITATION.cff and .zenodo.json: v1.11.0.

## What is left

What is left is as in v1.10.0.  Conjecture 1.6 follows from (G) and from
(C) for 25 to 30 centres, and neither is proved in general.  Proposition
2.41 narrows where a counterexample to (C) can be: most of its centres are
close to c, and from 27 centres on 24 of them are within 2.35 or less.
What (C) still needs is a per-direction localisation for 24 close centres
and a certificate that tells packings from fictitious configurations for
23, and the pair terms supply neither.
