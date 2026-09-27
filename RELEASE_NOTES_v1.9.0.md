# v1.9.0

Supplementary package for *The Sphere Packing Problem in Dimension 4 and the
Twenty-Four-Cell Conjecture*.  This release goes with the journal form of
the paper.  It adds the second order of the cell volume at the root system
(Section 2.9), an exact closed curve of contact cells of volume 8, a second
twenty-four-point code far from the root system, and the reduction of
Conjecture 1.6 to two explicit statements, (G) and (C) (Proposition 2.30),
and it proves (C) for every count of at least 31 centres within sqrt6
(Theorem 2.32).  Fifteen parts of the paper are now checked in Lean 4, eight
of them new.  Conjecture 1.6 remains open, and the paper says so: (G), and
(C) for 25 to 30 centres, are stated as the precise open statements, with
what a proof of each needs (Section 2.10).  The new code was developed, and its
computations run, with the assistance of Claude, an AI model made by
Anthropic.

## What the paper now proves

Everything proved in v1.8.0 still stands.  New:

- **Section 2.9, the second order at the root system.**
  - Lemma 2.19: the second variation of the volume of a 4-polytope under
    moving normals and support numbers; the volume is C^2 across changes of
    combinatorics.
  - Proposition 2.20: the Hessian H of the cell volume at the root system in
    closed form over the 96 triangles of the 24-cell; on the push-outs it is
    2 sum_edges eta_i eta_j - 4 sum eta_i^2.
  - Proposition 2.21: on the first-order packing cone, H >= -(sum delta_i)^2,
    with equality on the ray of one centre pushed out, by an exact
    copositivity certificate H + cc^T = P + B^T N B.  So the second-order
    model vol - 8 >= (2/3) S - (1/2) S^2 holds, and stays positive up to
    S = 4/3.
  - Lemma 2.23 and Proposition 2.24: pure push-outs integrated exactly, by
    Brunn-Minkowski.
  - Remark 2.22 says what the second order leaves: tilting is not monotone,
    the third-order remainder is not bounded, and the part of V(Y) outside
    the inversion hull needs a sharper bound than 13958 Theta^4 (it is below
    1e-6 along both rays up to 0.197).
  - Every volume is recomputed independently by intersecting halfspaces
    (multi_cap/second_order/independent_check.py): H against second
    differences (6e-8), the cone rows against exact distances (2e-9), the
    one-centre formula (1e-14), Lemma 2.23 on 400 push patterns.  Two
    statements of the draft these results came from are corrected: the
    eigenvalue -6.194 of H has multiplicity 8, and the push-out block is
    Adj - 4I.
- **Proposition 2.25 (the tilt block).**  The Hessian of the contact-cell
  volume in the directions alone is positive semidefinite: an exact LDL^T
  with 15 zero pivots (the 6 rotations and 9 further directions, none a
  strain), least positive eigenvalue 1/12.
- **Proposition 2.26 (a closed curve of contact cells of volume 8).**  Turn
  one A2 hexagon of roots by theta in its plane and tilt the other eighteen
  towards the orthogonal plane, with sin psi = 4C/(4C^2 + 3) and
  C = cos(pi/6 - theta).  The contact cell has volume exactly 8 for every
  theta in [0, pi/3], and the curve closes up at the root system.
  - For 0 < theta < pi/3 the directions are not a rotation of the root
    system: their largest inner product is at least 4C^2/(4C^2 + 3) > 1/2.
  - The proof integrates the hexagonal slices, with
    J1 = sqrt3 (4C^2 + 3)/(3C) and J2 = sqrt3 (16C^4 - 8C^2 + 9)/(12C^2),
    and reduces vol - 8 to the perfect square (4C - sigma(4C^2 + 3))^2.
  - So the root system is not a strict local minimum of the contact-cell
    volume, and no higher-order positivity can prove vol(Q_w) >= 8 near it.
- **Proposition 2.28 (a second twenty-four-point code).**  24 points of S^3
  with rational coordinates, pairwise inner products at most
  1/2 + 849/50000 = 0.51698, and d(W) >= 57/500 from the root system.  So
  statement (i) of "What is left" is false from slack 0.01698 on, well
  before 0.0374, where 25-point codes appear.
- **Proposition 2.30 (the two statements that remain).**  Conjecture 1.6
  follows from
  - (G): vol(V(Y) ∩ K(Y)) >= 8 for every packing set Y of exactly 24
    centres with 2 <= |y| < sqrt6 and T(Y) <= 8, with equality only at a
    copy of sqrt2 D4; and
  - (C): T(Y) > 8 for every packing set Y of at least 25 centres with
    2 <= |y| < sqrt6,
  where T(Y) is the right side of Lemma 2.15.  The proof combines Lemma
  2.15, the inversion hull (Proposition 2.13), Theorem 2.17 and Lemma 2.11.
  Remark 2.29 explains why the condition T(Y) <= 8 in (G) cannot be
  dropped.
- **Theorem 2.32 ((C) from thirty-one centres on).**  If at least 31 centres
  lie within sqrt6 of c, the union of the caps of Lemma 2.15 has volume
  below 3.0854, so vol(V_c) >= T(Y) > 8.0179.
  - The proof is a two-point kernel certificate in the style of Delsarte,
    with the distance from c as a continuous label: Chebyshev polynomials of
    the distance (degree 4) times the zonal polynomials U_k of S^3 (degree
    12).  If the kernel lies below the pair term on every pair the packing
    allows, the union is at most |Y| m + t/2 with m = -0.0728 < 0, so one
    certificate settles every count from 31 on.
  - Lemma 2.31 gives the pair term in closed form and proves that it
    decreases in the distances and increases in the inner product, so that
    its value at one corner bounds it on a box.
  - The certificate is proved in exact and ball arithmetic: exact LDL^T
    for positivity, the Schur complement made exact, the pair inequality by
    a Bernstein branch and bound (191 740 boxes, 54 248 values of the pair
    term in Arb), the per-point bracket on 356 intervals, and
    31 m + t/2 = 3.08537 < 9 pi^2/8 - 8 = 3.10330; under two minutes
    (multi_cap/radial_count_check.py).
- **Corollary 2.34 (no room beside a near root system).**  If 24 centres
  lie within rho of c and their directions are each within
  1/sqrt2 - a(sqrt6, rho) of a copy of the normalised roots (0.0930 for
  rho = 2.0161), no other centre lies within sqrt6; with Theorem 7.76, 24
  centres within 2/sqrt(1 - 4e-26) leave no room for another within sqrt6.
- **Theorem 2.36 ((G) along the push-outs).**  If the 24 centres lie on the
  rays of the normalised roots at distances 2 + delta_i and T(Y) <= 8, then
  every delta_i < 0.1971 (T on the one-centre ray exceeds 8 there, in ball
  arithmetic) and vol(V(Y) cap K(Y)) >= 8 + sum(delta_i)/5: the gain of
  Lemma 2.23 against a loss of at most 768 sum eta_i^4, the part of V(Y)
  outside the inversion hull lying in 24 vertex pyramids of volume 2 h_v^4.
  Corollary 2.37: the conjecture holds whenever the 24 centres within sqrt6
  lie on the rays of a root system (multi_cap/pushout_check.py).
- **Section 2.10, what (G) and (C) still need** (Remarks 2.33, 2.35 and
  2.38): why kernels on pairs stop at 31 (the same method bounds the kissing
  number of R^4 only by 25), a table of the best two-point bound against the
  largest union found for every count (Table 2), the two cases of (C) from
  25 to 30 (24 centres close to c, which needs the shape half of Theorem
  7.79 with a tolerance of 0.093 per direction; 23 centres close to c, which
  needs a certificate spent on triples), and the localisation and the volume
  bound that (G) needs, with the numbers of the gap.

## In Lean 4

Eight new parts, fifteen in all (Table 5 of the paper):

- lean/D4SecondOrder.lean: H, the cone and c built from the integral root
  system; the certificate of Proposition 2.21 by exact LDL^T; the push-out
  block; the identity of Proposition 2.20; the equality on the one-centre
  rays; the expansions of Lemma 2.23; the tilt block of Proposition 2.25.
- lean/D4HexagonLoop.lean: Proposition 2.26 in
  Q(sqrt3)[c, s]/(c^2 + s^2 - 1): J1 and J2 in closed form,
  3 J1^2 - 4 sqrt3 J2 = 32, the six slice conditions, the Gram matrix at
  theta = 0, and the perfect square.
- lean/D4SecondCode.lean: Proposition 2.28, the 24 rational points, their
  norms, the slack 849/50000 and the inner product 57/250 away from
  -1, -1/2, 0, 1/2, 1.
- lean/D4NearContact.lean: the exact arithmetic of epsilon_0 = 4e-26 in
  Theorem 2.18.
- lean/D4Rigidity.lean: Lemma 7.50 from the integral roots, the ranks, the
  11/16 diagonal and the ratio sqrt(96/11).
- lean/D4Cap.lean: the extremal cap theorem, the 25 vertices of Q,
  vol Q = 25/3, the three ranges of the case |S| = 1 and the 303 boxes.
- lean/certificate/D4Labelled*.lean and D4Omega*.lean: Theorem 2.17, all
  three regions (419 913, 39 551 and 9 545 boxes), with omega, A_* and the
  constants computed inside Lean in outward-rounded 256-bit interval
  arithmetic.
- lean/cardinality/: both certificates of Theorem 7.79 from their entries
  alone (degree 8: 441 intervals and 181 869 boxes; degree 10: 2 957
  intervals, about four hours).

The twelve standalone files pass `run_all.sh`
(lean/runs/run_all_2026-09-27_v1.9.0.log); the three Lake projects build.
No proof uses `sorry`.

## The paper in journal form

- The abstract (94 words) states only what is proved; it now includes the
  count from thirty-one on.
- The introduction, Section 2.8 and the conclusion are rewritten around
  Proposition 2.30.
- The reproducibility sections are replaced by one section, "Computer
  verification" (Section 23), which says what kind of computation each
  claim rests on and lists the formal part (Table 5), and by an appendix,
  "Methods of computation" (Appendix F).  File names of the package no
  longer appear in the text; the data and code availability statement
  points to the archive.
- New figures: Figure 11 (the second order), Figure 12 (the hexagon loop:
  the plane of the turned hexagon, a slice, the slices in three dimensions,
  and the tilt and largest inner product along the curve) and Figure 13
  (the spectra of the two codes, the deepest hole and the least contact
  cell against the slack, and the cases of Proposition 2.30), Figure 14
  (Theorem 2.32: the two-point bound against the count, the kernel of the
  certificate against the pair term, and what settles each count), and
  Table 2.
- The collision test of the TikZ figures is stricter: no two nodes overlap,
  pgfplots tick labels included, and no line, curve, marker or outline
  passes under any text.  All ten TikZ figures pass it
  (paper/figures_new/tikz_overlap_check.log); Figures 4, 6, 8 and 31 are
  redrawn to do so.
- "What is left" (end of Section 2.8) lists four statements, now with the
  reason each is needed; along the ray of one centre pushed out the slack
  reaches delta/4 = 0.049, beyond any statement about the slack alone.
- New reference: Hales and McLaughlin, The dodecahedral conjecture,
  J. Amer. Math. Soc. 23 (2010), 299-344, doi:10.1090/S0894-0347-09-00647-X,
  for why the three-dimensional analogue of the local statement fails; and
  Odlyzko and Sloane, J. Combin. Theory Ser. A 26 (1979), 210-214,
  doi:10.1016/0097-3165(79)90074-8, for the two-point bound on the kissing
  number.
- 177 pages; no overfull boxes and no undefined references.
- Section 20 and 28 other section, subsection and appendix headings are
  retitled; the sections and their numbers are unchanged, and no theorem,
  lemma or proposition number that existed in v1.8.0 changes.

## What is measured, not proved

- The second level at slack 0, degrees (10, 12): 24.9423, in 3.9 hours and
  11.8 GB on four cores and 16 GB (level2/).
- Without the packing constraint the contact-cell volume of 24 unit
  directions goes down to 7.96553, at inner products up to 0.579
  (multi_cap/second_order/contact_cell_scan.py).
- On codes of slack s, the least contact-cell volume found is 8 (the root
  system and the curves of Proposition 2.26) at s = 0.005, 0.01 and 0.02,
  and 7.99802 and 7.99288 at s = 0.03 and 0.035
  (contact_cell_constrained.py).
- (G) without the condition T <= 8 fails: a packing with
  vol(V(Y) ∩ K(Y)) = 4.750 and T(Y) = 10.82.  With the condition, none of 40
  random configurations enters the region T <= 8, and the 20 of 64 starts
  inside it that stay there all end at the root system with volume 8
  (cell_hull_search.py).
- The deepest hole of a 24-point code of slack s, the least
  max_i <theta, w_i>: 0.70711, 0.69280, 0.67817, 0.66322, 0.64409, 0.63238,
  0.58327 at s = 0, 0.004, 0.008, 0.012, 0.017, 0.02, 0.03.  A further
  centre within sqrt6 beside 24 centres within 2.0161 needs 0.6141, so
  numerically there is room for it only from about s = 0.025 on, while
  Theorem 7.79 gives s <= 0.008 (hole24.py).
- A labelled three-point certificate that would exclude that further centre
  does not exist at degree 7: with the 24 directions at slack 0.008 and the
  further one at inner product at most 0.6141, the best certificate has
  value 0 (multi_cap/cap_probe.py, sampled constraints,
  runs/cap_probe_d7.log).  The paper records this beside the earlier probe
  with 24 contacts.
- The best two-point bound on the union of the caps for each count
  (multi_cap/radial_count_sdp.py scan, runs/radial_count_scan.log): 3.330,
  3.360, 3.307, 3.253, 3.199, 3.144 for 25 to 30 centres, all above 3.10330;
  3.080 at 31, falling to 0.174 at 49.  So kernels on pairs stop at 31.
- The deepest hole of the second code of Proposition 2.28 is at inner
  product 0.659, above 0.6141 (floating point, 3000 local minimisations).

## New and changed files

- multi_cap/second_order/: tilt_block.py, tilt_quartic.py, hexagon_loop.py,
  contact_valley.py, contact_cell_scan.py, contact_cell_constrained.py,
  code24_second.py, code24_exact.py (with code24_exact.txt),
  cell_hull_search.py, hole24.py, with logs in runs/; see the README there.
- multi_cap/pushout_check.py (Theorem 2.36, log runs/pushout_check.log) and
  multi_cap/count_core_survey.py (Remark 2.35, log runs/count_core_survey.log).
- multi_cap/cap_probe.py (log runs/cap_probe_d7.log): the labelled
  three-point probe with an inner slack.
- multi_cap/radial_count_sdp.py (finds the certificate of Theorem 2.32, and
  with "scan" the bounds of Table 2), multi_cap/radial_count_check.py (the
  proof), multi_cap/radial_certificates/radial_31.json, and logs
  runs/radial_count_sdp_31.log, runs/radial_count_check_31.log,
  runs/radial_count_scan.log.
- lean/: the eight new parts above, with logs in lean/runs/.
- paper/figures_new/tikz_hexloop.tex, tikz_codes.tex and tikz_count31.tex
  (Figures 12, 13 and 14), with count31_data.py and its two data files.
- CITATION.cff and .zenodo.json: v1.9.0.

## What is left

Conjecture 1.6 follows from (G) and from (C) for 25 to 30 centres; neither
is proved in general, and (G) is proved along the push-outs of the root
system (Theorem 2.36).  (G) is a lower bound for the volume of a polytope over sets of
24 centres near the root system; the numerical searches find nothing below
8.  (C) from 25 to 30 needs a certificate that sees more than pairs; the
least value of T found for those counts is 8.264115, at 25.  Section 2.10
says what each needs and why the methods of the paper stop short of it.
Theorem 1.8 is unchanged.
