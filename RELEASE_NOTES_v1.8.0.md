# v1.8.0

Supplementary package for *The Sphere Packing Problem in Dimension 4 and the
Twenty-Four-Cell Conjecture*.  This release proves that the kissing number in
dimension four is stable under a small slack, sharpens the isolation radius
of the root system and the explicit neighbourhood of the near-contact
theorem, and says exactly why the robust reading of the LLM24 kernel cannot
reach moderate slack.  Conjecture 1.6 remains open, and the paper says so.
The new code was developed, and its computations run, with the assistance of
Claude, an AI model made by Anthropic.

## What the paper now proves

Everything proved in v1.7.0 still stands.  New or sharper:

- **Theorem 7.79 (the kissing number is stable).**
  - No 25 points of S^3 have pairwise inner products at most 1/2 + 0.008.
  - So at most 24 other centres of a unit-ball packing of R^4 lie within
    2/sqrt(1 - 0.016) = 2.0161 of any centre.
  - The proof is a Bachoc-Vallentin certificate of degree 10 on the enlarged
    domain [-1, 1/2 + 0.008], verified in exact rational and outward-rounded
    interval arithmetic (2.7 million boxes, under an hour).
  - Degree-8 certificates prove the same at slack 0.005 and 0.0065 in about
    80 seconds each.
- **Lemma 7.50 and Theorem 7.51, sharper.**
  - Every pair of T carries weight 11/16 in the image of the rigidity
    operator, since the Weyl group of F4 is transitive on the edges of the
    24-cell.
  - So ||g||_2 <= (sqrt 11/4) ||g||_1 on the image, and the isolation
    radius of the root system rises from 2/sqrt(577) = 0.08326 to
    2/sqrt(397) = 0.10038.
  - No constant of this argument can pass 0.244.  A single displaced
    direction has ||Lambda tau||_1 / ||tau||_2 = sqrt(96/11).
- **Theorem 2.18, second statement.**
  - The explicit neighbourhood grows from sum delta_i <= 4e-5 to 4e-3,
    a factor of 100.
  - Estimate (a) is sharper.  The positive and the negative parts of
    Lambda tau are bounded separately, which gives ||eps|| <= 2.38 (1 + delta) S
    in place of (384/71) (1 + delta) S.
  - Step (b) is new.  The derivative of the volume along the linear path is
    written exactly, facet by facet.  Each facet is compared, face by face,
    with the regular octahedron of the 24-cell.  The losses are summed over
    the 8-regular edge graph of the 24-cell, so every loss is of second
    order with no factor of the largest perturbation.
  - The bracket is at least 0.1 at 4e-3 and positive up to 4.5e-3.
  - epsilon_0 = 4e-26 is unchanged.
- **Remark 7.78 (the ceiling of Theorem 7.76).**
  - Even with the triple and quadruple terms set to zero, the robust reading
    pins the inner products only for kappa < 5.2e-7 delta^2.
  - So it cannot pass 4.7e-14 at the tolerance of Lemma 7.77, nor 9.2e-9
    under any rounding lemma.
  - A form of the classification at slack 1e-2 needs a kernel computed on
    the enlarged domain.

## Checked against the authors' own code

- level2/ writes the 490 zonal matrices of zonal/ in the file format of the
  code of de Laat, Leijenhorst and de Muinck Keizer, so that their
  second-level programme runs without their own construction of them (three
  days and 128 GB by their README).
- Their evaluate_zonal_matrix, reading these files, agrees exactly with
  zonal/zonal.py on eight entries up to |lambda| = 14.
- Their compute_PS, run from scratch one signature at a time, writes entries
  identical to ours, coefficient by coefficient: all 88 entries of the
  thirteen signatures run (|lambda| <= 10, lambda_2 <= 4).  Until now the
  zonal matrices had been checked only against independent facts.

## The second order at the root system (Section 2.9)

- Lemma 2.19: the second variation of the volume of a 4-polytope under
  moving normals and support numbers, and that the volume is C^2 (C^3 in
  fact) across changes of combinatorics.
- Proposition 2.20: the Hessian H of the cell volume at the root system, in
  closed form over the 96 triangles of the 24-cell; on the push-outs it is
  2 sum_edges eta_i eta_j - 4 sum eta_i^2.
- Proposition 2.21: on the first-order packing cone, H >= -(sum delta_i)^2,
  with equality on the ray of one centre pushed out.  Proved by an exact
  copositivity certificate H + cc^T = P + B^T N B, checked in exact
  rationals in three seconds.  So the tilts cost nothing more at second
  order than the pushes, where step (b) of Theorem 2.18 charges them
  73 S^2 to 98 S^2, and the second-order model vol - 8 >= (2/3) S - (1/2) S^2
  stays positive up to S = 4/3.
- Lemma 2.23 and Proposition 2.24: pure push-outs integrated exactly, by
  Brunn-Minkowski.
- Remark 2.22: what is still missing for statement (ii): tilting is not
  monotone (untilt.py), the third-order remainder is not bounded (on 300
  tilted packings vol - 8 - (2/3) S + (1/2) S^2 >= 0 holds in every case,
  remainder_test.py, floating point), and the part of V(Y) outside the
  inversion hull needs a sharper bound than 13958 Theta^4 (it is below 1e-6
  along both rays up to 0.197, multi_cap/hull_along_rays.py).
- The list of what is left (end of Section 2.8), corrected, now has
  four items: along the ray of one centre pushed out the slack reaches
  delta/4 = 0.049, beyond any statement about the slack alone (25-point codes
  exist at 0.0374), so twenty-three contacts with one near-contact need a
  classification of their own.
- Checked independently (multi_cap/second_order/independent_check.py, every
  volume recomputed by intersecting halfspaces): H against second
  differences (6e-8), the cone rows against exact distances (2e-9), the
  spectrum, the one-centre formula (1e-14), Lemma 2.23 on 400 push patterns.
  Two statements of the draft are corrected: the eigenvalue -6.194 has
  multiplicity 8, and the push-out block is Adj - 4I.
- In Lean: lean/D4SecondOrder.lean builds H, the cone and c from the
  integral root system inside the proof assistant, reads only the 43 orbit
  values of the certificate, and proves the certificate of Proposition 2.21
  (exact LDL^T, 30 zero pivots with zero rows), the push-out block and the
  identity of Proposition 2.20, the equality on the one-centre rays, and the
  expansions of Lemma 2.23 (grind).  About thirty seconds.
- In Lean, five more parts, thirteen in all:
  - lean/cardinality/, a new Lake project that checks the certificates of
    Theorem 7.79 (degree 8 at s = 0.0065, degree 10 at s = 0.008) from their
    entries alone in exact dyadic arithmetic: positivity by exact LDL^T, the
    bound, the inequality at |C| = 25, and both branch and bounds (441 and
    2 957 intervals; 181 869 boxes for degree 8, a run of several hours for
    degree 10).
  - lean/certificate/D4Labelled*.lean: the regions II_s and II_f of Theorem
    2.17 (39 399 and 9 545 boxes), with the shares, the packing bounds and the
    Gamma envelopes computed in Lean; only the omega tables, the slab constant
    and the cell integrals come from interval arithmetic.
  - lean/D4NearContact.lean: the exact arithmetic of epsilon_0 = 4e-26 in
    Theorem 2.18: -q >= 2.62e-4 for the quotient of the two-point polynomial,
    f >= 1.12e-8 off the windows, the constants of Lemma 7.77, the assembly of
    E(24, kappa) and E(25, kappa) from the bounds B_3, B_4, and the facet
    bracket; twelve theorems, twenty seconds.
  - lean/D4Rigidity.lean: Lemma 7.50 from the integral roots, the polynomial
    annihilating 4N, the ranks over Q (so the multiplicities 30, 29, 8, 21, 8),
    the 11/16 diagonal of the projector onto the image, and the ratio
    sqrt(96/11) of one displaced direction; the proof of the lemma said this
    was done in Lean's kernel, where only its input, the 96 tight pairs, was.
  - lean/D4Cap.lean: the extremal cap theorem: the 25 vertices of Q by exact
    enumeration, vol(24-cell) = 8 and vol(Q) = 25/3 by exact triangulation
    (the script had them in floating point), the cross-polytope enclosure,
    the fourth derivative of g, the three ranges of the case |S| = 1, and the
    303 boxes of the case |S| >= 2 with the same largest value 0.99755.
  - lean/certificate/D4Omega*.lean: the tables of omega, A_* and the
    constants of Theorem 2.17 computed inside Lean, in outward-rounded 256-bit
    interval arithmetic, from simplified closed forms,
    omega'(u) = (pi/16)(3u - 1)^2 / ((1 + u)^2 sqrt(1 - u^2)) and
    omega(u) = 4 pi [(9/32) arctan((t* - tau)/(1 + t* tau)) - (4 tau^3 - 24 tau + 11 sqrt2)/96]
    (new; multi_cap/omega_closed_form.py checks them symbolically and against
    the old form).  A_* > 8 - 0.0928555703, s(D) > 8 - A_*, the slab constant
    and the 256 cell integrals are confirmed, the covering bound exceeds 8
    for every m <= 22 and (9/8) pi^2 - 22 S(2) > 8.046, and the branch and
    bounds of regions II_s and II_f are rerun with the Lean tables (39 551 and
    9 545 boxes): they now take nothing from outside Lean but the
    certificate.  D4OmegaRegionI.lean states region I with them.
  Closed definitions that run a check now live in the Main modules
  (lean/certificate/D4CertMain.lean holds verifyDomain), since a precompiled
  module evaluates its closed definitions when it is loaded; D4CertStat now
  takes one run instead of two.

## What is measured, not proved

- The three-point bound on A(4, 1/2 + s) rises about 97 per unit of s from
  24.13 at s = 0.  So certificates of the kind behind Theorem 7.79 stop near
  s = 0.009 (sampled programmes, Figure 29(a)).
- 25 points of S^3 with inner products at most 0.53743 exist, a minimal
  angle of 57.49 degrees (code25_search.py).
- Letting the truncation radius of Lemma 2.15 grow with the distances moves
  the pair-term crossing along the root system pushed out evenly only from
  0.1551 to 0.1473 (d4_independent_check.py).
- The right side of Lemma 2.15 for every count M of centres within sqrt 6
  (count_survey.py, Figure 9):
  - the least value found is 8.264115 at M = 25, rising to 10.895 at M = 43;
  - all 196 local minima found for M >= 25 lie at or above 8.264;
  - every one of them has centres on the sphere of radius sqrt 6, whose caps
    are empty;
  - no packing of 44 to 49 centres was found in 300 random starts each.
- The second level at slack 0 on a machine with four cores and 16 GB
  (level2/, Figure 10):
  - degrees (4, 6): 32;
  - degrees (8, 10): 26.0000, in 42 minutes and 3.5 GB;
  - degrees (10, 12): 24.9423, in 3.9 hours and 11.8 GB, the most that
    fits in 16 GB;
  - the certificate's degrees (14, 16), where the bound is 24, need a larger
    machine for the programme itself.
- level2/las2_margin.jl is the programme with a margin in the two-point
  constraint that a classification of 24 points at positive slack needs;
  tested at degrees (4, 6).

## Corrections

- The pair bound at the root system truncated at r_* is 7.906940, not
  7.907070, in the paragraph after the second-order proposition.  The value
  7.906940 was already printed in Section 2.8.

## New and changed files

- New in multi_cap:
  - certify_cardinality.py, the proof of Theorem 7.79.  It takes the exact
    decimal threshold and runs its boxes to the least double at or above it.
  - cardinality_sdp.py, the search and the sweep.
  - cardinality_certificates/, three certificates.
  - robust_ceiling.py and code25_search.py.
  - facet_bounds_probe.py, a numerical check of the new step (b).
  - Logs in runs/.
- Changed in multi_cap:
  - rigidity_spectrum.py adds the exact 11/16 check and sqrt(96/11).
  - explicit_eps0.py uses the new constant of estimate (a) and the new
    step (b).
- New: multi_cap/second_order/ (Section 2.9; README.md there), with logs in
  its runs/.
- New: level2/ (above) and multi_cap/count_survey.py, with logs in
  level2/runs and multi_cap/runs/count_survey*.log.
- New: paper/figures_new/tikz_overlap_check.py, the collision test of the
  TikZ figures, run on all six of them in two modes (log:
  tikz_overlap_check.log).  No two nodes overlap, pgfplots tick labels
  included, and no line, curve, marker or outline passes under any text;
  a label may sit on a fill or a smooth shading.  Figures 4, 6, 8 and 28
  had labels touching lines or balls; they are redrawn.
- New: independent_verification/rebuilt_from_text/d4_independent_check.py,
  a recomputation of the paper's numerical claims from the statements
  alone.
- Paper:
  - New: Theorem 7.79, Remark 7.78 and Figure 29
    (figures_new/fig_cardinality.py); Section 2.9 with Figure 11; Figures 9
    and 10 (figures_new/tikz_counts.tex, tikz_level2.tex and
    tikz_secondorder.tex, TikZ, each passing
    figures_new/tikz_overlap_check.py).
  - The abstract is rewritten.
  - Revised: Lemma 7.50, Theorem 7.51, Theorem 2.18 (a new proof of its
    step (b)), "What is left", the
    abstract, the introduction, the code index and the data availability
    statement.
  - Figure 19(d) and Figure 28 redrawn; later figures renumbered.
  - No theorem, lemma or proposition number that existed in v1.7.0
    changes.

## What remains open

The case left by v1.7.0: a centre with at least twenty-four other centres
within sqrt 6, one of them at a distance between 2 + epsilon_0 and sqrt 6.
"What is left" in Section 2.8 now lists four statements that would close it:

1. the shape of 24 directions at slack of order 1e-2 (proved at 2e-26;
   its cardinality half, Theorem 7.79, at 0.008);
2. the volume for 24 centres near the root system between Sum delta_i = 4e-3
   and 0.155: the inequality vol V(Y) - 8 >= (2/3) S - (1/2) S^2, certified
   at second order in Section 2.9, proved for pure push-outs and true on
   every packing tested, together with a bound for the part of V(Y) outside
   the inversion hull.  Part (c) of Theorem 2.18 bounds that part by
   13958 Theta^4; along both rays it is below 1e-6 up to 0.197
   (hull_along_rays.py), so the bound, not the hull, is what is weak;
3. twenty-three contacts with one near-contact, a classification of codes
   of its own (the list of what is left in Section 2.8);
4. certificates for 24 centres away from the root system and for each
   count from 25, which need 1 and 3.

Conjecture 1.6 is not proved, and Theorem 1.8 is unchanged.
