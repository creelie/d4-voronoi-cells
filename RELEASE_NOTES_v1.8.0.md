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

## What is measured, not proved

- The three-point bound on A(4, 1/2 + s) rises about 97 per unit of s from
  24.13 at s = 0.  So certificates of the kind behind Theorem 7.79 stop near
  s = 0.009 (sampled programmes, Figure 26(a)).
- 25 points of S^3 with inner products at most 0.53743 exist, a minimal
  angle of 57.49 degrees (code25_search.py).
- Letting the truncation radius of Lemma 2.15 grow with the distances moves
  the pair-term crossing along the root system pushed out evenly only from
  0.1551 to 0.1473 (d4_independent_check.py).

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
- New: independent_verification/rebuilt_from_text/d4_independent_check.py,
  a recomputation of the paper's numerical claims from the statements
  alone.
- Paper:
  - New: Theorem 7.79, Remark 7.78 and Figure 26
    (figures_new/fig_cardinality.py).
  - Revised: Lemma 7.50, Theorem 7.51, Theorem 2.18 (a new proof of its
    step (b)), "What is left", the
    abstract, the introduction, the code index and the data availability
    statement.
  - Figure 16(d) and Figure 25 redrawn; later figures renumbered by one.
  - No theorem, lemma or proposition number that existed in v1.7.0
    changes.

## What remains open

The case left by v1.7.0: a centre with at least twenty-four other centres
within sqrt 6, one of them at a distance between 2 + epsilon_0 and sqrt 6.
Theorem 7.79 supplies the cardinality half of a quantitative form of the
twenty-four-point classification.  The shape half remains open: twenty-four
directions with slack of order 1e-2 lie near a root system.  Remark 7.78
shows it needs a new second-level kernel with a margin.
