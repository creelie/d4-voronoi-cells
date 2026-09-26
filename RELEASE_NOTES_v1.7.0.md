# v1.7.0

Supplementary package for *The Sphere Packing Problem in Dimension 4 and the
Twenty-Four-Cell Conjecture*.  This release makes the constant epsilon_0 of
the near-contact theorem explicit: epsilon_0 = 4e-26.  Conjecture 1.6 remains
open, and the paper says so.  The new code was developed, and its
computations run, with the assistance of Claude, an AI model made by
Anthropic.

## What the paper now proves

Everything proved in v1.6.0 still stands.  New:

- **Theorem 2.18, with numbers.**
  - If no centre lies at a distance from c strictly between
    2 + 4e-26 and sqrt 6, then vol(V_c) >= 8, with equality only at the
    24-cell.
  - The same holds whenever exactly 24 centres lie within sqrt 6, their
    directions are within 1/48 of the normalised roots (root sum of squares
    after the best rotation and matching), and the sum of the distances
    delta_i = d_i - 2 is at most 4e-5.  This neighbourhood uses no
    certificate.
- **Theorem 7.76 (twenty-four points, approximately).**
  - A set of directions whose pairwise inner products are at most
    1/2 + 2e-26 has at most 24 points.
  - When it has 24, each inner product lies within 3e-4 of
    {-1, -1/2, 0, 1/2}, and the set lies within 0.0113 of a copy of the
    normalised roots.

  The proof reads the certificate robustly.  In its sum-of-squares
  identities only the terms weighted by products of the factors
  (u + 1)(1/2 - u) can turn negative when an inner product u exceeds 1/2,
  and their size is bounded from the deposited data.
- **Lemma 7.77.**  If the Gram matrix of 24 unit vectors is within 3e-4 of
  a matrix g with entries in {-1, -1/2, 0, 1/2}, then g is the Gram matrix
  of the root system.  The reason: every 5 x 5 minor of the integral matrix
  2g is an integer of absolute value below 1, hence 0.  A Davis-Kahan
  estimate then bounds the distance from the root system.

Why epsilon_0 is so small: the bound for the sign-changing terms is about
1.2e14 kappa, and it has to stay below the value 2.9e-12 of the two-point
polynomial at distance 3e-4 from its double zeros.  The explicit
neighbourhood of the second statement is not small in this way.

## What is left

Conjecture 1.6 is left at a centre with at least 24 other centres within
sqrt 6, at least one of them at a distance between 2 + epsilon_0 and
sqrt 6.  Two pieces of work would close it.

- **25 to 49 centres.**  At most 49 fit within sqrt 6 (v1.6.0).  Numerically
  the truncated volume stays above 8.26 there.
- **Exactly 24 centres, at intermediate distances.**  This needs a lower
  bound for the part of the cell outside B(sqrt(3/2)), strong enough to
  reach the explicit neighbourhood of Theorem 2.18.  The size of the gap is
  measured by `truncated_search.py rays` (floating point):
  - the pair terms reach 8 only at sum delta_i = 0.155 when all 24 centres
    are pushed out evenly, and at delta = 0.197 when one centre is pushed
    out;
  - the explicit neighbourhood ends at sum delta_i = 4e-5.

  That is a factor of about four thousand.

## The manuscript (paper/)

- **Theorem 2.18** is restated with epsilon_0 = 4e-26 and the explicit
  neighbourhood.  Its proof is rewritten.
  - Step 1 uses Theorems 2.17 and 7.76 in place of compactness.
  - Step 2 gives the facet, moment and hull estimates with explicit
    constants.
- **Section 7.13 is new**, with Theorem 7.76, Lemma 7.77 and a footnote on
  the size of the constants.  It sits at the end of Section 7, so no
  theorem, lemma or section number changes.
- **One new figure.**  Figure 25 (TikZ, `figures_new/tikz_eps0.tex`) shows
  the chain of explicit bounds from kappa to the volume, and, to scale,
  the lower bound of sigma_2 near its double zero against the error level
  E.  Figures 25 to 32 of v1.6.0 become 26 to 33.
- **Other revisions.**
  - The abstract, Theorem 1.7(c) and the remark after it.
  - The opening of Section 2.8, "What is left" and the conclusion.
  - The index of notation and the code index.
  - The data and code availability statement names only release v1.7.0 and
    the archive DOI 10.5281/zenodo.22766562, with no version history; the
    bibliography entry of the archive likewise.
- **Build.**  166 pages, with no undefined references and no overfull or
  underfull boxes.

## The package

- `multi_cap/explicit_eps0.py` is new.  From the deposited certificate
  (`third_party/llm24-certificate`, folder `proofs/4_24`) it does the
  following, in sixteen seconds:
  - bounds the sign-changing sum-of-squares terms in ball arithmetic;
  - divides sigma_2 exactly by its zeros and bounds the quotient below;
  - checks the minor bounds of Lemma 7.77;
  - evaluates the estimates of Theorem 2.18 in exact arithmetic.

  Log: `multi_cap/runs/explicit_eps0.log`.
- `multi_cap/truncated_search.py` gains the mode `rays` (log
  `runs/truncated_search_rays.log`), quoted in "What is left".
- `multi_cap/three_point_probes.py` is new (exploration; log
  `runs/three_point_probes.log`).  Three-point certificates cannot force the
  pair inner products of a 24-point code into windows about the D4 values,
  and at degree 7 they cannot exclude 24 contacts plus one centre at sqrt 6.
  So the missing statement needs the second level of the hierarchy.
- `zonal/level2_numeric.py` and `zonal/level2_sampled.py` are new
  (exploration; logs in `zonal/runs`).  The first evaluates the second-level
  kernel of the certificate in floating point and reproduces it exactly
  (24, -1, sigma_2 to 2e-14, zeros on the root system).  The second samples
  the second-level programme: with 8000 quadruples at slack 0 it gives 19.42
  in place of 24, and it uses about 10 GB.  A faithful new second-level
  certificate is out of reach of this machine.
- `README.md`, `independent_verification/REPORT.md` (Addendum: v1.7.0),
  `CITATION.cff` and `.zenodo.json` are updated to v1.7.0.  The descriptions
  in `CITATION.cff` and `.zenodo.json`, which Zenodo shows on the record,
  describe the package as it stands, without a version history.

## Release

Publish the tag `v1.7.0` on `main`.  Zenodo files it under the concept DOI
10.5281/zenodo.22766562, which the paper cites.
