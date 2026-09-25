# v1.6.0

Supplementary package for *The Sphere Packing Problem in Dimension 4 and the
Twenty-Four-Cell Conjecture*.  This release proves the local bound for every
centre with at most twenty-three other centres within sqrt 6, whatever their
distances, and narrows Conjecture 1.6 accordingly.  It contains the changes
of v1.5.0 (merged, not archived as a release of its own).  The new code was
developed, and its computations run, with the assistance of Claude, an AI
model made by Anthropic.

## What the paper now proves

Everything proved in v1.5.0 still stands.  New:

- **Theorem 2.17 (at most twenty-three centres within sqrt 6).**  If at most
  23 centres other than c lie within sqrt 6 of c, then vol(V_c) > 8,
  whatever their distances from c and from one another.  With all 23 at
  distance 2 this is Theorem 7.40; the new point is that they need not
  touch c.

  The proof carries the three-point certificate of Theorem 7.73 from
  contacts to centres at mixed distances.  Inside B(sqrt(3/2)) the cell is
  an exact function of the pair terms (Lemma 2.15).  A centre that moves
  out from distance 2 to d hands back the part S(2) - S(d) of its cap, and
  it takes from each pair term only the part of the moving section that
  the other cap covers.  The certificate is allowed to spend 1/11 of what is
  handed back in each triple of centres: each centre lies in 231 = 21 * 11
  triples and each pair in 21.  That pays for every pair of directions
  closer than 60 degrees, which the packing allows only to centres that
  have moved out.  Three steps:

  1. a centre beyond D = 2.1648 settles the case by itself, since
     s(D) = 0.0928820 > 8 - A_* = 0.0928555;
  2. the triple sum of the new quantity Q_ijk is bounded below by the
     certificate lemma, which holds for any 23 unit vectors;
  3. Q_ijk >= Q_0(u, v, t) + Gamma_1 + Gamma_2 + Gamma_3, a function of the
     three inner products and three one-centre terms.  This is checked on
     three ranges of the largest inner product t:
     - t <= 1/2 is the inequality (C) of Theorem 7.73;
     - 1/2 < t <= 0.51 gains at least 0.1056 per unit of t, a new branch
       and bound over 30 051 boxes;
     - 0.51 <= t <= a_D uses the Gamma_i tabulated in ball arithmetic, a new
       branch and bound over 6 482 boxes.
- **Conjecture 1.6 is narrowed.**  What is left is a centre with at least
  24 other centres within sqrt 6, one of them at a distance between
  2 + epsilon_0 and sqrt 6.
- **Why 24 is different.**  At the root system the right side of Lemma 2.15
  is 7.906940 < 8.  The 24-cell reaches out to radius sqrt 2, and
  B(sqrt(3/2)) misses 0.093 of it, in 24 corners.  So with 24 or more
  centres the pair terms alone cannot decide.  Figure 8 draws the corners.
- **Numerical evidence (not proof).**  `truncated_search.py` minimises the
  right side of Lemma 2.15 over configurations of exactly M centres within
  sqrt 6, from twelve starts each.
  - M = 24: every run ends at the root system, 7.906940.
  - M = 25, 26 and 27: every local minimum found lies above 8.26.

  If that is the whole truth, the pair terms fail only at 24 centres near
  the root system, and 25 or more could be settled as in Theorem 2.17, with
  a certificate for each count.

## The manuscript (paper/)

- **Section 2.8.**
  - Theorem 2.17 with its full proof replaces Corollary 2.17 (twenty-three
    contacts), which it contains.  The corollary's covering-radius case now
    follows from Theorem 2.17 and Lemma 2.12; its deletion case is
    Corollary 2.14.  Numbers from Theorem 2.18 on are unchanged.
  - The paragraph "Why the count stops at twenty-three" is new.
  - Step 1 of the proof of Theorem 2.18 now uses Theorem 2.17.
  - The opening of the section and "What is left, and the evidence" are
    revised.
- **Three new figures.**  Figures after Figure 5 are renumbered by three.
  - Figure 6, TikZ (`figures_new/tikz_labelled.tex`): the share fr of a
    section cut by a second hyperplane, in three dimensions, and the three
    ranges of t with what settles each.
  - Figure 7 (`figures_new/fig_labelled.py`, collision test passed): a
    plane section through two centres at the packing boundary, and the
    bound along two lines past t = 1/2.
  - Figure 8, TikZ (`figures_new/tikz_stop.tex`): the section of the 24-cell
    by a coordinate plane, the corners that B(sqrt(3/2)) misses, and one
    corner magnified.
- **Three new footnotes** in Section 2.8: where the weight 1/11 comes from,
  why the certificate alone does not reach past t = 1/2, and the geometry
  of the missing corners.
- **Other revisions.**
  - Theorem 1.7 of the introduction now has three cases, (a) at most 23
    centres within sqrt 6.
  - The abstract and the conclusion are revised to the new claims.
  - The index of notation gains the symbols of Sections 2.7 and 2.8.
  - The code index and the data and code availability statement name the
    new programs.
  - The bibliography entry of the archive cites v1.6.0.
- **Duplicates.**  A repeated attribution sentence in the introduction and
  the conclusion is reworded.  A scan finds no sentence of more than 90
  characters twice in the text.
- **Build.**  The manuscript builds to 162 pages, with no undefined
  references and no overfull or underfull boxes.

## The package

- `multi_cap/labelled_certificate_check.py` is new: the verification of
  Theorem 2.17.  It does the following:
  - rebuilds the certificate and B exactly, and rechecks positivity;
  - evaluates the constants in ball arithmetic (python-flint);
  - reruns the branch and bound of (C);
  - runs the two new branch and bounds.

  It takes about six minutes.  Log: `multi_cap/runs/labelled_certificate_check.log`.
- `multi_cap/truncated_search.py` is new: the floating-point search quoted
  above.  Logs: `multi_cap/runs/truncated_search_M24.log` to `_M27.log`.
- `lean/D4Closure.lean` gains two theorems:
  - `amax_tangent`, the identity behind the tangent-plane bound of the
    packing condition, over every commutative ring (`grind`);
  - `triple_counts`, the counts 21 and 231 = 21 * 11 (`decide`, no axioms).

  `lean/run_all.sh` passes; log `lean/runs/run_all_2026-09-25_v1.6.0.log`.
- `independent_verification/logs/doi_audit_2026-09-25.md` is new: every
  DOI of the bibliography checked against the publisher's or indexer's
  page.
  - 49 of 51 are confirmed.
  - BBB26 is confirmed except for its version suffix `.v11`.
  - LLM24data is consistent with the archive downloaded through it, but
    could not be re-resolved here.

  No DOI points to a different work.
- `README.md`, `lean/README.md`, `independent_verification/REPORT.md`,
  `CITATION.cff` and `.zenodo.json` are updated to v1.6.0.

## Release

Publish the tag `v1.6.0` on `main`.  Zenodo files it under the concept DOI
10.5281/zenodo.22766562, which the paper cites.  The version DOI of v1.4.0 is
10.5281/zenodo.22940045.
