# v1.5.0

*v1.5.0 was merged but not archived as a release of its own; its changes
are part of v1.6.0, where Theorem 2.17 (at most twenty-three centres within
sqrt 6) replaces Corollary 2.17 below and the figures after Figure 5 are
renumbered by three.  See RELEASE_NOTES_v1.6.0.md.*

Supplementary package for *The Sphere Packing Problem in Dimension 4 and the
Twenty-Four-Cell Conjecture*.  This release narrows the case that v1.4.0
left open and adds a second machine-verification layer, in Lean 4, for the
exact content of the new section.  The new code was developed, like that of
v1.4.0, with the assistance of Claude, an AI model made by Anthropic.

## What the paper now proves

Everything proved in v1.4.0 still stands:

- the bound for every contact configuration;
- the bound for the Voronoi cell of a packing whose centres within 2 sqrt 2
  all touch it, or whose distances pass the criterion Phi > 8.

The new Section 2.8 ("Further reductions, and the conjecture near the
contact regime") proves vol(V_c) >= 8, with equality only at the 24-cell,
in four further cases.  The introduction collects them as Theorem 1.7.

- **(a) At most 22 other centres within sqrt 6** (Corollary 2.16).  The
  underlying fact is Lemma 2.15: for every packing, no point of the ball
  B(sqrt(3/2)) is closer to three centres than to the centre.  This makes
  the inclusion-exclusion identity with pair terms exact, and it gives a
  second distance criterion, Psi.
- **(b) Contacts that contain most of a root system** (Corollary 2.14).
  This means all but at most two roots of a copy of R, or, more generally,
  one orthogonal pair of roots at each vertex of its 24-cell.  The proof
  uses Proposition 2.13, the inversion hull: every centre outside a set Y
  leaves the convex hull of 0 and the points 4y/|y|^2 whole.
- **(c) 23 contacts** (Corollary 2.17) that are a deletion of a root
  system, or that have covering radius at most arccos(sqrt 6 / 4) =
  52.24 degrees.  This uses Lemma 2.12: every non-contact centre is at
  distance at least 4 cos(covering radius of the contacts).
- **(d) The contact regime** (Theorem 2.18).  There is an absolute
  epsilon_0 > 0 such that, if no centre lies at a distance between
  2 + epsilon_0 and sqrt 6, then vol(V_c) >= 8, with equality only at the
  24-cell.  Centres beyond sqrt 6 are unrestricted.  So **the D_4
  configuration is a strict local minimum of the cell volume among all
  packings**.  The proof has two steps:
  - a compactness argument, which uses Corollary 2.16 and Theorems 7.25
    and 7.73;
  - a local analysis at the root system, which uses the rigidity spectrum
    of Lemma 7.50 extended to packing slack (||eps|| <= 6 sum delta_i), the
    first-order expansion vol >= 8 + (2/3) sum delta_i - C (sum delta_i)^2,
    and the inversion hull.

  epsilon_0 is not explicit.

**What is still open (Conjecture 1.6):** a centre with at least 23 other
centres within sqrt 6, one of them at a distance between 2 + epsilon_0 and
sqrt 6, whose distances fail both Phi and Psi.  The twenty-four-cell
conjecture for all packings, and with it Delta_4 = pi^2/16, is not proved
unconditionally.  Closing that case would take an explicit epsilon_0, which
needs a quantitative form of Theorem 7.25, and a proof for the
intermediate distances.

## The manuscript (paper/)

- **Section 2.8** is new: Lemmas 2.12 and 2.15, Proposition 2.13,
  Corollaries 2.14, 2.16 and 2.17, and Theorem 2.18.
- **Theorem 1.7** is new in the introduction.  The density bound is now
  Theorem 1.8.
- **Other text** is updated: the abstract, Section 6, the conclusion, the
  data and code availability statement, and the Reproducibility section,
  which now counts eight Lean parts.
- **Two new figures:**
  - Figure 4 (TikZ): the tetrahedron in the ball of radius sqrt(3/2), and
    a hole of the contacts in a plane section.
  - Figure 5 (matplotlib): a vertex of the 24-cell with its three
    orthogonal root pairs, in 3D; the three-dimensional analogue of the
    inversion hull; and the first-order law near the root system.
- **Remark 7.36** now states only what the supplement lets a reader check.
  The floating-point measures (5.9424 on 2170 triples as deposited, 5.9851
  on 522 when rerun) lie below tau by more than 0.18.  That is strong
  evidence, not a proof, that no degree-6 three-point certificate exists.
  The rational measure the remark used to cite is not in the supplement,
  and nothing in the paper depends on the conclusion.
- **Conjecture 1.4:** its heading now says where it is proved (Corollary
  7.31), so the only conjecture left open in the paper is Conjecture 1.6.
- **Bibliography:** all 57 entries are cited, every citation is defined,
  and each of the 51 DOI links matches the DOI printed beside it
  (`independent_verification/check_bibliography.py`).  With `--online` the
  script also compares each DOI's registered title and year with the
  entry; that needs access to doi.org, api.crossref.org and
  api.datacite.org.
- **Title page:** the long provenance note is removed, and Deep
  Bhattacharjee (dagger) and Priyabrata Mandal (double dagger) are both
  corresponding authors.  The earlier preprint BBB26 is now mentioned in
  one sentence of the introduction.

## The package

- `multi_cap/closure_lemmas.py`: the exact content of Section 2.8, in
  sympy, integer and ball arithmetic.  Log: `multi_cap/runs/closure_lemmas.log`.
- `multi_cap/near_contact_probe.py`: floating-point exploration near the
  root system.  Log: `multi_cap/runs/near_contact_probe_200_seed5.log`.
- `lean/D4Closure.lean`: the same exact content in Lean 4, with no Mathlib
  and no `sorry`.
  - The four identities are proved over every commutative ring with
    `grind`.
  - The root-pair combinatorics is proved with `decide +kernel`, with no
    `native_decide`.
  - `lean/run_all.sh` now checks six files; the log is
    `lean/runs/run_all_2026-09-25.log`.
- `paper/figures_new/fig_closure.py` and `paper/figures_new/tikz_closure.tex`,
  with `paper/figures/fig_closure.png` and `paper/figures/fig_tikz_closure.png`.
- `CITATION.cff` and `.zenodo.json` carry version v1.5.0.

## Release

Publish the tag `v1.5.0` on `main`.  Zenodo files it under the concept DOI
10.5281/zenodo.22766562, which the paper cites.  The previous version DOI,
of v1.4.0, is 10.5281/zenodo.22940045.
