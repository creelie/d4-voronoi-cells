# v1.12.0

Supplementary package for *Voronoi Cells of Four-Dimensional Unit-Ball
Packings* (formerly *The Sphere Packing Problem in Dimension 4 and the
Twenty-Four-Cell Conjecture*).  The paper proves that a Voronoi cell of a
unit-ball packing of R^4 has volume at least 8, with equality only for the
twenty-four-cell of D_4, when at most 23 or at least 31 centres lie within
sqrt 6 of its centre, or none lies strictly between 2 + 4e-26 and sqrt 6.
The other cells reduce to two explicit statements, (G) and (C) for 25 to 30
centres, which together would imply that D_4 is the densest packing in
R^4.  Both remain open in this release.

## The paper

- New title and a short abstract that states exactly what is proved.
- The results on neighbours that do not touch the centre form Section 21,
  after the contact theory they depend on.
- Results of Sections 15 to 18 obtained by finite differences in double
  precision are stated as numerical observations; only exact computations
  are propositions.  The numbering is unchanged.
- One voice throughout: hedging, filler and informal phrasing removed.
- 182 pages; the source is ASCII and compiles cleanly.

## New exact results, checked in Lean

- Corollary 21.27, parts (ii) and (iii): 24 centres within
  rho = 2/sqrt(1 - 0.016) leave no room for a further centre within sqrt 6
  when their directions lie within sqrt 6 h = 0.22797 of a rotated root
  system in root-sum-square, or when the squared deviations of their
  pairwise inner products from those of the roots sum to at most 0.3059.
  `lean/D4HoleBudget.lean` checks the identities over every commutative
  ring and the rational bounds and the 200-interval table by kernel
  computation.
- Proposition 21.28: the same conclusion from the design defects alone,
  (103 S_1 + 191 S_2 + 280 S_3 + 230 S_4 + 159 S_5)/963 < 0.2678, with
  S_k the sum of G_k over all pairs of directions; no labelling is needed.
  `lean/D4DesignBudget.lean` checks the polynomial, its sign on
  [-1, 0.614039] (29 intervals), the bound and the 5-design identities of
  the roots.
- Proposition 18.5 and the congruence of the configurations grown from
  roots 5 and 12 with A_18: `lean/D4A18.lean`.

## What the new evidence says about (C)

Each of the three forms above would settle the case of (C) with 24 centres
within 2.0161, given the corresponding statement about 24-point codes of
slack 0.008.  Local search suggests that the statements hold with room:
the worst codes found reach 0.138 in the sum of squared deviations,
against 0.3059, and 0.0905 in the design defect, against 0.2678.  Proving
either is the missing piece.

- The second level of the Lasserre hierarchy at the degrees (14, 16) and
  slack 0.008 now runs on a 16 GB machine (`level2/fast/`): 86 iterations
  of about a quarter of an hour, 10.3 GB.  Its bound on the number of
  points is 24.555 (floating point, relative gap 1.3e-5), and its two-point
  polynomial charges almost nothing to any single pair.  So the plain
  programme cannot feed any of the three forms.
- The programme with exactly 24 points and the design defect as its
  objective does not either: its moment side reads 1.44 (floating point, a
  run stopped at relative infeasibility 3.6e-4), with S_5 near 6.4.  For
  its moments, no test polynomial of degree up to 11 excludes a further
  centre.
- Remark 21.29 records these values as evidence only; no certificate of
  the second level is rounded and checked.

The case with 23 close centres and (G) are unchanged: both need a
localisation that sees distances as well as directions, which no known
method provides (Remarks 21.29 and 21.37).

## New and changed files

- `lean/`: D4HoleBudget.lean, D4DesignBudget.lean, D4A18.lean and their
  generators; `run_all.sh` checks fifteen files, log in `lean/runs/`.
- `level2/fast/`: the chunked build, the triple-double Schur complement,
  the patched solver with checkpoints and warm starts (its licence
  included), the programme with exactly 24 points, and the analysis
  scripts `analyze_cert.jl`, `diag_fc_primal.jl` and `hole_test.py`, with
  all run logs in `runs/`.
- `level2/EXTERNAL_RUN.md`: what a larger machine would still add.
- `multi_cap/stability/`: the searches behind the tolerances above.
- `hessian_multidir/a18_other_starts.py`: the growth from roots 5 and 12.
- The file-by-file maps moved into the README files of `cap_certificate/`,
  `arc2_w1v2/`, `hessian_multidir/` and `multi_cap/`; release notes moved
  into `release_notes/`; the two legacy figure-composition scripts removed.
- `paper/D4.tex`, `paper/D4.pdf` (182 pages).
- `README.md`, `CITATION.cff` and `.zenodo.json`: the new title and v1.12.0.
