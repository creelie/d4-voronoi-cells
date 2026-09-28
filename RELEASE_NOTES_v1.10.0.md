# v1.10.0

Supplementary package for *The Sphere Packing Problem in Dimension 4 and the
Twenty-Four-Cell Conjecture*.  This release proves three new results on the
two statements to which the paper reduces Conjecture 1.6: the second-level
computation with a margin that earlier releases proposed for a large
machine cannot help at any size; twenty-four of the centres lie within
2.444 wherever (C) could fail from 25 to 33 centres; and in (G) at most two
of the twenty-four centres lie beyond 2.1, and at most one beyond 2.15.  It
also completes the table of two-point bounds for every count of centres,
corrects the values quoted for the three-point bound on the kissing number,
and adds seven references.
Conjecture 1.6 remains open.  (G), and (C) for 25 to 30 centres, are still
the two statements that would close it, and the paper says what each still
needs.  The new code was developed, and its computations run, with the
assistance of Claude, an AI model made by Anthropic.

## What the paper now proves

Everything proved in v1.9.0 still stands.  New:

- **Proposition 2.19 (the root system spends the margin).**  The
  second-level programme with a margin mu w(u) in its two-point constraint,
  w(u) = (u + 1)(u + 1/2)^2 u^2 (1/2 + s - u), yields
  w(u) <= (N - 24)/mu for every inner product of a 24-point code of slack s.
  The normalised roots of D4 are such a code, and w = 3s/8 on each of their
  96 pairs at 1/2, so (N - 24)/mu >= 36 s whatever the degree.  At s = 0.008
  that is 0.288, while max w < 0.019779 (ball arithmetic on 2 * 10^5
  subintervals), so the programme excludes no inner product.  The crossing is
  near s = 5.2e-4, far below the slack 0.008 of Theorem 7.79.  The
  (14, 16) run that level2/EXTERNAL_RUN.md planned for a 192 GB machine
  would therefore give nothing towards (C) or (G), and it is withdrawn.
  Remark 2.20 records why a weight that also vanishes at 1/2 cannot do
  better: at slack 0.008 there are 24-point codes with an inner product of
  0.1231, at least 0.123 from every inner product of the root system (floating
  point).
- **Lemma 2.38 (twenty-four centres within 2.444).**  If 25 to 33
  centres lie within sqrt 6 and T(Y) <= 8, then at least 24 of them lie
  within 2.444, and the 24 closest satisfy T(Z) <= 8 + (M - 24) S(d_25).
  The proof combines the subset inequality T(X) <= T(Y) + sum S(|y|) with
  Theorem 2.17, whose Steps 1 and 2 give T(X) >= 8 + 2.6491e-5 for any
  23 centres, against at most 10 S(2.444) = 2.5325e-5 for the others; the
  constants are in ball arithmetic at 200 bits
  (multi_cap/twentyfour_close_check.py).  At thirty centres, Remark 2.40
  now reads T <= 8.00368 for the 24 closest, in place of 8.005.
- **Proposition 2.43 (where the twenty-four centres lie).**  If exactly 24
  centres lie within sqrt 6 and T(Y) <= 8, then at most two lie beyond 2.1
  and at most one beyond 2.15.  These are the case certificates of Theorem
  2.39 run at M = 24 with one split each; their exact values, 3.07883 and
  3.09251, lie below 9 pi^2/8 - 8 = 3.10330.  The caps alone allow three and
  two.  At the next counts (at most one beyond 2.1, none beyond 2.2 or 2.25)
  the programme of the same degree already exceeds the level, at 3.1248,
  3.1422 and 3.1306 in floating point.
- **Table 2 covers every count from 25 to 49.**  The best two-point
  bound on the union of the caps for each count (degree 14 in the angle and 5
  in the distance), the largest union found and the least T found; the
  counts 33, 35, 36, 38, 39, 41, 42, 44, 45, 47 and 48 are new, and Figure
  15(a) is redrawn with all of them.  The bound falls below 3.10330 from
  31 on and reaches 0.174 at 49.

## Corrections

- The paper credited Bachoc and Vallentin with a three-point bound of 24.10
  on the kissing number of R^4.  Their paper gives 24.5797; later
  computations lowered it to 24.0569 at degree 16 (Machado and Oliveira) and
  24.0472 at degree 20 (Leijenhorst).  All four places are corrected, and the
  value 24.13 of this project's own computation is marked as the bound at
  degree 10.
- Remark 2.37: the "room of at least 0.35" for a three-point route in
  the case of 23 close centres exists only over real configurations.
  Relaxations admit a root system with one centre pushed out and a
  fictitious centre near sqrt 6 whose union of caps reaches 9 pi^2/8 - 8,
  and neither the kernel of Theorem 2.34 nor a direction-only three-point
  kernel excludes it.  The remark now says so and points to Lemma
  2.38.
- The footnote near Figure 11, Remark 2.37 and level2/README.md
  no longer present the second level with a margin as the route to (C);
  level2/EXTERNAL_RUN.md and the header of level2/las2_margin.jl carry the
  floor 36 s.

## The paper

- Abstract rewritten (60 words): it states exactly what is proved.
- New figures: Figure 9 (the weight of the margin programme against the
  floor 36 s, and max w against 36 s as functions of the slack), Figure
  16 (the cap volume S(d) with the radii of Lemma 2.38, and the
  caps of the far centres against the room of Theorem 2.17), and Figure 17 (the counts of
  Proposition 2.43 against the cap budget, and the value of each certificate).  All
  TikZ figures pass the collision test (paper/figures_new/tikz_overlap_check.py):
  no two labels overlap and no line, curve or marker passes under any text.
- New references, each DOI taken from the reference lists of the papers
  themselves: Altschuler and Pérez-Garrido (arXiv:1301.4884), Böröczky and
  Glazyrin (arXiv:1711.06012), Cohn, de Laat and Salmon (arXiv:2206.15373),
  Dostert, de Laat and Moustrou (SIAM J. Optim. 31 (2021),
  doi:10.1137/20M1351692), Leijenhorst's thesis (TU Delft 2025,
  doi:10.4233/uuid:91af805a-376c-4ef8-aec5-e6ce08ae20a7), Machado and
  Oliveira (Experiment. Math. 27 (2018), doi:10.1080/10586458.2017.1286273)
  and Mittelmann and Vallentin (Experiment. Math. 19 (2010),
  doi:10.1080/10586458.2010.10129070).  CVXPY, Clarabel and Gorin and López,
  listed before but never cited, are now cited where they are used.
- 182 pages; no overfull boxes, no undefined references, no missing glyphs,
  all fonts embedded.

## DOIs

paper/tools/verify_dois.py resolves every DOI of the bibliography through
Crossref or DataCite and compares the registered title with the cited one.
The build container cannot reach doi.org, Crossref or DataCite directly.  Of
the 59 DOIs, 27 appear verbatim in the reference lists of the papers in the
project, 10 more were resolved through the Crossref API with matching titles,
the five arXiv DOIs follow from the arXiv identifiers, and two are the
authors' own (Preprints.org and Zenodo).  The remaining 15 were confirmed by hand
against the publishers' records.  Run the script once more on a machine with
network access before submission.

## New and changed files

- level2/margin_floor_check.py: Proposition 2.19, exactly and in ball
  arithmetic.
- multi_cap/twentyfour_close_check.py: the constants of Lemma 2.38.
- multi_cap/radial_certificates/case_24_r2.1.json, case_24_r2.15.json and
  their specs, with runs/radial_case_{sdp,check}_24_r2.1.log and
  runs/radial_case_{sdp,check}_24_r2.15.log: Proposition 2.43.
- multi_cap/radial_case_sdp.py: a spec key "antipodal" that refines the
  sampled pairs along u = -1, where the caps are disjoint.  Without it the
  programme can return a kernel slightly positive at some antipodal pairs,
  which the exact branch and bound then refuses.
- multi_cap/runs/radial_count_scan.log: the eleven new counts of Table 2.
- paper/tools/verify_dois.py: the DOI check.
- paper/figures_new/: tikz_margin.tex with margin_data.py, tikz_lemma.tex
  with lemma_data.py, tikz_radial.tex with radial_data.py, and their data files; tikz_count31.tex,
  tikz_codes.tex and tikz_hexloop.tex redrawn.
- Comments of the scripts, the Lean files and the README index renumbered to
  the paper.
- CITATION.cff and .zenodo.json: v1.10.0.

## What is left

Conjecture 1.6 follows from (G) and from (C) for 25 to 30 centres; neither
is proved in general.  The work of this release makes precise what is
missing, and that it is mathematics and not computation:

- (C) for 25 to 30 centres with 24 of them close to c needs a per-direction
  localisation: each of 24 directions of slack 0.008 within 0.093 of a
  rotated root system (Corollary 2.36).  Floating-point searches find at most
  0.064, so the statement looks true, but no known bound sees single
  directions, and the second level with a margin provably cannot.
- (C) with only 23 centres within 2.0161 needs a certificate that tells
  real configurations from the fictitious ones above.
- (G) needs a localisation that uses the distances through T, and a volume
  bound in the neighbourhood it gives.  The pair terms place all but two of the 24 centres
  within 2.1 and all but one within 2.15 (Proposition 2.43), but at this
  degree they cannot exclude one centre far out, and they say nothing about
  the directions.
