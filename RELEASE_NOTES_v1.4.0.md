# v1.4.0

Supplementary package for *The Sphere Packing Problem in Dimension 4 and the
Twenty-Four-Cell Conjecture*, after an independent re-verification of every
computation of the paper, with the corrections it led to in the paper and in
the package.

The re-verification, the two new programs of `multi_cap/` and the corrections
were developed with the assistance of Claude, an AI model made by Anthropic.

## What the paper now claims

- **Contact configurations: proved in full.**  For every contact
  configuration, meaning the directions to neighbours at distance exactly 2,
  vol >= 8, with equality only at the root system.  This is the contact count
  argument: at most 22 contacts (Theorem 1.3), exactly 23 (Theorem 7.40) and
  exactly 24 (Theorem 7.25), assembled in Corollary 7.31.  Every certificate
  it uses reproduces.
- **Voronoi cells of packings: proved under a stated hypothesis.**
  Theorem 1.5 ("General local bound") proves vol(V_c) >= 8, with equality
  only at the 24-cell, in two cases.  The first is when every centre within
  2 sqrt 2 of c touches it.  The second is when the distances d_i of those
  centres satisfy Phi(d_1, d_2, ...) > 8, where Phi is the covering criterion
  of Proposition 2.10.  Phi depends on the distances alone, not on the
  directions.  Theorem 1.7 gives the density bound pi^2/16 for periodic
  packings whose centres satisfy that hypothesis.
- **Open: crowded near-contacts (Conjecture 1.6).**  The case left open is a
  centre that has at least 23 neighbours within 2 sqrt 2, some of which do
  not touch it, and whose distances fail the criterion.  The density bound
  Delta_4 = pi^2/16 for all packings follows once this case is settled.  A
  numerical search of it finds nothing below 8; that is evidence, not proof.

v1.3.0 asserted the general bound without that hypothesis.  It did so through
Lemma 2.7 (shell localisation) and Lemma 2.9 (radial reduction), and both are
false as stated there:

- A centre at distance 2.9, beyond 2 sqrt 2, cuts the cell of a deletion of a
  root, from 25/3 to 1328453/160000 = 8.3028.
- Pulling non-contact neighbours in to distance 2 need not give a packing.  A
  packing of 49 balls has an all-contact corner of volume 6.594, while its
  true cell has volume 19.216.

In every example the true cell is above 8.  The details are in
`independent_verification/proof_gaps/`.

## The manuscript (paper/)

- **Section 2.4 (Lemma 2.7).**  Shell localisation is corrected.  A centre
  at distance at least 2R does not cut B(c, R), and a far centre can still
  cut the cell outside B(c, sqrt 2).  The far-centre example is given.
- **Section 2.6 (Lemma 2.9).**  Monotonicity in the contact radii is
  corrected.  The all-contact corner is a lower bound, and it is a contact
  configuration only when every active neighbour already touches the
  centre.  The 49-ball example is given.
- **New Section 2.7, "Neighbours that do not touch".**  It contains:
  - Proposition 2.10, the distance criterion;
  - Lemma 2.11, which uses the covering radius of 45 degrees of the root
    directions: 24 contacts leave no other centre within 2 sqrt 2, and the
    cell is the 24-cell;
  - the proof of Theorem 1.5;
  - Table 1 (thresholds);
  - Figure 3 (new, `fig_noncontact`);
  - the paragraph on the numerical search of the open case.
- **Theorem 1.1** (single deviation) is restated for the cell cut out by the
  active neighbours, with the hypothesis its proof uses: the deviating
  direction replaces a root that no active neighbour takes.
- **Other sections.**  The following are revised to the claims above: the
  title-page note, the abstract, the introduction, Theorems 1.5 and 1.7,
  new Conjecture 1.6, Section 6, the polar form, Corollary 7.31, the
  conclusion ("Three things remain open"), and the data and code
  availability statement.
- **Remark 7.36** now says that its rational measure is not in the
  supplement, and records the floating-point measures that are.  The
  deposited one has value 5.9424 on 2170 triples; a rerun gives 5.9851 on
  522 triples.
- **Bibliography**: the Zenodo entry cites v1.4.0 under the concept DOI.
- **Appendix "Computational methodology"**: new subsections on the two new
  scripts and on the re-verification.

## The package

- `independent_verification/`: the report (`REPORT.md`), the logs of every
  run, the proof-gap examples, and `gegenbauer_check.py`,
  `compare_figures.py` and `numscan.py`.
- `multi_cap/shell_reduction.py`: the distance criterion in ball arithmetic
  (python-flint).  Log: `multi_cap/runs/shell_reduction.log`.
- `multi_cap/shell_neighbour_search.py`: the numerical search of the open
  case.  Logs: `multi_cap/runs/shell_neighbour_search_*.log`.
- `paper/figures_new/fig_noncontact.py` and `paper/figures/fig_noncontact.png`.
- `zonal/run_sos4.sh` runs all fifty sum-of-squares blocks of constraint 4
  (v1.3.0 omitted five).  It checks that its groups cover every block,
  takes its paths as arguments, and exits 1 on a failure.
- Documentation fixes: `lean/README.md`, `README.md`,
  `third_party/llm24-certificate/README.md`, and a comment in
  `lean/D4Stress.lean`.
- `CITATION.cff` (Mandal's affiliation) and `.zenodo.json` carry version
  v1.4.0.

## Release

Published on 24 September 2026 as the tag `v1.4.0` on `main`
(https://github.com/creelie/d4-voronoi-cells/releases/tag/v1.4.0) and
archived by Zenodo under the concept DOI 10.5281/zenodo.22766562, which
always resolves to the newest version and is the identifier the paper cites.
The version DOI of v1.4.0 is recorded in `README.md` once Zenodo shows it.
