# Independent verification of the D4 package (v1.3.0), 24 September 2026

Scope: the manuscript `paper/D4.tex` ("The Sphere Packing Problem in
Dimension 4 and the Twenty-Four-Cell Conjecture") and everything in this
repository that supports it.  The repository was re-run from a fresh clone.
The Lean files were built with the pinned toolchain.  Every script was run in
a scratch copy, so that nothing it writes could overwrite shipped data.  Its
output was compared with the logs shipped in `multi_cap/runs/`,
`zonal/runs/` and `lean/runs/`.  The paper was compiled and checked against
the package.

## Verdict

**The computations hold up.**  Every machine check the paper relies on was
reproduced, and none failed:

- all seven Lean verifications, with the axiom reports the paper states;
- the exact and interval-arithmetic certificates for 23 contacts
  (Theorem 7.73) and for the single-deviation cap inequality (Section 13);
- the independent checks of the de Laat–Leijenhorst–de Muinck Keizer
  certificate for 24 contacts, including the zonal-matrix construction and
  the four constraint identities;
- the three independent 600-cell enumerations.

Details, and the steps still pending when this was written, are below.

**The proof does not, as written, establish the main theorems.**  The volume
results are proved for *contact configurations*: sets of directions pairwise
at least 60 degrees apart, at distance 2.  Theorems 1.5 and 1.6 are about the
Voronoi cell of an arbitrary packing.  The only bridge between the two is
Lemma 2.7 (shell localisation) and Lemma 2.9 (radial reduction), and both are
false as stated and used.  `proof_gaps/README.md` gives concrete, checked
counterexamples to the two lemmas, not to the theorems.

- **Radial reduction (Lemma 2.9).** Pulling every neighbour closer than
  2 sqrt 2 in to distance 2 does not in general produce a packing or a
  contact configuration.  One valid 49-ball packing gives a corner cell of
  volume 6.594 < 8.
- **Shell localisation (Lemma 2.7).** Neighbours at distance 2 sqrt 2 or
  more can cut the cell.  One cuts the deletion cell from 25/3 to 8.3028.

In every example the true cell still has volume at least 8.  A numerical
search found nothing below 8 in the single-deviation setting.  The theorems
may well be true.  The missing step, though, is where non-contact neighbours
enter the 24-cell conjecture, and it is a real gap, not a typo.  It should be
resolved before the paper claims the conjecture.

**The package had a handful of defects**, fixed in this commit (section
"Changes made").  The most important: `zonal/run_sos4.sh` omitted 5 of the 50
sum-of-squares blocks of the four-point constraint, and hard-coded the
authors' paths.  The shipped result is consistent, because those blocks came
from an unlogged earlier run, but a fresh run of the script as shipped could
not have reproduced it.

## Environment

| | |
| --- | --- |
| machine | Linux 6.18 x86_64, 4 cores, 15 GB, no swap |
| Lean | 4.34.0-rc2 (commit 6a10ac8c), the pinned toolchain, from the GitHub release |
| Python | 3.11.15; numpy 2.4.6, scipy 1.17.1, sympy 1.14.0, mpmath 1.3.0, python-flint 0.9.0, cvxpy 1.9.3, clarabel 0.11.1, matplotlib 3.11.2 (3.10.9 for the figure comparison), Pillow 12.3.0 |
| C | gcc 13.3.0, GMP 6.3.0 |
| TeX | TeX Live 2023 (pdfTeX 1.40.25), latexmk 4.83 |

The Lean release server and doi.org/Crossref are blocked by this machine's
network policy.  The toolchain came from the identical GitHub release asset.
The DOIs could not be re-resolved online (see "Bibliography").

## 1. Lean (all seven verifications)

| part | result | time | axioms |
| --- | --- | --- | --- |
| `D4Stress.lean` (17 theorems) | accepted | | none |
| `D4Meet.lean` (23 theorems) | accepted | | none |
| `D4Certificate.lean` (11 theorems) | accepted | | propext, Classical.choice, Quot.sound |
| `D4InnerProducts.lean` (7 theorems) | accepted | | propext, Classical.choice, Quot.sound |
| `D4RootLattices.lean` (7 theorems) | accepted | | propext + native_decide |
| (the five files, `run_all.sh`) | | 1 min 27 s | |
| `cell600/` (16 theorems) | `lake build` succeeds | 2 min 29 s | 13 kernel theorems: none (`cells_rootlike`: propext); 3 enumeration theorems: propext, Quot.sound + native_decide |
| `certificate/` `domain_ok` | `lake build` succeeds; `verifyDomain = true` | 1479 s for `D4CertMain` | propext, Classical.choice, Quot.sound + native_decide |

* The axiom report of `run_all.sh` is identical, line for line, to the
  shipped `lean/runs/run_all_2026-09-21.log`.
* `monomialCount` evaluates to 449, as the Python expansion says.
  `D4CertStat` prints `(0, 421881)`: status 0 after 421881 boxes, the count
  `lean/README.md` states.  The whole `certificate/` build took 44 min wall
  clock on a loaded machine.
* No `sorry`, `admit`, user `axiom`, `implemented_by`, `@[extern]` or `unsafe`
  anywhere in the Lean sources.
* The statements were read against the paper.  They say what the README says
  they say.  In particular, the LDL^T positivity test, the Sturm-sequence
  argument (with q(-1), q(1/2) != 0 proved) and the 600-cell search (the
  clique-cover pruning and the exactly-once enumeration) are sound.  So is the
  branch and bound of `D4CertDomain`: its domain discards, the omega terms
  dropped on boxes below 1/3 (omega >= 0), and the second-order Taylor form
  with inherited Hessian bounds.  What `native_decide` does not certify is
  stated in `lean/README.md`: the omega tables computed in Python, and the
  Taylor and monotonicity mathematics.
* Regenerated from their sources, `D4Certificate.lean`,
  `D4InnerProducts.lean`, `cell600/D4Cell600Enum.lean` and
  `certificate/D4CertData.lean` are **byte-identical** to the committed files.
  `llm24_p2.txt`, rewritten by the fresh run of `llm24_certificate_check.py`
  from the deposited data, is also byte-identical.  So the whole chain, from
  the 4TU archive through p_2 to the kernel-checked Lean file, reproduces
  exactly.

## 2. The certificate for 24 contacts (Theorem 7.25, Proposition 7.29)

* Archive: the two parts in `third_party/llm24-certificate/` join to
  `LasserreSphericalCodes.zip`, MD5 `02acd5270f7b3fa799abdeb5291706fd`
  (matches), 331 files under `proofs/4_24`.
* `llm24_certificate_check.py` (steps 1, 2, 4, 6, 7): 10 of 10 checks pass
  in 641 s.  The numbers the paper quotes all come out: 127 blocks (60 + 2 + 15
  + 50), total dimension 3726, largest 350, least Cholesky pivot
  1.38 x 10^-15 at 256 bits, objective exactly 24, p_2 of degree 16 vanishing
  exactly at -1, -1/2, 0, 1/2.  The output matches the shipped log except for
  its last two lines, which the current script words differently; the shipped
  log was made before that wording changed.
* `p2_zeroset_check.py` (Fractions, no FLINT): zero set {-1, -1/2, 0, 1/2},
  multiplicities 1, 2, 2, 1, q < 0 on [-1, 1/2].
* `D4InnerProducts.lean`: see section 1.
* `zonal/` (steps 3 and 5): ZONAL-PENDING

## 3. The certificate for 23 contacts (Theorem 7.73)

* `certificate_check.py`: PASS.  Exact LDL^T of the nine matrices, the bound
  0.0929 > 8 - A_* = 0.092855570294 in interval arithmetic, 449 monomials, and
  the branch and bound to level 44 with every box verified.  The output is
  identical to `multi_cap/runs/certificate_check_d8.log` apart from timings.
* `D4Certificate.lean` and `certificate/` (`domain_ok`): see section 1.  The
  Lean branch and bound is a second, independent implementation.

## 4. Other exact results

| script | claim | result |
| --- | --- | --- |
| `cap_certificate/cap_inequality_certificate.py` | Section 13, 21 checks | 21 of 21 PASS |
| `cell600_exact.py`, `cell600_enum.c`, Lean `cell600/` | Proposition 7.63: 0, 5, 115 | all three agree; Python and C logs identical to shipped |
| `octahedral48_exact.py` | 48 unit quaternions | identical to shipped log |
| `root_lattices_rank4.py` | census 1, 3, 23, 393; maxima 2, 6, 12, 24 | identical to shipped log; agrees with `D4RootLattices.lean` |
| `rigidity24.py` | constant stress on 96 pairs | identical to shipped log |
| `root_deletions_exact.py`, `rigidity23.py`, `rigidity_spectrum.py`, `three_point_reduction.py`, `covering_bound.py`, `extendability.py`, the `m24_*` scripts, `staged_*` | as in the README | all exit 0 with no FAIL; `three_point_reduction.py` identical to shipped log |

All 42 scripts the paper names (batch 1) ran to completion, exited 0 and
printed no FAIL, Traceback or MISMATCH.  Three of the `m24_*` scripts rewrite
floating-point SDP outputs in `continuation_out/` (`m24_primal_d6.npy`,
`m24_dir_d10/12/14.npy`).  The regenerated files differ from the committed
ones, as solver output from a different cvxpy/Clarabel version will.  The
paper treats them as exploration, and `m24_exact_reduction.py` passes on both
the committed and the regenerated primal measure.  See also finding F4.

## 5. The rest of the package (batch 2)

BATCH2-PENDING

## 6. Figures

All ten figure scripts in `paper/figures_new/` run and pass their label
collision tests.  With matplotlib 3.10.9, the version the shipped figures
were made with, every one of the ten PNGs is **pixel-identical** to the shipped
file.  The bytes differ only in the PNG encoding, which depends on the local
Pillow/zlib.  Under matplotlib 3.11.2 the renderings differ slightly.

## 7. The manuscript

* `latexmk -pdf D4.tex` on TeX Live 2023: 147 pages, no errors, no LaTeX
  warnings, no undefined or multiply defined references or citations, no
  overfull boxes, no `??` in the output.
* The text of the rebuilt PDF is identical to the shipped `paper/D4.pdf`
  except for the date (`\today`), so the shipped PDF is current.
* All 68 files the paper names exist in the package.
* All 189 references of the form "Theorem 7.73", "Section 7.5.3", "(7.33)"
  in the READMEs, Lean sources, scripts and logs were resolved against the
  compiled paper.  All point at an item of the right type, and a scan of the
  titles found them pointing at the right result.  There was one stale number
  (fixed): `D4Stress.lean` called the equilibrium relation (7.29); it is (7.33).
* Quoted numbers: NUMSCAN-PENDING

## 8. Bibliography

doi.org and Crossref are blocked here, so the 50 DOIs could not be
re-resolved online.  Checked from knowledge, every DOI of the classical
references carries the right venue, volume and pages: Viazovska, CKMRV,
Musin 2008 and 2018, Bachoc–Vallentin, Cohn–Elkies, Cohn–Kumar (both),
Cohn–Zhao, CJKT, Hales (1997, 2005, 2017), Lasserre, Parrilo, de Laat–
Vallentin, DGS, Schoenberg, Blichfeldt, Korkine–Zolotareff, Fejes Tóth,
Schütte–van der Waerden (both), Gorin–López, Voronoi, Curry–Schoenberg,
Farouki, Roth–Whiteley, Quickhull, NumPy, SciPy, SymPy and Arb.  The 2026
Clarabel entry could not be checked.

## Findings

F1 (proof, open): the radial reduction.  See `proof_gaps/README.md`.

F2 (proof, open): shell localisation.  See `proof_gaps/README.md`.

F3 (package, fixed): `zonal/run_sos4.sh` omitted blocks 0, 1, 3, 4 and 5 of
the 50 sum-of-squares blocks of the four-point constraint.  Its own comment
counts seven large blocks and gives four of them groups.  It also hard-coded
`/home/claude/aud` and `/tmp/claude-0/`.  The shipped log shows the running
total already holding 74606 terms before the first listed group, so those
blocks came from an earlier, unlogged run.  The script now takes the data
folder, `ps.txt` and a work directory as arguments, groups all 50 blocks, and
refuses to run if the groups do not cover 0..49 exactly once.  (A first
version of that guard used the variable name `GROUPS`, which bash reserves;
it is now `SOS_GROUPS`.)

F4 (paper, open): Remark 7.36 says a rational measure on 632 of 670 grid
triples, with value 5.9826875 and an exact LDL^T, *proves* that no degree-6
certificate exists.  Neither the measure nor any script that builds or checks
it is in the package.  The shipped `m24_primal_d6.npy` is a different,
floating-point measure: 522 triples, value 5.942410.  A fresh run of
`m24_primal_sdp.py` gives 2170 triples and 5.985098.  The Remark is not on
the path of the main theorems, but as it stands its "proof" cannot be
checked from the package.

F5 (documentation, fixed): `third_party/llm24-certificate/README.md` listed
three archive parts (`-00, -01, -02`); there are two.  It also gave the size
as "234 MB" beside the archive name.  The zip is 152065368 bytes.  The 234 is
the unpacked data as `du -h` reports it (243787660 bytes).  Both lines are now
exact.  The paper's "234 megabytes in 331 files under proofs/4_24" is that
unpacked figure, and is left as it is.

F6 (documentation, fixed): the README's directory entry and usage section for
`lean/` predated `D4RootLattices.lean` ("the first four", four `lean`
commands).  `lean/README.md` said "Thirty-seven theorems across the two
files"; `D4Stress.lean` and `D4Meet.lean` have 17 + 23 = 40.  The README said
`three_point_reduction.py` takes "about twenty minutes"; it takes seconds (its
Monte Carlo is seeded), as the authors' own `named_scripts_run.log` records.

F7 (metadata, fixed): `CITATION.cff` lacked P. Mandal's affiliation, which the
paper and `.zenodo.json` give.

F8 (documentation, noted): `RELEASE_NOTES_v1.3.0.md` says every figure script
"reproduces the shipped PNG byte for byte".  That holds pixel for pixel with
matplotlib 3.10.9, not byte for byte in general.  The same notes say D4.pdf
"is delivered separately (it is 21 MB)"; it is in `paper/` and is 18 MB.  The
release notes are a record of v1.3.0 and were not edited.

F9 (code, noted): `gcc -Wall -Wextra` on `zonal/psker.c` gives one
`-Wmaybe-uninitialized` warning in `canon()`.  It is a false positive: every
byte of `t` is written before the `memcpy`.

## Changes made

| file | change |
| --- | --- |
| `zonal/run_sos4.sh` | F3: portable arguments; all 50 blocks; coverage guard |
| `lean/D4Stress.lean` | comment: equation (7.29) -> (7.33); rechecked with Lean |
| `lean/README.md` | "Thirty-seven" -> "Forty" |
| `README.md` | `D4RootLattices.lean` in the directory entry and usage; runtime of `three_point_reduction.py` |
| `third_party/llm24-certificate/README.md` | part list and sizes (F5) |
| `CITATION.cff` | Mandal's affiliation (F7) |
| `independent_verification/` | this report, the proof-gap examples, and the logs of every run |

The manuscript `paper/D4.tex` was not changed.  F1, F2 and F4 are for the
authors to resolve.
