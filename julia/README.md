# D4Checks: the computer checks of the paper, in Julia

This project re-implements in Julia the computer checks on which proofs of the paper
depend (section "Computer verification" and appendix "Methods of computation"), with the
same arithmetic guarantees as the original programmes:

* **exact** steps use `Rational{BigInt}` (LDL^T factorisations, polynomial identities,
  Sturm sequences, values of bounds, count-vector maxima);
* **interval** steps use either Arb balls (`Arblib.jl`, 128 to 256 bits) or
  outward-rounded `Float64` intervals (`prevfloat`/`nextfloat` after every operation, the
  same operations in the same order as the original), and a sign is decided only when the
  enclosure excludes zero;
* floating point is used only where the original uses it and where it decides nothing
  (filters before an Arb evaluation, the search for thresholds that are then proved).

Nothing here is a new proof: each module is a second implementation of an existing
check, and each compares its box counts and values with the log of the original run
(`data/logs/`).  What is and is not covered is stated in the two tables below; the
second table matters as much as the first.

## Running

    julia --project=. -e 'using Pkg; Pkg.instantiate()'       # Arblib, NPZ, JSON, SHA, CodecZlib
    julia --project=. run_all.jl --list                       # the checks and their tiers
    julia --project=. run_all.jl --tier quick                 # about 2 minutes
    julia --project=. run_all.jl --tier standard --fetch-llm24 /some/dir     # about 42 minutes
    julia --project=. run_all.jl --tier full --llm24 /some/dir/LasserreSphericalCodes/proofs/4_24
    julia --project=. run_all.jl --tier full --only 'combo30$'               # one check

Julia 1.11 was used (1.10 should do).  One thread is enough; nothing is parallel.
`run_all.jl` writes one log per check to `logs/` and a table to `logs/SUMMARY.txt`, and
exits with status 0 only if no check failed.

The data set of de Laat, Leijenhorst and de Muinck Keizer (146 MB) is not in `data/`.
`--fetch-llm24 DIR` downloads `LasserreSphericalCodes.zip` from 4TU.ResearchData
(DOI 10.4121/74ce1c25-6fca-4680-8a36-e9c18e7e9594), checks

    SHA-256  2f806fac9f5f8644c53f9ad1590d920d346b61549d5d8860729ae6125a7568cc
    MD5      02acd5270f7b3fa799abdeb5291706fd   (the checksum published by the repository)
    size     152065368 bytes

and unpacks it (needs `unzip`).  Without `--llm24` or `--fetch-llm24` the two checks that
read it are reported as SKIPPED.

## What is ported

"Same counts" means that the Julia run reproduces the box and interval counts of the
original log exactly; where it does not, the difference and its cause are stated.

| Julia (`src/`) | original programme, log | supports | arithmetic | agreement with the original log |
|---|---|---|---|---|
| `anc_check_exact.jl` | `anc/check_exact.py` | the explicit constants of the annular and exterior arguments | exact rationals; rational-endpoint intervals for sqrt, atan, pi | all values |
| `anc_annular.jl` | `anc/annular_constants.py`, `cap_comparison.py`, `cap_profiles_26.py` (`anc/*.log`) | the annular constants; the cap comparison; the 17 845 count profiles at 26 centres (12 378 excluded) | exact rationals and rational-endpoint intervals | all values and counts |
| `anc_exterior.jl` | `anc/exterior_test.py` (`anc/exterior-computation.log`) | the exterior computation | exact rationals and rational-endpoint intervals | all values |
| `cap_certificate.jl` | `cap_certificate/cap_inequality_certificate.py` | thm:cap-inequality, rem:cap-boxes: the 25 vertices of Q, vol Q = 25/3, the three ranges of the case \|S\| = 1, the 303 boxes | exact rationals, Sturm sequences | 303 boxes, largest corner value 0.99755050 (steps A, B, D, E, F, F'; the two floating-point cross-checks C and G are not ported) |
| `constants.jl` | `multi_cap/twentyfour_close_check.py` | lem:twentyfour-close | Arb, 200 bits | all inequalities; room 2.6492e-5 |
| `constants.jl` | `multi_cap/pushout_check.py` (`runs/pushout_check.log`) | thm:G-pushout | Arb, 128 bits; exact rationals | T = 8.000015233 and 7.906940156 (the qhull comparison, floating point and not part of the proof, is not ported) |
| `covering.jl` | `multi_cap/covering_bound.py` | thm:covering-bound, cor:m22, tab:covering-bound | certified bisection and Arb (the original is floating point) | 8.046376, 7.916728, 7.798989, density 0.632749 (the Monte-Carlo and quadrature cross-checks are not ported) |
| `threepoint.jl`, `run_certificate` | `multi_cap/certificate_check.py` (`runs/certificate_check_d8.log`) | thm:certificate (at most 23 centres) | exact LDL^T and exact bound; omega tables in Arb; branch and bound in outward-rounded Float64 | bound identical; 313 822 boxes against 313 780 (see "Differences") |
| `threepoint.jl`, `run_labelled` | `multi_cap/labelled_certificate_check.py` (`runs/labelled_certificate_check.log`) | thm:m23-labelled, regions I, II_s, II_f | the same, with the Gamma envelopes in Arb | II_f identical (6 482 / 1 242); II_s 30 050 / 3 410 against 30 051 / 3 408 |
| `cardinality.jl` | `multi_cap/certify_cardinality.py` (`runs/certify_cardinality_*.log`, `gap_closure/density/card26_d10_certify.log`) | thm:kissing-stable (N = 25), thm:twenty-six (N = 26) | exact LDL^T, exact B, outward-rounded Float64 branch and bounds, exact final inequality | same counts on the degree-8 certificates (197 intervals, 115 532 boxes; 221 and 150 225) |
| `radial.jl`, `run_count` | `multi_cap/radial_count_check.py` (`runs/radial_count_check_31.log`) | thm:count31 | exact LDL^T; Bernstein subdivision with an a-priori rounding allowance; Pi at box corners in Arb | same counts: 191 740 closed, 180 310 set aside, 54 248 Arb values, 356 intervals |
| `radial.jl`, `run_case` | `multi_cap/radial_case_check.py` (`runs/radial_case_check_*.log`, `gap_closure/density/check_M.log`) | thm:count30 (the tree, `case_30.json`), prop:C-radial and prop:G-radial (the single splits `case_M_rR.json`), prop:levels (`density_M.json`) | as above, plus the tree check and the dynamic programme over count vectors in exact rationals | same counts on every file run (see `logs/SUMMARY.txt`) |
| `combo.jl` | `multi_cap/combo30_check.py`, `combo_case_check.py`, `combo_mixtype_check.py` (`runs/combo*_check.log`) | thm:count30 fifth case, thm:count29, prop:count28-few, prop:count28-mixed, prop:count27-few, prop:levels at M = 24 and 26 | exact LDL^T; pairs by tensor Bernstein bounds with Pi in Arb; bins as in `radial.jl`; triples by Taylor forms in outward-rounded Float64; the bound over all count vectors in exact rationals | `combo24_level`: every count identical (3 361 pair boxes, 285 intervals, 104 028 triple boxes, bound 3.253872); `combo30`: thresholds and bound 3.087203 identical, see "What was run" |
| `llm24.jl`, `run_check` | `multi_cap/llm24_certificate_check.py` (`runs/llm24_certificate_check.log`) | prop:verified, thm:twentyfour-points: steps 1, 2, 4, 6, 7 of the verification of the deposited certificate | exact rationals; Cholesky in Arb at 256 bits; Sturm sequence over Q | 127 blocks, total dimension 3726, least pivot 1.38e-15, 81 blocks also by exact LDL^T, objective 24, p_2 equal coefficient by coefficient to the Python output, multiplicities 1, 2, 2, 1 |
| `llm24.jl`, `run_eps0` | `multi_cap/explicit_eps0.py` (`runs/explicit_eps0.log`) | thm:near-contact, thm:m24-robust: the explicit eps_0 | Arb at 256 bits; exact rationals on the grid of step 1e-5 | all values (q_min 0.000262008, kappa* = 2e-26, the brackets) |
| `gnbhd.jl`, `run` | `anc/gnbhd.py`, `anc/gnbhd_bins.py` of the TeX source (`data/logs/gnbhd/log_thmBC.txt`, `log_thmC_48.txt`) | thm:G-nbhd, app:G-nbhd, tab:G-nbhd-bins: lem:T-near-roots and thm:one-large-push at e_0 = 1/48, the first case and the 29 bins | Arb, 200 bits; grid values exact rationals | S(2), Pi_0, T_root, L, the two bounds of lem:T-near-roots, the 16 quoted constants, m = 0.029208 in the first case, and every row (interval, y, m) of the table, which prints m rounded down; 62 checks (see "Differences", 5) |
| `gnbhd.jl`, `run_hexagon` | the statements after thm:hexagon-rigidity (no original programme) | thm:hexagon-rigidity, cor:hexagon-close: F(1/48, 17/20000) < 1/48, F(7/100, 9/5000) < 7/100, 2/sqrt(325) = 0.110940..., 2/sqrt(1-2s) = 2.0017021... | exact rationals (the two square roots removed by squaring twice); Arb as a cross-check | all statements; 15 checks |
| `profile_levels.jl` | `anc/profile_levels.py` of the TeX source (`data/logs/density/density_levels.out`, `density_levels_tight.out`, `combo26_tight_check.log`) | the profile-level bound on U(Y) for 24 to 30 centres and the density bound it gives | exact rationals rounded up to multiples of 1e-13 (integer numerators, overflow-checked); t = z^T A_0^{-1} z by exact elimination; S(e_b) and the density in Arb at 200 bits | all values: the seven S(e_b), the five (t/2, m_0, m_1), the profile counts, the seven maxima with and without n_0 <= 24, max G = 3.2888647945 at (22,0,3,0,0,0,0), density < 0.631495029884106 < 0.631496; 21 checks. With `--tight-log` (conditional, see below): 3.2826739, density < 0.630995132110131 < 0.630996; 20 checks |

## What is not ported

| original | supports | why not, and what covers it |
|---|---|---|
| `zonal/*.py`, `zonal/psker.c` (steps 3 and 5 of the verification of the cited certificate: the zonal matrices and the four polynomial identities; about five hours) | prop:verified | Not ported.  **The Julia code therefore does not re-verify the cited certificate completely**: it repeats the five steps of `llm24_certificate_check.py` and not the polynomial identities. |
| the fast branch and bound `multi_cap/pair_bb_fast.py` | thm:count30 (`case_30.json`, 185 million boxes closed) | `radial.jl` is the reference algorithm (one box at a time, exact corners); at about 5 000 boxes per second the tree would take about 20 hours.  It was not run. |
| the Bernstein sign certificates of the arcs (`arc1_v1w1/`, `arc2_w1v2/`, `core/`, `hessian_multidir/`, `swap_configs/`) | sec:arc1-cert, sec:arc2-cert (exact, in Z[sqrt 2] and Z[sqrt 5]) | Not ported (not in the brief of this port; a large body of SymPy code). |
| `multi_cap/closure_lemmas.py`, `shell_reduction.py`, `root_meet.py`, `rigidity23.py`, `rigidity24.py`, `rigidity_spectrum.py`, `cell600_exact.py` with `cell600_enum.c`, `root_lattices_rank4.py`, `root_deletions_exact.py`, `octahedral48_exact.py`, `second_order/`, `m24_*.py`, `verification/refcell_verify.py` | sec:closure, sec:noncontact, prop:meet22, prop:meet21, prop:deletion-rigid, lem:rigidity-spectrum, prop:cell600, lem:root-lattice, sec:second-order | Not ported.  Most of these finite computations are also checked in the Lean development (tab:formal of the paper); `shell_reduction.py` is not. |
| the options `PRECHECK` (as `bb=false`), `ONLY` (as `only=`) are ported; `CHECKPOINT`, `REFINE_DIVERSE`, `REFINE_GRID` of `combo_case_check.py` are not | | Bookkeeping and floating-point search aids; they decide nothing. |
| all semidefinite programmes (`*_sdp.py`, `combo_gen*.py`, ...) and the floating-point surveys | the search for certificates | By design: a certificate is found in floating point and only its exact check enters a proof. |

### The density bound: `profile_levels` and `profile_levels_tight`

`profile_levels` uses the constants of the typed certificate at 26 centres that are built
into `anc/profile_levels.py` (and checks that the log `data/logs/density/recheck_combo26_full.log`
implies them); it gives max G = 3.2888648 and density < 0.631496.  `profile_levels_tight`,
the bound stated in the paper, replaces the pair and triple thresholds by those of the run
of `combo_case_check.py` with pair margin 2e-7, whose log is
`data/logs/density/combo26_tight_check.log` (one uninterrupted run of the unmodified
programme, PASS, 3 383 s on three cores); it gives max G = 3.2826739 and density < 0.630996.
The Julia code reads the thresholds from the log and checks that it is one full PASS
naming all fifteen kinds; the thresholds themselves are proved by that run, and again by
the full-tier entry `combo26_tight` (ported, not run here).

## What was run, and what was only ported

Every entry of `run_all.jl --list` in the tiers `quick` and `standard` was run in one
pass and passes: 54 entries, 557 checks, 2 488 s (about 41 minutes) on one core, of which
the cited certificate takes 8 minutes; the table is `logs/SUMMARY.txt`.

Of the tier `full`, four entries were run and pass: `labelled_certificate_full` (117 s),
`case_27_r2.1` (1 087 011 boxes, 272 s), `case_30_r2.05` (2 888 432 boxes, 964 s) and
`cardinality_26_d10` (1 287 272 boxes, 692 s); see `logs/FULL_TIER_RUNS.txt`.
`cardinality_d10_s008` was started and stopped by a time limit before it ended.  The others are ported and were run through reduced
smoke tests only:

* `combo.jl` was validated in full on `combo24_level` (every count identical to the
  Python log), on two of the nine triple inequalities of `combo30` (FFF 202 959 boxes,
  BBF 500 251 boxes: identical), and without branch and bound on `combo30` (all 15
  thresholds agree with the Python log to its seven printed digits; bound 3.087203) and
  on the mixed-typing certificate `combo27_k14f2g2amt_d10` (48 of 49 thresholds agree;
  see "Differences").  The options `refine_triples`, `refine_pairs`, `skip_covered`,
  `region`, `extra_cases` and `mono` are ported from the Python line by line and were
  **not** exercised on data.
* The long runs of the original campaign (about a hundred certificates for 27 and 28
  centres, one to six hours each on four cores) were not repeated.

Commands for the long runs, with what the original logs give (`data/logs/`):

    # thm:count30, fifth case: 47 646 pair boxes, 4 131 intervals, 8 275 056 triple boxes, bound 3.087203
    julia --project=. run_all.jl --tier full --only 'combo30$'            # about 50 minutes
    # prop:levels at M = 26: 9 996 pair boxes, 6 323 intervals, 7 512 629 triple boxes, bound 3.329050 < 3.3352
    julia --project=. run_all.jl --tier full --only 'combo26_level'
    # thm:count29: 21 895 pair boxes, 17 762 intervals, 16 301 757 triple boxes, bound 3.101462
    julia --project=. run_all.jl --tier full --only 'combo29'
    # thm:kissing-stable at degree 10: 1 479 intervals, 2 642 578 boxes;  thm:twenty-six: 790 intervals, 1 287 272 boxes
    julia --project=. run_all.jl --tier full --only 'cardinality_(d10_s008|26_d10)'
    # thm:m23-labelled with region I: 419 913 boxes in the Lean count of the paper
    julia --project=. run_all.jl --tier full --only 'labelled_certificate_full'
    # the tree of thm:count30: 4 leaves, 184 976 789 boxes closed in all (about 20 hours here)
    julia --project=. run_all.jl --tier full --only 'case_30$'

A case of the campaign is run from Julia directly, for instance prop:count28-few at most
13 centres within 2.0161 (`runs/combo28_lo13_check.log`: 3 279 356 pair boxes, 3 698
intervals, 15 805 222 triple boxes, bound 3.100609 at (13, 0, 11, 0, 0, 0, 4)):

    include("src/D4Checks.jl"); using .D4Checks
    D4Checks.Combo.run("data/combo/combo28_lo13_d8.npz", 8; case="data/combo/cases/case28_lo13.json",
                       mg2=5e-6, mg3=5e-7, mgm=1e-6, border=true, triple_extra=Dict("FFF" => 1e-5))

The keyword arguments are the options of `combo_case_check.py` (see the header of
`src/combo.jl`); the options of each original run are in `multi_cap/README.md` and
`gap_closure/README.md` of the original package.  A long run can be split, as with `ONLY`
in the Python: one run with `triples=false` (positivity, pairs, bins, bound) and runs with
`only=[...]` that together name every triple kind; the thresholds do not depend on the
split.  If a run raises a threshold after meeting a point above it (see below), the
bound must be taken with the raised threshold.

## Differences from the original runs

1. **Tables of omega** (`threepoint.jl`).  The Python builds its tables of omega, omega'
   and the bound on omega'' with mpmath intervals at 53 bits; the Julia uses Arb series.
   Both are rigorous enclosures, but the bound M2 on omega'' differs (0.291014 here,
   0.289712 there in one table), so a few boxes are split differently: 313 822 boxes
   against 313 780 for thm:certificate, 30 050 / 3 410 against 30 051 / 3 408 in region
   II_s.  The exact parts (positivity, the bound B) are identical.
2. **Rounding of the solver's matrices** (`combo.jl`).  The Python clips eigenvalues with
   numpy's `eigh` before rounding to dyadic rationals; the Julia uses `LinearAlgebra.eigen`.
   The rounded matrices can differ in their last bits.  Each is proved positive by exact
   LDL^T, which is all the argument needs, and on the certificates run the box counts
   came out identical.
3. **Thresholds** (`combo.jl`).  The Python sets each threshold c2, c3 to "largest value
   found + margin", the largest value coming from samples of numpy's generator refined
   by SLSQP.  The Julia samples the same grids plus points of its own generator and
   refines by a pattern search that follows the boundary of the Gram domain.  A
   threshold is valid once the branch and bound proves it, whatever its origin.  On
   `combo24_level` and `combo30` all thresholds agree with the Python log to the seven
   digits it prints.  On `combo27_k14f2g2amt_d10` one of 49 (triple A1A1A2) is lower by
   3.9e-8, at a margin of 2e-8: the Julia search missed a peak that the Python found.  In
   that situation the branch and bound meets a point above the threshold; the Julia
   then climbs from that point, raises the threshold and proves the inequality again
   (`retries`, tested by lowering a threshold on purpose), where the Python stops and is
   rerun with `TRIPLE_EXTRA`.
4. **Speed.**  The Julia code is single-threaded and unoptimised; the pair branch and
   bound of `radial.jl` runs at about 5 000 boxes per second, the triple branch and
   bound of `combo.jl` at about 3 500 boxes per second at 425 monomials.

5. **Rounding in tab:G-nbhd-bins.**  The table prints each certified margin m rounded
   down to four decimals (for instance 0.176582 is printed 0.1765).  `gnbhd.jl` checks
   that every printed value is below the certified lower bound by less than one unit of
   the last digit, and checks the two bounds the proof uses, m >= 0.147 on the bins and
   m >= 0.0292 in the first case.

## Layout

    Project.toml, Manifest.toml   the environment
    run_all.jl                    driver; tiers quick / standard / full
    src/                          one module per group of checks (D4Checks.jl includes them)
    data/certificates/            certificate_d8.npz (thm:certificate), llm24_p2.txt
    data/cardinality_certificates/, data/radial_certificates/, data/combo/   the certificates and case files
    data/logs/                    the logs of the original runs, read for the expected values
    logs/                         written by run_all.jl
