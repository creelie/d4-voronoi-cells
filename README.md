# Supplementary code for Voronoi Cells of Four-Dimensional Unit-Ball Packings

Programs, certificates, run logs and Lean 4 checks for

> D. Bhattacharjee, U. Bhattacharya, P. Mandal and S. Bhattacharya,
> *Voronoi Cells of Four-Dimensional Unit-Ball Packings*, arXiv:2609.25120.

Every exact value the paper asserts is printed in the paper, so its
derivations can be followed without running anything. The one computation
the paper cannot print is the verification of the certificate of
Section 7.5.3, whose inputs are the deposited data in `third_party/` and
whose run logs are in `multi_cap/runs` and `zonal/runs`. The package lets every derivation, that one included, be
repeated independently.

## What the paper proves

Let V_c be the Voronoi cell of a centre c of a packing of unit balls in
R^4. The paper proves vol(V_c) >= 8, with equality only for the regular
24-cell of the D_4 root lattice, in each of the following cases:

- every centre that can cut the cell touches c (contact configurations,
  Sections 7 to 20, with the classification of twenty-four contacts proved
  in Section 7.5);
- at most twenty-three other centres lie within sqrt 6 of c (Theorem 21.8);
- at least thirty-one lie within sqrt 6 (Theorem 21.25);
- no centre lies strictly between 2 + 4e-26 and sqrt 6 (Theorem 21.9 and
  Section 21.1).

For the remaining centres, Proposition 21.23 reduces the bound to two
explicit statements: (G), a volume bound when exactly twenty-four centres
lie within sqrt 6, and (C), a bound on the caps when twenty-five to thirty
do. Together they imply that D_4 gives the densest packing of R^4. Neither
is proved; Section 21.3 records what is proved about them and where each
method stops.

## Citation

The package is archived on Zenodo under the concept DOI
10.5281/zenodo.22766562, which resolves to the latest release;
`CITATION.cff` and `.zenodo.json` carry the metadata. It is maintained at
https://github.com/creelie/d4-voronoi-cells and can also be obtained from
Deep Bhattacharjee, itsdeep@live.com.

## Layout

    paper/              the manuscript (D4.tex, D4.pdf), its figures and
                          the scripts that draw them (figures_new/)
    cap_certificate/    the cap inequality of Section 13, end to end
    arc1_v1w1/          first-arc certificates (Sections 9 and 10)
    arc2_w1v2/          second-arc certificates and the symmetry of the
                          fundamental triangle (Sections 8, 11 and 12)
    hessian_multidir/   joint-Hessian computations (Sections 15 to 18)
    swap_configs/       swap configurations at finite angle (Section 20)
    multi_cap/          the contact cases of Section 7 and the non-contact
                          cases of Section 21: covering bound, twenty-three
                          contacts, the certificate of de Laat, Leijenhorst
                          and de Muinck Keizer, stability of the kissing
                          number, the distance criterion, the labelled
                          certificate, the counts and radial bounds for (C)
                          and (G); second_order/ (Section 21.2) and
                          stability/ (the corollary "No room beside a near
                          root system") have their own README files
    zonal/              the zonal matrices and the four polynomial
                          identities of the certificate of Section 7.5.3,
                          in exact arithmetic (see zonal/README.md)
    level2/             the second-level programme of de Laat, Leijenhorst
                          and de Muinck Keizer, set up from zonal/ and solved
                          on the enlarged domain [-1, 1/2 + s] (see
                          level2/README.md)
    lean/               Lean 4 checks, no Mathlib: thirteen files checked
                          by run_all.sh and three Lake projects (see
                          lean/README.md and Table 6 of the paper)
    third_party/        the certificate data set of de Laat, Leijenhorst
                          and de Muinck Keizer, with its licence and checksum
    independent_verification/
                        an independent re-verification of the package, with
                          its report and logs
    verification/       cross-checks of the reference cell
    core/               arithmetic and Hessian routines shared by the scripts
    data/               cached intermediate results (.pkl) of the longer
                          symbolic derivations
    misc/               utility and diagnostic scripts
    release_notes/      the notes of each release

The README files of cap_certificate/, arc2_w1v2/, hessian_multidir/ and
multi_cap/ describe their scripts one by one, with the results of the paper
each one supports.

## Exact and numerical

Every step the paper uses as a proof is computed in exact rational
arithmetic (fractions.Fraction, sympy.Rational, Lean's Rat), in exact
dyadic arithmetic, or in outward-rounded interval and ball arithmetic
(mpmath intervals, python-flint's arb). High-precision floating point is
used only for numerical estimates, which the paper labels as such and never
uses as proofs. In particular no script in hessian_multidir/ proves the
multi-direction positivity conjecture: each result there is either an exact
statement about a specific finite configuration or numerical evidence
reported as such.

## Requirements

    Python >= 3.9
    sympy >= 1.10
    mpmath >= 1.2
    numpy and scipy (the double-precision scans, cap_inequality_certificate.py,
      vertex_degeneracy_check.py and the scripts of multi_cap/stability/)
    cvxpy >= 1.9 with an SDP solver such as Clarabel (gram_sos_lib.py,
      quartic_fit_and_check.py and the sampled programmes of multi_cap/)
    python-flint >= 0.9 (llm24_certificate_check.py, shell_reduction.py,
      closure_lemmas.py, labelled_certificate_check.py, explicit_eps0.py and
      the radial checks of multi_cap/)
    Julia 1.10 with the package of de Laat, Leijenhorst and de Muinck Keizer
      (level2/ only)
    Lean 4, the version pinned in lean/lean-toolchain (lean/ only)

## Usage

Each script runs from any working directory, for example

    python cap_certificate/cap_inequality_certificate.py
    python multi_cap/covering_bound.py
    python multi_cap/certificate_check.py
    python multi_cap/certify_cardinality.py multi_cap/cardinality_certificates/cert_d8_t0.50650.npz 1e-6 1e-5 1e-5
    python multi_cap/llm24_certificate_check.py /path/to/LasserreSphericalCodes/proofs/4_24
    python multi_cap/shell_reduction.py
    python multi_cap/closure_lemmas.py
    python multi_cap/labelled_certificate_check.py
    python multi_cap/radial_count_check.py multi_cap/radial_certificates/radial_31.json
    python multi_cap/pushout_check.py
    python multi_cap/twentyfour_close_check.py
    python verification/refcell_verify.py

with the arguments each README lists. The Lean files are checked with

    cd lean && sh run_all.sh
    (cd lean/cell600 && lake build)
    (cd lean/certificate && lake build)
    (cd lean/cardinality && lake build)

and the generated Lean files are regenerated with the gen_*.py scripts
beside them (lean/README.md lists them). The C enumeration of the 600-cell
runs as

    (cd multi_cap && cc -O2 -o cell600_enum cell600_enum.c && ./cell600_enum)

after cell600_exact.py has written cell600_graph.txt. Scripts that read a
cached intermediate result look for it in data/ relative to the package
root.

## Running times

Most scripts finish in seconds. The slowest are the exact region
derivations in arc1_v1w1/ (regionA_derive.py, region6_derive.py,
regionC_derive.py, regionD_derive.py and regionE_derive.py), whose symbolic
cancellation can take an hour or more, the directional sweeps in
hessian_multidir/, about eight minutes each, and the formal check of the
degree-10 certificate of Theorem 7.79 in lean/cardinality/, about four
hours on one core. The cached files in data/ spare the consuming scripts
the wait.

A few symbolic steps do not terminate, and the paper reports that they do
not. Those scripts stop at an explicit wall-clock bound, which an
environment variable overrides, and report what they established:

    arc2_w1v2/bandBcurve_certificate.py   BANDB_BUDGET_SECONDS (240) and
                                          BANDB_PER_SIMPLEX_SECONDS (90);
                                          the moving-volume step does not
                                          complete (Section 12.2), and
                                          Proposition 8.5 makes the direct
                                          certificate unnecessary
    misc/sumtest.py                       SUMTEST_STAGE_SECONDS (120)
    swap_configs/swap_exact_volume.py     SWAPVOL_STAGE_SECONDS (300)
    swap_configs/swap_exact_volume3.py    SWAPVOL_STAGE_SECONDS (300)
    arc1_v1w1/bisect_crossover.py         BISECT_ITERATIONS (40), about
                                          three minutes

## Licence

See LICENSE. The data set in third_party/ keeps its own licence (MIT),
reproduced beside it.
