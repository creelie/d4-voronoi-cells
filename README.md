# The Sphere Packing Problem in Dimension 4 and the Twenty-Four-Cell Conjecture: supplementary code

Certificates, the programmes that check them, run logs and Lean 4 checks for

> D. Bhattacharjee, U. Bhattacharya, P. Mandal and S. Bhattacharya,
> *The Sphere Packing Problem in Dimension 4 and the Twenty-Four-Cell
> Conjecture*, arXiv:2609.25120.

Every exact value that a proof in the paper uses is stated in the paper. The
objects it cannot print, the matrices of the semidefinite certificates and the
box lists of the branch and bounds, are here, with the programmes that check
them. The paper's section "Positive kernels and certificates" proves the
positivity theorems that turn each of these certificates into a proof, states
what each check must establish, and lists the certificates with the size of
each check. The region-by-region computations for one deviating contact and
the computations for several deviating contacts are in two appendices of the
paper; no proof uses them.

## What the paper proves

Let V_c be the Voronoi cell of a centre c of a packing of unit balls in
R^4. The paper proves vol(V_c) >= 8, with equality only for the regular
24-cell of the D_4 root lattice, in each of the following cases:

- every centre within 2 sqrt 2 of c touches c (contact configurations), with
  the classification of twenty-four contacts proved from the certificate of
  de Laat, Leijenhorst and de Muinck Keizer;
- at most twenty-three other centres lie within sqrt 6 of c;
- at least twenty-nine lie within sqrt 6;
- exactly twenty-eight lie within sqrt 6, and at most fourteen of them
  within 2.0161, or at most fifteen with at most two beyond 2.35;
- exactly twenty-four lie within sqrt 6, on the rays of a root system, or
  the contacts contain enough of a root system;
- no centre lies strictly between 2 + 4e-26 and sqrt 6.

With one deviating contact direction the bound has an elementary proof, by
the extremal cap of a cross-polytope. With no assumption on the centres, the
same lower bounds for the cells give density at most 0.63477 for every
packing of unit balls in R^4; this is below the three-point bound 0.63611 of
Cohn, de Laat and Salmon, and above pi^2/16 = 0.61685.

For the remaining centres the paper reduces the bound to two explicit
statements: (G), a volume bound when exactly twenty-four centres lie within
sqrt 6, and (C), a bound on the caps when twenty-five to twenty-eight do. Together
they imply that D_4 gives the densest packing of R^4. **Neither is proved**;
(C) is proved at twenty-eight only when at most fourteen centres lie within
2.0161, or fifteen with at most two beyond 2.35.
The paper records what is proved about them and where each method stops;
`gap_closure/` holds the floating-point explorations against them, none of
which is a certificate, and the logs of the exact checks behind the density
bound of the paper, whose certificates are in multi_cap/radial_certificates/.

## Citation

The release v2.0.0 accompanies version 3 of arXiv:2609.25120. The concept
DOI 10.5281/zenodo.22766562, which the paper cites, covers every version
and resolves to the latest release; the earlier package v1.12.0 is
10.5281/zenodo.23076993. `CITATION.cff` and `.zenodo.json` carry the
metadata. The package is maintained at
https://github.com/creelie/d4-voronoi-cells and can also be obtained from
Deep Bhattacharjee, itsdeep@live.com.

## Layout

    paper/              the manuscript (D4.tex, D4.pdf), its figures and the
                          scripts that draw them (figures_new/)
    cap_certificate/    the cap inequality of the section "Positivity for
                          every deviation direction", end to end
    multi_cap/          the contact cases and the non-contact cases: covering
                          bound, twenty-three contacts, the certificate of de
                          Laat, Leijenhorst and de Muinck Keizer, stability of
                          the kissing number, the distance criterion, the
                          labelled certificate, the counts and radial bounds
                          for (C) and (G), and the push-outs; second_order/
                          and stability/ have their own README files
    zonal/              the zonal matrices and the four polynomial identities
                          of the certificate of the twenty-four-point case, in
                          exact arithmetic (see zonal/README.md)
    level2/             the second-level programme of de Laat, Leijenhorst and
                          de Muinck Keizer, set up from zonal/ and solved on
                          the enlarged domain [-1, 1/2 + s] (see level2/README.md)
    gap_closure/        floating-point explorations of (G) and (C), and the
                          logs of the density bound (see gap_closure/README.md)
    lean/               Lean 4 checks, no Mathlib (see lean/README.md and the
                          table "The formally verified parts of the paper")
    third_party/        the certificate data set of de Laat, Leijenhorst and
                          de Muinck Keizer, with its licence and checksum
    independent_verification/
                        an independent re-verification of the package, with
                          its report and logs
    verification/       cross-checks of the reference cell
    core/               arithmetic routines shared by the scripts
    data/               cached intermediate results (.pkl)
    misc/               utility and diagnostic scripts
    release_notes/      the notes of each release

Four directories hold the elementary computations on deviations from the
root system: the region certificates of the two arcs of the fundamental
triangle (arc1_v1w1/, arc2_w1v2/), the joint Hessian at several deviating
directions and the dense configurations A_18 and A_20 (hessian_multidir/),
and the swap paths (swap_configs/). The single-deviation theorem itself is
proved by the cap inequality (cap_certificate/), and the multi-direction case
through the classification of twenty-four contacts; these four directories
give the independent checks and the numerical evidence that the paper reports
beside those proofs.

## Exact and numerical

Every step the paper uses as a proof is computed in exact rational
arithmetic (fractions.Fraction, sympy.Rational, Lean's Rat), in exact
dyadic arithmetic, or in outward-rounded interval and ball arithmetic
(mpmath intervals, python-flint's arb). High-precision floating point is
used only for numerical estimates, which the paper labels as such and never
uses as proofs. The floating-point programmes of gap_closure/ and the sampled programmes
of multi_cap/ only find candidates; a candidate enters a proof only after it
is rounded to rationals and checked exactly or in interval arithmetic.

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
    (cd lean/count && lake build)

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
degree-10 certificate of the theorem "The kissing number is stable" in
lean/cardinality/, several hours
on one core. The cached files in data/ spare the consuming scripts
the wait.

A few symbolic steps do not terminate, and the paper reports that they do
not. Those scripts stop at an explicit wall-clock bound, which an
environment variable overrides, and report what they established:

    arc2_w1v2/bandBcurve_certificate.py   BANDB_BUDGET_SECONDS (240) and
                                          BANDB_PER_SIMPLEX_SECONDS (90);
                                          the moving-volume step does not
                                          complete, as the section on the
                                          second arc reports
    misc/sumtest.py                       SUMTEST_STAGE_SECONDS (120)
    swap_configs/swap_exact_volume.py     SWAPVOL_STAGE_SECONDS (300)
    swap_configs/swap_exact_volume3.py    SWAPVOL_STAGE_SECONDS (300)
    arc1_v1w1/bisect_crossover.py         BISECT_ITERATIONS (40), about
                                          three minutes

## Licence

See LICENSE. The data set in third_party/ keeps its own licence (MIT),
reproduced beside it.
