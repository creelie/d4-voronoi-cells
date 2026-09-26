# Supplementary code package (Code_4)

Python scripts supporting the computations in:

  "The Sphere Packing Problem in Dimension 4 and the
   Twenty-Four-Cell Conjecture"
  Deep Bhattacharjee, Ushashi Bhattacharya, Priyabrata Mandal,
  Shounak Bhattacharya

Every exact value the paper asserts is printed in the paper itself, so
the derivations can be followed without running anything; the one
computation that the paper cannot print is the verification of the
certificate of Section 7.5.3, whose inputs are the deposited data in
third_party/llm24-certificate and whose run logs are in multi_cap/runs and
zonal/runs. The package is here so that every derivation, that one
included, can be repeated independently.

## Availability

The package is maintained at https://github.com/creelie/d4-voronoi-cells
and archived on Zenodo through that repository (CITATION.cff and
.zenodo.json carry the metadata of the archive). It is also obtainable
from Deep Bhattacharjee <itsdeep@live.com>.

Archive: DOI 10.5281/zenodo.22766562; cite it.  The release that
corresponds to the paper is v1.7.0.

## What is new in v1.7.0

The constant epsilon_0 of Theorem 2.18 is explicit: epsilon_0 = 4e-26.  If
no centre lies at a distance from c strictly between 2 + epsilon_0 and
sqrt 6, then vol(V_c) >= 8, with equality only at the 24-cell.  Before,
epsilon_0 came from a compactness argument.  Three new pieces make it a
number:

- Theorem 7.76 (twenty-four points, approximately).  A set of directions
  with pairwise inner products at most 1/2 + 2e-26 has at most 24 points,
  and 24 of them lie within 0.0113 (root sum of squares) of a copy of the
  normalised roots.  The proof reads the certificate of de Laat,
  Leijenhorst and de Muinck Keizer robustly: in its sum-of-squares
  identities only the terms whose weight contains a factor
  (u + 1)(1/2 - u) with u > 1/2 can turn negative, and their size is
  bounded from the deposited data.
- Lemma 7.77.  If a Gram matrix of 24 unit vectors is within 3e-4 of a
  matrix g with entries in {-1, -1/2, 0, 1/2}, then 2g is integral and
  every 5 x 5 minor of it is an integer of absolute value below 1, so
  zero; hence g is the Gram matrix of the root system, and a Davis-Kahan
  estimate bounds the distance to it.
- Theorem 2.18, second statement.  The root system is a strict local
  minimum in an explicit neighbourhood that uses no certificate: exactly
  24 centres within sqrt 6, directions within 1/48 of the normalised roots
  (root sum of squares), and sum of delta_i <= 4e-5.  The constants of
  the facet and hull estimates are explicit.

What is left of Conjecture 1.6 is unchanged in shape: a centre with at
least 24 other centres within sqrt 6, one of them at a distance between
2 + epsilon_0 and sqrt 6.  The new code was developed, and its
computations run, with the assistance of Claude, an AI model made by
Anthropic.

  paper/              Theorem 2.18 restated with epsilon_0 = 4e-26 and
                        the explicit neighbourhood, its proof rewritten
                        with explicit constants; Section 7.13
                        (Theorem 7.76, Lemma 7.77) is new; the abstract,
                        Theorem 1.7, the opening of Section 2.8, "What is
                        left", the conclusion, the notation, the code index
                        and the data availability statement are revised.
                        One new figure: Figure 25
                        (figures_new/tikz_eps0.tex, TikZ), the chain of
                        explicit bounds and, to scale, sigma_2 against the
                        error level; later figures are renumbered by one.
  multi_cap/truncated_search.py rays
                      Where the pair terms reach 8 near the root system:
                        sum of delta_i = 0.155 pushed out evenly, and
                        delta = 0.197 for one centre, against the explicit
                        neighbourhood 4e-5 of Theorem 2.18.  Log:
                        runs/truncated_search_rays.log.
  multi_cap/explicit_eps0.py
                      The constants of Theorems 7.76 and 2.18: the sizes
                        of the sign-changing sum-of-squares terms of the
                        certificate in ball arithmetic, the exact division
                        of sigma_2 by its zeros, the minor bounds of Lemma
                        7.77 and the facet estimates in exact arithmetic.
                        Sixteen seconds.  Log: runs/explicit_eps0.log.

## What is new in v1.6.0

Theorem 2.17: a centre with at most 23 other centres within sqrt 6 has a
Voronoi cell of volume greater than 8, whatever the distances of those
centres.  The three-point certificate of Theorem 7.73 is carried from
contacts to centres at mixed distances: a centre that moves out from
distance 2 hands back part of its cap, and one eleventh of that, spent in
each triple of centres, pays for every pair of directions closer than
60 degrees that the packing then allows.  Conjecture 1.6 is narrowed to a
centre with at least 24 other centres within sqrt 6, one of them at an
intermediate distance; at 24 the pair terms cannot decide, since at the
root system they give 7.906940 < 8.  The new code was developed, and its
computations run, with the assistance of Claude, an AI model made by
Anthropic.

  paper/              Theorem 2.17 with its full proof replaces Corollary
                        2.17 (twenty-three contacts), which it contains;
                        Theorem 1.7 of the introduction, the abstract, the
                        opening of Section 2.8, "What is left" and the
                        conclusion are revised to it.  Three new figures:
                        Figure 6 (figures_new/tikz_labelled.tex, TikZ),
                        Figure 7 (figures_new/fig_labelled.py) and
                        Figure 8 (figures_new/tikz_stop.tex, TikZ); later
                        figures are renumbered by three.  New footnotes in
                        the proof, and new rows in the index of notation.
  multi_cap/labelled_certificate_check.py
                      The verification of Theorem 2.17: the constants in
                        ball arithmetic, a rerun of the branch and bound of
                        (C), and two new branch and bounds (30 051 and
                        6 482 boxes).  Log: runs/labelled_certificate_check.log
                        (about six minutes).
  multi_cap/count_bound_sqrt6.py
                      At most 49 other centres lie within sqrt 6 of a
                        centre: their directions have inner products
                        below 2/3, and an exact Delsarte polynomial of
                        degree 13 (Sturm's theorem, sympy) gives
                        N <= 49.577.  Log: runs/count_bound_sqrt6.log.
  multi_cap/truncated_search.py
                      Floating-point minimisation of the truncated volume
                        with exactly M centres within sqrt 6 (evidence, not
                        proof).  Logs: runs/truncated_search_M*.log.
  lean/D4Closure.lean Two more theorems: the identity behind the packing
                        bound a(d1, d2) <= 1/2 + (d1 + d2 - 4)/4, and the
                        counts 21 and 231 = 21 * 11 of the proof.
  independent_verification/logs/doi_audit_2026-09-25.md
                      Every DOI of the bibliography checked against the
                        publisher's or indexer's page.

## What is new in v1.5.0

The case left open by v1.4.0 (Conjecture 1.6, a centre crowded by
near-contacts) is narrowed in a new Section 2.8 of the paper, with a second
machine-verification layer in Lean for its exact content.  The new code was
developed, like that of v1.4.0, with the assistance of Claude, an AI model
made by Anthropic.

  paper/              Section 2.8, "Further reductions, and the conjecture
                        near the contact regime": the hole lemma, the
                        inversion hull (every centre outside a set Y leaves
                        the convex hull of 0 and the points 4y/|y|^2 whole),
                        contacts that contain most of a root system, no
                        three centres cutting the same point of
                        B(sqrt(3/2)), at most 22 centres within sqrt 6, the
                        case of 23 contacts, and Theorem 2.18: there is
                        epsilon_0 > 0 such that if no centre lies at a
                        distance between 2 + epsilon_0 and sqrt 6, the cell
                        has volume at least 8, with equality only at the
                        24-cell -- the D_4 configuration is a strict local
                        minimum among all packings.  Summarised in the
                        introduction as Theorem 1.7.  Two new figures
                        (fig_closure, drawn by figures_new/fig_closure.py;
                        fig_tikz_closure, drawn in TikZ by
                        figures_new/tikz_closure.tex).  The long note on the
                        title page is removed; Deep Bhattacharjee and
                        Priyabrata Mandal are both corresponding authors.
                        What v1.5.0 left open: a centre with at least 23
                        centres within sqrt 6, one of them at a distance
                        between 2 + epsilon_0 and sqrt 6 (narrowed to 24
                        in v1.6.0).
  multi_cap/closure_lemmas.py
                      The exact content of Section 2.8, in sympy, integer
                        and ball arithmetic.
  multi_cap/near_contact_probe.py
                      Exploration near the root system (not proof).
  lean/D4Closure.lean The same exact content in Lean 4: the identities over
                        every commutative ring (grind), the root-pair
                        combinatorics by decide +kernel.  run_all.sh now
                        checks six files.

## What is new in v1.4.0

An independent re-verification of every computation of the paper, and the
corrections it led to, in the paper and in the package.  The re-verification,
the two new programs of multi_cap/ below and the corrections were developed
with the assistance of Claude, an AI model made by Anthropic.

  independent_verification/
                      A re-run of the whole package from a fresh copy: the
                        seven Lean verifications, both lake builds, the
                        certificates for 23 and 24 contacts (steps 3 and 5
                        of the certificate of Section 7.5.3 included, all
                        fifty sum-of-squares blocks), every script against
                        its shipped log, the figures and the manuscript.
                        REPORT.md lists what was found and corrected; logs/
                        holds every run; proof_gaps/ the two counterexamples
                        below.
  paper/              The reduction from the Voronoi cell of a packing to the
                        cell of its contacts is corrected (Section 2.4, 2.6
                        and the new Section 2.7, "Neighbours that do not
                        touch").  Shell localisation and radial reduction, as
                        stated in v1.3.0, fail when a neighbour does not
                        touch the centre: a centre at distance 2.9 cuts the
                        cell of a deletion (8.3028 instead of 25/3), and the
                        all-contact corner of a packing of 49 balls has volume
                        6.594 while its true cell has 19.216.  The paper now
                        proves the bound for contact configurations in full,
                        and for the Voronoi cell of a packing whenever the
                        centres within 2 sqrt 2 all touch it or their
                        distances pass a covering criterion (Proposition
                        2.10, Theorem 1.5); the case of a centre crowded by
                        near-contacts is stated as Conjecture 1.6, with the
                        evidence of a numerical search.  Theorem 1.1 is
                        restated for the cell of the active neighbours with
                        the hypothesis its proof uses.  Title-page note,
                        abstract, introduction, Section 6, the conclusion,
                        the data and code availability statement and one
                        new figure (fig_noncontact) follow.
  multi_cap/shell_reduction.py
                      The distance criterion in ball arithmetic: covering
                        radius 45 degrees of the root directions, Phi above
                        8.044 for every m <= 22 contacts, the thresholds of
                        Table 1, the two counterexamples.
  multi_cap/shell_neighbour_search.py
                      The numerical search of the open case (exploration,
                        not proof); its runs are in multi_cap/runs.
  zonal/run_sos4.sh   Runs all fifty sum-of-squares blocks of constraint 4,
                        refuses to report success unless every block ran,
                        and exits 1 on a failure.

## What is new in v1.3.0

  paper/              The manuscript (D4.tex, with D4.pdf beside it in a
                        tagged release) and its figures.  The
                        twenty-four-contact case is now proved in the text:
                        the bound of the second level of the hierarchy with
                        its equality case, the positivity of Gram kernels,
                        and the finite properties of the deposited certificate
                        (Section "Twenty-four contacts, and what is left").
                        No statement of an unrefereed source is used.
  paper/figures_new/  The matplotlib scripts of the data figures, each with a
                        label collision test that it has to pass before it
                        writes (figstyle.py); the ray-traced panels of the
                        composites are reused unchanged (oldpanels.py, old/).
                        make_p2_pickle.py builds p2.pkl, the exact
                        coefficients of the two-point polynomial, from
                        multi_cap/llm24_out/llm24_p2.txt.
  third_party/llm24-certificate/
                      The certificate data set of de Laat, Leijenhorst and
                        de Muinck Keizer (MIT licence), redistributed with
                        its licence notice, a fetch script that checks the
                        MD5 recorded in the paper, and verify.sh, which runs
                        the seven-step verification.  See its README.
  lean/D4RootLattices.lean
                      A seventh Lean verification: the finite content of the
                        lemma from the four inner products to the D_4 root
                        system (see lean/README.md); lean/run_all.sh checks
                        the five standalone files.
  multi_cap/p2_zeroset_check.py
                      A standalone exact-arithmetic check (Fractions, Sturm
                        sequences) of the zero set of the two-point
                        polynomial, independent of python-flint.

## Directory Structure

  paper/              The manuscript, its figures and their scripts
  independent_verification/
                      The re-verification of v1.4.0 (see above)
  third_party/        The redistributed certificate data set (see above)
  cap_certificate/    The cap inequality, end to end (Section 13)
  multi_cap/          The multi-cap, polar and boundary reformulations of
                        the 23-point case, the covering bound, the
                        extendability criterion, the rigidity of the
                        deletion configuration, the failure of both
                        local decompositions, the root-system meets,
                        the codes of the 600-cell, the continuation
                        in the slack, the cell inside a ball and the
                        three-point certificate that settles the case
                        (Section 7.2 to VII L)
  data/               Cached intermediate results (.pkl) for the longer
                        symbolic derivations, so that the consuming
                        scripts can be run without first re-running the
                        step that produced them
  core/               Core arithmetic and Hessian routines
  hessian_multidir/   Joint Hessian computations (Sections 15-19)
  arc1_v1w1/          First-arc certificate scripts (Section 10)
  arc2_w1v2/          Second-arc certificate scripts (Section 12)
  swap_configs/       Swap-configuration finite-angle proofs (Section 20)
  verification/       Independent cross-checking scripts
  zonal/              The two steps of the LLM24 certificate that the
                        authors' own implementation needs 128 GB for:
                        the zonal matrices and the four polynomial
                        identities that use them (see zonal/README.md)
  misc/               Utility and diagnostic scripts
  lean/               Eight Lean 4 verifications: the equilibrium stress of
                        Proposition 7.48 and of the corollary after it, with
                        the surrounding root-system combinatorics
                        (D4Stress.lean); the finite half of
                        Propositions 7.55 and 7.56 on how much of a root
                        system a contact configuration can hold
                        (D4Meet.lean); the exact half of the certificate
                        of Theorem 7.73, positivity of its matrices and
                        the value of its bound over the rationals
                        (D4Certificate.lean, generated from the
                        certificate by gen_certificate_lean.py); the
                        zeros of the two-point polynomial of the
                        certificate of de Laat, Leijenhorst and de Muinck
                        Keizer, the last step of its verification
                        (D4InnerProducts.lean, generated from
                        multi_cap/llm24_out/llm24_p2.txt by
                        gen_innerproducts_lean.py); the finite content of
                        the lemma from the four inner products to the D_4
                        root system, the census of norm-2 Gram matrices in
                        rank at most 4 (D4RootLattices.lean); and, in
                        the Lake project cell600/, the enumeration behind
                        Proposition 7.63, that every 23-point code of
                        minimal angle 60 degrees among the vertices of the
                        600-cell is an inscribed 24-cell less a vertex;
                        and, in the Lake project certificate/, the interval
                        branch and bound of Theorem 7.73 re-done in exact
                        dyadic arithmetic (D4CertDomain.lean), the
                        polynomial expanded inside Lean and only the
                        tables of bounds for omega and its derivatives
                        taken from certificate_check.py (D4CertData.lean,
                        generated by gen_data.py), the theorem settled by
                        native_decide (D4CertMain.lean).
                        D4Closure.lean (v1.5.0, extended in v1.6.0)
                        proves the identities of Section 2.8, those of
                        Theorem 2.17 included, over every commutative ring
                        and its root-pair combinatorics and counts by
                        kernel computation.
                        The six files need no Mathlib and no lakefile:
                        run "lean D4Stress.lean", "lean D4Meet.lean",
                        "lean D4Certificate.lean",
                        "lean D4InnerProducts.lean",
                        "lean D4RootLattices.lean" and
                        "lean D4Closure.lean" (or "sh run_all.sh")
                        with the toolchain pinned in lean-toolchain.
                        The two projects are built with
                        "lake build" inside lean/cell600 and
                        lean/certificate (no Mathlib either). See
                        lean/README.md

## hessian_multidir/ -- file-by-file map to paper sections

  extend_hessian.py, c2_constant_test.py, c2_hessian_relation_check2.py,
  triality_search.py, triality_search2.py, triality_exact.py,
  chain_length.py, m3_hessian.py
                              Earlier exploratory and chain-configuration
                              Hessian scans (Sections 15, 18).

  joint_hessian_closed_form.py
                              Derivation of the m=2 adjacent-pair cross-Hessian
                              closed form used in Section 15. Recorded
                              in-file and in the paper (Section 15, "A first
                              step toward a first-principles derivation for
                              m=2") as obtained by numerical fitting, not a
                              first-principles symbolic derivation.

  vertex_degeneracy_check.py
                              Exact vertex-enumeration check referenced in
                              Section 15: confirms two adjacent facets of the
                              24-cell share exactly 3 of their 6 vertices, and
                              examines the vertex-degeneracy reduction (which
                              4 of 10 candidate 3-subsets of a vertex's other
                              active facets remain jointly feasible after one
                              facet is perturbed).

  multidir_exact_zero_hessian_hp.py
                              High-precision (35-45 digit mpmath) three-step-size
                              scaling argument establishing the exact quartic
                              degeneracy at the m=18 configuration A_18
                              (Section 17, "exact singularity"): confirms the
                              near-null subspace is genuinely 4-dimensional
                              (eigenvalues shrink by a clean factor of ~4 per
                              halving of the step size, the signature of an
                              exact zero rather than a small positive floor).

  multidir_nullspace_broad_sample.py, nullspace_broad_sample_A18.log
                              Broadened directional sampling at A_18 (Section
                              16): 36 directions spanning the full 4D near-null
                              subspace (4 basis eigenvectors, 6 pairwise sums,
                              6 pairwise differences, 20 random directions),
                              quartic-coefficient estimator evaluated at two
                              independent step sizes. Log file is the actual
                              run output quoted in Theorem
                              "Broadened sampling at A_18".

  multidir_nullspace_broad_sample_m20.py
                              Same broadened-sampling method applied to the
                              independent m=20 configuration A_20 (Section 18):
                              confirms a 5-dimensional near-null subspace and
                              samples it with 30 directions, all strictly
                              positive.

  quartic_sample_dense.py, quartic_dense_A18_data.tsv,
  quartic_dense_A18_run.log
                              A further-overdetermined directional dataset at
                              A_18 (166 directions, 4.7x overdetermined vs the
                              35 unknowns of a general quartic form in 4
                              variables), superseding the 36-direction dataset
                              above for the purpose of the least-squares fit
                              below. Resumable by design (checks its own
                              output file and only computes missing
                              directions) since this dataset was originally
                              collected across several cloud-sandbox container
                              restarts.

  gram_sos_lib.py             Standard Gram-matrix/semidefinite-programming
                              method (via cvxpy) for testing whether a quartic
                              form in n variables admits a sum-of-squares
                              certificate. Self-validates against two textbook
                              quartics when run directly (`python
                              gram_sos_lib.py`): a perfect square (must find a
                              certificate) and the Choi--Lam polynomial (must
                              NOT find one -- the classical 1977 example of a
                              non-negative quaternary quartic with no SOS
                              representation).

  quartic_fit_and_check.py    Fits the 35 coefficients of the general quartic
                              form at A_18 by least squares from
                              quartic_dense_A18_data.tsv, then runs
                              gram_sos_lib.py's SDP check on the fitted
                              quartic (Section 18, "A Gram-matrix sum-of-
                              squares check"). Reports the fit residual and
                              the SDP result plainly; states in its own output
                              exactly what the result does and does not
                              establish about Conjecture 1.4.

  multidir_chain_hessian_extended.py, hp_volume.py
                              Shared dependencies (root-system construction,
                              tangent-space bases, and the high-precision
                              volume routine) needed by every script above;
                              duplicated here from core/ so this directory is
                              runnable on its own.

  a18_stabilizer_group.py    Identifies A_18's own exact symmetry. Within the full 384-element
                              signed-coordinate-permutation automorphism
                              group of the D4 root system (exact integer
                              arithmetic throughout -- every element is a
                              0,+-1 monomial matrix), finds the subgroup
                              stabilising A_18 as a SET of 18 root indices
                              has order exactly 48 (a hyperoctahedral
                              B_3-type subgroup: signed permutations of 3
                              coordinates, with one coordinate's identity
                              and sign left untouched). This is Proposition
                              "a18-stabiliser", Proposition 18.5 of the paper.

  a18_symmetry_representation_check.py
                              Verifies the order-48 group is not a
                              coincidence of index labels but a genuine
                              symmetry of the actual dynamics at A_18: (1)
                              builds the induced 54-dimensional
                              representation of the group on the joint
                              tangent space and confirms it commutes with
                              the finite-difference joint Hessian to
                              ~4e-6 (consistent with the h=0.02 step size
                              used, not a real discrepancy); (2) confirms
                              the quartic-coefficient estimator a4 is
                              invariant under the group's action on the
                              confirmed 4-dimensional near-null subspace
                              (a4(g.v) agrees with a4(v) to ~1e-5-8e-5
                              across 6 sampled non-identity elements); (3)
                              computes the character of the group's action
                              on that 4D subspace and the resulting
                              commutant dimension (exactly 2), showing the
                              near-null subspace splits into exactly two
                              non-isomorphic irreducible pieces under this
                              group (their exact dimensions are not
                              identified here). Does NOT prove positivity
                              of the quartic form; narrows what an eventual
                              exact derivation's ansatz would need to cover.

  a18_exact_character_table.py
                              Builds the complete, exact
                              (integer/rational arithmetic, no floating
                              point) character table of the order-48
                              stabiliser group -- all 10 conjugacy classes,
                              all 10 irreducible characters constructed
                              directly from the group's natural
                              representations (triv, two independent sign
                              characters and their product, the defining
                              3-dim representation std3 and its three
                              twists, and a 2-dim representation pulled
                              back from S3 and its twist), verified pairwise
                              orthonormal. Decomposes the EXACT character of
                              the 54-dimensional joint tangent
                              representation (via the closed form
                              chi_54(g) = n_fix(g)*(trace_4x4(g)-1)) into
                              these 10 irreducibles, confirming all
                              multiplicities are non-negative integers
                              summing to 54 -- a self-checking computation,
                              since a first attempt at embedding the
                              abstract group back into 4x4 (assuming the
                              wrong fixed coordinate) gave impossible
                              non-integer multiplicities, caught by exactly
                              this check.

  a18_nullspace_irrep_match.py
                              Recomputes the 54x54 joint
                              Hessian and its 4-dimensional near-null
                              subspace, tags all 48 stabiliser elements
                              (not a sample) by their EXACT conjugacy class,
                              and compares the per-class near-null
                              character against the two candidate
                              irreducible-pair decompositions that matched
                              the earlier (aggregate-only) trace
                              distribution. Finds perfect per-class
                              constancy and an exact, unambiguous match to
                              eps_perm (+) (std3 (x) det) -- resolving an
                              ambiguity the aggregate data alone could not.
                              This is Proposition "a18-nullspace-irrep",
                              Proposition 18.6 of the paper.

  a18_invariant_quartic_basis.py
                              Proves, via exact power-sum
                              (Newton's-identity) character formulas
                              (integer/rational arithmetic throughout, no
                              floating point), that the space of quartic
                              forms invariant under the full order-48 group
                              has dimension exactly 5 (not the general 35
                              for an unconstrained quaternary quartic), and
                              constructs an explicit basis
                              {x0^4, x0^2*S2, x0*P3, S4, S22} independently
                              via exact Reynolds-operator projection,
                              confirming the two methods agree. This is
                              Proposition "a18-invariant-dim" (with its
                              proof), Proposition 18.7 of the paper.

  a18_fit_and_sos_check.py
                              End-to-end numerical
                              pipeline building on the three scripts above.
                              Constructs an orthonormal basis of the actual
                              54-dimensional ambient tangent space realising
                              eps_perm(+)(std3 (x) det) concretely (Reynolds
                              projectors onto the two isotypic pieces, then
                              an exact intertwiner found via SVD nullspace
                              of a Sylvester-type equation), evaluates the
                              quartic-coefficient estimator at 7 directions
                              isolating the 5 basis invariants using the
                              project's high-precision (mpmath, 40-digit)
                              volume routine at two independent step sizes
                              (sigma=0.02, 0.01), fits c1..c5 by least
                              squares (both step sizes separately and
                              combined), and re-runs the Gram-matrix SOS
                              check on all three fits. Reports plainly that
                              an initial attempt using plain double
                              precision gave fits disagreeing by 100-700%
                              between step sizes -- discarded, not reported
                              as reliable, exactly the failure mode
                              Section 16 (musin-degeneracy) already
                              documents for this configuration. This is
                              Proposition "a18-sos-reduced", Proposition 18.8
                              of the paper. Does NOT prove positivity: the
                              fitted coefficients remain finite-difference
                              numerical estimates, not an exact symbolic
                              derivation, and the SOS margin, though larger
                              and more stable than the unconstrained
                              35-coefficient fit, remains thin.

## arc2_w1v2/ -- the symmetry of the fundamental triangle

  third_edge_symmetry.py     Resolves the open question left in the paper's
                              Section 8.3 ("Where the worst
                              direction sits"): identifies and
                              verifies, both by exact sympy/Q(sqrt2)
                              computer algebra and as a numerical sanity
                              check, the isometry S = diag(1,1,1,-1)
                              mapping the fundamental triangle's third edge
                              (v1-v2, off both named arcs) exactly onto the
                              already-certified first arc (v1-w1). This is
                              Lemma "third-edge-symmetry" and Proposition
                              "third-edge-cert" in Section 8.4 of the paper.

  full_triangle_symmetry.py  The bigger result this led to: the group
                              generated by S (above) and H (the Hadamard
                              matrix already used in misc/confirm_H_full_arc_symmetry.py)
                              has order 6, fixes u0 throughout, and realises
                              ALL SIX permutations of {v1,w1,v2} -- i.e. it
                              is the full symmetric group S_3. In particular
                              M=S*H maps the SECOND arc's own parametrisation
                              exactly onto the first arc's (proved
                              symbolically for symbolic t, not just sampled
                              points, then cross-checked numerically across
                              every combinatorial sub-part of the second
                              arc, including the B-bounded band). This is
                              Lemma "s3-symmetry", Proposition
                              "arc2-is-arc1", and Corollary
                              "restated-conjecture-holds" in Section 8.5 of
                              the paper: it proves Conjecture
                              (Direction-of-Deviation Positivity, RESTATED)
                              in full on the fundamental triangle's
                              boundary, superseding the direct
                              region-by-region attempts on the second arc
                              elsewhere in this directory (region_top_*.py,
                              bandC_certificate.py,
                              bandBcurve_certificate.py -- kept in the
                              package as an independent, partial
                              cross-check and as a record of how
                              this connection was found; the
                              bandBcurve_certificate.py background run was
                              stopped once this result was in hand).

                              Neither script touches the fundamental
                              triangle's two-dimensional interior. That is
                              settled instead by cap_certificate/, which
                              covers the entire transverse 2-sphere at once
                              and does not use the triangle at all.

## arc2_w1v2/ -- the closed-form cap bound

  cap_reformulation_check.py The supporting script for
                              the paper's strongest single-deviation
                              result. The paper's Section "A closed-form
                              cap bound" proves, with no floating-point
                              step anywhere, that the defect equals
                              1/3 minus the volume of the cap cut from a
                              fixed 25-vertex polytope Q by the deviated
                              half-space; that this cap is bounded in
                              closed form by an exact integration over
                              the cone K with apex 2*u0 over the
                              octahedral facet,
                                 cap <= (1/3)(2 - sec t)^4
                                        / prod_i (1 - tan^2 t * a_i^2);
                              and that the resulting inequality is
                              strictly positive for EVERY deviation
                              direction on the full transverse sphere at
                              every tilt below arccos(c0) = 42.1759...
                              degrees, c0 the real root of
                              14c^3 - 18c^2 + 7c - 1. Those are exact
                              proofs, not computations. This script
                              exists only so that a reader can confirm
                              the objects are what the paper says they
                              are: it checks vol(Q) = 25/3, that
                              Q above the facet plane is exactly the
                              pyramid of volume 1/3, that the facet is a
                              regular octahedron of volume 4/3 with the
                              orthogonal roots as vertices, that the
                              closed form agrees with a direct polytope
                              volume, that the lower bound really is a
                              lower bound at random (tilt, direction)
                              pairs, that the cubic threshold is right,
                              and that the pyramid's cap is empty above a
                              right angle. Its final item -- the
                              numerical maximum 0.058874... of the
                              24-cell's own cap, a factor 5.66 below the
                              1/3 that would be needed -- is reported in
                              the paper as evidence about the UNCLOSED
                              range (tilt >= 42.18 degrees) and is
                              labelled as such both there and in this
                              script's own output.

## arc2_w1v2/ -- the interior investigated directly

  interior_concavity_probe.py, interior_concavity_convergence_test.py
                              Tests, and rules out, the most natural
                              candidate for closing the triangle's INTERIOR
                              for free now that its boundary is proved:
                              whether defect(theta, .) is concave along
                              every geodesic through the triangle (which
                              would bound any interior point below by the
                              boundary values). The probe script finds small
                              positive second differences (apparent
                              non-concavity) at several sampled points; the
                              convergence-test script runs one of them down
                              properly, via Richardson-style step-size
                              shrinkage across 8 orders of magnitude
                              (h=0.08 to h~3e-4): d2(h)/h^2 converges
                              cleanly to +0.3056, confirming a genuine
                              positive second derivative (local convexity),
                              not grid noise and not a combinatorial
                              breakpoint (vertex count confirmed constant
                              at 37 throughout the tested window). This is
                              Section 8.6's falsified-shortcut result: a
                              real (negative) finding, not an unattempted
                              guess.

  interior_region_top_stability.py, interior_combinatorial_map.py
                              Maps the triangle's own combinatorial
                              structure directly. The stability script
                              shows the boundary's simple 32-vertex "region
                              top" type does NOT persist into the interior
                              (most interior points give 34, 37, or 40
                              vertices instead). The map script grids a
                              full triangle at theta=1.45 (861 points) and
                              finds (at least) four distinct region types
                              -- 32 (25%, edge-adjacent), 34 (17%), 37 (53%,
                              the plurality "generic" type), and 40 (5%,
                              a central island) -- with the partition
                              visibly respecting the S_3 symmetry of
                              full_triangle_symmetry.py (in particular the
                              reflection fixing v1 and swapping w1,v2).
                              This is Section 8.6's structural map: a real
                              scope reduction (one-sixth of the triangle
                              suffices, by the S_3 action, to know the
                              whole interior's combinatorial structure) but
                              not a certificate -- no volume formula or
                              positivity proof is attempted for any of the
                              four region types.

## cap_certificate/ -- the main theorem, in one script

  cap_inequality_certificate.py
                              Reproduces the whole of Section 13 end to
                              end, in the order the section proves it:

                                (A) exact rational vertex enumeration of
                                    Q, giving 25 vertices and vol(Q)=25/3,
                                    with the pyramid of volume 1/3
                                    recovered as the difference from the
                                    24-cell;
                                (B) the exact check that every vertex of Q
                                    has l^1 norm at most sqrt(2) in the
                                    orthonormal root frame, which is the
                                    enclosure Q inside B, together with the
                                    fact that the bound is attained;
                                (C) the closed form for the cap of the
                                    cross-polytope as a divided difference
                                    of a(2 sqrt a - 1)_+^4, checked against
                                    directly computed polytope volumes at
                                    60 random directions;
                                (D) the symbolic identity for the fourth
                                    derivative, 3(20t^2-3)/(2t^5) at
                                    a = t^2, and the vanishing of the third
                                    derivative at a = 1/4, which together
                                    give the monotonicity of the divided
                                    difference in each node;
                                (E) exact real-root isolation for the three
                                    polynomial inequalities covering the
                                    case of at most one node above 1/4;
                                (F) the adaptive subdivision covering the
                                    remaining case: 303 boxes, corner
                                    bounds entirely in fractions.Fraction,
                                    largest value 0.99755050;
                                (G) an independent check of the conclusion:
                                    directly enumerated Voronoi cell
                                    volumes at 300 random (tilt, direction)
                                    pairs, all with nonnegative defect.

                              Steps (A), (B), (D), (E) and (F) are exact;
                              (C) and (G) are double-precision
                              cross-checks that confirm the exact work and
                              are not relied on by it. The script prints
                              PASS or FAIL for each of its 21 checks and
                              exits nonzero if any fails. Runtime under
                              ten seconds.

## multi_cap/ -- the 23-point case, reformulated and settled

  closure_lemmas.py
                              Supports Section 2.8 (Lemmas 2.12 and 2.15,
                              Proposition 2.13, Corollaries 2.14 and 2.16,
                              and the angle quoted after Theorem 2.17).
                              Symbolic checks of the identities,
                              integer enumeration of the orthogonal root
                              pairs at the vertices of the 24-cell and of
                              every set of at most four deleted roots, and
                              ball arithmetic for S(2), the bound
                              9 pi^2/8 - 22 S(2) > 8.046, the thresholds of
                              Psi beside those of Phi, and arccos(sqrt6/4).
                              A few seconds.  Log: runs/closure_lemmas.log.
                              The same exact content is in
                              lean/D4Closure.lean.

  labelled_certificate_check.py
                              Supports Theorem 2.17 (at most 23 centres
                              within sqrt 6).  Rebuilds the certificate of
                              Theorem 7.73 and its bound exactly, rechecks
                              the positivity of its matrices, evaluates the
                              constants s(D), a_D, fr(1, 1, 1/2), r, kappa
                              and c in ball arithmetic (python-flint), reruns
                              the branch and bound of (C), and runs the two
                              branch and bounds of the ranges II_s and II_f
                              of t, with the tables of the Gamma_i in ball
                              arithmetic.  About six minutes.  Log:
                              runs/labelled_certificate_check.log.  Option
                              --skip-region-1 leaves out the rerun of (C).

  explicit_eps0.py DATA      Supports Theorems 7.76 and 2.18 (the explicit
                              epsilon_0 = 4e-26).  DATA is the folder
                              proofs/4_24 of the certificate in
                              ../third_party/llm24-certificate.  Bounds the
                              sum-of-squares terms of the certificate whose
                              weight can change sign (ball arithmetic,
                              python-flint), divides sigma_2 exactly by its
                              zeros and bounds the quotient below, checks
                              the minor bounds of Lemma 7.77, and evaluates
                              the facet and hull estimates of Theorem 2.18
                              in exact arithmetic.  Sixteen seconds.  Log:
                              runs/explicit_eps0.log.

  truncated_search.py M [starts] [seed]
  truncated_search.py rays
                              Supports "What is left" after Theorem 2.18.
                              Floating point, exploration: minimises the
                              right side of Lemma 2.15 over configurations
                              of exactly M centres within sqrt 6.  Logs:
                              runs/truncated_search_M24.log to _M27.log.
                              With "rays" it prints where that right side
                              reaches 8 on two families through the root
                              system (sum of delta_i = 0.155 pushed out
                              evenly, delta = 0.197 for one centre).  Log:
                              runs/truncated_search_rays.log.

  near_contact_probe.py
                              Supports the remarks after Theorem 2.18.
                              Floating point, exploration: the inversion
                              hull beside the volume, and 200 volume
                              minimisations near the root system.  Log:
                              runs/near_contact_probe_200_seed5.log.

  shell_reduction.py
                              Supports Section 2.7 (Proposition 2.10,
                              Lemma 2.11, Theorem 1.5, Table 1). The
                              covering bound with every centre within
                              2 sqrt 2 counted by the cap it cuts, as a
                              function Phi of the distances alone, stopped
                              at 45 degrees so that no centre beyond
                              2 sqrt 2 enters.  Checks in ball arithmetic
                              (python-flint): the covering radius of the
                              root directions is 45 degrees; Phi > 8.044
                              for m <= 22 contacts; the thresholds of
                              Table 1; the two counterexamples to the
                              uncorrected Lemmas 2.7 and 2.9.  About 50 s.
                              Log: runs/shell_reduction.log.

  shell_neighbour_search.py
                              Supports the paragraph "The open case,
                              searched" of Section 2.7 and Figure 3(c).
                              Minimises the exact volume of the cell, with
                              its exact gradient, over configurations of
                              22 to 26 centres under the packing
                              constraints (SLSQP); optionally one
                              neighbour held at distance 2 + delta.
                              Floating point: exploration, not proof.
                              Logs: runs/shell_neighbour_search_*.log.

  root_deletions_exact.py
                              Supports Propositions 7.55 and 7.56.
                              Exact integer vertex enumeration of the cell
                              left when j = 1, 2, 3 roots are removed from
                              D_4. Scaling the roots to integer vectors of
                              squared length 2 makes every vertex of the
                              cell rational, of the form m/e with m
                              integral and e a positive integer, so the
                              circumradius test is the integer comparison
                              |m|^2 <= 2 e^2 and the test for a direction
                              that can be added is |m|^2 = 2 e^2.
                              Reports, over all 24, 276 and 2024 cases:
                              circumradius exactly 2 with the removed roots
                              as the only attaining directions, except at
                              the 96 triples pairwise at 60 degrees, where
                              the circumradius is sqrt6.
                              No floating point anywhere.

  cell600_exact.py
                              Establishes Proposition 7.63. The 120
                              vertices of the 600-cell with doubled
                              coordinates in Z[phi], every inner product
                              computed exactly in Z[phi]; two vertices are
                              closer than 60 degrees exactly when joined by
                              an edge (36 degrees), so the subsets with all
                              inner products at most 1/2 are the
                              independent sets of a 12-regular graph. The
                              script finds the 25 inscribed 24-cells as the
                              24-cliques of the graph of root-system
                              angles, then enumerates completely the
                              independent sets through a fixed vertex:
                              none of size 25, five of size 24 (the cells
                              through that vertex), and 115 of size 23,
                              every one inside a cell, so 115 x 120 / 23 =
                              600 = 25 x 24 in all. Writes the graph and
                              the cell list to cell600_graph.txt. Runtime
                              about ninety seconds; no floating point
                              anywhere. Run with --sat for an additional
                              SAT cross-check through python-sat (slow).

  cell600_enum.c
                              Independent C implementation of the size-23
                              enumeration of cell600_exact.py, on 64-bit
                              bitsets, reading cell600_graph.txt:
                                cc -O2 -o cell600_enum cell600_enum.c
                                ./cell600_enum
                              Prints 115, 600 and 0 and PASS in under two
                              seconds. The Lean project lean/cell600/ is
                              the third implementation.

  cell600.py
                              The sampling that preceded the enumeration:
                              two million maximal independent sets of the
                              same graph drawn by greedy growth along
                              random orders, sizes reported (10 to 22 and
                              24, never 23). Kept for the record;
                              superseded by cell600_exact.py.

  inradius_search.py
                              Supports Section 7.10. Maximises the inradius
                              g(W) of conv(W), the cosine of the covering
                              radius, over 23-point configurations with all
                              inner products at most 1/2, by a trust-region
                              sequential linear programme with analytic
                              gradients: the facet offsets of the convex
                              hull and the pairwise inner products are
                              linearised in a tangent move of bounded size,
                              the LP maximising the smallest linearised
                              offset is solved (HiGHS through scipy), and
                              the move is kept only if the true inradius,
                              read from the facet equations of a fresh
                              hull, improves without any contact constraint
                              being violated. Three kinds of start in
                              rotation: random points brought to
                              feasibility through a relaxed-then-tightened
                              schedule; perturbed deletions; twenty roots
                              plus three random directions. Reports every
                              feasible endpoint's g and the multiset of its
                              inner products, and flags any g > 1/2.
                              Exploration, double precision throughout.
                                python inradius_search.py [starts] [seed]

  octahedral48_exact.py
                              The 48 unit quaternions of the binary
                              octahedral group (two 24-cells in dual
                              position), coordinates doubled in Z[sqrt2];
                              enumerates every subset of size 23 and 24
                              with inner products at most 1/2 and finds
                              the two root systems and their 48 deletions
                              only. Exact; a few seconds.

  three_point_reduction.py
                              Supports Section 7.11 and the r_* forms of
                              Propositions 7.23, 7.41 and 7.44. Part 1:
                              the pair-only truncated-volume bound at the
                              limits r_23 and r_* (bracket, weight of a
                              60-degree pair, value at a deletion, pairs
                              needed for 8). Part 2: the three-point bound
                              at r_4 = arcsin sqrt(3/8), with the triple
                              cap measure by Monte Carlo, at a deletion
                              (8.140848) and at the root system
                              (7.968684). Part 3: the pair-angle linear
                              relaxation of the bound for both weights
                              (0.04157 of 0.08327 at r_23; 0.07380 of
                              0.09286 at r_*) and the minimising measure.
                              Under a minute (the Monte Carlo sample is seeded).

  truncated_volume.py
                              Supports Section 7.11. Evaluates the volume
                              of the cell inside a ball of radius R by a
                              fixed quasi-random quadrature on S^3
                              (deterministic, 400000 points), checks it
                              against the exact volumes at the root system
                              and a deletion, then minimises it over
                              23-point configurations with inner products
                              at most 1/2 + delta on a decreasing schedule
                              of delta, by the trust-region sequential LP
                              of inradius_search.py with the analytic
                              gradient of the quadrature. Prints, per
                              level, the least truncated volume found, the
                              exact volume of that configuration, and the
                              number of endpoints. Exploration.
                                python truncated_volume.py [starts] [seed] [points] [R]
                              (default R = sqrt(8/5); R = 1.224744871391589
                              for sqrt(3/2)). About two hours per run.

  three_point_sdp.py
                              Supports Section 7.12: the three-point
                              (Bachoc-Vallentin) relaxation of the pair
                              inequality of Theorem 7.69, and the
                              certificate of Theorem 7.73. Builds the
                              Gegenbauer polynomials of S^3, the matrices
                              Y_k for n = 4 (Legendre polynomials, in the
                              Chebyshev basis T_i(u) T_j(v)) and their
                              symmetrisation, imposes the condition (C) of
                              Lemma 7.72 on a grid of admissible triples
                              and solves the semidefinite programme with
                              cvxpy and Clarabel. First mode: maximise the
                              bound, check on about 1.15 million further
                              triples, add the worst 3000 and repeat;
                              reports the bound reduced by the largest
                              violation (degree 6: 0.09011, degree 8:
                              0.09523, against the target 0.09286; the
                              two-point programme alongside: 0.0712).
                              Second mode (fifth argument): fix the bound
                              and maximise the least slack; writes
                              continuation_out/certificate_d8.npz (and a
                              plain-text copy, certificate_d8.txt, with
                              every number printed exactly). Needs cvxpy
                              and Clarabel (pip install cvxpy clarabel).
                              One to three minutes per solve, a quarter
                              of an hour for five rounds.
                                python three_point_sdp.py 8 30 CLARABEL 5
                                python three_point_sdp.py 8 30 CLARABEL 5 0.0929

  certificate_check.py
                              The proof of Theorem 7.73: verifies the
                              certificate in exact rational and interval
                              arithmetic, sharing no code with the solver.
                              Step 1, exact LDL^T of the nine matrices and
                              f_k >= 0; step 2, the bound computed exactly
                              against 8 - A_* from the closed form of A_*
                              in mpmath interval arithmetic; step 3, the
                              polynomial P expanded exactly (449 monomials,
                              degree 16, symmetric), omega and its two
                              derivatives tabulated from closed forms in
                              interval arithmetic, and the inequality (C)
                              verified on the ordered admissible domain by
                              an interval branch and bound with the
                              second-order Taylor form on each box (313780
                              boxes verified, 13517 discarded, 44 levels,
                              about six minutes). Stops at the first
                              failure with the box or the counterexample.
                              Needs numpy, mpmath, sympy.
                                python certificate_check.py 8 1e-6

  root_lattices_rank4.py
                              Supports Lemma 7.26, the combinatorial half
                              of the twenty-four-point classification.
                              Enumerates every positive definite Gram
                              matrix with 2 on the diagonal and -1, 0, 1
                              off it of order 1 to 4 (1, 3, 23, 393 of
                              them), counts the integer solutions of
                              x^T G x = 2 within the rigorous box
                              |x_i| <= 15 by exact integer arithmetic,
                              and reports the largest counts 2, 6, 12, 24,
                              the count 24 occurring only at determinant
                              4 with the neighbour profile of D_4 (next:
                              20 at determinant 5, A_4). About a minute.
                              Log: runs/root_lattices_rank4.log.
                                python root_lattices_rank4.py

  llm24_certificate_check.py
                              Supports Theorem 7.25, Proposition 7.29 and
                              Section 7.5.3: an independent verification,
                              sharing no code with the authors' Julia
                              package, of the certificate of de Laat,
                              Leijenhorst and de Muinck Keizer (data set
                              doi:10.4121/74ce1c25-6fca-4680-8a36-e9c18e7e9594,
                              LasserreSphericalCodes.zip, md5
                              02acd5270f7b3fa799abdeb5291706fd, redistributed
                              in third_party/llm24-certificate; pass the
                              path of its proofs/4_24 folder). Of the seven steps of
                              their verification it repeats five: the
                              reading and format of the data (127 blocks
                              of total dimension 3726, the block of every
                              signature having exactly as many rows as
                              there are admissible tuples); positive
                              definiteness of all 127 blocks by Cholesky
                              in ball arithmetic at 256 bits (python-flint
                              arb, every pivot a ball inside the positive
                              reals; least pivot about 1.38e-15), the 81
                              blocks of size at most 16 also by exact
                              rational LDL^T; the structure of all 125
                              sum-of-squares prefactors (each a
                              nonnegative multiple of a domain weight);
                              the objective K(empty, empty) = 24 exactly;
                              and the two-point polynomial p_2, degree
                              16, computed exactly from the data and
                              shown by a Sturm sequence over Q to vanish
                              on [-1, 1/2] exactly at -1, -1/2, 0, 1/2
                              (multiplicities 1, 2, 2, 1). Steps 3 and 5,
                              the construction of the zonal matrices and
                              the polynomial identities that use them,
                              are done in zonal/ (see zonal/README.md).
                              Writes
                              llm24_out/llm24_p2.txt, from which
                              lean/gen_innerproducts_lean.py generates
                              lean/D4InnerProducts.lean. Needs
                              python-flint. About eight minutes.
                              Log: runs/llm24_certificate_check.log.
                                python llm24_certificate_check.py /path/to/LasserreSphericalCodes/proofs/4_24

  run_llm24_full_verification.sh
                              The computation this package does not
                              contain: on a machine with 128 GB of
                              memory and 8 cores, installs Julia 1.10,
                              runs the authors' complete verification
                              (zonal matrices and the polynomial
                              identities included, about three days),
                              then llm24_certificate_check.py on the
                              same data, and keeps both logs. Needs
                              LasserreSphericalCodes.zip in the current
                              directory.
                                bash run_llm24_full_verification.sh

  symmetric_search.py
                              Supports Section 7.10 ("Configurations with a
                              symmetry"). For every rotation type of order
                              n <= 12 (angles 2 pi a/n, 2 pi b/n in two
                              orthogonal planes) and the two improper
                              involutions, lists every orbit structure
                              with 23 points in all (73 structures) and
                              runs the slack continuation inside the
                              class, the orbit representatives being the
                              parameters, from a number of random starts;
                              reports per structure the feasible
                              endpoints at slack 0, the best inradius,
                              the best at slack 0.01 and whether the best
                              endpoint is a deletion. Exploration; about
                              an hour for 16 starts per structure.
                                python symmetric_search.py [starts] [seed]

  runs/                       The recorded output of the runs quoted in
                              the paper: cell600_exact.log and
                              cell600_enum.log (the enumeration, Python
                              and C), three passes of
                              slack_continuation.py (seeds 1, 2, 3 with
                              40, 80 and 60 fresh starts per level; the
                              table in Section 7.10 takes the largest
                              inradius at each level over the three), and
                              inradius_search_300_seed7.log (300 direct
                              maximisations at slack 0: 94 feasible
                              endpoints, all deletions), and
                              symmetric_search_16_seed1.log,
                              three_point_reduction.log, the two
                              truncated_volume runs (truncated_volume_r4.log
                              and truncated_volume_r3.log), the relaxation
                              runs three_point_sdp_d6.log and
                              three_point_sdp_d8.log, the certificate run
                              three_point_sdp_certificate_d8.log, its
                              verification certificate_check_d8.log, the
                              rejected altered certificate
                              certificate_check_altered.log (f_0 raised by
                              1e-5), the Lean check
                              D4Certificate_lean.log, the build log of
                              the Lean branch and bound
                              lean_certificate_build.log and its counter
                              run lean_certificate_stat.log (421881
                              boxes), the root-lattice census
                              root_lattices_rank4.log, the exact
                              48-point enumeration octahedral48_exact.log,
                              rigidity24.log, and the re-verification of
                              the certificate of de Laat, Leijenhorst and
                              de Muinck Keizer, llm24_certificate_check.log,
                              with its Lean file's run
                              D4InnerProducts_lean.log. The directory
                              continuation_out/ holds the best
                              configuration of every level of the third
                              pass as a 23 x 4 array (.npy). The runs used
                              numpy 2.4, scipy 1.17 (HiGHS), cvxpy 1.9,
                              Clarabel 0.11, mpmath 1.3, sympy 1.14,
                              python-flint 0.9 and Python 3.11; the endpoints of the searches
                              are local optima and a different platform
                              may reach different ones, whereas the
                              certificate is a fixed file and its
                              verification is deterministic.

  slack_continuation.py
                              Supports Section 7.10 and its table. For
                              delta on a schedule from 0.05 down to 0,
                              estimates h(delta), the largest inradius over
                              23-point configurations with inner products
                              at most 1/2 + delta, by the optimiser of
                              inradius_search.py started from the best
                              configurations of the previous level
                              (restored to the tighter constraint) and from
                              fresh random starts. Prints one line per
                              level: h(delta), the covering radius of the
                              best configuration, the number of feasible
                              endpoints and how many are deletions, and
                              the largest inner products of the best.
                              About a quarter of an hour per pass.
                                python slack_continuation.py [fresh] [seed]

  root_meet.py
                              Supports Propositions 7.55, 7.56 and
                              Remark 7.58, in five parts:

                                (i)   no pair of removed roots destroys a
                                      whole couple of complementary
                                      supports;
                                (ii)  512 of the 2024 triples do, in eight
                                      support patterns, four stars and four
                                      triangles;
                                (iii) the symmetry group of the root
                                      system, generated here and of order
                                      1152, is transitive on the 96 triples
                                      that are pairwise at 60 degrees;
                                (iv)  for the standard star, every
                                      direction of the doorway has first
                                      coordinate at least 1/sqrt2 and the
                                      other three nonnegative, and no two
                                      of them are 60 degrees apart unless
                                      both are removed roots;
                                (v)   one rung lower, at four removed
                                      roots, only four of the 31 orbits
                                      open a doorway at all, and in each
                                      the least achievable largest pairwise
                                      inner product among three directions
                                      kept an angle away from every root
                                      rises above 1/2 as soon as that angle
                                      is positive.

                              Parts (i) to (iii) are exact; parts (iv) and
                              (v) are sampling and minimisation, and the
                              output says so.


  multi_cap_reformulation.py
                              Supports Proposition 7.5, Open Problem 7.6
                              and Remark 7.3.
                              Checks, in order:

                                (A) Q_D is bounded for every packing-valid
                                    D of size 1 to 5, exhaustively: 24,
                                    276, 2024, 10626 and 42504 sets, no
                                    failures;
                                (B) at size 6, exactly 24 of the 134596
                                    packing-valid sets fail, and the
                                    complement of the first of them lies
                                    in a closed half-space;
                                (C) vol(V) = vol(Q_D) - vol(union of the
                                    caps that the tilted half-spaces cut
                                    from Q_D), against direct polytope
                                    volumes;
                                (D) vol(Q_D) - 8 equals |D|/3 exactly when
                                    D has no adjacent pair, and exceeds it
                                    otherwise;
                                (E) the term-by-term reduction fails: in a
                                    random search vol(E_j) comes out both
                                    above and below its single-deviation
                                    value, the shortfall exceeding 0.1
                                    against defects of order 0.3, so the
                                    single-deviation theorem gives no
                                    bound on the individual terms.

                              (A), (B) and (D) are exhaustive over the
                              stated ranges; (C) and (E) are random
                              searches and say so in the output. Eleven
                              checks, about four and a half minutes.

  polar_surface_reformulation.py
                              Supports Section 7.3: Proposition 7.8,
                              Corollary 7.9, Lemma 7.10, Lemma 7.11 and
                              Remarks 7.12 and 7.13. Eighteen checks,
                              grouped here as follows:

                                (a) the cell is the polar dual of the
                                    convex hull of the contact
                                    directions: at the root configuration
                                    the two volumes are 2 and 8;
                                (b) the cell depends on the configuration
                                    only through that hull, tested by
                                    adjoining hull points to random
                                    configurations;
                                (c) vol(V) = (1/4) * total facet 3-volume,
                                    exactly 32/4 = 8 at the roots, with
                                    all 24 facets of 3-volume 4/3;
                                (d) the same identity on configurations
                                    obtained by dropping four roots,
                                    perturbing and separating again;
                                (e) rho(g) = sqrt((1-g)/(1+g)) decreasing,
                                    rho(1/2) = 1/sqrt(3);
                                (f) every facet contains the 3-ball of
                                    radius 1/sqrt(3), by direct inradius
                                    computation;
                                (g) the facet-local bound at 24 contacts,
                                    8*pi/(3*sqrt(3)) = 4.8368, and the
                                    inscribed-ball bound pi^2/2 = 4.9348,
                                    both short of 8;
                                (h) the tight structure at a root facet,
                                    in exact rational arithmetic: eight
                                    tight neighbours, projected Gram
                                    values in {-1, -1/3, 1/3}, cutting out
                                    the regular octahedron of volume 4/3;
                                (i) enlarging the configuration never
                                    increases the cell volume;
                               (i') the square-antiprism configuration of
                                    Remark 7.12: its nine directions in
                                    R^4 are packing-valid, its facet has
                                    3-volume 16 sqrt(2) - 64/3 =
                                    1.2940836646 against the octahedron's
                                    4/3, and the ceiling for any
                                    facet-local bound is therefore
                                    96 sqrt(2) - 128 = 7.7645, below both
                                    8 and the covering bound;
                                (j) the polar volume integral reproduces
                                    8, the Jensen bound returns 7.7351,
                                    and the Mahler-type product is 16
                                    against a conjectural 32/3.

                              (h) is exact; (j) uses Monte Carlo integration
                              over S^3 at four million samples and says so
                              in its output; the rest is double-precision
                              polytope arithmetic on exactly specified
                              configurations. The pairwise-repulsion step
                              that generates test configurations is a
                              sampling device and enters no argument.
                              Runtime under a minute.

  covering_bound.py           Supports Section 7.4: Lemma 7.15,
                              Theorem 7.16, Corollary 7.17,
                              Proposition 7.18 and Table I, which
                              together are Theorem 1.3 of the paper (the
                              local bound at any centre with at most 22
                              contacts). Nineteen checks, grouped here
                              as follows:

                                (a) the cap-area formula
                                    C(r) = pi(2r - sin 2r), against its
                                    value at r = pi and against direct
                                    sampling at four radii;
                                (b) the radial identity
                                    vol = (1/4) int sec^4(delta) at the
                                    root configuration, returning 8 and
                                    recovering the covering radius 45 deg;
                                (c) the layer-cake rewriting used in the
                                    proof, against the direct form;
                                (d) the same identity on random
                                    configurations;
                                (e) the closed form (pi m / 3) tan^3 r_m
                                    against numerical quadrature of the
                                    same estimate, agreeing to 1e-14;
                                (f) the whole of Table I for
                                    m = 5 .. 24, and the monotonicity in m
                                    that the proof of Corollary 7.17 uses;
                                (g) that 22 is exactly the largest m at
                                    which the bound reaches 8;
                                (h) no violation of the bound at the root
                                    configuration, at 100 of its subsets
                                    of sizes 20 to 23, or at 50 random
                                    packing-valid configurations;
                                (i) the constants quoted in the text:
                                    8.046376 at m = 22, 7.798989 at
                                    m = 24, 7.916728 at m = 23, and the
                                    implied density 0.632749;
                                (j) Proposition 7.18: the closed form for
                                    phi'(s), the convexity of phi, and
                                    the fact that equal Voronoi cell
                                    areas reproduce the global bound and
                                    minimise the per-cell sum (against
                                    1500 random area splittings);
                                (k) Proposition 7.21: every spherical
                                    Voronoi cell has circumradius at
                                    least arccos sqrt(5/8) = 37.7612
                                    degrees, checked against the
                                    configurations directly, together
                                    with the fact that it moves the
                                    volume bound only in the fifth
                                    decimal;
                                (l) the split of the shortfall at m = 24
                                    into about 0.057 of overlap and about
                                    0.140 of truncated tail, and the
                                    total overlap of about 2.0 from the
                                    96 pairs of the root system at 60
                                    degrees.

                              (b), (c) and (d) are Monte Carlo and say so
                              in the output; the table, the closed form
                              and the polytope volumes are not.
                              Runtime under a minute.

  saturation_search.py        Supports Section 7.7: Proposition 7.46 and
                              Remark 7.47. Takes the covering radius as
                              the objective from the start, which is the
                              quantity the open case is about. It uses
                              that g(W) is the inradius of conv(W) about
                              the origin, so it comes off a convex hull
                              exactly. Four checks:

                                (a) the 24 roots give g = 1/sqrt(2);
                                (b) all 24 deletions give g = 1/2 exactly,
                                    with zero spread;
                                (c) no 23-point code is a single cyclic
                                    orbit: 528 rotation types, each
                                    decided by a one-variable linear
                                    programme, none feasible, smallest
                                    violation 0.0740002839 at 60 digits,
                                    attained at (a,b) = (8,17);
                                (d) direct minimisation of max_I |z_I|
                                    from perturbed deletions, random
                                    starts and partial root systems.

                              Check (c) is a proof; check (d) is
                              exploration and the output says so. Over 202
                              starts it reached 12 contact configurations,
                              every one of them a deletion. Runtime about
                              three minutes.

  covering_multiplicity.py    Supports Section 7.7: Proposition 7.43
                              and Proposition 7.44. Two parts:

                                (a) the Cauchy-Schwarz bound on the total
                                    overlap of the 60-degree caps of a
                                    saturated configuration, 155.172054
                                    at m = 23, against 162.566121 at the
                                    deletion configuration, a slack of
                                    4.5 per cent;
                                (b) a linear programme over all measures
                                    on [60, 180] degrees of mass 253
                                    subject to Bonferroni at 24 radii and
                                    Gegenbauer positivity at 12 degrees,
                                    whose minimum of the pair weight is
                                    0.041573648 against the 0.083272432
                                    that Proposition 7.41 needs.

                              Part (b) is the sharp statement about the
                              route: nothing reading only the pair angles
                              gets past half way. Runtime about four
                              minutes.

  pair_budget.py              Supports Section 7.7: Proposition 7.41
                              and Remark 7.42. Proposition 7.23 evaluates
                              the pairwise estimate at the deletion
                              configuration and gets 7.997885, three
                              pairs short of 8. The deletion
                              configuration is extendable, so it is not
                              one of the configurations Problem 7.33 is
                              about, and this script computes what the
                              same estimate would need from one that is.
                              Five checks:

                                (a) the per-pair weight w(gamma) and its
                                    collapse away from 60 degrees;
                                (b) the 7.997885 of the deletion
                                    configuration recovered from its 88
                                    tight pairs;
                                (c) that 91 pairs at 60 degrees carry the
                                    estimate to 8.000652 while 90 give
                                    7.999729;
                                (d) the elementary degree bound of 10 on
                                    S^2, hence at most 115 tight pairs on
                                    23 directions;
                                (e) that 91 lies strictly between 88 and
                                    115.

                              The quadratures are adaptive; nothing here
                              is Monte Carlo. Runtime about two minutes.

  spherical_code_23.py        Supports Section 7.9: Proposition 7.59,
                              Theorem 7.60 and Remark 7.62. Contact
                              configurations of m directions are
                              spherical codes of m
                              points on S^3 of minimal angle at least 60
                              degrees, and the script measures how much
                              room such a code has. Riesz continuation,
                              minimising sum |w_i - w_j|^{-s} over the
                              sphere with s = 4, 16, 64, 256, 1024,
                              followed by a direct softmax reduction of
                              the largest inner product. Four sizes:

                                m = 24, where the answer is 1/2;
                                m = 22, where a code with room to spare
                                    exists;
                                m = 25, where no contact configuration
                                    exists;
                                m = 23, from 1500 random starts.

                              At m = 23 it reports the smallest largest
                              inner product reached, the number of
                              starts that fell below 1/2, and, for every
                              outcome within 5e-3 of feasibility, the
                              multiset of Gram entries, which is an
                              O(4) invariant and therefore separates
                              configurations not related by an isometry.
                              Floating point throughout. This is
                              exploration and the output says so; it is
                              not a proof that no saturated 23-point
                              configuration exists. Runtime about forty
                              minutes.

  extendability.py            Supports Section 7.5: Theorem 7.25,
                              Proposition 7.38 and Corollary 7.39. Six
                              checks:

                                (a) the root configuration has
                                    circumradius sqrt(2), g = 1/sqrt(2),
                                    cell volume 8, and is saturated, as
                                    a 24-point kissing configuration must
                                    be;
                                (b) the root system minus one root has
                                    circumradius exactly 2, g exactly
                                    1/2 and cell volume 25/3, and the
                                    deleted root is what extends it, so
                                    it sits on the boundary of the
                                    criterion;
                                (c) the three forms of the criterion
                                    (admissible extra direction,
                                    circumradius at least 2, covering
                                    radius at least 60 degrees) agree on
                                    random configurations;
                                (d) the root configuration has Gram
                                    values in {-1,-1/2,0,1/2}, which is
                                    what Lemma 5.1 of de Laat,
                                    Leijenhorst and de Muinck Keizer
                                    asserts of every 24-point
                                    configuration;
                                (e) the numerical chain of Corollary
                                    7.26;
                                (f) a search for a saturated 23-point
                                    configuration, over the 24 one-root
                                    deletions and random starts, which
                                    found none: the largest g attained
                                    was 1/2 itself, and the random
                                    starts reached no 23-point contact
                                    configuration at all. This is
                                    exploration, labelled as such in the
                                    output, and is not evidence that none
                                    exists.

                              Runtime about 3 minutes.

  second_order_estimate.py    Supports Lemma 7.22 and Proposition 7.23:
                              the covering estimate with every pairwise
                              overlap put back. Four checks, the
                              quadratures in 30-digit arithmetic:

                                (a) three contact directions have
                                    circumradius at least
                                    arccos sqrt(2/3) = 35.2644 degrees,
                                    by the Rayleigh-quotient argument and
                                    against 40000 random admissible
                                    triples;
                                (b) r_23 = 34.6106 and r_24 = 34.0987 both
                                    lie below that, so below r_m the caps
                                    meet only in pairs and
                                    inclusion-exclusion is exact;
                                (c) the closed form for the lens measure
                                    against direct sampling on S^3;
                                (d) the resulting bound: 7.997885 at
                                    m = 23 and 7.858738 at m = 24, both
                                    below 8. At m = 23 only 0.002115 of
                                    the target is left unaccounted for,
                                    and no argument built from cap
                                    measures and pairwise intersections
                                    can recover it.

                              Runtime a few seconds.

  local_cell_obstruction.py   Supports Remark 7.19: the shape deficit
                              that Proposition 7.18 leaves cannot be
                              collected one cell at a time. Five checks:

                                (a) at the root configuration the mean of
                                    sec^4 over a spherical Voronoi cell
                                    is 16/pi^2 = 1.6211389, so the
                                    per-cell inequality is sharp;
                                (b) nine directions of R^3 with pairwise
                                    inner products at most 1/3 exist, so
                                    nine contacts at 60 degrees are
                                    admissible; the cell of such a
                                    direction has mean 1.5732, and 1.5954
                                    when the nine are moved out to 61
                                    degrees, where every inner product is
                                    strictly below 1/2;
                                (c) on a 22-point configuration there are
                                    directions whose nearest contact lies
                                    outside the Delaunay cell containing
                                    them, missing every vertex of that
                                    cell by more than 0.1 in cosine: the
                                    Delaunay decomposition does not
                                    localise the integrand at all;
                                (d) the closed form
                                    4 E[(1+M)^-4] / E[(1+q)^-2], with
                                    E[(1+M)^-4] = 1/5 exactly, gives
                                    1.543643 at the regular simplex of
                                    edge 60 degrees, 4.78 per cent short
                                    of the target, against direct
                                    sampling at edge 60 and 62 degrees;
                                (e) the Delaunay cells of D4 are 24
                                    congruent spherical octahedra of
                                    circumradius 45 degrees, and no four
                                    roots are pairwise at 60 degrees, so
                                    the regular simplex occurs nowhere in
                                    D4 though it does occur as a Delaunay
                                    cell of a contact configuration.

                              The integrals over S^3 are Monte Carlo and
                              the output says so; the margins are
                              percentages, not last digits.
                              Runtime about 30 seconds.

  rigidity23.py               Supports Proposition 7.48, the
                              infinitesimal rigidity of the deletion
                              configuration. Five checks, all in exact
                              integer arithmetic on the unnormalised
                              roots, with no floating-point step:

                                (a) deleting one root leaves 23
                                    directions and 88 tight pairs, the
                                    four values of <a_i,a_0> occurring
                                    1, 8, 6 and 8 times;
                               (a') the six pair types listed in the
                                    proof are the only ones that occur
                                    among the 88;
                                (b) the stress weights y_ij, which are
                                    1, 2 or 3 according to the pair
                                    (<a_i,a_0>, <a_j,a_0>), are strictly
                                    positive on every tight pair;
                                (c) they satisfy the equilibrium relation
                                    sum_j y_ij a_j + mu_i a_i = 0 for
                                    every i, the 23 x 4 residual matrix
                                    being exactly zero;
                                (d) the space of motions holding all 88
                                    pairs at equality has dimension
                                    exactly 6, and the infinitesimal
                                    rotations already span it.

                              Together (b), (c) and (d) say that every
                              first-order motion of the deletion
                              configuration through contact
                              configurations is a rotation.
                              Runtime well under a second.

  rigidity24.py               Supports the corollary after Proposition
                              7.42, the same statement for the root
                              system itself with nothing deleted. Four
                              checks in exact integer arithmetic: 24
                              directions and 96 tight pairs, every
                              direction in 8 of them; the constant
                              stress (weight 1 on every tight pair, -4
                              on the diagonal) is positive on the pairs;
                              it is in equilibrium, the eight
                              neighbours of every root summing to four
                              times the root; and the space of motions
                              holding all 96 pairs at equality has
                              dimension exactly 6, spanned by the
                              rotations. Runtime a few seconds.

## Scripts that run under a time bound

A few of the symbolic derivations do not terminate, and that
non-termination is itself something the paper reports rather than a
defect to be worked around.  Those scripts run under an explicit,
overridable wall-clock bound, report what they established, and exit
cleanly rather than running without end:

  arc2_w1v2/bandBcurve_certificate.py   BANDB_BUDGET_SECONDS (240) and
                                        BANDB_PER_SIMPLEX_SECONDS (90).
                                        The moving-volume step does not
                                        complete; Section 12.2 says so,
                                        and Proposition 8.5 makes the
                                        direct certificate unnecessary.

  misc/sumtest.py                       SUMTEST_STAGE_SECONDS (120).
                                        Compares four ways of collapsing
                                        the same sum; the two that cancel
                                        the whole sum at once do not
                                        settle, which is the point of the
                                        experiment.

  swap_configs/swap_exact_volume.py     SWAPVOL_STAGE_SECONDS (300).
  swap_configs/swap_exact_volume3.py    Same. The raw signed sum is what
                                        the later stages consume, so a
                                        simplify() that does not settle
                                        costs nothing.

  arc1_v1w1/bisect_crossover.py         BISECT_ITERATIONS (40). Not a
                                        non-termination: forty exact
                                        bisections simply take about
                                        three minutes.

Everything else runs to completion, and most of it finishes in seconds.
The slowest are the exact region derivations in arc1_v1w1/
(regionA_derive.py, region6_derive.py, regionC_derive.py,
regionD_derive.py and regionE_derive.py), whose symbolic cancellation
step is dominated by a handful of hard terms and can run for an hour or
more, and the directional sweeps in hessian_multidir/, at roughly seven
to eight minutes each. The cached .pkl files in data/ are there so
that the consuming scripts need not wait for them.

## Requirements

  Python >= 3.9
  sympy >= 1.10
  mpmath >= 1.2
  numpy, scipy (vertex_degeneracy_check.py,
    cap_inequality_certificate.py, and the double-precision scans)
  cvxpy >= 1.9, with a working SDP solver (e.g. CLARABEL, bundled with
    cvxpy's default install) -- gram_sos_lib.py and quartic_fit_and_check.py
    only
  python-flint >= 0.9 -- llm24_certificate_check.py, shell_reduction.py,
    closure_lemmas.py and labelled_certificate_check.py only

## Usage

Each script can be run directly, from any working directory:
  python cap_certificate/cap_inequality_certificate.py
  python multi_cap/multi_cap_reformulation.py
  python multi_cap/polar_surface_reformulation.py
  python multi_cap/covering_bound.py
  python multi_cap/extendability.py
  python multi_cap/spherical_code_23.py
  python multi_cap/pair_budget.py
  python multi_cap/covering_multiplicity.py
  python multi_cap/saturation_search.py
  python multi_cap/rigidity23.py
  python multi_cap/rigidity24.py
  python multi_cap/llm24_certificate_check.py /path/to/LasserreSphericalCodes/proofs/4_24
  python multi_cap/local_cell_obstruction.py
  python multi_cap/second_order_estimate.py
  python multi_cap/root_deletions_exact.py
  python multi_cap/root_meet.py
  python multi_cap/cell600_exact.py
  python multi_cap/inradius_search.py 300 1
  python multi_cap/slack_continuation.py 40 1
  python multi_cap/shell_reduction.py
  python multi_cap/closure_lemmas.py
  python multi_cap/near_contact_probe.py 200 5
  python multi_cap/shell_neighbour_search.py 240 1 3
  python multi_cap/shell_neighbour_search.py 60 7 3 0.05
  python verification/refcell_verify.py
  python arc1_v1w1/regionA_boxcover.py
  python hessian_multidir/multidir_nullspace_broad_sample.py
  etc.

The Lean developments are checked with
  lean lean/D4Stress.lean
  lean lean/D4Meet.lean
  lean lean/D4Certificate.lean
  lean lean/D4InnerProducts.lean
  lean lean/D4RootLattices.lean
  lean lean/D4Closure.lean
  (cd lean/cell600 && lake build)
  (cd lean/certificate && lake build)
and D4Certificate.lean, D4InnerProducts.lean and
lean/certificate/D4CertData.lean are regenerated with
  (cd lean && python3 gen_certificate_lean.py)
  (cd lean && python3 gen_innerproducts_lean.py)
  (cd lean/certificate && python3 gen_data.py)
and the C enumeration with
  (cd multi_cap && cc -O2 -o cell600_enum cell600_enum.c && ./cell600_enum)
after cell600_exact.py has written cell600_graph.txt.

Scripts that consume a cached intermediate result look for it in data/
relative to the package root, so the working directory does not matter.
Where a script is one stage of a pipeline, the matching *_derive.py
produces what *_map.py and *_boxcover.py consume; running the derive
step refreshes the cached file in place.

All computations use exact rational arithmetic (fractions.Fraction or
sympy.Rational) wherever a proof is claimed; mpmath high-precision
floating-point is used only for numerical estimates, which are labelled
as such in the paper and are never treated as proofs. In particular, no
script in hessian_multidir/ constitutes a proof of the
multi-direction positivity conjecture: every result there is either an exact
demonstration of a specific finite fact (e.g. the shared-vertex count in
vertex_degeneracy_check.py) or numerical evidence explicitly reported as
such (the near-null-subspace scaling arguments and the directional
sampling), exactly as the paper itself describes them.
