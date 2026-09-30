# hessian_multidir/

The joint-Hessian computations of Sections 15 to 19 of the paper, file by file.



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
