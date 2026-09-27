# The second order at the root system

The scripts of Section 2.9 of the paper (Lemma 2.19 to Remark 2.25, Figure 11).
Run them from this directory.

| file | what it does |
| --- | --- |
| `rational_model.py` | the form H of the second variation, the 120 rows of the first-order packing cone and c = 2 1_eta, in exact rationals from the integral root system |
| `verify_cone_certificate.py` | the proof of m = -1: H + c c^T = P + B^T N B with N >= 0 and P positive semidefinite, by an exact LDL^T (three seconds, no floating point) |
| `exact_certificate.pkl`, `.txt` | the certificate: N constant on 43 orbits of the signed-permutation group, 35 values positive |
| `cone_min.py` | m in floating point: the best feasible point over 62 starts, and the Shor relaxation with the products of the cone rows, both -1 |
| `certify_m.py`, `exact_cert.py`, `exact_cert2.py`, `exact_cert3.py` | how the certificate was found and made exact |
| `hessian.py` | the form in floating point (writes `H.npy`), its spectrum, and a finite-difference check |
| `exact_push.py` | the one-centre formula, Lemma 2.23 and Proposition 2.24 on examples |
| `untilt.py` | tilting is not monotone: 103 of 237 tilted packings below their untilted versions, by up to 3.6e-5 |
| `remainder_test.py` | 300 packings with tilts on the edges of the cone: vol - 8 - (2/3) S + (1/2) S^2 >= 0 on all; writes `runs/remainder_rows.dat` |
| `independent_check.py` | recomputes every volume by intersecting halfspaces: H against second differences (6e-8), the cone rows against exact distances (2e-9), the spectrum, the one-centre formula (1e-14), Lemma 2.23 on 400 push patterns |

Logs are in `runs/`.  What is proved here is the second order (the certificate)
and the push-out statements (Lemma 2.23, Proposition 2.24, by the proofs in the
paper); the remainder and the non-monotonicity are floating-point experiments.

Two statements of the draft these files came with are corrected in the paper:
the eigenvalue -6.194 of H has multiplicity 8, not 6, and on the push-outs H has
matrix Adj - 4I, not 2 Adj - 4I (its form is 2 sum over edges eta_i eta_j - 4 sum eta_i^2).
