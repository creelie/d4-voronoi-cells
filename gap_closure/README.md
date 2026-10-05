# gap_closure: floating-point explorations of the open cases

This directory holds the computations run against the two statements that the
paper leaves open (Section "The statements (G) and (C)"): statement (C) for
25 to 30 centres within sqrt6, and statement (G). **None of them yields a
certificate.** They are floating-point semidefinite programmes on sampled
constraints, kept so that the negative results quoted in the paper and in the
pull request can be reproduced. No result here is used in any proof, except
the exact checks logged in `density/check_M.log` (see the last section).

The level that (C) asks for is `9 pi^2/8 - 8 = 3.10330`: a programme proves a
case only if its bound, after the correction for the constraints that the
samples miss, is below that level.

## C30: the residual case at thirty centres

The case left by `thm:count30`: 24 centres within 2.25 (22 within 2.05, 23
within 2.15) and 6 beyond 2.4.

| script | kernel | result (log) |
|---|---|---|
| `combo30.py` | distance-labelled two-point kernel plus a three-point kernel typed by distance class (A = [2, 2.05], B = (2.05, 2.25], F = [2.4, sqrt6)) | degree 0: 3.1212 (`combo30_d0.log`); degree 2: 3.1211 (`combo30_d2.log`); degree 6, three rounds: raw 3.0991, corrected 3.1160 (`combo30_d6b.log`) |
| `combo30p.py` | the same, with the sample set pruned so that more rounds fit in memory | degree 6: raw 3.0995, corrected 3.1103 (`combo30p_d6.log`); degree 8, two rounds: raw 2.8946 then 2.9338, corrected 4.4901 then 3.8157 (`combo30p_d8.log`, `combo30p_d8b.log`); resumed with 3000 pair and 1200 triple samples kept per kind (`combo30q_d8.log`): raw 2.9349 and 2.9388, corrected 3.3005 and 3.1528 in its first two rounds, the excess now coming from the triples |
| `typed_delsarte.py`, `code_feasible.py`, `code_feasible2.py`, `code_feas_gen.py` | typed linear programming bound and local search for two-shell codes with points in holes | exploration only (`slackhole_d6.log`) |
| `case_value.py` | value of the two-point case programme of the paper for given count constraints | used to choose splits |
| `multi_cap/combo30_check.py` | the exact check that a certificate of `combo30p.py` has to pass: exact positivity, the pair inequalities by Bernstein branch and bound with Pi in Arb, the bins, the typed triple inequalities by Taylor branch and bound, and the bound over the five count vectors | its exact polynomials agree with the solver's to 2e-16 on test points; no certificate has been submitted to it yet |

The raw value is the optimum on the samples; the corrected value adds, for each
kind of pair and triple, the largest violation found between the samples times
the number of such pairs or triples. At degree 6 the raw value rises with each
round of added samples and the corrected value stays above 3.10330.

## CM: twenty-five to twenty-nine centres

`combo_gen2.py` generalises `combo30.py` to an arbitrary case file
(`case29_all.json`: 29 centres, types A, B, F, seven bins and the count
constraints proved by `prop:C-radial`), and bounds the maximum over all count
vectors by linear-programming duality instead of enumerating them
(`combo_gen.py` enumerates them and runs out of memory). With the two-point
kernel alone (degree 0 in the three-point part) the value is 3.1967
(`c29all2_d0.log`), far above the level.

## musin: where the two-point kernel fails

`dual_where.py` reads off the optimal dual of the two-point programme of
`thm:count31` (`multi_cap/radial_count_sdp.py`): the fictitious configuration
that the kernel cannot exclude. At 25 centres it puts all 25 at distance 2,
with pair inner products near 0.5, -0.15, -0.35 and -0.85; at 30 it puts 23.6
at distance 2 and 6.4 near sqrt6 (`dual_where.log`). The kernel is blind in the
way Delsarte's bound for the kissing number of R^4 (25.56) is blind.

`musin_scan.py` adds Musin's relaxation (Ann. of Math. 168, 2008): the pair
inequality is dropped for nearly antipodal pairs, and their excess is bounded
centre by centre by the number of centres that fit in a cap about the
antipode. At 30 centres the bound stays at 3.1441 for every threshold from
0.95 to 0.6 (`musin_scan_30.log`): the dual has no surplus near the antipode.

## localisation: the room for a budgeted localisation

Both open statements reduce to one statement (L): 24 centres with
T(Y) <= 8 + eta have directions within an explicit root-sum-square distance
rho(eta) of a root system. The residual case at thirty needs rho < 0.206 at
eta = 0.00368 (`rem:thirty-left`, `cor:no-room`).

`budget_far.py` maximises the distance to the nearest root system under
T(Y) <= 8 + eta and the packing conditions, from root systems pushed out at
random (`budget_far.log`). The largest distance found is 0.0971 at
eta = 0.00368 (19 feasible ends of 24 starts) and 0.1291 at eta = 0.0334.
These are local searches near the root system, not bounds.

`room_table.py` gives, for one further centre at distance r, the budget
eta = S(r) and the distance that `cor:no-room` (ii) needs, with a(sqrt6, rho)
replaced by a(r, rho24) (`room_table.log`). `budget_far.py eta starts seed rmax`
holds the 24 centres within rmax. At the worst budget, eta = S(2.0161) = 0.127,
the searches reach 0.339, 0.378 and 0.378 for rmax = 2.25, 2.35 and 2.444
(`budget_far_r.log`), against 0.348, 0.288 and 0.231 needed: for a close
further centre the root-sum-square route does not suffice, and the hole form of
(L) is needed.

`hole_room.py` holds the 25th centre at a given distance r and minimises T
over 25 centres (`hole_room.log`): the least T found is 8.290, 8.287, 8.277,
8.259, 8.276, 8.275 and 8.264 at r = 2.0161, 2.05, 2.1, 2.2, 2.3, 2.4 and
2.449, so the hole form of (L) has a margin of about 0.26 at every distance.

`link_census.py` lists, point by point, the near neighbours (inner product at
least 0.4), the distance of the link from a cube and the Voronoi cell on S^3
(`link_census.log`). At the root system and at the ends of `budget_far.py`
every point has 8 near neighbours with a cube-like link; the second code has
none (6, 7 and 9 near neighbours).

`s3_cell.py` shows that a cell-by-cell volume bound on S^3 cannot prove (L):
with 8 neighbours at 60 degrees whose link is the square antiprism, the
Voronoi cell of a direction has volume at most about 0.817, below the
pi^2/12 = 0.8225 of the root system (`s3_cell.log`, Monte Carlo). A proof of
(L) along the lines of dimension three therefore has to move volume between
neighbouring cells.

## density: a bound for every packing

Section "A density bound from the cells alone" of the paper (`prop:levels`,
`thm:density-cells`) bounds the union of caps U(Y) for every count M from 24
to 30, with no case split except at M = 26, where the case N(2.03) = 26 is
empty by `thm:twenty-six`; so every cell has volume at least
9 pi^2/8 - max L_M = 7.7532 and every packing of unit balls in R^4 has density
at most 0.63649. This is weaker than the three-point bound 0.63611 of Cohn,
de Laat and Salmon, and far from pi^2/16 = 0.61685.

| file | what it is |
|---|---|
| `levels.py` | floating point: the least level L that a certificate can reach for M centres, by bisection; used to choose the levels |
| `cert_M.log` | the certificate search, `multi_cap/radial_case_sdp.py M --level L`, writing `multi_cap/radial_certificates/density_M.json`; for M = 26, `cert_26_card.log`, with the spec `density_spec_26_card.json` (split at N(2.03), level 3.3503) |
| `check_M.log` | the exact check, `multi_cap/radial_case_check.py`: positivity by exact LDL^T, K <= Pi by branch and bound with Pi in Arb, the bin bounds, and the bound with no assumption on the count vector. All seven pass: L_M = 3.3352, 3.3274, 3.3501, 3.3090, 3.2587, 3.2080, 3.1538 for M = 24 to 30 (`check_26_level3353.log` is the earlier check of M = 26 with no split, 3.3524) |
| `case26_all.json`, `c26all_d2.log` | floating point: a two-point kernel plus a typed three-point kernel at M = 26, value 3.34979; not enough |

A bound below 0.63611 needs L_26 <= 3.34549. The route tried here splits the
26-centre case by counts at 2.0161 and 2.03 (`split_level.py`, floating point):

| leaf | counts | level | what removes its complement |
|---|---|---|---|
| A | N(2.0161) <= 23, N(2.03) <= 25 | `split26_A.log` | N(2.03) >= 26: at most 25 points of S^3 with inner products at most 0.51468 (`card26_d10_*.log`) |
| B | N(2.0161) = 24, N(2.03) <= 24 | `split26_B.log` | N(2.03) >= 25: the typed exclusion at t2 = 0.5114 below |

With N(2.03) <= 25 left out, leaf B is too high: {N(2.0161) = 24, N(r2) <= 24}
gives 3.34717 at r2 = 2.0248 and 3.34708 at r2 = 2.025. Within leaf B the
optimum sits where a centre lies in (2.03, 2.05]: with N(2.05) <= 24 the level
is 3.33516 (`split26_B1.log`), with N(2.05) >= 25 it is 3.34471
(`split26_B2.log`), so an exclusion at r2 = 2.05 (t2 = 0.5163) would bring the
level to 3.34336, that of leaf A.

A cardinality bound alone does not move the two-point level: adding
N(2.0793) <= 24 at M = 25 (no 25 points of S^3 with inner products at most
0.5374, true numerically, the best 25-point code found having 0.537429,
`multi_cap/runs/code25_best.txt`) leaves the level at 3.33003
(`ladder25_t5374.log`).

The typed exclusion: 24 directions with inner products at most t1 = 0.508
(centres within 2.0161) leave no further direction with inner product at most
t2 with all of them (`multi_cap/typed_cardinality_sdp.py`; a corrected Z below
0 excludes the code in floating point, `multi_cap/typed_cardinality_check.py`
is the exact check). A centre within r2 of the origin has t2 = a(2.0161, r2):
0.5101 for r2 = 2.0248, 0.51135 for 2.03, 0.5163 for 2.05, 0.6141 for sqrt6.

| log | degree, t2 | result |
|---|---|---|
| `typed_d10_5101.log` | 10, 0.5101 | sampled Z +0.0255: no certificate of degree 10 (the sampled programme is a relaxation) |
| `typed_d14.log` | 14, 0.5114 | sampled -1.95, corrected +6.60 after one round (old sampling) |
| `typed2_d10_control.log` | 10, 0.508 | control, 25 points at 0.508, which `thm:kissing-stable` excludes: corrected Z -1.22 at round 2, so the refinement converges |
| `typed2_d12_5101.log` | 12, 0.5101 | corrected Z -1.26 at round 4: excluded in floating point |
| `typed2_d14_5101.log` | 14, 0.5101 | corrected Z +9.58 after two rounds |
| `typed2_d12_5114.log` | 12, 0.5114 | the exclusion leaf B needs: corrected Z 25.9, 16.7, 21.2, 0.506 in rounds 1 to 4 (sampled -0.97 to -0.42); the run was stopped by the memory limit in round 5 |
| `typed2_d12_5163.log` | 12, 0.5163 | sampled Z +0.011: no certificate of degree 12 |
| `typed2_d12_6141.log`, `typed2_d14_6141.log` | 12 and 14, 0.6141 | sampled Z +0.167 and +0.018: no certificate; an exclusion at sqrt6 would settle (C) when 24 centres lie within 2.0161 |
| `typed2_d12_071_control.log`, `typed2_d14_071_control.log` | 12 and 14, 0.71 | control: the deep holes of the 24-cell are at inner product 0.7071, so no certificate exists; sampled Z +0.951 and +0.192 |

Degrees 16 and 18 at t2 = 0.6141 are run on a larger machine
(`bigmachine/run_typed.py`): degree 16 gives sampled Z +0.0015 in each of
three rounds, so no certificate of degree 16 either. The sampled values at
sqrt6 fall about tenfold per degree (+1.87, +0.167, +0.018, +0.0015 at
degrees 10 to 16); those of the 0.71 control, where no certificate can
exist, fall too (+5.56, +0.951, +0.192 at degrees 10 to 14), so the fall
alone does not point towards a certificate.
