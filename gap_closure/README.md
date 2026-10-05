# gap_closure: floating-point explorations of the open cases

This directory holds the computations run against the two statements that the
paper leaves open (Section "The statements (G) and (C)"): statement (C) for
25 to 30 centres within sqrt6, and statement (G). **None of them yields a
certificate.** They are floating-point semidefinite programmes on sampled
constraints, kept so that the negative results quoted in the paper and in the
pull request can be reproduced. No result here is used in any proof.

The level that (C) asks for is `9 pi^2/8 - 8 = 3.10330`: a programme proves a
case only if its bound, after the correction for the constraints that the
samples miss, is below that level.

## C30: the residual case at thirty centres

The case left by `thm:count30`: 24 centres within 2.25 (22 within 2.05, 23
within 2.15) and 6 beyond 2.4.

| script | kernel | result (log) |
|---|---|---|
| `combo30.py` | distance-labelled two-point kernel plus a three-point kernel typed by distance class (A = [2, 2.05], B = (2.05, 2.25], F = [2.4, sqrt6)) | degree 0: 3.1212 (`combo30_d0.log`); degree 2: 3.1211 (`combo30_d2.log`); degree 6, three rounds: raw 3.0991, corrected 3.1160 (`combo30_d6b.log`) |
| `combo30p.py` | the same, with the sample set pruned so that more rounds fit in memory | degree 6: raw 3.0995, corrected 3.1103 (`combo30p_d6.log`); degree 8, two rounds: raw 2.8946 then 2.9338, corrected 4.4901 then 3.8157 (`combo30p_d8.log`, `combo30p_d8b.log`); the third round needs more than 9 GB of memory |
| `typed_delsarte.py`, `code_feasible.py`, `code_feasible2.py`, `code_feas_gen.py` | typed linear programming bound and local search for two-shell codes with points in holes | exploration only (`slackhole_d6.log`) |
| `case_value.py` | value of the two-point case programme of the paper for given count constraints | used to choose splits |

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
