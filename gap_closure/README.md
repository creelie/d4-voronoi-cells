# gap_closure: floating-point explorations of the open cases

This directory holds the computations run against the two statements that the
paper leaves open (Section "The statements (G) and (C)"): statement (C) for 25
to 27 centres within sqrt6 and statement (G), and the cases they settled,
the residual case at thirty centres, the whole case at twenty-nine, and
twenty-eight: first with at most fourteen centres within 2.0161 (fifteen
with at most three beyond 2.35, sixteen with none), then all but exactly
twenty-three within 2.2, all of them within 2.0161, with four beyond 2.35,
and last that case, with its close-close-far triples counted, and
twenty-seven with at most ten, then exactly eleven, twelve, thirteen and
fourteen, centres within 2.0161. They are floating-point semidefinite programmes on sampled
constraints. **Three of them yield certificates**, and convex combinations of
several yield thirty-one more: the round-3 certificate of `C30/combo30p.py`, which
`multi_cap/combo30_check.py` proves exactly and which is the fifth case of
`thm:count30`, and the certificates of the runs `c29r4` and `c28lo13b` of
`CM/combo_direct.py`, which `multi_cap/combo_case_check.py` proves exactly and
which are `thm:count29` and `prop:count28-few`; and four combinations of
certificates of `CM/combo_direct.py`, one of fourteen found by `CM/c28_hull.py`
for at most twenty-one within 2.0161, and one of nine for exactly twenty-two,
one of six for exactly twenty-four and one of eight for exactly twenty-three
with at most three beyond 2.35, found by `CM/c28_hull_multi.py`, and four
for single vectors with twenty-three within 2.0161, four beyond 2.35 and the
last centre in (2.05, 2.1], (2.1, 2.15], (2.15, 2.2] or (2.2, 2.35], found by
`CM/c28_hull_case.py` (the last one checked with `REGION_TRIPLES`, which
counts its close-close-far triples), which `multi_cap/combo_case_check.py`
proves exactly and which together are `prop:count28-mixed`, (C) at
twenty-eight; and thirty-three combinations for twenty-seven centres, one of six
certificates for at most ten within 2.0161 and at most two beyond 2.35, one
of seven for at most ten within 2.0161 and exactly three beyond 2.35, one of
six for exactly eleven within 2.0161 and at most two beyond 2.35, one of
seven for exactly eleven and exactly three, one of seven for exactly twelve
and at most two, one of nine for exactly twelve and exactly three, one of
nine for exactly thirteen and at most two and one of ten for exactly
thirteen and exactly three, found by `CM/c28_hull_multi.py` with
`CASEM=27`, and two for exactly fourteen within 2.0161 and exactly three
beyond 2.35, one of nine for the 219 vectors with a centre in (2.1, 2.35]
(`CM/case27_k14f3g2.json`) and one of seven for the 11 others
(`CM/case27_k14f3g1.json`), found by `CM/c28_hull_case.py`, the second with
the middle type narrowed to (2.05, 2.1] (`CM/case27_k14f3g1_narrow.json`) and
checked with thresholds per combination of bins, one of thirteen for the 39
vectors with exactly fourteen within 2.0161, at most two beyond 2.35 and none
in (2.1, 2.35] (`CM/case27_k14f2g1.json`), found by `CM/c28_hull_case.py` on
the narrowed case `CM/case27_k14f2g1_narrow.json` and checked in the same
way, and one of eleven for the 192 vectors with exactly fifteen within
2.0161, exactly three beyond 2.35 and a centre in (2.1, 2.35]
(`CM/case27_k15f3g2.json`), found by `CM/c28_hull_case.py`, which
`multi_cap/combo_case_check.py` proves exactly, and one of twelve for the 127
vectors with exactly fourteen within 2.0161, at most two beyond 2.35, none in
(2.1, 2.2] and one in (2.2, 2.35] (`CM/case27_k14f2g2a_bins.json`), found by
`CM/c27_hull_mt.py` with each certificate typing the bin (2.2, 2.35] as B or
as F, which `multi_cap/combo_mixtype_check.py` proves exactly, and one of
ten for the 1868 vectors with exactly fourteen within 2.0161, at most two
beyond 2.35 and a centre in (2.1, 2.2] (`CM/case27_k14f2g2b.json`), found by
`CM/c27_hull_sub.py` with `SUB=2,3,4,5` (a threshold of its own for each
combination of bins with a middle bin), which `multi_cap/combo_case_check.py`
proves exactly, and one of seven for the 15 vectors with exactly sixteen
within 2.0161, exactly three beyond 2.35, none in (2.1, 2.2] and one in
(2.2, 2.35] (`CM/case27_k16f3g2a.json`), found by `CM/c27_hull_sub.py` with
`SUB=2,5`, which `multi_cap/combo_case_check.py` proves exactly
(`CM/check_c27k16f3g2a.log`: U <= 3.093588 at (16, 4, 3, 0, 0, 1, 3), 120.0
million triple boxes), and one of five for the 9 vectors with exactly
fifteen within 2.0161, exactly three beyond 2.35, none in (2.1, 2.35] and
at most eight in (2.0161, 2.05] (`CM/case27_k15f3g1lo8_bins.json`), found
by `CM/c27_hull_mt.py` under five typings of the bins, 0.93 of the weight
typing the bin (2.05, 2.1] as A, which `multi_cap/combo_mixtype_check.py`
proves exactly (`CM/check_c27k15f3g1lo8mt.log`: U <= 3.098238 at (15, 8, 1,
0, 0, 0, 3), 208.9 million triple boxes), and one of nine for the 8 vectors
with exactly sixteen within 2.0161, exactly three beyond 2.35, none in
(2.1, 2.35] and at most seven in (2.0161, 2.05]
(`CM/case27_k16f3g1lo7_bins.json`), found by `CM/c27_hull_mt.py` under five
typings of the bins, 0.94 of the weight typing the bin (2.05, 2.1] as A,
which `multi_cap/combo_mixtype_check.py` proves exactly
(`CM/check_c27k16f3g1lo7mt.log`: U <= 3.101848 at (16, 7, 1, 0, 0, 0, 3),
292.0 million triple boxes), and one of six for the 150 vectors with
exactly sixteen within 2.0161, exactly three beyond 2.35 and a centre in
(2.1, 2.2] (`CM/case27_k16f3g2b_bins.json`), found by `CM/c27_hull_mt.py`
under three typings of the bins, 0.81 of the weight typing the bin
(2.05, 2.1] as A, which `multi_cap/combo_mixtype_check.py` proves exactly
(`CM/check_c27k16f3g2bmt.log`: U <= 3.101991 at (16, 7, 0, 1, 0, 0, 3),
630.6 million triple boxes), and one of seven for the 36 vectors with
exactly fifteen within 2.0161, at most two beyond 2.35 and none in
(2.1, 2.35] (`CM/case27_k15f2g1_bins.json`), found by `CM/c27_hull_mt.py`
under five typings of the bins, 0.52 of the weight typing the bin
(2.05, 2.1] as A, which `multi_cap/combo_mixtype_check.py` proves exactly
(`CM/check_c27k15f2g1mt.log`: U <= 3.102529 at (15, 0, 11, 0, 0, 0, 1),
256.4 million triple boxes), and two of eight for the 30 vectors each
with exactly fifteen within 2.0161, at most two beyond 2.35, at least one in
(2.05, 2.1] and exactly one in (2.1, 2.35], lying in (2.1, 2.15] or in
(2.15, 2.2] (`CM/case27_k15f2g2bH4_bins.json`,
`CM/case27_k15f2g2bH5_bins.json`), found by `CM/c27_hull_mt.py` under five
typings of the bins, 0.67 of the weight typing the bin (2.05, 2.1] as A,
which `multi_cap/combo_mixtype_check.py` proves exactly
(`CM/check_c27k15f2g2bH4mt.log`: U <= 3.099453 at (15, 8, 1, 1, 0, 0, 2),
511.1 million triple boxes; `CM/check_c27k15f2g2bH5mt.log`: U <= 3.099138
at (15, 8, 1, 0, 1, 0, 2), 500.8 million triple boxes), and one of eight
for the 6 vectors with exactly fifteen within 2.0161, at most two beyond
2.35, none in (2.05, 2.1] and exactly one in (2.1, 2.35], lying in
(2.1, 2.2] (`CM/case27_k15f2g2bH0_bins.json`), found by `CM/c27_hull_mt.py`
under the same typings with the peaks of a dense search added to the
samples, which `multi_cap/combo_mixtype_check.py` proves exactly
(`CM/check_c27k15f2g2bH0mt2.log`: U <= 3.102737 at (15, 9, 0, 0, 1, 0, 2),
436.5 million triple boxes), and one of six for the 23 vectors with
exactly sixteen within 2.0161, at most one beyond 2.35 and none in
(2.1, 2.35] (`CM/case27_k16f2G1a_bins.json`), found by `CM/c27_hull_mt.py`
under five typings, 0.50 of the weight typing the bin (2.05, 2.1] as A,
which `multi_cap/combo_mixtype_check.py` proves exactly
(`CM/check_c27k16f2G1amt.log`: U <= 3.102925 at (16, 0, 10, 0, 0, 0, 1),
254.0 million triple boxes), and one of six for the 8 vectors with exactly
sixteen within 2.0161, exactly two beyond 2.35, none in (2.1, 2.35] and at
least one in each of (2.0161, 2.05] and (2.05, 2.1]
(`CM/case27_k16f2G1bm_bins.json`), found by `CM/c27_hull_mt.py` under five
typings, 0.66 of the weight typing the bin (2.05, 2.1] as A, which
`multi_cap/combo_mixtype_check.py` proves exactly
(`CM/check_c27k16f2G1bmmt.log`: U <= 3.102229 at (16, 8, 1, 0, 0, 0, 2),
337.7 million triple boxes), one of six for the 103 vectors with exactly
fifteen within 2.0161, at most two beyond 2.35, none in (2.1, 2.2] and at
least one in each of (2.05, 2.1] and (2.2, 2.35]
(`CM/case27_k15f2g2aY_bins.json`), found by `CM/c27_hull_mt.py` under five
typings with the dense-search peaks of a first mix fed back, 0.68 of the
weight typing the bin (2.05, 2.1] as A and 0.45 typing (2.2, 2.35] as F,
which `multi_cap/combo_mixtype_check.py` proves exactly
(`CM/check_c27k15f2g2aYmt2.log`: U <= 3.103081 at (15, 8, 1, 0, 0, 1, 2),
653.7 million triple boxes), one of four for the 6 vectors with exactly
fifteen within 2.0161, at most two beyond 2.35, none in (2.05, 2.2] and at
least three in (2.2, 2.35] (`CM/case27_k15f2g2aX3_bins.json`), found by
`CM/c27_hull_mt.py` under five typings, 0.14 of the weight typing the bin
(2.2, 2.35] as F, which `multi_cap/combo_mixtype_check.py` proves exactly
(`CM/check_c27k15f2g2aX3mt.log`: U <= 3.102002 at (15, 7, 0, 0, 0, 3, 2),
383.2 million triple boxes), one of eight for the 36 vectors with
exactly sixteen within 2.0161, at most two beyond 2.35, exactly one in
(2.1, 2.2], none in (2.2, 2.35] and at least four in (2.05, 2.1]
(`CM/case27_k16f2Hn3_bins.json`), found by `CM/c27_hull_mt.py` under five
typings, 0.65 of the weight typing the bin (2.05, 2.1] as A and 0.13 typing
(2.1, 2.2] as F, which `multi_cap/combo_mixtype_check.py` proves exactly
(`CM/check_c27k16f2Hn3mt.log`: U <= 3.096243 at (16, 4, 4, 0, 1, 0, 2),
610.1 million triple boxes), one of eight for the 12 vectors with
exactly sixteen within 2.0161, at most two beyond 2.35, exactly one in
(2.1, 2.2], none in (2.2, 2.35] and two or three in (2.05, 2.1]
(`CM/case27_k16f2Hm2_bins.json`), found by `CM/c27_hull_mt.py` under five
typings for the 18 vectors with one to three in (2.05, 2.1], 0.73 of the
weight typing the bin (2.05, 2.1] as A and 0.14 typing (2.1, 2.2] as F,
which `multi_cap/combo_mixtype_check.py` proves exactly
(`CM/check_c27k16f2Hm2mt.log`: U <= 3.101346 at (16, 6, 2, 1, 0, 0, 2),
707.0 million triple boxes), one of seven for the 67 vectors with
exactly sixteen within 2.0161, at most two beyond 2.35, none in (2.1, 2.2],
at least one in (2.2, 2.35] and at least three in (2.05, 2.1]
(`CM/case27_k16f2G2an3_bins.json`), found by `CM/c27_hull_mt.py` under five
typings for the 79 vectors with at least two in (2.05, 2.1], 0.67 of the
weight typing the bin (2.05, 2.1] as A and 0.34 typing (2.2, 2.35] as F,
which `multi_cap/combo_mixtype_check.py` proves exactly
(`CM/check_c27k16f2G2an3mt.log`: U <= 3.102140 at (16, 5, 3, 0, 0, 1, 2),
616.0 million triple boxes), and four of eight, ten, five and six for the
single vectors (15, 11, 0, 0, 0, 1, 0), (15, 10, 0, 0, 0, 1, 1),
(15, 10, 0, 0, 0, 2, 0) and (15, 9, 0, 0, 0, 2, 1), with exactly fifteen
within 2.0161, at most one beyond 2.35, none in (2.05, 2.2] and one or two
in (2.2, 2.35] (`CM/case27_k15f2g2aX1v*_bins.json`), found by
`CM/c27_hull_mt.py` for one vector at a time under five typings, the last
with the dense-search peaks of a first mix fed back, which
`multi_cap/combo_mixtype_check.py` proves exactly
(`CM/check_c27k15f2g2aX1v1110mt.log`: U <= 3.098938, 102.5 million triple
boxes; `CM/check_c27k15f2g2aX1v1011mt.log`: U <= 3.100562, 225.5 million;
`CM/check_c27k15f2g2aX1v1020mt.log`: U <= 3.096137, 118.0 million;
`CM/check_c27k15f2g2aX1v921mt2.log`: U <= 3.103026, 368.0 million);
together they are
`prop:count27-few`, (C) at 27 with at most
14 within 2.0161, or with exactly 15 or exactly 16 within 2.0161, exactly 3
beyond 2.35 and at least one in (2.05, 2.35], or with exactly 15 within
2.0161, at most 2 beyond 2.35 and either none in (2.2, 2.35] and at most
one in (2.1, 2.2], or none in (2.1, 2.2], at least one in (2.2, 2.35] and,
if none in (2.05, 2.1], at least 3 in (2.2, 2.35] or at most 1 beyond 2.35, or with exactly 16 within 2.0161, at most 2 beyond 2.35 and
either none in (2.1, 2.35] and, if 2 lie beyond 2.35, at least one in each
of (2.0161, 2.05] and (2.05, 2.1], or exactly one in (2.1, 2.2], none in
(2.2, 2.35] and at least 2 in (2.05, 2.1], or none in (2.1, 2.2], at least
one in (2.2, 2.35] and at least 3 in (2.05, 2.1]. The rest are kept so that the
values quoted in the paper and in the pull requests can be reproduced, and no
result of theirs is used in any proof, except the exact checks logged in
`density/check_M.log` (see the last section).

The level that (C) asks for is `9 pi^2/8 - 8 = 3.10330`: a programme proves a
case only if its bound, after the correction for the constraints that the
samples miss, is below that level.

## C30: the residual case at thirty centres (proved)

The fifth case of `thm:count30`: 24 centres within 2.25 (22 within 2.05, 23
within 2.15) and 6 beyond 2.4. The two-point programmes of the tree do not
close it (`rem:thirty-pairs`); the kernel on typed triples does.

| script | kernel | result (log) |
|---|---|---|
| `combo30.py` | distance-labelled two-point kernel plus a three-point kernel typed by distance class (A = [2, 2.05], B = (2.05, 2.25], F = [2.4, sqrt6)) | degree 0: 3.1212 (`combo30_d0.log`); degree 2: 3.1211 (`combo30_d2.log`); degree 6, three rounds: raw 3.0991, corrected 3.1160 (`combo30_d6b.log`) |
| `combo30p.py` | the same, with the sample set pruned so that more rounds fit in memory | degree 6: raw 3.0995, corrected 3.1103 (`combo30p_d6.log`); degree 8, two rounds: raw 2.8946 then 2.9338, corrected 4.4901 then 3.8157 (`combo30p_d8.log`, `combo30p_d8b.log`); resumed with 3000 pair and 1200 triple samples kept per kind (`combo30q_d8.log`): raw 2.9349, 2.9388, 2.9414, 2.9422 and corrected 3.3005, 3.1528, 3.0062, 2.9971 in four rounds; the round-3 certificate `combo30q_d8_r3.npz` is the one proved (copied to `multi_cap/radial_certificates/combo30_d8.npz`) |
| `typed_delsarte.py`, `code_feasible.py`, `code_feasible2.py`, `code_feas_gen.py` | typed linear programming bound and local search for two-shell codes with points in holes | exploration only (`slackhole_d6.log`) |
| `case_value.py` | value of the two-point case programme of the paper for given count constraints | used to choose splits |
| `multi_cap/combo30_check.py` | the exact check that a certificate of `combo30p.py` has to pass: exact positivity, the pair inequalities by Bernstein branch and bound with Pi in Arb, the bins, the typed triple inequalities by Taylor branch and bound, and the bound over the five count vectors | **PASS** on `combo30q_d8_r3.npz` with margins 2e-5, 5e-6, 2e-6: 47 646 pair boxes, 4131 bin intervals, 8 275 056 triple boxes, largest bound 3.087203 at counts (22, 1, 1, 0, 6) < 3.10330 (`C30/combo30_check_r3.log`, 6891 s) |

The raw value is the optimum on the samples; the corrected value adds, for each
kind of pair and triple, the largest violation found between the samples times
the number of such pairs or triples. At degree 6 the raw value rises with each
round of added samples and the corrected value stays above 3.10330.

## CM: twenty-five to twenty-nine centres (twenty-nine proved)

`combo_gen2.py` generalises `combo30.py` to an arbitrary case file
(`case29_all.json`: 29 centres, types A, B, F, seven bins and the count
constraints proved by `prop:C-radial`), and bounds the maximum over all count
vectors by linear-programming duality instead of enumerating them
(`combo_gen.py` enumerates them and runs out of memory). With the two-point
kernel alone (degree 0 in the three-point part) the value is 3.1967
(`c29all2_d0.log`), far above the level.

`combo_direct.py` runs the same programme through Clarabel's own interface
instead of cvxpy, in about half the memory (environment: `DSM` the
factorisation, `THREADS`, `TAG`). The case files `case25_all.json` to
`case28_all.json` carry the constraints of `prop:C-radial` at 25 to 28
centres (5330, 17 845, 42 927 and 74 550 count vectors), and `case28_t4.json`
splits the close centres of 28 into four types. At 29 and degree 8
(`c29all_d8.log`, continued from its samples in `c29x_d8.log`) the raw value
is 3.0287, 3.0568, 3.0592, 3.0605 in four rounds and the corrected value
comes down to 3.1026, but the corrected values of these early rounds were too
low (see the next paragraph): the precheck of the round-2 certificate of the
run with the compass search (`c29r_d8.log`) gives 3.196 (`pre_c29r_r2.log`).

The run that closes 29 continues that one. Resumed with `REFINE=100`
(`c29r3_d8.log`), its raw value is 3.06337 to 3.06413 in four rounds, with
corrected values 4.86, 3.09964, 5.78 and 3.227: pruning to the samples nearest
to violation dropped whole regions of pairs, which came back as large
violations a round later. With `KEEP_BASE=1`, which keeps the starting grid of
pair samples through every pruning, the resumed run `c29r4` (`c29r4_d8.log`)
gives raw 3.06421 and corrected 3.10016 in its first round. The precheck of
that certificate at margins 5e-6, 5e-7, 1e-6 is 3.099635
(`pre_c29r4_r1_small.log`), and the exact check
`multi_cap/combo_case_check.py case29_all.json
radial_certificates/combo29_d8.npz 8 5e-6 1e-6 1e-6` passes with 3.101462 <
3.10330 (`multi_cap/runs/combo29_check.log`, 104 minutes): this is
`thm:count29`, and `combo29_d8.npz` is `c29r4_d8_r1.npz`. Four types at 29
(`c29t4c_d8.log`) and degree 10 (`c29d10_d10.log` here, which diverged, and
`bigmachine/records/c29all_d10.log` on a 16-core machine, raw 2.9353,
3.0004, 3.0105, 3.0119 in four rounds) were not needed.

Below 29 the programme stays above the level on its samples, as the solver
returns it: at 28, 3.10668 with three types (`c28x_d8.log`, one round)
and 3.11169 with four (`c28t4_d8.log`, round 2); at 27, 3.1570 (round 1 on the
larger machine, `bigmachine/records/c27all_d8.log`); at 26, 3.30766 at
degree 6 (`c26s_d6.log`) and 3.2456 at degree 8 (`c26h_d8.log`). Only the
degree-6 solve at 26 ends `optimal`; the degree-8 solves end
`optimal_inaccurate`, and every degree-10 and degree-12 solve at 27 and 28
ends with `NumericalError` (its solutions still give the certificates of the
mixes below). Samples only lower the true optimum, so these values suggest,
but do not prove, that no certificate of that degree exists. At 29 degree
10 lowers the raw value by about 0.05 against degree 8, so 28 is run at
degree 10 next. Runs on a larger machine go through
`bigmachine/run_combo.py`.

Each round of `combo_gen2.py` (and so of `combo_direct.py`) looks for the
violations between the samples by random points and grids and, unless
`REFINE=0`, by a compass search from the worst of them and from the worst
samples (`REFINE` starts per kind, default 60), which only moves to admissible
points of larger value. On the round-2 certificate at 29 that search finds
violations two to four times larger than the random points for the triples
(for example 1.8e-4 against 8.0e-5 for ABF), which is why the corrected values
of the earlier rounds were too low. `multi_cap/combo_case_check.py` with
`PRECHECK=1` computes the thresholds of the exact check and the bound they
give, without the branch and bounds, to say whether a certificate is worth
the full check.
Each round adds the `NEWP` worst new pair samples (default 3000) and the
`NEWT` worst new triple samples (default 400) per kind, the triples drawn from
`TPROBE` random points (default 8000) and a grid, before the compass search.
With `BORDER=1` the exact check rounds the bordered matrix
[[A_0, z], [z^T, t]] as a whole instead of setting t = z^T A_0^{-1} z after
rounding A_0, which keeps t at the solver's value when A_0 is nearly singular
(at 28 and degree 10 the difference is 0.007 in the bound).
At 28 and degree 10 the rounds after the first stopped improving: the raw
value stayed at 3.0913 and 3.0911 while the triple constraints failed by about
1e-3 between the samples (corrected 6.48 and 6.45). The solver meets its own
samples to 5e-7; the failures lie on the boundary of the Gram domain, the
triples of directions in one 2-plane through the centre, which the grid and
the random points reach only thinly. With `TCOPLANAR=n` every kind of triple
carries a fixed n-by-n grid of such coplanar triples, equally spaced in the
two angles and kept through every pruning, and the probes add a grid of twice
that density (`c28m_d10.log`). Its first round gave 3.09166 on the samples;
the triple violations between them fell from about 1e-3 to between 4e-5 and
4.4e-4, and the corrected value from 6.45 to 3.957. The worst violations now
sit inside the Gram domain, at triples with three nearly equal inner products
near 0.31, so the second round adds those points; a bound at 28 needs the
triple violations below about 2e-6. Round 2 gave 3.09216 on the samples and
3.227 corrected, round 3 3.09435 and 4.247: the violations move from one kind
to another, and the worst of rounds 1 to 3 sit at triples fixed by the
exchange of two points of the same type (BBB at (0.3316, 0.3316, 0.3316),
ABB at (-0.6188, -0.6188, 0.2783), AAB at (-0.6524, 0.3879, 0.3879)), where a
symmetric kernel has its critical points and from which the pruning between
rounds drops the samples once they are satisfied. With `TSYM=n` every kind
with two points of one type carries a fixed n-by-n grid on that symmetric
slice (u13 = u23 or u12 = u13, and for three equal types also the diagonal),
kept through every pruning, and the probes add one of twice that density
(`c28s_d10.log`, continued from the round-3 certificate). Its first round gave
3.09379 on the samples, but with the symmetric triples held the solver moved
weight onto the pairs and used the gaps between the pair samples, thinned to
5000 per kind by the pruning: the pair violations reached 0.56 and the
corrected value 134.9. The second round put those pair points back and gave
3.09513 (corrected 29.4); the run was stopped in round 3. Over the five rounds
at degree 10 the value on the samples went 3.09166, 3.09216, 3.09435, 3.09379,
3.09513, while the violations between samples moved from kind to kind and never
fell below about 1e-4. The sampled value bounds the optimum over all admissible
pairs and triples from below, and the margins of the exact check (2e-5 on each
of the 378 pairs and 2e-6 on each of the 3276 triples, as at 29 and 30) cost
about 0.014, so the exact check would need that optimum below about 3.089,
which the samples already exceed. The whole case at 28 is out of reach of the
kernel at degree 10.

The case at 28 does not split usefully at the count that limits the bound
at 29. With at least 19 of the centres within 2.05 (`case28_hiA.json`), the
kernel of three-point degree 6 gives 3.23113 in one round, the value of the
whole case at that degree (`c28hiA_d6.log`). At degree 10 the round-3
certificate is flat in the counts as well (`count_effect.py`,
`c28m_r3_counts.log`): its largest value over the count vectors, 3.09436, is at
22 centres within 2.0161, one in (2.05, 2.1], one in (2.2, 2.35] and four
beyond 2.35, and asking for 10 to 21 centres within 2.05, 24 within 2.2, 25 to
27 within 2.35, or at most 20 to 22 within 2.0161 changes it by at most
1.6e-4. Sharper radial counts would not lower the bound at 28; the obstacle
is the arrangement of 22 nearly touching centres with six farther out.

The case does split at the count of centres within 2.0161. Replacing the
condition N(2.0161) <= 24 of `case28_all.json` by N(2.0161) <= k
(`case28_lo{k}.json`) lowers the kernel steadily; in one round of three-point
degree 6 (floating point, `c28lo{k}_d6.log`):

| k | 24 (whole case) | 21 | 20 | 19 | 15 | 13 | 8 |
|---|---|---|---|---|---|---|---|
| value | 3.23113 | 3.21554 | 3.21021 | 3.20482 | 3.18288 | 3.17164 | 3.14307 |

At degree 8 the value drops by about 0.11 (k = 19: 3.09003 in round 1,
`c28lo19_d8.log`; k = 20: 3.10234, `c28lo20_d8.log`, stopped when the smaller
limits came in). `c28lo13_d8.log` ran k = 13 at degree 8 for four rounds: 3.06625, 3.07708, 3.07872, 3.07913 on the samples, corrected at best 3.14343 (round 3; in round 4 the pair functions bulged between samples); `c28lo13_d8b.log` continues it with the starting pair grid kept through every pruning (`KEEP_BASE=1`): round 1 gave 3.07929 on the samples and 3.09797 corrected, and the precheck of that certificate (`CM/c28lo13b_d8_r1.npz`, margins 5e-6, 5e-7, 1e-6, `BORDER=1`, `pre_c28lo13b_r1_small.log`) gives 3.100569 at (13, 0, 11, 0, 0, 0, 4), below 3.10330. Its first exact check (`CM/check_c28lo13b_r1.log`, 3.5 hours) proved the positivity, the six pair inequalities, the bins and nine of the ten triple inequalities, and failed on the kind FFF: a box with value 7.693729e-04 above the threshold 7.693708e-04. That threshold is the largest sampled value refined by SLSQP from the 40 best samples, plus margin3, and the sampling stopped below a narrow peak near 7.7176e-04 (ascent from 600 starts). The rerun raises that one threshold by 1e-5 (`TRIPLE_EXTRA=FFF:1e-5` in `multi_cap/combo_case_check.py`), which adds 4e-5 to the bound at (13, 0, 11, 0, 0, 0, 4), where four triples are of kind FFF, and passes (`CM/check_c28lo13b_r2.log`, copied to `multi_cap/runs/combo28_lo13_check.log`): 15 805 222 triple boxes and the largest bound 3.100609 at (13, 0, 11, 0, 0, 0, 4), below 3.10330. This was the first form of `prop:count28-few`, and `c28lo13b_d8_r1.npz` is `multi_cap/radial_certificates/combo28_lo13_d8.npz`. Over all count vectors the certificate reaches no further: with its proved thresholds the bound is 3.105310 at k = 14 and 3.149929 at k = 24, between 0.0043 and 0.0047 more for each further close centre (`CM/c28_kscan.py`, `CM/c28_kscan.log`). It does reach two more cases once the centres beyond 2.35 are bounded, since the worst vectors have four of them: at most 14 within 2.0161 with at most 3 beyond 2.35 (3.102452 at (14, 0, 11, 0, 0, 0, 3), `CM/case28_k14f3.json`) and at most 15 within 2.0161 with none beyond (3.101488 at (15, 0, 13, 0, 0, 0, 0), `CM/case28_k15f0.json`); at 15 with one beyond it would be 3.103335, just above. The pair, bin and triple inequalities do not involve the counts, so the passing check proves these cases once the exact largest bound over their count vectors lies below the level; the `EXTRA_CASES` option of `multi_cap/combo_case_check.py` computes it (`multi_cap/runs/combo28_extra_cases.log`). Narrower bins then reach one centre further. Three centres in the bins b, b' and b'' have inner products at most a(rho_b, rho_b'), a(rho_b, rho_b'') and a(rho_b', rho_b''), with rho_b the upper end of bin b, because a increases in each argument; so the inequality on triples of their kind holds on that smaller set with a smaller constant, and the same goes for pairs. The options `REFINE_TRIPLES` and `REFINE_PAIRS` of `multi_cap/combo_case_check.py` compute such constants with the same margins and check them by branch and bound; with `ONLY=NONE` the kinds themselves are not checked again (the passing check above covers them, with the same thresholds). For the triples in the bins (1, 3, 7) and (1, 6, 7) the constant drops from 3.2152e-4 to 3.1880e-4, for (3, 3, 3) from 5.4355e-4 to 5.3549e-4, for (1, 1, 3) from 5.6233e-4 to 5.6168e-4, and for the pairs in (3, 6) and (6, 6) from 6.59702e-2 to 6.59505e-2 and 6.59697e-2; the bins (1, 7, 7) and (1, 1, 6) give no measurable gain. The check (`multi_cap/runs/combo28_refined_check.log`, 13 315 362 further boxes, 8 792 068 of them for the combinations that lower a constant) gives 3.102386 at (14, 0, 9, 0, 0, 1, 4) over the 62 865 count vectors with at most 14 within 2.0161 (`CM/case28_k14.json`), and 3.102823 at (15, 0, 8, 0, 0, 3, 2) over the 55 611 with at most 15 within 2.0161 and at most 2 beyond 2.35 (`CM/case28_k15f2.json`), both below 3.10330; at 15 with 3 beyond it would be about 3.10479 at (15, 0, 8, 0, 0, 2, 3), above. Smaller margins then reach one case further on each side. The pair inequalities take the margin 1e-6, except the kind AA, which keeps 5e-6 because it comes within 1e-5 of its largest value all over the face of its domain where the two centres touch (`PAIR_MARGIN=AA:5e-6` in `multi_cap/combo_case_check.py`); the triples take 5e-8 (FFF still 1e-5 more) and the bins 1e-7. The full check with these margins and the refined bins (`multi_cap/runs/combo28_squeeze_check.log`: 3 532 617 pair boxes, 8 511 bin intervals, 17 006 005 triple boxes and 14 396 003 refined boxes) gives 3.099739 at (14, 0, 9, 0, 0, 1, 4) over the 62 865 count vectors with at most 14 within 2.0161, 3.102197 at (15, 0, 8, 0, 0, 2, 3) over the 62 063 with at most 15 within 2.0161 and at most 3 beyond 2.35 (`CM/case28_k15f3.json`), and 3.102468 at (16, 0, 7, 0, 0, 5, 0) over the 27 423 with at most 16 within 2.0161 and none beyond 2.35 (`CM/case28_k16f0.json`), all below 3.10330 (the floating-point precheck is `multi_cap/runs/combo28_squeeze_precheck.log`). This is `prop:count28-few`. With its thresholds the bound is 3.104750 at (15, 0, 8, 0, 0, 1, 4), with 4 beyond 2.35, and 3.104002 at (16, 0, 7, 0, 0, 4, 1), with 1 beyond, both above, and it rises to 3.114620 at 17 close and 3.148022 at 24 (`CM/c28_open_cases.py`, `CM/c28_open_cases.log`, floating point); narrower bins do not close that gap, so the next step needs a certificate built for more close centres. `CM/c28_config_test.py` checks the implementation on actual sets of 28 directions and distances, where the bound B(Y) of the kernels must lie above the union of caps U(Y): on the 12-close witness U = 2.263412 and B = 2.593852, and on the root system with four centres in deep holes U = 3.196365 and B = 3.702009 (`CM/c28_config_test.log`). A further round from those samples (`c28lo13_d8c.log`) gave 3.07941 and 3.10486 corrected and was stopped for the check.

At 27 centres the same split at degree 6 (one round each, floating point;
`c27lo16_d6.log`, `c27lo18_d6.log`, and `c27k20s_d6.log`, `c27k22s_d6.log`,
`c27k24s_d6.log` with exactly k within 2.0161, `case27_k{k}.json`):

| k | 24 (= whole case, 3.28290) | 22 | 20 | at most 18 | at most 16 |
|---|---|---|---|---|---|
| value | 3.28289 | 3.27258 | 3.26203 | 3.25132 | 3.24048 |

The bound rises by about 0.005 per close centre, as at 28, and each value
lies about 0.052 above the one at 28 with the same k (3.23113 at 24, 3.21021
at 20). Exactly 24 within 2.0161 is the binding part, as it is at 25 and 26.

At degree 10 the whole case at 27 (`c27k10_d10.log`) gives 3.12320 on the
samples, nearly flat in the number k of centres within 2.0161 (3.1224 at
k = 0, 3.1232 at k = 22). Rounds aimed at parts of the case, with the
penalty `TRACE_REG=3e-4` and the starting grids kept (`KEEP_BASE=1`), do
better on their samples: at most 14 within 2.0161 (`case27_lo14.json`;
`c27lo14b_d10.log` from the samples of the whole case, `c27lo14c_d10.log`
from those of the first) gave 3.09976 and 3.10122 on the samples (3.23083
and 3.15011 corrected), and at most 10 (`case27_lo10.json`,
`c27lo10_d10.log`, from the samples of `c27lo14c`) 3.07963 (3.12474
corrected). Mixed with the older certificates (`CM/c28_hull_multi.py` with
`CASEM=27 LOWS=0,14,20,21,22,24`) the two rounds for at most 14 give 3.126974
over at most 14 within 2.0161 and 3.105624 over at most 10, where the worst
vectors have three centres beyond 2.35; over the vectors with at most two
beyond 2.35 (`FILTER="v[6]<=2"`) the same certificates give 3.102609, and
3.102616 with the dense-search maxima of `CM/c28_peak2.py`
(`CM/peak2_c27lo10f2.log`, `CM/c27_peakcorr.py`). Its exact check
(`CM/mix27_lo10f2.npz` on `case27_lo10f2.json`, margins 1e-6, 1e-7 and 1e-7,
precheck `CM/pre_c27lo10f2.log`) passes: 3.102938 at (10, 10, 5, 0, 0, 0, 2)
over 22 374 vectors, after 1 193 081 pair boxes and 99 741 630 triple boxes
(`multi_cap/runs/combo27_lo10f2_check.log`, 58 minutes). With the round for
at most 10 added, the mix over the 2530 vectors with at most 10 within 2.0161
and exactly 3 beyond 2.35 (`FILTER="v[6]==3"`) gives 3.098631, 3.098724 with
the dense-search maxima (`CM/peak2_c27lo10f3.log`), and its exact check
(`CM/mix27_lo10f3.npz` on `case27_lo10f3.json`, margins 2e-6, 5e-7 and
1e-7, precheck 3.100556, `CM/pre_c27lo10f3.log`) passes: 3.100556 at
(10, 10, 4, 0, 0, 0, 3), after 576 954 pair boxes and 105 075 559 triple boxes
(`multi_cap/runs/combo27_lo10f3_check.log`, 50 minutes). Together these two
settle 27 with at most 10 within 2.0161. A check
of the first mix over at most 9 within 2.0161 and any number beyond 2.35
(`CM/check_c27lo9.log`, precheck 3.100668) was stopped once the second mix
covered that set. With the type counts fixed, the bound moves with the bins
only through the bin constants, and the constant of [2, 2.0161] exceeds that
of (2.0161, 2.05] by 0.0053 to 0.0055 in these mixes, so each further centre
within 2.0161 costs about 0.005. The vector where the first mix was worst,
(10, 11, 3, 0, 0, 0, 3), looks unrealisable: the best of 200 penalised
minimisations of the overlaps and of the departures from the bins still
misses, with two centres at distance 1.9853 (`CM/c27_config_search.py`,
`CM/c27_config_search.log`; floating point, evidence only).

At exactly 11 within 2.0161 the same certificates, the round for at most 10
included, give 3.104043 over all the vectors, again worst at three beyond
2.35, and 3.101403 over those with at most two beyond 2.35 (3.101433 with the
dense-search maxima, `CM/peak2_c27k11f2.log`). The exact check of that
combination (`CM/mix27_k11f2.npz` on `case27_k11f2.json`, margins 1e-6, 1e-7
and 1e-7, precheck 3.101756, `CM/pre_c27k11f2.log`) passes: 3.101756 at
(11, 9, 5, 0, 0, 0, 2) over 2034 vectors, after 613 982 pair boxes and
105 397 437 triple boxes (`CM/check_c27k11f2.log`, copy
`multi_cap/runs/combo27_k11f2_check.log`, 47 minutes). With the same
certificates the bound is 3.120346 at exactly 14 within 2.0161, 3.141904 at
18, 3.157926 at 21 and 3.170496 at 24, each worst at three beyond 2.35
(floating point). A round aimed at exactly 11 within 2.0161
(`case27_k11.json`, `c27k11_d10.log`, from the samples of `c27lo10`) gave
3.08563 on its samples (3.11374 corrected). With it the mix over the 230
vectors with exactly 11 within 2.0161 and 3 beyond 2.35 gives 3.097257
(3.097926 with the dense-search maxima, `CM/peak2_c27k11f3.log`), and its
exact check (`CM/mix27_k11f3.npz` on `case27_k11f3.json`, margins 2e-6, 5e-7
and 1e-7, precheck 3.099786, `CM/pre_c27k11f3.log`) passes: 3.099786 at
(11, 8, 5, 0, 0, 0, 3), after 1 275 045 pair boxes and 119 805 055 triple
boxes (`CM/check_c27k11f3.log`, copy `multi_cap/runs/combo27_k11f3_check.log`,
66 minutes). The four checks together prove (C) at 27 with at most 11
within 2.0161.

The round for exactly 11 also lowers the bound further up: with it the mixes
give 3.101177 at exactly 12 within 2.0161 with at most two beyond 2.35 and
3.102693 with three (`case27_k12f2.json`, `case27_k12f3.json`), and 3.1137
at 14, 3.1354 at 18, 3.1515 at 21 and 3.1652 at 24 (floating point). The
first precheck at 12 with at most two beyond 2.35 (`CM/pre_c27k12f2.log`,
3.100839) left the AAF and AFF thresholds up to 2.25e-6 below the
dense-search maxima of `CM/peak2_c27k12f2.log` (a narrow AAF peak at inner
products (-0.106, 0.608, 0.608), on the edge of the range); with those two
raised (`TRIPLE_EXTRA=AAF:2.4e-6,AFF:1.2e-6`) the precheck gives 3.101565
(`CM/pre_c27k12f2x.log`), and its exact check passes: 3.101565 at
(12, 6, 7, 0, 0, 0, 2) over 2034 vectors, after 2 056 809 pair boxes and
125 144 848 triple boxes (`CM/check_c27k12f2.log`, copy
`multi_cap/runs/combo27_k12f2_check.log`, 75 minutes). So (C) also holds at
27 with exactly 12 within 2.0161 and at most two beyond 2.35.
A round aimed at exactly 12 (`case27_k12.json`, `c27k12_d10.log`, from the
samples of `c27k11`) gave 3.09128 on its samples (3.10975 corrected). With it
the mix over the vectors with exactly 12 within 2.0161 and 3 beyond 2.35
(`case27_k12f3.json`) gives 3.099786 (3.099910 with the dense-search maxima,
`CM/peak2_c27k12f3.log`). Its precheck (`CM/pre_c27k12f3.log`, 3.101307) left
the AAF threshold 7.3e-7 below the dense-search maximum, at the same narrow
peak as at 12 with at most two beyond 2.35; with it raised
(`TRIPLE_EXTRA=AAF:1.3e-6`) the exact check passes: 3.101859 at
(12, 6, 6, 0, 0, 0, 3) over 230 vectors, after 975 258 pair boxes and
134 115 496 triple boxes (`CM/check_c27k12f3.log`, copy
`multi_cap/runs/combo27_k12f3_check.log`, 66 minutes). The six checks
together prove (C) at 27 with at most 12 within 2.0161.

At exactly 13 within 2.0161 the mixes with the round for 12 give 3.103476
with at most two beyond 2.35 and 3.105368 with three (`case27_k13f2.json`,
`case27_k13f3.json`, floating point), above the level. A round aimed at
exactly 13 (`case27_k13.json`, `c27k13_d10.log`, from the samples of
`c27k12`) gave 3.09689 on its samples (3.13850 corrected). With it the mix
for at most two beyond 2.35 gives 3.102635 (3.102903 with the dense-search
maxima of `CM/peak2_c27k13f2.log`); its precheck (`CM/pre_c27k13f2.log`,
margins 1e-6, 1e-7 and 1e-7, every triple threshold 1e-7 above the
dense-search maximum) gives 3.103225, and its exact check passes: 3.103225 at
(13, 5, 7, 0, 0, 0, 2) over 2034 vectors, after 2 587 765 pair boxes and
143 913 003 triple boxes (`CM/check_c27k13f2.log`, copy
`multi_cap/runs/combo27_k13f2_check.log`, 91 minutes). So (C) also holds at
27 with exactly 13 within 2.0161 and at most two beyond 2.35. With three
beyond 2.35 the mix gives
3.103560 after six iterations, still above the level. A round aimed at
exactly 14 (`case27_k14.json`, `c27k14_d10.log`, from the samples of
`c27k13`) gave 3.10236 on its samples (3.10928 corrected). With it the mix
for 13 within 2.0161 and three beyond 2.35 gives 3.101190 in floating point
(margins of the check and climbed maxima included, worst at
(13, 4, 7, 0, 0, 0, 3)). Its precheck (`CM/pre_c27k13f3.log`, margins 1e-6,
1e-7 and 1e-7) gives 3.102153 but leaves the AAF threshold 3.6e-8 below the
dense-search maximum (`CM/peak2_c27k13f3.log`); with it raised
(`TRIPLE_EXTRA=AAF:2e-7`) the exact check passes: 3.102225 at
(13, 3, 8, 0, 0, 0, 3) over 230 vectors, after 2 493 357 pair boxes and
162 369 735 triple boxes (`CM/check_c27k13f3.log`, copy
`multi_cap/runs/combo27_k13f3_check.log`, 100 minutes). The eight checks
together are the part of `prop:count27-few` with at most 13 within 2.0161. At exactly 14 the mixes first gave 3.106334
with at most two beyond 2.35 and 3.107087 with three, above the level; with
the rounds described below they give 3.104966 with at most two beyond
2.35 (with the four rounds for at most two beyond 2.35; 3.104543 on the
samples of the mixing LP), still above the level (floating point, margins of
the check and climbed maxima included), and exactly three beyond 2.35 is
proved, as described after the rounds. The
optima of the rounds on their own samples rise by about 0.0055 per close
centre (3.08563, 3.09128, 3.09689 and 3.10236 for 11 to 14); if that trend
holds, a round for 15 lands near 3.108 on its samples. That round
(`case27_k15.json`, `c27k15_d10.log`, from the samples of `c27k14`, 3.0 hours)
gave 3.10766 on its samples (3.11850 corrected), above the level even before
any gap between samples, so degree 10 in this form very likely cannot settle
15 or more within 2.0161 at 27 (evidence, not a proof: the solve carries the
trace penalty). A second round for exactly 14 (`c27k14b_d10.log`, from the
samples of `c27k14`) gave 3.10240 on its samples (3.11028 corrected), no gain.
Rounds for exactly 14 with at most two beyond 2.35 alone (`case27_k14f2.json`,
from the samples of `c27k14b` and then each from the previous one) have more
room on their samples: 3.10067 (3.24496 corrected, `c27k14f2_d10.log`),
3.10083 (3.12683, `c27k14f2b_d10.log`), 3.10091 (3.10880,
`c27k14f2c_d10.log`) and 3.10094 (3.11619, `c27k14f2d_d10.log`), all reported
optimal_inaccurate by the solver. The mixes with them stall near 3.105, since
even on the samples of the mixing LP the best combination gives 3.1045,
above the level. Two rounds for exactly 14 within 2.0161 and exactly three
beyond 2.35 (`case27_k14f3.json`, `c27k14f3_d10.log` from the samples of
`c27k14b`, then `c27k14f3b_d10.log` from those of the first) gave 3.09438
(3.23804 corrected) and 3.09444 (3.13891). With the first, the mix for all
230 vectors with exactly three beyond 2.35 stopped at 3.103356 (with the
dense-search peaks fed back to the mixing programme), just above the level,
at vectors whose middle centres all lie within 2.1. Split by whether a centre lies in (2.1, 2.35],
the 219 vectors with one mix to 3.096158 and pass the exact check:
3.098199 at (14, 5, 4, 1, 0, 0, 3) after 1 255 354 pair boxes and 112 999 796
triple boxes (`CM/check_c27k14f3g2.log`, copy
`multi_cap/runs/combo27_k14f3g2_check.log`, 80 minutes, margins 2e-6, 5e-7
and 1e-7 with FFF raised 1e-6; dense peaks `CM/peak2_c27k14f3g2a.log`). The
11 others are mixed with the middle type narrowed to (2.05, 2.1]
(`case27_k14f3g1_narrow.json`, peaks on `case27_k14f3g1_proxy.json`), which
with the second round gives 3.101813 (precheck 3.101755,
`CM/pre_c27k14f3g1_b.log`), and checked on `case27_k14f3g1.json` with the
thresholds of every pair and triple kind with a middle centre taken per
combination of the bins 1, 2, 3 and 7 (`REFINE_PAIRS`, `REFINE_TRIPLES`,
`SKIP_COVERED=1`): 3.102269 at (14, 8, 2, 0, 0, 0, 3) after 2 486 047 pair
boxes and 159 075 380 triple boxes (`CM/check_c27k14f3g1.log`, copy
`multi_cap/runs/combo27_k14f3g1_check.log`, 158 minutes). The same split
for at most two beyond 2.35 gives 3.104588 for the vectors with a centre in
(2.1, 2.35] and, for the 39 others narrowed (`case27_k14f2g1_narrow.json`),
3.104173, then 3.103629 with a round built for that narrowed case
(`c27f2g1_d10.log`, from the samples of `c27k14f2c` inside the narrowed
ranges: 2.99574 on its samples, 13.27450 corrected); with its second round
(`c27f2g1b_d10.log`, 3.02158 on its samples, 3.60058 corrected) the mix
gives 3.102577 (precheck 3.102851) and the exact check passes: 3.102974 at
(14, 3, 8, 0, 0, 0, 2) after 1 980 579 pair boxes and 220 104 690 triple
boxes (`CM/check_c27k14f2g1.log`, copy `multi_cap/runs/combo27_k14f2g1_check.log`,
167 minutes, margins 3e-7, 5e-8 and 1e-7). The vectors with a centre in
(2.1, 2.35] are split again at 2.2. The 1868 with a centre in (2.1, 2.2]
(`case27_k14f2g2b.json`) mix to 3.102039 with each middle bin a type of its
own for the thresholds (`CM/c27_hull_sub.py`, `SUB=2,3,4,5`; precheck
3.102749, `CM/pre_c27k14f2g2b.log`), and their exact check passes:
3.102942 at (14, 5, 5, 0, 1, 0, 2) after 4 552 754 pair boxes and
1 712 578 763 triple boxes (margins 5e-7, 5e-8 and 1e-7, 10 certificates).
Container restarts stopped it twice, so it ran in parts: the pair pieces were
proved by the start of 03:43 on 2026-10-20 (`CM/check_c27k14f2g2b.log`) and
stored in the checkpoint `CM/check_c27k14f2g2b.ckpt`, a resume from 05:17
(`CM/check_c27k14f2g2b_r2.log`) stored the triple setups in
`CM/check_c27k14f2g2b.ckpt.setup`, and the resume from 06:02
(`CM/check_c27k14f2g2b_r3.log`, copy `multi_cap/runs/combo27_k14f2g2b_check.log`)
proved the 77 triple pieces and the bound and lists every piece with its
result. The 127 with every middle centre in
(2.05, 2.1] or (2.2, 2.35] (`case27_k14f2g2a.json`) stopped at 3.103583,
even with a round built for them (`c27f2g2b_d10.log`, 3.10055 on its
samples, 3.11126 corrected, from the samples of the round `c27f2g2_d10.log`
for all vectors with a centre in (2.1, 2.35]). Since the three-point
inequality of a certificate holds for every assignment of its types to the
points, each certificate may type each bin as any of A, B and F; mixing
(certificate, typing) columns (`CM/c27_hull_mt.py`, typings AABBBBF and
AABBBFF, the second giving (2.2, 2.35] the type F) gives 3.102001 (precheck
3.102597, `CM/pre_c27k14f2g2amt.log`), and its exact check
(`multi_cap/combo_mixtype_check.py` on the per-bin case
`case27_k14f2g2a_bins.json`, `CM/check_c27k14f2g2amt.log`, copy
`multi_cap/runs/combo27_k14f2g2amt_check.log`, 375 minutes, margins 1e-6,
1e-7 and 1e-7) passes: 3.103112 at (14, 6, 4, 0, 0, 1, 2) after 3 838 506
pair boxes and 598 800 571 triple boxes, 12 certificates, 0.38 of the weight
on the typing AABBBFF and 0.34 on the round `c27f2g2b`. Widening
the far type to [2.2, sqrt 6] (`case27_k14f2w.json`) gives 3.386, so the
certificates cannot be used outside their far range. At exactly 15 within
2.0161 and three beyond 2.35 the vectors with a centre in (2.1, 2.35]
(`case27_k15f3g2.json`, 192 vectors) mix to 3.098688, below the level, and
pass the exact check: 3.099309 at (15, 4, 4, 1, 0, 0, 3) after 2 459 509
pair boxes and 103 106 474 triple boxes (`CM/check_c27k15f3g2.log`, copy
`multi_cap/runs/combo27_k15f3g2_check.log`, 94 minutes, margins 1e-6, 1e-7
and 1e-7 with AFF and FFF raised 1e-6 after the dense search
`CM/peak2_c27k15f3g2a.log`). The 10 others give 3.108423 even narrowed;
three rounds built for them with A = [2, 2.0161] and B = (2.0161, 2.1]
(`case27_k15f3g1_t1.json`; `c27k15g1t_d10.log`, `c27k15g1t2_d10.log`,
`c27k15g1t3_d10.log`: 2.92918, 3.07254 and 3.08924 on their samples,
20.46, 3.97 and 3.33 corrected), mixed with free typing, leave the mixing
LP at 3.10494 on its samples, above the level, at (15, 8, 1, 0, 0, 0, 3)
(`CM/sub/mix_c27k15f3g1mt.log`, `CM/sub/mix_c27k15f3g1mt4.log`). A fourth
and a fifth round (`c27k15g1t4_d10.log` with the trace penalty raised to
1e-3, `c27k15g1t5_d10.log`) reach 3.09745 and 3.10009 on their samples but
3.18203 and 3.18672 between them, and the full mix with them still settles
at 3.1048 to 3.1054, so degree 10 stops here; a round of degree 12 from the
samples of the fifth (`c27k15g1d12.log`) needs more than 10 GB and was
killed for memory beside an exact check.
These ten vectors look unrealisable: with 15 within 2.0161, 9 in
(2.0161, 2.1] and 3 beyond 2.35 the best of 200 penalised minimisations
still leaves two centres at 1.98754 (`CM/c27_config_search_k15f3g1.log`),
and the best of 10000 tries at 1.98906 (`CM/c27_config_search_k15g1f3.log`),
and at (15, 8, 1, 0, 0, 0, 3) at 1.97948 (`CM/c27_config_search_k15f3g1w.log`),
while 24 within 2.08 and 3 beyond 2.35 are found within ten tries
(`CM/c27_config_search_24f3_208.log`) and 24 within 2.06 and 3 beyond 2.35
are not (best distance 1.99325, `CM/c27_config_search_24f3_206.log`). The
pair and triple functions of the best mix are flat (88 to 90 % of the
sampled triples within 1e-4 of each maximum, `CM/c27_triple_profile.py`,
`CM/c27_triple_profile.log`), so counting triples by region, as
`REGION_TRIPLES` does at twenty-eight, gains nothing there. Floating
point, evidence only.
Mixed one vector at a time (`CM/mix27_k15f3g1_pervector.txt`), the ten
give about 2.9787 + 0.014 a at (15, a, 9 - a, 0, 0, 0, 3) for a = 0..8
(3.090543 at a = 8) but 3.105044 at a = 9, where all 24 centres within 2.1
lie within 2.05, so that one vector is what stops the whole mix. The nine
with a <= 8 (`case27_k15f3g1lo8.json`) mix to 3.090624 under five typings,
0.93 of the weight typing the bin (2.05, 2.1] as A (`CM/mix27_k15f3g1lo8mt.log`),
and pass the exact check (above): 3.098238 at (15, 8, 1, 0, 0, 0, 3). For
(15, 9, 0, 0, 0, 0, 3) a round for the four vectors with 15 or 16 within
2.0161, 24 within 2.1 (at most one beyond 2.05) and 3 beyond 2.35
(`case27_k15a9.json`, `c27k15a9_d10.log`, from the samples of `c27k14f3b`)
reaches 3.10380 on its own samples, already above the level, and 5.55
corrected; with it the single-vector mix gives 3.104398
(`CM/mix27_k15f3g1a9.log`), so that vector stays open.
The same search finds packings in the other open parts. With at least 15
within 2.0161, 20 within 2.1, 21 within 2.15, 22 within 2.2 and 25 within
2.35 (the counts of prop:C-radial when at most 2 lie beyond 2.35) the first
600 tries of seed 1 leave two centres at 1.999815
(`CM/c27_config_search_k15f2mono.log`), but seed 2 finds a packing at try
4318 of 10000, at (15, 0, 5, 1, 1, 3, 2), with most radii at the ends of
their ranges and 28 distances equal to 2
(`CM/c27_config_search_k15f2mono_b.log`); 15 within 2.0161, 10 in
(2.0161, 2.35] and 2 beyond 2.35 with no other counts pack at once
(`CM/c27_config_search_k15f2free.log`). At least 16 within 2.0161 with the
counts 20, 21, 22 and 24 within 2.1, 2.15, 2.2 and 2.35 pack at try 276, at
(16, 0, 4, 1, 1, 2, 3) (`CM/c27_config_search_k16mono.log`). At least 17
within 2.0161 with these counts give no packing in 600 tries of seed 1 (best
distance 1.999775, `CM/c27_config_search_k17mono.log`) and 10000 of seed 2
(1.999301, `CM/c27_config_search_k17mono_b.log`). So 15 within 2.0161 with
at most 2 beyond 2.35 and 16 with 3 beyond hold real packings and need
certificates; 17 or more may not occur. Floating point, evidence only.
At 16 and 17 the vectors with a centre in (2.1, 2.35] give 3.106091 and
3.113100. At 16 within 2.0161 and 3 beyond 2.35 (`case27_k16f3g2.json`,
165 vectors) each middle bin a type of its own lowers the mix to 3.103447
at the climbed maxima, and split at 2.2 the 15 vectors with none in
(2.1, 2.2] mix to 3.087267 and are proved (above), while the 150 with a
centre in (2.1, 2.2] (`case27_k16f3g2b.json`) stop at 3.103284 on the
samples and 3.103497 at the climbed maxima, 0.0002 above the level. A round
for the 165 from the samples of `c27k14f3b` (`c27k16f3g2_d10.log`) reaches
3.00148 on its samples and 72.85 between them; a second round from its
samples (`c27k16f3g2b_d10.log`) was stopped before it solved, since mixed
typing does better: with most of the weight typing the bin (2.05, 2.1] as A
(typings AABBBBF, AABBBFF and AAABBBF, `CM/c27_hull_mt.py`) the 150 mix to
3.101422 at the climbed maxima, the precheck at margins 1e-6, 1e-7 and 1e-7
gives 3.101991 (`CM/pre_c27k16f3g2bmt.log`), and their exact check passes
at that value (`CM/check_c27k16f3g2bmt.log`, above). The 9 vectors with 16 within 2.0161,
3 beyond 2.35 and none in (2.1, 2.35] mix to 3.114348 with the bin
(2.05, 2.1] a type of its own, and look unrealisable too: with 16 within
2.0161, 8 in (2.0161, 2.1] and 3 beyond 2.35 the best of 10000 tries leaves
two centres at 1.98290 (`CM/c27_config_search_k16g1f3.log`; floating point). With the rounds for 13 and 14 added to the pool, one iteration of
the mix gives 3.119733 at 16 within 2.0161, 3.131219 at 18, 3.144641 at 21
and 3.163428 at 24, each worst at three beyond 2.35 (floating point, margins
of the check included; before these two rounds the values were 3.1243,
3.1349, 3.1513 and 3.1680). For exactly 15 within 2.0161, at most two
beyond 2.35 and a centre in (2.1, 2.2] (`case27_k15f2g2b.json`, 1650
vectors) a round from the samples of `c27f2g2` (`c27k15f2g2b_d10.log`)
reaches 3.10173 on its samples and 3.21195 corrected; mixed with 28 earlier
certificates, each middle bin a type of its own for the thresholds, it gives
3.103550 on the first samples of the mixing LP, but 3.106087 and 3.106412
once the climbed maxima are fed back (the pool alone, on fewer samples,
gave 3.105865). A second round from its samples (`c27k15f2g2b2_d10.log`),
first stopped to give the machine to a degree 12 round and then run again,
reaches 3.10177 on its samples and 3.12584 corrected; a third, from its
samples (`c27k15f2g2b3_d10.log`), was stopped before it solved, to leave
memory for the mixes below. The
pool gives 3.107342 over the 36 vectors with 15 within 2.0161, at most two
beyond 2.35 and none in (2.1, 2.35], and 3.107266 over the 115 with none in
(2.1, 2.2] and one in (2.2, 2.35] (floating point). The degree 12 round for
the ten vectors with 15 within 2.0161 and three beyond 2.35, from the
samples of the fifth round above, fits in 9.6 GB alone on the machine. At
trace penalty 3e-4 it stopped on InsufficientProgress, its inequalities
broken by up to 0.11 on the samples (`c27k15g1d12b.log`); at 1e-3 it
reaches 3.04362 on its samples, against 3.10009 for the best round of
degree 10, but 11.22727 corrected (`c27k15g1d12c.log`). The second round
of degree 12, from its samples, reaches 3.07205 on its samples and 4.06246
corrected (`c27k15g1d12d.log`), and the third 3.08024 and 3.35771
(`c27k15g1d12e.log`; the degree 10 chain gave 3.08924 and 3.33 at its third
round), and the fourth 3.08373 and 3.35796 (`c27k15g1d12f.log`): the
values between samples stopped falling, as the degree 10 chain did at
3.18, so the chain was stopped there. Mixed with the degree 10
certificates the four rounds lower the floor only from 3.104545 to
3.104383 on the samples (3.105111 at the climbed maxima). For 15 within
2.0161, at most 2 beyond 2.35 and none in (2.1, 2.35] (36 vectors) the
rounds above get no weight (mix 3.107655, 3.108208 at climbed maxima), and
a degree 10 round of its own on the same narrowed types
(`case27_k15f2g1_t1.json`, from the samples of the fifth degree 10 round)
stopped on a numerical error with the trace penalty 3e-4 and no usable
certificate (`c27k15f2g1t_d10.log`); the retry with 1e-3 reaches 2.98886
on its samples and 16.12364 between them (`c27k15f2g1tb_d10.log`), and
mixed with the pool it gets weight 0.001 and leaves the mix at 3.10746 on
the samples (3.10855 at the climbed maxima). A second round from its
samples (`c27k15f2g1t2_d10.log`) was stopped before it solved, to give
the machine to the round for 16 within 2.0161.
Mixed typing reaches further on these parts (`CM/c27_hull_mt.py`, floating
point, evidence only). Over the 1650 vectors with 15 within 2.0161, at most
2 beyond 2.35 and a centre in (2.1, 2.2], the typings AABBBBF, AABFFBF,
AABBBFF and AABFFFF (middle bins typed F) with the second round give
3.104466 at (15, 2, 7, 1, 0, 0, 2), 0.23 of the weight on that round; on
that mix only the 66 vectors with exactly one centre in (2.1, 2.35] are
above the level. The 1584 with two or more, mixed on their own, reach
3.102613 (`CM/mix27_k15f2g2bEmt.log`), but their dense search
(`CM/peak2_c27k15f2g2bEmt.log`) finds peaks with a far centre up to 6.8e-6
above the mix's thresholds, which lifts them to 3.104916 with 25 vectors
above the level; re-mixed with those peaks fed back as samples they reach
3.103303 (`CM/mix27_k15f2g2bEmt2.log`), 2e-6 below the level with no
margins. Per vector, at the margins of a check (1e-6, 1e-7, 1e-7),
only the 60 with exactly one centre in (2.1, 2.2] and one in (2.2, 2.35]
are above the level (3.103502); they (`CM/case27_k15f2g2bE11.json`) mix
on their own to 3.101942 (`CM/mix27_k15f2g2bE11mt.log`), and their exact
check runs after their dense search (`CM/peak2_c27k15f2g2bE11mt.log`,
`CM/check_c27k15f2g2bE11mt.log`). The other 1524
(`CM/case27_k15f2g2bEr.json`, which drops the 60 with the case key
`exclude`) peak at 3.101633 under that re-mix; after the dense search of the
re-mix (`CM/peak2_c27k15f2g2bEmt2.log`) their exact check started at margins
1e-6, 1e-7 and 1e-7, was stopped when the thresholds showed the peak of
B3B3B4 5e-8 above its threshold (that kind would fail), and runs again with
B3B3B4 raised by 2e-7 (`CM/check_c27k15f2g2bErmt2.log`; 3.102853 in
floating point with its own thresholds). The 66 split into the 6 with no
centre in (2.05, 2.1] (`case27_k15f2g2bH0.json`) and the 30 each with one in
(2.05, 2.1] and the middle centre in (2.1, 2.15] or in (2.15, 2.2]
(`case27_k15f2g2bH4.json`, `case27_k15f2g2bH5.json`). The last two mix to
3.096621 and 3.096048; their dense searches
(`CM/peak2_c27k15f2g2bH4mt.log`, `CM/peak2_c27k15f2g2bH5mt.log`) find peaks
up to 3.8e-6 above the thresholds, about 3.0995 and 3.0991 with them, the
margins (2e-6, 5e-7, 1e-6) and extra room, and their exact checks pass
(`CM/check_c27k15f2g2bH4mt.log`: U <= 3.099453 at (15, 8, 1, 1, 0, 0, 2),
741 minutes; `CM/check_c27k15f2g2bH5mt.log`: U <= 3.099138 at (15, 8, 1, 0,
1, 0, 2), 728 minutes; above). The 6 mix to 3.101844, 3.102614 with the
peaks their dense search finds; re-mixed with those peaks fed back they
reach 3.102414 (`CM/mix27_k15f2g2bH0mt2.log`), the dense search of the
re-mix (`CM/peak2_c27k15f2g2bH0mt2.log`) finds every peak within 8e-7 of
the thresholds, and the exact check at margins 5e-7, 5e-8 and 1e-7
(`CM/check_c27k15f2g2bH0mt2.log`) verified every inequality but one: the
triple kind A2B4F7, whose dense-search peak lies 2.6e-7 above its threshold
and which had no extra room, failed (value 1.532870e-3 above 1.532793e-3).
The rerun raises A2B4F7 by 4e-7, takes the other pieces from the
checkpoint and passes (U <= 3.102737 at (15, 9, 0, 0, 1, 0, 2), above), so
all 66 vectors with exactly one centre in (2.1, 2.35] are proved. Over the 36 with none in (2.1, 2.35] a pool
of 14 certificates under five typings, among them AAABBBF, gives 3.099707,
0.2 of the weight on `c27f2g1b` typed AAABBBF (`case27_k15f2g1_bins.json`);
its dense search (`CM/peak2_c27k15f2g1mt.log`) finds an A2B3B3 peak 1.6e-5
above the samples. Its exact check (`CM/check_c27k15f2g1mt.log`, margins
2e-6, 5e-7 and 1e-6) verifies every inequality, but with the extra room put
on ten triple thresholds (2e-5 on A2B3B3 alone) the bound is 3.103547 at
(15, 4, 6, 0, 0, 0, 2), above the level, so that run fails; with the extra
room on A2B3B3, A1A2B3, A2A2B3, B3B3B3 and A1A1B3 lowered to 1e-6, 1e-6,
5e-7, 5e-7 and 5e-7 the same thresholds give 3.102529 in floating point,
and a rerun that checks only those five kinds again, the rest taken from
the checkpoint, passes at that value (U <= 3.102529 at (15, 0, 11, 0, 0, 0,
1), above). Over the 115 with
none in (2.1, 2.2] and one in (2.2, 2.35] the same pool gives 3.106653,
and only the 12 with no centre in (2.05, 2.1] are above the level; with the
bin (2.2, 2.35] split at 2.25 the mixing programme for those 12 stays at
3.106459 on its samples. Per vector, under the mix of the 115, the 3 with
one centre in (2.2, 2.35] give 3.1065, the 3 with two 3.1049 and the 6
with three or more at most 3.10341. The 6 (`CM/case27_k15f2g2aX3.json`) mix
on their own to 3.098894 (`CM/mix27_k15f2g2aX3mt.log`), about 3.100984 with
the peaks of their dense search, and their exact check at margins 1e-6,
1e-7 and 1e-7 passes (`CM/check_c27k15f2g2aX3mt.log`: U <= 3.102002 at
(15, 7, 0, 0, 0, 3, 2)) after one triple kind, A1B6B6, failed at a narrow
peak its dense search had missed and was checked again with 6e-7 more room,
the rest taken from the checkpoint. The 3 with two (`CM/case27_k15f2g2aX2.json`)
mix to 3.102991, but to 3.104405 with the peaks of their dense search
(`CM/peak2_c27k15f2g2aX2mt.log`, A2F7F7 1.3e-5 above), and re-mixed with
those peaks fed back the mixing programme stays at 3.103657 on its samples,
above the level with no margins; they and the 3 with one centre in
(2.2, 2.35] have 23 to 26 centres within 2.05, next to
(15, 9, 0, 0, 0, 0, 3). Mixed one vector at a time under five typings
(`CM/case27_k15f2g2aX1v*.json`, `CM/mix27_k15f2g2aX1v*mt.log`), four of these six pass their exact checks, (15, 11, 0, 0, 0, 1, 0)
(`CM/check_c27k15f2g2aX1v1110mt.log`: U <= 3.098938, 102.5 million triple
boxes), (15, 10, 0, 0, 0, 2, 0) (`CM/check_c27k15f2g2aX1v1020mt.log`:
U <= 3.096137, 118.0 million) and (15, 10, 0, 0, 0, 1, 1)
(`CM/check_c27k15f2g2aX1v1011mt.log`: U <= 3.100562, 225.5 million, after
a first run failed on A2B6F7 and was resumed from its checkpoint with more
room there) and (15, 9, 0, 0, 0, 2, 1)
(`CM/check_c27k15f2g2aX1v921mt2.log`: U <= 3.103026, 368.0 million, mixed
again with the dense-search peaks of a first mix fed back, after a first
run failed on A1A1A1 and A2B6F7 and was resumed from its checkpoint with
more room there), and these four are in the paper;
(15, 9, 0, 0, 0, 1, 2), mixed on its own, stays at 3.106600 and is open,
and (15, 8, 0, 0, 0, 2, 2), mixed on its own, gives 3.103115 but 3.105325
with the peaks of its dense search, and mixed again with those peaks
3.103269 at small margins (`CM/mix27_k15f2g2aX1v822mt2.npz`), whose dense
search runs (`CM/peak2_c27k15f2g2aX1v822mt2.log`); it is open. The 103
others, mixed on their own, reach 3.101923 (`CM/mix27_k15f2g2aYmt.log`), but
their dense search (`CM/peak2_c27k15f2g2aYmt.log`) finds peaks with a far
centre up to 4.9e-6 above the thresholds (3.102708 with them); re-mixed with
those peaks fed back they reach 3.102697 (`CM/mix27_k15f2g2aYmt2.log`), the
dense search of the re-mix (`CM/peak2_c27k15f2g2aYmt2.log`) finds every peak
within 8.4e-7 of the thresholds, and the exact check at margins 5e-7, 5e-8
and 1e-7, with extra room on eight triple kinds (about 3.103077 in floating
point with the peaks, margins and extra room), passes after a container
restart, resumed from its checkpoint (`CM/check_c27k15f2g2aYmt2.log`:
U <= 3.103081 at (15, 8, 1, 0, 0, 1, 2), 7 585 532 pair boxes, 1367 bin
intervals, 653 692 639 triple boxes; copy `multi_cap/runs/combo27_k15f2g2aY_check.log`). At 16 within 2.0161, 3
beyond 2.35 and none in (2.1, 2.35], the 8 vectors with at most 7 in
(2.0161, 2.05] (`case27_k16f3g1lo7.json`) mix to 3.099519, their dense
search (`CM/peak2_c27k16f3g1lo7mt.log`) puts every peak within 1.0e-6 of
the samples, and their exact check passes
(`CM/check_c27k16f3g1lo7mt.log`: U <= 3.101848 at (16, 7, 1, 0, 0, 0, 3),
429 minutes); (16, 8, 0, 0, 0, 0, 3) is left, like (15, 9, 0, 0, 0, 0, 3).
At 16 within 2.0161 and at most 2 beyond 2.35 the parts follow those at
15. With none in (2.1, 2.35] the 33 vectors give 3.104063 together; the 23
with at most one beyond 2.35 (`CM/case27_k16f2G1a.json`) mix to 3.101015
and the 8 with two beyond 2.35 and neither all nine middle centres in
(2.0161, 2.05] nor all in (2.05, 2.1] (`CM/case27_k16f2G1bm.json`) to
3.101333; the 23 pass their exact check (`CM/check_c27k16f2G1amt.log`:
U <= 3.102925 at (16, 0, 10, 0, 0, 0, 1), 351 minutes, after a first run
whose extra room counted the gap of each peak twice and gave 3.104441),
and so do the 8 (`CM/check_c27k16f2G1bmmt.log`: U <= 3.102229 at
(16, 8, 1, 0, 0, 0, 2), resumed from its checkpoint after a restart of the
machine); (16, 0, 9, 0, 0, 0, 2) and
(16, 9, 0, 0, 0, 0, 2), with 25 centres within 2.1, give 3.104116 and
3.103832 and are open. With exactly one centre in (2.1, 2.2] and none in
(2.2, 2.35] the 60 vectors give 3.107876, and per vector only the 24 with
at most three in (2.05, 2.1] are above the level; the 36 others
(`CM/case27_k16f2Hn3.json`) mix to 3.092122 and pass their exact check
(`CM/check_c27k16f2Hn3mt.log`: U <= 3.096243 at (16, 4, 4, 0, 1, 0, 2),
946 240 pair boxes, 3354 bin intervals, 610 085 598 triple boxes, resumed
from its checkpoint after restarts of the machine; copy
`multi_cap/runs/combo27_k16f2Hn3_check.log`); a mix for the 18 with one to three in
(2.05, 2.1] gives 3.103291, and at most 3.099892 on the 12 with two or
three (`CM/case27_k16f2Hm2.json`), whose exact check passes
(`CM/check_c27k16f2Hm2mt.log`: U <= 3.101346 at (16, 6, 2, 1, 0, 0, 2),
1 256 629 pair boxes, 5353 bin intervals, 707 041 803 triple boxes, about
3.101345 expected with the peaks of its dense search
`CM/peak2_c27k16f2Hm2mt.log`, the margins and extra room; copy
`multi_cap/runs/combo27_k16f2Hm2_check.log`); the 6 with exactly one
(`CM/case27_k16f2Hn1.json`) mix to 3.102636, but to 3.104253 with the
peaks of their dense search (`CM/peak2_c27k16f2Hn1mt.log`) and the
margins; mixed again with those peaks they give 3.103299
(`CM/mix27_k16f2Hn1mt2.npz`), and with the peaks of the new mix about
3.1038 on the 2 with two beyond 2.35, which are open, and 3.103244 on the
4 with at most one (`CM/case27_k16f2Hn1f1.json`), whose exact check runs
(`CM/check_c27k16f2Hn1f1mt.log`). With none in (2.1, 2.2] and one or more in
(2.2, 2.35] the 103 vectors give 3.111974, each centre in (2.05, 2.1]
lowering the bound by about 0.003; a mix for the 79 with at least two in
(2.05, 2.1] gives 3.103361, and 3.101372 on the 67 with at least three
(`CM/case27_k16f2G2an3.json`), whose exact check passes
(`CM/check_c27k16f2G2an3mt.log`: U <= 3.102140 at (16, 5, 3, 0, 0, 1, 2),
3 590 588 pair boxes, 1600 bin intervals, 616 011 287 triple boxes, about
3.102149 expected with the peaks of its dense search
`CM/peak2_c27k16f2G2an3mt.log`, the margins and extra room; copy
`multi_cap/runs/combo27_k16f2G2an3_check.log`); a mix
for the 12 with exactly two in (2.05, 2.1] gives 3.103355, and 3.100707 on
the 9 with at least two in (2.2, 2.35] (`CM/case27_k16f2G2ae2b.json`),
whose exact check passes (`CM/check_c27k16f2G2ae2bmt.log`: U <= 3.102282
at (16, 5, 2, 0, 0, 2, 2), 3 612 007 pair boxes, 1604 bin intervals,
641 317 975 triple boxes, about 3.102292 expected with the peaks of
`CM/peak2_c27k16f2G2ae2bmt.log`, the margins and extra room; copy
`multi_cap/runs/combo27_k16f2G2ae2b_check.log`; not yet in the paper,
which takes these 9 together with the 3 below). Of the 3 with exactly
two in (2.05, 2.1] and one in (2.2, 2.35], (16, 8, 2, 0, 0, 1, 0), mixed
on its own (`CM/case27_k16f2v8210.json`), passes its exact check
(`CM/check_c27k16f2v8210mt.log`: U <= 3.099725, 230.4 million triple
boxes; not yet in the paper), (16, 7, 2, 0, 0, 1, 1), mixed on its own
(`CM/case27_k16f2v7211.json`), gives 3.098756 (3.099854 with the peaks of
`CM/peak2_c27k16f2v7211mt.log`) and its exact check runs
(`CM/check_c27k16f2v7211mt.log`), and (16, 6, 2, 0, 0, 1, 2), mixed on its
own (`CM/case27_k16f2v6212.json`), gives 3.100933 (3.102064 with the peaks of
`CM/peak2_c27k16f2v6212mt.log`) and its exact check runs
(`CM/check_c27k16f2v6212mt.log`). The 1372 vectors with a centre in (2.1, 2.2] and at least two in
(2.1, 2.35] mix to 3.107964 with 225 above the level (scratchpad scout, 9
certificates under four typings, `CM/mix27_k16f2Emt.log`); under that mix
the 1109 with at least two in (2.1, 2.2] and at least three in (2.1, 2.35]
(`CM/case27_k16f2Ex.json`) give 3.102171 at (16, 8, 0, 0, 2, 1, 0), but
3.104914 at (16, 0, 6, 1, 1, 1, 2), with 38 above the level, with the
peaks of their dense search (`CM/peak2_c27k16f2Exmt.log`), so they are
mixed again with those peaks; the 81 with exactly two in (2.1, 2.2] and
none in (2.2, 2.35] (up to 3.1042 under that mix,
`CM/case27_k16f2E2.json`), mixed on their own, give 3.101190 at
(16, 7, 0, 1, 1, 0, 2) and their dense search runs
(`CM/peak2_c27k16f2E2mt.log`); the 144 with exactly one in (2.1, 2.2] and
one to three in (2.2, 2.35] (up to 3.1078 under that mix) are being mixed
on their own and stand at 3.106889 after three rounds. All of these values are floating
point.

At twenty-five centres with at most 14 within 2.0161 (`case25_lo14.json`:
eight bins, the far range split at 2.444, the counts of row 25 of
`prop:C-radial` and 24 within 2.444 by `lem:twentyfour-close`) the pool of
27 and 28 mixes to 3.457, and a round of degree 10 built for the case
(`c25lo14_d10.log`, from the samples of `c27k10` round 2, `TRACE_REG=3e-4`)
gives 3.15848 on its samples and 3.23999 corrected. All 24 type-count
vectors sit at 3.1585, and with fewer allowed within 2.0161 the bound of the
same certificate drops by 0.00575 per centre, to 3.0779 with none
(`CM/combo_perkey.py`, `CM/c25lo14_perkey.log`), so this form of degree 10
gets below the level only with at most four within 2.0161. Floating point,
evidence only.

At k = 15 the kernel of three-point degree 10 has more room on the samples.
`c28lo15_d10.log`, from the samples of round 2 of `c28s_d10.log`, gave 3.05709
on the samples and 3.92182 corrected, and its continuation `c28lo15b_d10.log`
gave 3.05883 and 3.15940, then 3.05807 and 3.40545; every round ended in a
numerical error of the solver. `c28lo15c_d10.log` continued it with a penalty
on the traces of the three-point blocks in the objective (`TRACE_REG=3e-4` in
`CM/combo_direct.py`), to damp the swings between samples, and was stopped in
its first round, slowed by the exact check beside it, to give its place to
degree 8 at k = 16 (`c28lo16_d8.log`, from the samples of the certificate of
`prop:count28-few`). Its first round gave 3.09593 on the samples (3.14856
corrected), 0.0165 above the certificate at k = 13; with the thresholds of
`multi_cap/runs/combo28_squeeze_precheck.log` that certificate bounds its case
0.0154 above its value on the samples, so a certificate for all of k = 16
would land near 3.111, and the run was stopped. `c28k15f4_d8.log` aimed instead
at the one case with 15 close that the check at smaller margins leaves, 15
within 2.0161 and 4 beyond 2.35 (`case28_k15f4.json`), where the certificate
of `prop:count28-few` gives 3.104750: its first round gave 3.08978 on the
samples. With the thresholds its solver found on its own samples, the
certificate of `prop:count28-few` already gives 3.07924 at k <= 13 (the
programme's 3.07929), 3.09042 on that case and 3.09607 at k <= 16
(`CM/c28_sampled_eval.py`, `CM/c28_sampled_eval.log`), so a certificate built
for either case gains at most about 6e-4 on the samples, less than the
1.45e-3 by which the level is missed. At degree 8 the kernel is saturated near
15 close centres. `c28lo16_d10.log` returned to degree 10 with the trace
penalty, at k = 16, from the samples of `c28lo15b_d10.log`; it was killed for
memory in its first round, with no round finished, when the mixing programme
below was started beside it.

A convex combination of certificates is again a certificate: the positivity
conditions are convex, and every inequality is linear in the kernel, with its
threshold taken for the combination. `CM/c28_hull.py` takes the certificates
at 28 with the same types (degrees 6, 8 and 10; a certificate of lower degree
enters with its three-point blocks embedded in the layout of degree 10, which
changes none of its functions) and finds the convex weights that minimise the
bound over a chosen set of count vectors. With only the weights unknown the
programme is a linear programme, so it takes about 3 million samples, and it
adds the local maxima of each combination, found by compass search, as
further samples until its thresholds stop moving. The certificates of degree
10 stopped on numerical errors of the solver and have eigenvalues near -1e-6
in their blocks. The exact check sets negative eigenvalues to zero before it
rounds, and for these certificates that moves the triple polynomials by up to
1e-5, so the programme projects every certificate onto the positive
semidefinite matrices first (without that step its first values, 3.08622 for
15 close with 4 beyond 2.35 and 3.09929 over the whole case, were not what
the check would see; the precheck of that mix gave 3.11286). With the
projection, in floating point and with the margins of the check, the best mix
gives 3.13159 over the whole case, 3.10618 at k = 17, 3.10076 over every
vector with k <= 16, 3.09916 over 15 close with 4 beyond and 16 close with 1
beyond, and 3.10014, 3.10039 and 3.10074 at 16 close with 2, 3 and 4 beyond.
The two mixes below the level that cover the open cases at k <= 16 are
`multi_cap/radial_certificates/combo28_mix1516_d10.npz` (15 close with 4
beyond, 16 close with at most 1 beyond) and
`multi_cap/radial_certificates/combo28_mix16_d10.npz` (k <= 16). Their
prechecks at margins 2e-6 (AA 5e-6), 2e-7 and 5e-7 give 3.098860 on
`CM/case28_k15f4.json` and 3.098060 on `CM/case28_k16f1.json`
(`multi_cap/runs/combo28_mix1516_precheck.log`), and 3.101264 at (16, 7, 1,
0, 0, 0, 4) on `CM/case28_lo16.json` (`multi_cap/runs/combo28_mix16_precheck.log`),
all below 3.10330. The precheck's threshold for the kind AAA lies 2e-7 below
a local maximum that the mixing programme found (and for AAB in the second
mix it has only 1e-8 to spare), so the exact checks raise those thresholds
(`TRIPLE_EXTRA=AAA:4e-7`, and `AAB:2e-7` in the second). The exact checks of
both mixes at three-point degree 10 started on 8 October; in both, every pair
inequality and the bins passed, and both stopped while their triples were
being checked when the container was restarted (`multi_cap/runs/combo28_mix1516_check.log`).
For k = 17, no mix of these certificates reaches the level even on its
samples (3.10447 there), and at degree 8 a certificate built for one case
gains only about 6e-4 on the samples. Certificates of degree 10 built for one
value of k do reach it once mixed. `CM/c28k17_d10.log` builds one for exactly
17 close (`CM/case28_k17.json`): round 1 gives 3.07522 on its samples, but its
triple inequalities exceed their constants between the samples by up to
1.4e-4 (corrected 3.31061), and round 2 was killed for memory. Mixed with the
earlier certificates it gives 3.101579 at k = 17
(`multi_cap/radial_certificates/combo28_mix17_d10.npz`, precheck 3.101230 in
`multi_cap/runs/combo28_mix17_precheck.log`). `CM/c28k18_d10.log` builds one
for exactly 18 close from the samples of the 17-close run, which carry its
refinements; the gaps are far smaller (3.08058 on the samples, corrected
3.13951), and the mix gives 3.100714 at k = 18, 3.10532 at 19, 3.10982 at 20
and 3.12489 at 24 (3.13365 before). One mix then covers every vector with
k <= 18 at 3.100721 (`multi_cap/radial_certificates/combo28_lo18_d10.npz`,
precheck 3.100554 at (18, 3, 3, 0, 0, 0, 4) in
`multi_cap/runs/combo28_lo18_precheck.log`, with the case file
`CM/case28_lo18.json`). The exact checks of the mixes for k <= 16 and k = 17
(restarted on 12 October) were stopped by hand in their triples
(`multi_cap/runs/combo28_mix16_check.log`, `multi_cap/runs/combo28_mix17_check.log`),
since the check of this one covers both; it runs with `TRIPLE_EXTRA=AAA:5e-7`
and the interval evaluation in C (`POLY_EVAL=c`, `multi_cap/ipoly.c`, bit for
bit the numpy one), in `multi_cap/runs/combo28_lo18_checkc.log`; the same
check in numpy (`multi_cap/runs/combo28_lo18_check.log`) closed the same boxes
at its first two levels about 22 times slower and was stopped. The check
**passes** in 52 minutes: 3.101381 at (18, 6, 3, 0, 0, 0, 1) over the 70 719
vectors with k <= 18, after 181 820 pair boxes and 106 993 576 triple boxes.
The chain continues: `CM/c28k19_d10.log` builds a certificate for exactly 19 close from
the samples of the 18-close run (3.08563 on its samples, corrected 3.10827),
and `CM/c28k20_d10.log` one for 20 close from its samples, each certificate
starting from the samples of the one before, towards one mix for the whole
case. Mixed with the 19-close certificate the pool gives 3.096464 over k <= 19,
3.100958 at k = 20 and 3.116546 at k = 24 (floating point, climbed maxima,
margins of the check). The room at k <= 19 allows the triple margin 1.5e-6
(precheck 3.100907, `multi_cap/runs/combo28_lo19_precheck.log`), and the exact
check of that mix (`multi_cap/radial_certificates/combo28_lo19_d10.npz`)
**passes** in 47 minutes (`multi_cap/runs/combo28_lo19_check.log`): 3.100907 at
(19, 2, 7, 0, 0, 0, 0) over the 71 995 vectors with k <= 19, after 112 492 pair
boxes and 134 464 632 triple boxes. `CM/c28k20_d10.log` then builds a
certificate for exactly 20 close from the 19-close samples (3.09051 on its
samples, corrected 3.11096); mixed in, it gives 3.098129 over k <= 20, 3.102751
at k = 21, 3.107033 at 22 and 3.115317 at 24. The k <= 20 mix
(`multi_cap/radial_certificates/combo28_lo20_d10.npz`, precheck 3.099289 at the
triple margin 5e-7) **passes** its exact check in 54 minutes
(`multi_cap/runs/combo28_lo20_check.log`): 3.099289 at (20, 3, 4, 0, 0, 0, 1)
over the 72 996 vectors with k <= 20, after 168 390 pair boxes and 157 796 046
triple boxes. A jump from the 20-close samples straight to 24 close failed
numerically (`CM/c28k24_d10.log`, bound -23804), so the chain goes on one
step at a time: `CM/c28k21_d10.log` builds a certificate for exactly 21 close
from the 20-close samples (3.09515 on its samples, corrected 3.19358). A mix
over 38 certificates with it (`multi_cap/radial_certificates/combo28_lo21_d10.npz`,
fourteen nonzero weights, among them three certificates at 26 centres) gives
3.101466 over k <= 21, 3.105672 at k = 22, 3.105741 over k <= 22 and
3.114126 at 24 (floating point, climbed maxima, margins of the check); its
precheck at the triple margin 5e-7 gives 3.103041
(`multi_cap/runs/combo28_lo21_precheck.log`), and the exact check **passes**
in 54 minutes (`multi_cap/runs/combo28_lo21_check.log`): 3.103041 at
(21, 4, 2, 0, 0, 0, 1) over the 73 722 vectors with k <= 21, after 120 628
pair boxes and 168 808 634 triple boxes. This is the first half of
`prop:count28-mixed`, which contains the earlier checks at k <= 18, 19 and 20.
`CM/c28k22_d10.log` builds the next certificate, for exactly 22 close, from the
21-close samples (3.09890 on its samples, corrected 3.13908), and
`CM/c28k22b_d10.log` (two rounds), `CM/c28k22c_d10.log` and
`CM/c28k22d_d10.log` refine it from the samples of the round before (corrected
3.11819, 3.10373, 3.10201, 3.10563; `CM/c28k22e_d10.log` was stopped). Mixes
over all vectors with k <= 22 stay above the level at the earlier margins
(3.104638 with the round-4 certificate); at margins 2e-7 and 2e-8 their
precheck gives 3.103286, too close to the level, and that exact check was
stopped (`CM/check_c28lo22c.log`). Two changes close k = 22. First, the
margins: the pair inequality on two A centres closes with the margin 2e-7 in
place of 5e-6 (about 190 000 boxes for the round-4 mix), and the triple kind
AAA closes with 5e-8 in 19.7 million boxes, against 13.4 million at 5e-7 for
the k <= 21 mix. Second, k <= 21
is proved, so the mix only has to cover the 451 vectors with exactly 22 close,
and there its weights can favour the 22-close certificates:
`CM/c28_hull_multi.py` (the multi-target form of `CM/c28_hull.py`, with the
reported margins set by `MG_AA`, `MG2`, `MG3`) gives 3.101887 at
(22, 0, 2, 0, 0, 0, 4) with margins 5e-7 for the pairs and 5e-8 for the
triples, against 3.103487 over k <= 22 (`CM/c28_hull_k22d.log`; the
combination is `CM/mix28_k22d.npz`, copied to
`multi_cap/radial_certificates/combo28_k22_d10.npz`: nine nonzero weights,
four of them on the 22-close rounds, one on a certificate at 26 centres). Its
first exact check at margins 1e-6 (pairs, AA included), 2e-7 (triples) and
1e-7 (bins) fails on the triple kind AAF (`CM/check_c28k22d.log`): that
inequality has a narrow peak which the 650 000 samples and the local ascent
from their 40 best miss. A denser search (`CM/c28_peak.py`, 15 million samples
and 300 starts, `CM/peak_k22d.log`) puts its largest value at 2.326572e-03,
3.87e-7 above the first run's value, and the rerun with the AAF threshold
raised by 4e-7 (`TRIPLE_EXTRA=AAF:4e-7`) **passes** in 64 minutes
(`CM/check_c28k22d2.log`, copied to `multi_cap/runs/combo28_k22_check.log`):
3.102969 at (22, 0, 2, 0, 0, 0, 4) over the 451 vectors with k = 22, after
64 474 pair boxes, 2339 bin intervals and 207 350 728 triple boxes. With the
check at k <= 21 this is the first part of `prop:count28-mixed`, at most 22
within 2.0161. `CM/c28k23b_d10.log` builds a certificate for exactly 23 close
from the round-4 22-close samples (3.09585 on its samples, corrected 4.16711),
and `CM/c28k23c_d10.log`, `CM/c28k23d_d10.log` and `CM/c28k23e_d10.log` refine
it from the samples of the round before (3.09733, 3.09776, 3.09792 on their
samples, corrected 3.26595, 3.12889, 3.12978; `CM/c28k23d_try1_d10.log` is a
first run of round 3, killed when the container ran out of memory).
`CM/c28k23f_d10.log` is a fifth round with the trace penalty `TRACE_REG=1e-3`
in place of 3e-4 (3.10027 on its samples, corrected 3.17106: worse on its own,
but it lowers the mixes). The mixes with these rounds (`CM/c28_hull_multi.py`, margins
5e-7 and 5e-8, each triple threshold raised to the largest value that
`CM/c28_peak2.py` finds by dense sampling and local ascent; the climbing of the
mixer alone misses peaks on the faces where two centres touch) give 3.105094,
3.104861 and 3.104693 for exactly 23 close after rounds 2, 3 and 4, and
3.104594 and 3.103447 for exactly 24 after rounds 3 and 4: the 23-close rounds
lower the 24-close mix more than the 23-close one. Most of what separates the
24-close mix from the level is the margins, at the 2600 AAA, 325 AAB and 325
AAF triples and 325 AA pairs of its worst vector. At margins 1e-7 for the pairs
and bins and 1e-8 for the triples its precheck gives 3.103128
(`CM/pre_c28k24e.log`), and the exact check of the mix (`CM/mix28_k24e.npz`,
copied to `multi_cap/radial_certificates/combo28_k24_d10.npz`: six nonzero
weights, 0.90 of them on the 23-close rounds 3 and 4) at margins 1e-7, 3e-8
and 1e-7, with the thresholds raised by 1e-6 for the triple kinds ABB, AFF,
BBB, BBF, BFF, FFF and by 5e-7 for ABF, which have no triples, or 26, at the
worst vectors, **passes** in 44 minutes (`CM/check_c28k24e.log`, copied to
`multi_cap/runs/combo28_k24_check.log`): 3.103206 at (24, 2, 0, 0, 0, 1, 1)
over the 126 vectors with k = 24, after 29 367 pair boxes, 7156 bin intervals
and 117 675 808 triple boxes. This completes `prop:count28-mixed`, at most 22
or exactly 24 within 2.0161; since at most 24 lie within 2.0161
(`thm:kissing-stable`), exactly 23 close was then the only open case at 28.
There, at the same small margins, the round-4 mix exceeds the level only at
three of the five vectors with four centres beyond 2.35 (at most 3.10437) and
stays at most 3.10298 at every other vector; mixing for the worst vector alone
gives the same value. So exactly 23 is split by the number beyond 2.35.
`CM/case28_k23fle3.json` holds the 246 vectors with at most 3 beyond 2.35. The
mix with round 5 (`CM/mix28_k23fle3f.npz`, copied to
`multi_cap/radial_certificates/combo28_k23fle3_d10.npz`: eight nonzero weights,
0.62 on two 22-close certificates and 0.35 on the 23-close rounds 2 to 5) gives
3.102605 after the dense correction and 3.102649 in its precheck
(`CM/pre_c28k23fle3f.log`), and its exact check at margins 5e-7, 1e-7 and
1e-7, with the thresholds raised by 5e-7 for AAF (whose peak, on the face where
the AF inner product is largest, the precheck missed by 2e-7) and by 1e-6 for
FFF (at most one triple at these vectors), **passes** in 55 minutes
(`CM/check_c28k23fle3f.log`, copied to `multi_cap/runs/combo28_k23fle3_check.log`):
3.102775 at (23, 0, 0, 0, 0, 4, 1), after 37 565 pair boxes, 2719 bin intervals
and 162 022 640 triple boxes. So at 28 only exactly 23 within 2.0161 with 4
beyond 2.35 stays open. `CM/case28_k23f4.json` is the case of those five
vectors; the mixes give 3.10408 there at the small margins and about 3.1040
with no margins, above the level, so only better certificates can close it.
`CM/c28k23f4a_d10.log` is a certificate targeted at them from the samples of
round 5 (3.09712 on its samples, corrected 3.23109). `CM/case28_k23f4c.json`
is the same five vectors with the close type narrowed to A = [2, 2.0161] and
B = [2.0161, 2.35]: at these vectors every A centre lies within 2.0161, so the
inequalities with an A centre need to hold only there, which cuts off the AAF
peak at the AF inner product 0.6176; `CM/c28k23f4c_d10.log` is a round on it
from the samples of c28k23f4a. It ends DualInfeasible after 4384 s: a case
with a single type-count vector leaves the programme unbounded, and a round on
such a case needs a second, proved vector beside it (`CM/case28_k23f4h.json`
adds 24 close with 4 beyond 2.35; `CM/c28k23f4h_d10.log`: 3.07461 on its
samples, 5.09653 corrected).

The five vectors are then taken one at a time. (23, 1, 0, 0, 0, 0, 4) needs no
new certificate: the thresholds of the passing check of `case28_k23fle3.json`
give 3.102361 there (`CM/pre_c28k23f4m_extra.log`, `EXTRA_CASES`). For each
of the other four, `CM/c28_hull_case.py` mixes the pool for that vector alone
on a case whose types are narrowed to the bins it occupies
(`CM/case28_k23f4g{2,3,4,5}_narrow.json`: A = [2, 2.0161], B = the bin of the
last centre, F), with 400 climbing starts per kind and the dense-search peaks
of an earlier mix as seeds (`PEAKS`); the dense search itself
(`CM/c28_peak2.py`) needs contiguous bins and runs on a proxy case with the
same upper ends (`CM/case28_k23f4g*_proxy.json`), since the triple domains
depend only on the upper ends. `CM/remap_bins.py` recomputes the bin
constants of the mix on the usual seven bins, and `combo_case_check.py` checks
it on `CM/case28_k23f4n{2,3,4,5}.json` with the thresholds refined to the
occupied bins and `SKIP_COVERED=1`. Dense-corrected at margins 1e-7/1e-7/1e-8
the mixes give 3.103249 for (2.05, 2.1], 3.101123 for (2.1, 2.15], 3.101210
for (2.15, 2.2] and 3.103356 for (2.2, 2.35], and the first three
**pass** their exact checks (`CM/check_c28k23f4n2p.log`,
`CM/check_c28k23f4n3g.log`, `CM/check_c28k23f4n4g.log`): 3.103281, 3.101860
and 3.102128. So at 28 only (23, 0, 0, 0, 0, 1, 4) stays open: 23 within
2.0161, one in (2.2, 2.35] and 4 beyond 2.35. There the mix gives 3.103285
with no margins at all, too close to the level for a check; splitting its
bin at 2.3 gains nothing, since the radial function rises to 2.35 and is flat
from 2.28 on, and a round for that vector alone (`CM/case28_k23f4g5r.json`,
from the samples of c28k23f4h, with the proved 24 + 4 vector beside it)
reached 3.07957 on its samples but 3.66384 between them
(`CM/c28k23f4g5r_d10.log`), as first rounds do. The radial function of the mixes falls by about
0.0022 from 2 to 2.0161, so splitting the close bin at 2.002
(`CM/case28_k23f4n5s.json`) lowers the bound of every vector with a close
centre beyond 2.002 by about 3e-4 per centre, while 23 centres within 2.002
have smaller inner products (the AF bound 0.6126 instead of 0.6140, on which
the AAF peak lies); in floating point this takes the (2.2, 2.35] mix from
3.103356 to 3.103251 (dense-corrected, margins 1e-7/1e-7/1e-8). A round
with the close type narrowed to [2, 2.002] (`CM/case28_k23f4g5x2.json`)
ended DualInfeasible, and a mix with its weights chosen for close centres
within 2.002 gives 3.103234 with the split, hardly better. The precheck of
the split case gives 3.103176, but two of its triple thresholds sat up to
5.4e-8 below dense-search peaks, so the exact check
(`CM/check_c28k23f4n5s.log`, from 12:43 UTC on 2026-10-13) raises each
refined triple threshold above its peak (`TRIPLE_EXTRA`), which should have
given about 3.10328. It **failed** at 14:21 on the triples of bins 1-1-8 (two
centres within 2.002 and one beyond 2.35), at a box centre 7e-10 above the
threshold, after eleven other inequalities had passed. The true maximum of
that inequality on its domain is 2.7474807e-4, at the interior point
(0.42875, -0.81298, -0.81298), a narrow ridge that falls by 3e-6 within 0.01;
the dense search and the precheck had climbed only the peak on the face
where a close and the far direction are at their largest inner product
(2.7466e-4 after the split). The ridge lies inside the domain of every split,
so the split gains nothing at a = 23 and the route is dead. The check now has
`REFINE_GRID=<n>`, which also climbs from every local maximum of an n^3 grid
(n = 121 finds the ridge, n = 41 does not), and with it the precheck of the
unsplit vector (`CM/pre_c28k23f4n5u_grid.log`) gives 3.103360 at margins
1e-7/1e-8/1e-7, the value quoted before. Round 2 of the single-vector round
(`CM/c28k23f4g5r2_d10.log`, 14:40 to 16:08) reached 3.08341 on its samples
and 3.27305 between them. A count of the close-close-far triples does what
the split could not. On this mix the AAF inequality comes within 1.5e-7 of
its maximum only where both close directions have inner product near the
largest value 0.614 with the far one, or both near -0.813. With one of the
two inner products in [-0.7944, 0.58] its maximum is 2.7460032e-4, against
2.7474849e-4 on the whole domain (121-grid local maxima and SLSQP, scratchpad
search). With one at most -0.7944 and the other at least 0.58 it is
2.7420025e-4. Two exact facts bound how many close directions can sit at
either end for a given far direction y, among directions with pairwise inner
products at most c0 = 0.507954. At most 3 have inner product at most -0.7944:
K of them have a sum S with |S|^2 <= K + K(K-1) c0 and <S, -y> >= 0.7944 K. At
most 9 have inner product at least 0.58: their projections to the link
sphere of y have pairwise inner products at most
(c0 - 0.58^2)/(1 - 0.61403^2), so caps on S^2 about them are disjoint. So at
most 36 of the 253 close pairs per far centre take the full threshold, and
the bound drops by 4 * 217 * 1.48e-7, to 3.10323.
`combo_case_check.py` has this as `REGION_TRIPLES=<b>-<b>-<p>:<s1>:<s2>`. The
two new thresholds are proved by branch and bound on sub-boxes, the symmetry
of the polynomial is checked coefficient by coefficient, and the bound takes
the worst split of each pivot's pairs. The exact check
(`CM/check_c28k23f4n5r.log`, with the options of the grid precheck and
`REGION_TRIPLES=1-1-7:-0.7944:0.58`, from 16:09 UTC on 2026-10-13)
**passes** at 17:52: 3.103231, with thresholds 2.746103e-4 (middle) and
2.742084e-4 (mixed) against 2.747585e-4 on the bins, after 21 697 pair boxes,
2713 bin intervals and 200 740 860 triple boxes (copied to
`multi_cap/runs/combo28_k23f4n5_check.log`; certificate
`multi_cap/radial_certificates/combo28_k23f4n5_d10.npz`). With it
`prop:count28-mixed` is (C) at twenty-eight, and (C) is proved from
twenty-eight centres on. The values on the samples rise by about 0.005 per step
(3.07522, 3.08058, 3.08563, 3.09051, 3.09515 for 17 to 21 close), and each mix
lies about 0.006 above the sample value of its last certificate, which is why
the chain was not continued to 24 close: the 24-close case closed with the
23-close certificates instead. Extra certificates at 26,
27 and 29 centres change the mixes by less than 0.001
(`CM/c29d10_d10_r1.npz`, a failed run with coefficients near 1e12, has to be
left out of any mix).
A certificate for some k reduces (C) at 28 to a statement about directions
alone: no 28 centres satisfying the radial counts of 28 (21 within 2.1, 22
within 2.15, 23 within 2.2, 24 within 2.35) have k + 1 of them within 2.0161.
Since a(r, s) increases in both distances, it suffices to rule out the
directions with every centre at the largest distance its count allows.
`CM/code_feas_dist.py` searches for such directions by local minimisation of
the squared violations (`CM/c28_direction_search.log`): with 2000 starts it
fits 9 centres within 2.0161 and none of 10 to 20 (least largest violation
1.6e-3 at 10, 2.8e-2 at 20). Random starts miss rare arrangements: moving one
centre of a 9-centre code in from 2.1 to 2.0161 in small steps
(`CM/c28_continue.py`, `CM/c28_chain.py`) gives a code with 10, and a search
that makes the 21 distances within 2.1 variables and rewards those at 2.0161
(`CM/c28_maxclose.py`) finds 12 in 2000 starts and never 13
(`CM/c28_close12_witness.txt`, slack 1e-5; from it `CM/c28_extend.py` misses
13 by 4.2e-3, and 600 steps of basin hopping, `CM/c28_hop.py`, stay at 12). At 27 the same search finds 18. So the limit k at 28 has to be
at least 12, and at 27 at least 18. A miss in
a local search is evidence, not a proof, and the statement is of the same kind
as the 24-plus-one exclusion below, which three-point bounds do not settle: the
typed three-point bound on 20 centres within 2.0161, one at 2.1, three at 2.35
and four at sqrt6 stays at 0 at degree 6, and so does the one on 23 centres
within 2.0161, one at 2.35 and four at sqrt6 at degrees 8 and 10 (`CM/c28_typed_dirs.py`).

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
rho(eta) of a root system. Through `cor:no-room`, the residual case at thirty
would have needed rho < 0.206 at eta = 0.00368; the kernel on typed triples
settled that case instead (`thm:count30`), and (L) remains the route for
twenty-five to twenty-seven centres, for the cases of twenty-eight that
`prop:count28-mixed` leaves (exactly twenty-three within 2.0161), and for
(G).

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
to 30. At 25 and 27 to 30 the bound comes from two-point certificates with no
split; at M = 24 and M = 26 it comes from the combined kernel of `CM/`
(labelled pairs and typed triples), run against a level with
`multi_cap/combo_case_check.py`. At M = 26, degree 8, with LEVEL=3.3352 on
`CM/case26_all.json` and `multi_cap/radial_certificates/combo26_level_d8.npz`
(the round-3 certificate of `CM/c26h_d8.log`), it passes with
U(Y) <= 3.329050 (`multi_cap/runs/combo26_level_check.log`). At M = 24,
degree 6, with LEVEL=3.3291 on `CM/case24_all.json` and
`multi_cap/radial_certificates/combo24_level_d6.npz` (the round-1 certificate
of `CM/c24lvl_d6.log`, 3.21901 on its samples), it passes with
U(Y) <= 3.253872 (`multi_cap/runs/combo24_level_check.log`); dropping the
pair terms gives only 24 S(2) = 3.335124. At M = 25 the same kernel at degree
6 stays at 3.313 on its samples, and its round-1 certificate prechecks at
3.3768, above the two-point level 3.3274 (`CM/c25lvl_d6.log`). Its largest
value is at the bin counts (17, 0, 4, 1, 1, 2, 0), with exactly the seventeen
centres within 2.05 that `prop:C-radial` proves, as the level at 26 is
largest at its fifteen. So every cell has
volume at least 9 pi^2/8 - max L_M = 7.7742, the largest level being
L_26 = 3.3291, and every packing of unit balls in R^4 has density at most
0.63477. This is below
the three-point bound 0.63611 of Cohn, de Laat and Salmon, and above
pi^2/16 = 0.61685; it is a density bound, not statement (G) or (C).

| file | what it is |
|---|---|
| `levels.py` | floating point: the least level L that a certificate can reach for M centres, by bisection; used to choose the levels |
| `cert_M.log` | the certificate search, `multi_cap/radial_case_sdp.py M --level L`, writing `multi_cap/radial_certificates/density_M.json`; for M = 26, `cert_26_card.log`, with the spec `density_spec_26_card.json` (split at N(2.03), level 3.3503) |
| `check_M.log` | the exact check, `multi_cap/radial_case_check.py`: positivity by exact LDL^T, K <= Pi by branch and bound with Pi in Arb, the bin bounds, and the bound with no assumption on the count vector. All seven pass: two-point levels 3.3352, 3.3274, 3.3501, 3.3090, 3.2587, 3.2080, 3.1538 for M = 24 to 30 (`check_26_level3353.log` is the earlier check of M = 26 with no split, 3.3524); the paper takes 3.3291 at M = 26 from the combined kernel above |
| `case26_all.json`, `c26all_d2.log` | floating point: an early two-point kernel plus a typed three-point kernel at M = 26 with these coarse types and no counts, value 3.34979; superseded by `CM/case26_all.json` with the counts of `prop:C-radial`, which certifies 3.3291 |

Before the combined kernel, a bound below 0.63611 was sought with the pair
kernel alone, which needs L_26 <= 3.34549. The route tried there splits the
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

More further directions do not help either. `typed_card_m.py` is the same
programme with m further directions, at inner product at most t2 with the 24
and at most t3 among themselves (two centres within sqrt6 at distance at
least 2 have inner product at most 2/3). For 27 centres with 24 within 2.0161
(m = 3, t2 = 0.6141, t3 = 2/3) degree 10 gives sampled Z +3.33
(`typed_m24_3_d10.log`): the three-point bound does not exclude 24 centres
within 2.0161 at 27 centres either.

Counting the near neighbours of the further direction does not help.
Directions of the 24 at inner product above t1 with it lie between 52.1 and
59.5 degrees from it and at least 70.3 degrees apart on its link sphere, so
by the Tammes problem there are at most 9 of them. `typed_card_hole.py` relaxes
only those: with exactly mp of them, the pairs and triples through the further
direction split by how many of their directions are near neighbours, each
kind with its own constant. At degree 10 (`typed_hole_d10.log`) the sampled Z
is +0.245, +0.534, +0.753, +0.950 and +1.123 for mp = 1, 3, 5, 7 and 9,
against -2.5 for mp = 0 (25 directions at 0.508): one near neighbour already
defeats the three-point bound.

A four-point kernel does not help at the degrees that fit here.
`C30/typed4pt.py` adds to the typed three-point kernel of `C30/typed3pt.py` a
kernel that fixes two of the directions, e1 and e2, expands every other
direction x in Chebyshev polynomials of (<e1, x>, <e2, x>) (and, with
`WITH_S=1`, of <e1, e2>), and multiplies by the zonal functions of the plane
orthogonal to e1 and e2; summed over the code it is a sum over pairs, triples
and quadruples, and `typed4pt.py selftest` checks that grouping against the
direct sum over the ordered poles (difference 1e-10). On 25 directions with
inner products at most t the three-point part of degree 8 alone excludes the
code up to t = 0.507 (corrected value -0.0096 there) and gives 0 from 0.508
on; at degree 4 it gives 0 already at 0.5. Four-point parts of degree 3 and 5
at t = 0.508 (over three-point degree 8), and of degree 3 (with or without
<e1, e2>) and 5 at t = 0.5 (over three-point degree 4), all settle at the
value 0 once the samples are refined. The negative values of their first
rounds, down to -31, came from gaps between the quadruple samples (0.035 to
0.053 in round 1), and fell to 0 as the gaps closed. On the 24-plus-one
exclusion at sqrt6 the four-point part of degree 3 gives 0 as well. Degree 5
with <e1, e2> needs more memory than this container has. Commands and output
are in `C30/typed4pt_runs.log`.

Local bounds on the number of neighbours do not help either.
`CM/local_counts.py` adds to the programme of `CM/combo_direct.py`, for each
pair of types (s, t) and each band u >= alpha of inner products of
directions, a bound kappa(s, t, alpha) on the number of centres of type t in
that band around a centre of type s: projected to the link sphere S^2 they are
separated, and the solved cases of the Tammes problem (up to 14 points; the
bound of Fejes Toth beyond) limit how many fit, for example 9 A-neighbours of
an A centre with u >= 0.45 and 7 A-neighbours of an F centre with u >= 0.6.
Each pair inequality is relaxed by multipliers on its bands, and the bound
pays n_s kappa for each. On the whole case at 28 and degree 6 the dual of the
programme describes the relaxed configuration: 24 centres of type A at
distance 2 with about 8.5 neighbours each at u = 0.5, which the count 9
allows, and four of type F at sqrt6, each with about 7.5 A-neighbours at the
largest inner product 0.612. Bands from 0.2 to 0.55 leave the optimum at
3.23114 (3.23113 without them); bands at 0.6 and 0.61 bind for the F centres
and lower it to 3.22891, while the relaxed configuration moves half a
neighbour to just below the band (`CM/local_counts_c28.log`). The near-D4
arrangements that defeat the kernels are locally possible; what excludes
them is global.
