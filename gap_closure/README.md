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
| `combo30p.py` | the same, with the sample set pruned so that more rounds fit in memory | degree 6: raw 3.0995, corrected 3.1103 (`combo30p_d6.log`); degree 8: `combo30p_d8.log` |
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
