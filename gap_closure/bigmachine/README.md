# bigmachine: runs that need more memory than the cloud container

The cloud container has 15 GB of memory and 4 cores. The programmes here are
run on a larger machine and write their logs next to themselves.

## run_typed.py: twenty-four close centres and one more

The typed three-point programme of `multi_cap/typed_cardinality_sdp.py`: 24
points of S^3 with pairwise inner products at most t1 = 0.508 (centres within
2.0161) and a further point with inner product at most t2 with each of them.
A corrected Z below 0 says that no such code exists.

| t2 | radius | what an exclusion gives |
|---|---|---|
| 0.5101 | 2.025 | the count {24 within 2.0161, at most 24 within 2.025} at 26 centres, level 3.34491, below the three-point density bound |
| 0.5114 | 2.03 | the same with 2.03, level 3.34468 |

At degree 10 the programme does not exclude the code: at t2 = 0.5114 the
sampled Z is +0.066 and the corrected Z +0.179 (`gap_closure/density/typed_sweep_b.log`).
Degrees 16, 18 and 20 are run here, with t2 = 0.5101 and 0.5114 (t2 = 0.508
is 25 centres within 2.0161, which thm:kissing-stable already excludes).

    python -m pip install --user numpy scipy cvxpy clarabel
    python gap_closure/bigmachine/run_typed.py            # degrees 16 18 20
    python gap_closure/bigmachine/run_typed.py 22 --t2 0.5101     # other choices

The log `run_typed.log` records the machine, the package versions and one line
per round. Floating point on sampled constraints: a corrected Z below 0 is what
an exact check would then have to confirm, and no result here is used in a
proof until it has one.
