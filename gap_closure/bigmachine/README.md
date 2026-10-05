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
| 0.6141 | sqrt6 | statement (C) whenever 24 centres lie within 2.0161: no further centre fits within sqrt6 |
| 0.51135 | 2.03 | the count {24 within 2.0161, at most 24 within 2.03} at 26 centres, one of the two leaves of the density split in `../README.md` |

At degree 14 in the container the sampled Z at t2 = 0.6141 is +0.018
(`../density/typed2_d14_6141.log`); degrees 16 and 18 are run here.

    python -m pip install --user numpy scipy cvxpy clarabel
    python gap_closure/bigmachine/run_typed.py 16 18 --t2 0.6141

The log `run_typed.log` records the machine, the package versions and one line
per round. Floating point on sampled constraints: a corrected Z below 0 is what
an exact check would then have to confirm, and no result here is used in a
proof until it has one.
