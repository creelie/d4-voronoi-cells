# The second level on an ordinary machine

This directory lets the code of

> D. de Laat, N. M. Leijenhorst, W. H. H. de Muinck Keizer,
> *Optimality and uniqueness of the D4 root system*, arXiv:2404.18794,
> data: 4TU.ResearchData, doi:10.4121/74ce1c25-6fca-4680-8a36-e9c18e7e9594

set up and solve its second-level programme without the construction of the
zonal matrices that its README puts at three days and 128 GB.  The zonal
matrices computed in `../zonal` (about two hours, under 400 MB) are written in
the file format that code reads, and it takes them as its own.

It also compares the two constructions directly, which `../zonal` could not
do: entries computed from scratch by the authors' `compute_PS` agree with ours
coefficient by coefficient.

## The files

| file | what it does |
| --- | --- |
| `zonalstore_4.tar.gz` | the 490 entries of P(S) for n = 4, in the authors' cache format (150 KB) |
| `pkl2cache.py` | writes that cache from `ps.txt.reduced.pkl`, the reduced output of `../zonal`; the tarball is its output |
| `setup_cache.sh` | installs the cache into a copy of the authors' package |
| `las2_slack.jl` | the second-level bound on the enlarged domain [-1, 1/2 + s], for several s |
| `xcheck_zonal.jl` | evaluates entries of Z_lambda with the authors' `evaluate_zonal_matrix` from the installed cache |
| `their_entries.jl` | runs the authors' `compute_PS` one signature at a time, in a clean folder |
| `compare_entries.py` | compares the entries it writes with ours, in exact rationals |
| `runs/` | the logs |

## Running it

Install Julia 1.10 and the authors' package as its README says
(`LasserreSphericalCodes.zip` from the 4TU record, then `Pkg.instantiate()`),
and from its folder:

```
sh /path/to/level2/setup_cache.sh .
julia --project=. -t 4 /path/to/level2/las2_slack.jl 8 10 128 B.txt 0 1//200
```

`setup_cache.sh` also writes empty placeholders for the integrand files
`cache/zonalstore/pol-2-*.txt`; `compute_PS` looks for those before it looks
for the entries themselves, and with them present it prepares nothing.

## The two constructions agree

`xcheck_zonal.jl` evaluates eight entries of Z_lambda, with |lambda| up to 14,
at the inner products (1/3, -2/7, 1/5, -1/4, 2/9, 3/11) through the authors'
code reading our cache, and `xcheck_zonal.py` finds them equal to the values
of `../zonal/zonal.py`, including one reached through the authors' exchange of
the two indices (`runs/xcheck_zonal.log`).  So the format is read as intended.

The entries themselves can be compared as well.  The authors' construction
needs its 128 GB for all signatures at once; run one signature at a time it
fits in a few gigabytes for the smaller ones.  `compare_entries.py` reads the
files written by their `compute_PS`, run from scratch in a clean folder, and
compares them with ours.  It was run for thirteen signatures, (2, 0), (3, 1),
(4, 0), (4, 2), (5, 1), (5, 3), (6, 2), (6, 4), (7, 1), (8, 0), (8, 2), (9, 1)
and (10, 0), and every one of their 88 entries is identical, coefficient by
coefficient (`runs/compare_entries.log`).  The signatures with larger lambda_2
cost more: (10, 2) passed 3 GB, and was stopped there.

## What the programme gives here

At slack 0, on a machine with four cores and 16 GB:

| degrees (d1, delta) | bound | iterations | time | memory |
| --- | --- | --- | --- | --- |
| (4, 6) | 32 | | 100 s in all | |
| (8, 10) | 26.0000 | 68 | 42 min on one core | 3.5 GB |

At (10, 12) an iteration takes about 220 s on four cores and the programme
holds 11.5 GB, eight times the time and three times the memory of (8, 10).
The certificate of the authors uses (14, 16), where the bound is 24.  The
logs are in `runs/`.
