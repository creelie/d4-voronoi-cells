# The second-level runs this machine cannot do

What was run here (`README.md`): the second-level programme of de Laat,
Leijenhorst and de Muinck Keizer with our zonal matrices, at degrees up to
(d1, delta) = (10, 12), 24.9423 at slack 0 in 3.9 hours and 11.8 GB.  What
was avoided here is the authors' construction of the zonal matrices (their
128 GB step); the programme itself at the degrees of their certificate,
(14, 16), where the bound is 24, has not been run on this machine.

## What it would take

Measured on four cores:

| degrees | memory | time |
| --- | --- | --- |
| (8, 10) | 3.5 GB | 42 min on one core, 68 iterations |
| (10, 12) | 11.8 GB | 3.9 h on four cores, 62 iterations |

From (8, 10) to (10, 12) the memory grew by a factor of about 3.4 and the time
of an iteration by about 8.  If the next steps behave alike:

| degrees | memory (estimate) | time on 4 cores (estimate) |
| --- | --- | --- |
| (12, 14) | about 40 GB | about 1.5 days |
| (14, 16) | about 135 GB | about 12 days (fewer with more cores) |

These are extrapolations from two points and could be off by a factor of two
either way.  The (14, 16) figure is consistent with what the authors report
for their own run.

## The runs, in order

On a machine with at least 64 GB for (12, 14), or 192 GB for (14, 16):

1. Install Julia 1.10 and the authors' package (`LasserreSphericalCodes.zip`
   from the 4TU record, then `Pkg.instantiate()` in its folder), and install
   our zonal cache into it:

       sh /path/to/level2/setup_cache.sh .

2. Reproduce the certificate's value at slack 0 (a check that everything is
   in place; it should print 24 up to the solver's tolerance at (14, 16)):

       julia --project=. -t 8 /path/to/level2/las2_slack.jl 14 16 128 B.txt 0

3. The bound on the enlarged domain, which decides whether the cardinality
   half of Theorem 7.79 extends (it does wherever B < 25):

       julia --project=. -t 8 /path/to/level2/las2_slack.jl 14 16 128 B.txt 1//500 1//200 1//100

4. The margin, which is what statement (i) of "What is left" needs: with N a
   little above the bound B(s) of step 3,

       julia --project=. -t 8 /path/to/level2/las2_margin.jl 14 16 128 1//200 N margin.txt

   It maximises mu with A_2K({x, y}) + SOS_2(u) + mu w(u) = 0,
   w(u) = (u + 1)(u + 1/2)^2 u^2 (1/2 + s - u).  Every inner product of a
   24-point code with slack s then has w(u) <= (N - 24)/mu, so it lies within
   about (8 (N - 24)/mu)^(1/2) of -1/2 or 0, and near -1 or in [1/2, 1/2 + s].

## What the answer can and cannot give

If (N - 24)/mu is small enough that the window is below about 0.01, the
inner products of every 24-point code at slack s lie within 0.01 of the root
system's values, and Lemma 7.77 (the robust root-lattice step) with Theorem
7.51 (the rigidity radius) places the code near the root system.  That is
statement (i) at slack s, which would have to be made rigorous (an exact
rational certificate, as for the authors' own run).  It would not by itself
prove Conjecture 1.6: statements (ii), (iii) and (iv) of "What is left"
remain.  And the window is of order the square root of N - 24: if N - 24 grows
linearly in s, as it does at the first level, a slack of 0.005 gives a window
near 0.1, too wide.
