# Twenty-seven centres: how the parts of prop:count27-few were found

This file keeps the record of the mixes, rounds and checks behind
prop:count27-few (pull request #40); the pull request description points
here, as the record outgrew its length limit. Values are floating point
(evidence only) unless an exact check is named; only a PASS of
`multi_cap/combo_mixtype_check.py` or `multi_cap/combo_case_check.py` proves
anything.

With the type counts fixed, the bound depends on the bins only through the bin
constants, and the constant of [2, 2.0161] exceeds that of (2.0161, 2.05] by
0.0053 to 0.0055, so each further centre within 2.0161 costs about 0.005 and
each count needs its own round. With the round for exactly 14 the mixes gave
3.106334 at exactly 14 within 2.0161 with n_7 <= 2 and 3.107087 with n_7 = 3;
with the later rounds (a second for exactly 14, c27k14b_d10.log, 3.10240 on
its samples; one for 15; and four for exactly 14 with n_7 <= 2 alone,
c27k14f2{,b,c,d}_d10.log, 3.10067, 3.10083, 3.10091 and 3.10094 on their
samples but 3.24496, 3.12683, 3.10880 and 3.11619 between them) they give
3.104966 with n_7 <= 2 and 3.105798 with n_7 = 3, both still above the level.
The n_7 <= 2 mixes stall there: on the samples of the mixing LP (random
samples plus climbed maxima) the best combination gives 3.104543. Split in the
same way, n_1 = 14 with n_7 <= 2 gives 3.10459 over the vectors with a centre
in (2.1, 2.35], and 3.10417 over the 39 others with type b narrowed
(gap_closure/CM/case27_k14f2g1_narrow.json); a round built for those 39 on the
narrowed range (c27f2g1_d10.log, 2.99574 on its samples, 13.27 between them)
brings that to 3.10363, and with a second such round (c27f2g1b_d10.log) the
mix passes (above). The vectors with a centre in (2.1, 2.35] are split again
at 2.2: the 1868 with a centre in (2.1, 2.2] mix to 3.102039 with each middle
bin a type of its own for the thresholds and pass their exact check (above);
the 127 with every middle centre in (2.05, 2.1] or (2.2, 2.35] stopped at
3.103583 even with a round built for them (c27f2g2b_d10.log), and mix to
3.102001 with mixed typing (precheck 3.102597), and pass their exact check
(above). Widening the far range to [2.2, sqrt 6] to absorb bins 5 and 6 gives
3.386 and is dead. With n_1 = 15 and n_7 = 3 the 192 vectors with a centre in
(2.1, 2.35] are proved (above); the other 10 give 3.10842 (narrowed), and five
rounds built for them (c27k15g1t{,2,3,4,5}_d10.log on case27_k15f3g1_t1.json;
2.92918, 3.07254, 3.08924, 3.09745 and 3.10009 on their samples, but still
3.18 to 3.19 between them) mixed with free typing leave the mix at 3.1048 to
3.1054, above the level, at (15, 8, 1, 0, 0, 0, 3), so degree 10 stops there;
degree-12 rounds from the samples of the fifth (c27k15g1d12{,b,c,d,e,f}.log)
stopped falling between samples as the degree-10 chain did, and mixed with the
degree-10 certificates they give 3.104383 on the samples and 3.105111 at the
climbed maxima, still above the level (Remark 7.41 now says so). The nine of
the ten with a centre in (2.05, 2.1] are now proved (above); for (15, 9, 0, 0,
0, 0, 3) a round for the four vectors with 15 or 16 close, 24 within 2.1 (at
most one beyond 2.05) and 3 far (c27k15a9_d10.log, from the samples of
c27k14f3b) reaches 3.10380 on its own samples, above the level, and 5.55
corrected, and with it that vector alone mixes to 3.104398
(gap_closure/CM/mix27_k15f3g1a9.log), still above. With n_1 = 15, n_7 <= 2 and
a centre in (2.1, 2.2] (1650 vectors) a round built for them
(c27k15f2g2b_d10.log) reaches 3.10173 on its samples but 3.21195 between them,
and the mix with it gives 3.103550 on the first samples of the mixing LP but
3.106412 once the climbed maxima are fed back; the other n_1 = 15, n_7 <= 2
vectors give 3.107342 (none in (2.1, 2.35], 36 vectors) and 3.107266 (none in
(2.1, 2.2], one in (2.2, 2.35], 115 vectors; 3.107203 with mixed typing) on
the pool; for the 36 the rounds for n_1 = 15, n_7 = 3 get no weight
(3.107655), and a degree-10 round of their own on the narrowed types
(gap_closure/CM/case27_k15f2g1_t1.json) stopped on a numerical error with
trace penalty 3e-4 and no usable certificate (c27k15f2g1t_d10.log); the retry
with 1e-3 reaches 2.98886 on its samples and 16.12364 between them
(c27k15f2g1tb_d10.log), and the mix with it stays at 3.10746 on the samples
(3.10855 at the climbed maxima, weight 0.001 on it); a second round from its
samples (c27k15f2g1t2_d10.log) was stopped before it solved, to give the
machine to the round for 16 close. 16 and 17 close with n_7 = 3 and a centre
in (2.1, 2.35] give 3.10609 and 3.11310. At 16 close with n_7 = 3, each middle
bin a type of its own lowers the 165 vectors with a centre in (2.1, 2.35] to
3.103447 at the climbed maxima; split at 2.2, the 15 with none in (2.1, 2.2]
are proved (above) and the 150 with one in (2.1, 2.2] stop at 3.103284 on the
samples and 3.103497 at the climbed maxima, 0.0002 above the level; a round
for the 165 (c27k16f3g2_d10.log, from the samples of c27k14f3b) reaches
3.00148 on its samples and 72.85 between them, and a second round from its
samples (c27k16f3g2b_d10.log) was stopped before it solved, since mixed typing
with the bin (2.05, 2.1] typed A (AAABBBF) brings the 150 to 3.101422 at the
climbed maxima (precheck 3.101991, gap_closure/CM/pre_c27k16f3g2bmt.log), and
they are now proved (above). The 8 with 16 close, 3 far, none in (2.1, 2.35]
and at most 7 in (2.0161, 2.05] mix to 3.099519 and are now proved (above);
(16, 8, 0, 0, 0, 0, 3) is left, like (15, 9, 0, 0, 0, 0, 3). Over the 1650
with 15 close, n_7 <= 2 and a centre in (2.1, 2.2], typing middle bins as F
with the second round (c27k15f2g2b2_d10.log, 3.10177 on its samples, 3.12584
corrected) gives 3.104466; on that mix only the 66 with exactly one centre in
(2.1, 2.35] are above the level. The 1584 others mix to 3.102613, but their
dense search finds peaks with a far centre up to 6.8e-6 above the mix's
thresholds (3.104916 with them, 25 vectors above the level), and re-mixed with
those peaks fed back they reach 3.103303
(gap_closure/CM/mix27_k15f2g2bEmt2.log), 2e-6 below the level with no margins;
per vector at check margins (1e-6/1e-7/1e-7) only the 60 with exactly one
centre in (2.1, 2.2] and one in (2.2, 2.35] are above the level (3.103502), so
they get a mix of their own (3.101942, gap_closure/CM/mix27_k15f2g2bE11mt.npz;
exact check gap_closure/CM/check_c27k15f2g2bE11mt.log, 3.102924 in floating
point with its own thresholds; it failed one triple kind, A1A2B6, by 2.1e-8 at
a narrow peak its dense search missed, and resumes from its checkpoint with
2e-7 more room there), and the other 1524
(gap_closure/CM/case27_k15f2g2bEr.json; case files may now drop vectors with
the key exclude) peak at 3.101633; their exact check at margins 1e-6/1e-7/1e-7
was stopped when its thresholds showed the B3B3B4 peak 5e-8 above its
threshold, and runs again with B3B3B4 raised by 2e-7
(gap_closure/CM/check_c27k15f2g2bErmt2.log, 3.102853 in floating point with
its own thresholds). These two checks, with 83 and 57 kinds of triples, were
paused for memory at 17:58 on 2026-10-26 and resume from their checkpoints
when fewer checks run. The 66 split into the 6 with no centre in (2.05, 2.1]
(3.101844, 3.102614 with the peaks of their dense search; re-mixed with those
peaks to 3.102414, whose own dense search finds every peak within 8e-7 of the
thresholds; exact check at margins 5e-7/5e-8/1e-7,
gap_closure/CM/check_c27k15f2g2bH0mt2.log, verified every inequality but
A2B4F7, whose dense-search peak lies 2.6e-7 above its threshold and which had
no extra room; the rerun raises it by 4e-7, takes the rest from the checkpoint
and passes, above) and the 30 each with one centre in (2.05, 2.1] and the
middle centre in (2.1, 2.15] or (2.15, 2.2] (3.096621 and 3.096048, about
3.0995 and 3.0991 with the dense-search peaks, margins 2e-6/5e-7/1e-6 and
extra room; both now proved, above). The 36 with 15 close, n_7 <= 2 and none
in (2.1, 2.35] mix to 3.099707 with a pool of 14 certificates and the typing
AAABBBF; their exact check first failed on its bound (3.103547, from the extra
room on the triple thresholds) and with that room lowered on five kinds passes
at 3.102529 (above); of the 115 with none in (2.1, 2.2] and one in (2.2,
2.35], the 103 with a centre in (2.05, 2.1] mix to 3.101923, 3.102708 with the
far-centre peaks their dense search finds; re-mixed with those peaks fed back
they reach 3.102697, the dense search of the re-mix finds every peak within
8.4e-7 of the thresholds, and their exact check at margins 5e-7/5e-8/1e-7
(about 3.103077 in floating point with peaks, margins and extra room) passes
at 3.103081 (above); the 12 others stay at 3.1066, and splitting (2.2, 2.35]
at 2.25 leaves their mixing programme at 3.106459 (floating point). Per vector
the 3 with one centre in (2.2, 2.35] give 3.1065, the 3 with two 3.1049 and
the 6 with three or more at most 3.10341; the 6 mix on their own to 3.098894
and are now proved (above); the 3 with two mix to 3.102991 but 3.104405 with
their dense-search peaks, and re-mixed with those peaks the mixing programme
stays at 3.103657, above the level with no margins, so they and the 3 with one
stay open as sets. Mixed one vector at a time under five typings
(gap_closure/CM/case27_k15f2g2aX1v*.json), four of these six pass exact checks
on 2026-10-26: (15, 11, 0, 0, 0, 1, 0) at 3.098938 (102.5 million triple
boxes, check_c27k15f2g2aX1v1110mt.log), (15, 10, 0, 0, 0, 2, 0) at 3.096137
(118.0 million, check_c27k15f2g2aX1v1020mt.log) and (15, 10, 0, 0, 0, 1, 1) at
3.100562 (225.5 million, check_c27k15f2g2aX1v1011mt.log, after a first run
failed on A2B6F7 and was resumed with 3e-7 more room there) and (15, 9, 0, 0,
0, 2, 1) at 3.103026 (368.0 million, check_c27k15f2g2aX1v921mt2.log, mixed
again with the peaks of a first mix fed back, after a first run failed on
A1A1A1 and A2B6F7 and was resumed with more room there), and these four are
now proved (above); (15, 9, 0, 0, 0, 1, 2) alone stays at 3.106600 and is
open, and (15, 8, 0, 0, 0, 2, 2) alone gives 3.103115, 3.105325 with its
dense-search peaks and 3.103269 re-mixed with them, but 3.103774 with the
peaks of the re-mix even with no margins, and a third re-mix with both sets of
peaks stalls at 3.103548, so it is open. The 9 vectors with 16 close, n_7 = 3
and none in (2.1, 2.35] mix to 3.114348 with the bin (2.05, 2.1] a type of its
own (floating point); with that bin typed A the 8 of them with at most 7 in
(2.0161, 2.05] are proved (above). With 16 close and n_7 <= 2 the parts follow
those at 15 (floating point): with none in (2.1, 2.35] the 23 with n_7 <= 1
mix to 3.101015 and the 8 with n_7 = 2 other than (16, 0, 9, 0, 0, 0, 2) and
(16, 9, 0, 0, 0, 0, 2) to 3.101333; the 23 and the 8 are proved (above), the
two left out giving 3.104116 and 3.103832; mixed on its own, (16, 0, 9, 0, 0,
0, 2) gives 3.082699 and is proved (above), while (16, 9, 0, 0, 0, 0, 2) alone
stays at 3.103729 and is open; with exactly one centre in (2.1, 2.2] and none
in (2.2, 2.35] the 60 give 3.107876, only the 24 with at most three in (2.05,
2.1] above the level, and the 36 others mix to 3.092122 and are now proved
(above); a mix for the 18 with one to three in (2.05, 2.1] gives 3.103291, at
most 3.099892 on the 12 with two or three (case27_k16f2Hm2.json), 3.101345
with their dense-search peaks (peak2_c27k16f2Hm2mt.log), margins
1e-6/1e-7/1e-7 and extra room, and are now proved (above); the 6 with one
there mix to 3.102636 (case27_k16f2Hn1.json) but 3.104253 with their
dense-search peaks (peak2_c27k16f2Hn1mt.log) and the margins, and re-mixed
with those peaks to 3.103299 (mix27_k16f2Hn1mt2.npz), about 3.1038 with the
peaks of the new mix on the 2 with n_7 = 2, which are open (mixed one at a
time 3.104478 and 3.1078), and 3.103244 on the 4 with n_7 <= 1
(case27_k16f2Hn1f1.json), whose exact check runs (check_c27k16f2Hn1f1mt.log);
with none in (2.1, 2.2] and one or more in (2.2, 2.35] the 103 give 3.111974,
and the 67 with at least three in (2.05, 2.1] mix to 3.101372, 3.102149 with
their dense-search peaks (peak2_c27k16f2G2an3mt.log), margins 1e-6/1e-7/1e-7
and extra room, and are now proved (above); a mix for the 12 with exactly two
in (2.05, 2.1] gives 3.103355, 3.100707 on the 9 with at least two in (2.2,
2.35] (case27_k16f2G2ae2b.json), 3.102292 with their dense-search peaks
(peak2_c27k16f2G2ae2bmt.log), margins and extra room, and their exact check
passes at U <= 3.102282 at (16, 5, 2, 0, 0, 2, 2) (3 612 007 pair boxes, 1604
bin intervals, 641 317 975 triple boxes, PASS at 13:19 on 2026-10-26;
check_c27k16f2G2ae2bmt.log, multi_cap/runs/combo27_k16f2G2ae2b_check.log,
certificate multi_cap/radial_certificates/combo27_k16f2G2ae2b_d10.npz; not yet
in the paper, which takes these 9 together with the 3 below), while the 3 with
one in (2.2, 2.35] stay at 3.1032 together; of these, (16, 8, 2, 0, 0, 1, 0)
mixed on its own passes its exact check at 3.099725 (230.4 million triple
boxes, check_c27k16f2v8210mt.log; not yet in the paper), and the other two are
mixed one at a time ((16, 7, 2, 0, 0, 1, 1) to 3.098756, 3.099854 with its
dense-search peaks, exact check running, check_c27k16f2v7211mt.log; (16, 6, 2,
0, 0, 1, 2) to 3.100933, 3.102064 with its dense-search peaks, exact check
running, check_c27k16f2v6212mt.log); of the 24 with at most one in (2.05,
2.1], a mix for the 12 with three or more in (2.2, 2.35] stalls at 3.104560,
under which the 6 with four or more give 3.100833, 3.101610 with their
dense-search peaks (exact check running, check_c27k16f2G2al4mt.log), and the 6
with exactly three are mixed one at a time; of the 12 with one or two there,
the 6 with one in (2.05, 2.1] are mixed one at a time ((16, 9, 1, 0, 0, 1, 0)
to 3.100745 and (16, 8, 1, 0, 0, 2, 0) to 3.098354, exact checks running,
check_c27k16f2v9110mt.log and check_c27k16f2v8120mt.log; (16, 8, 1, 0, 0, 1,
1) to 3.101428, dense search running), and the 6 with none there are open
((16, 10, 0, 0, 0, 1, 0) alone 3.103160); the 1372 vectors with a centre in
(2.1, 2.2] and at least two in (2.1, 2.35] mix to 3.107964 with 225 above the
level; under that mix the 1109 with at least two in (2.1, 2.2] and at least
three in (2.1, 2.35] give 3.102171 (case27_k16f2Ex.json) but 3.104914 with
their dense-search peaks (peak2_c27k16f2Exmt.log), 38 above the level, and a
re-mix of all 1109 ran out of memory, so they are split: under that mix with
those peaks the 887 with at least three in (2.1, 2.2] give 3.101331 and the
150 with exactly two there and at least two in (2.2, 2.35] 3.102821, and both
exact checks run (check_c27k16f2Ex3mt.log, check_c27k16f2Ex2bmt.log), while
the 72 with exactly two there and one in (2.2, 2.35], which hold all 38 above
the level, are re-mixed with the peaks; the 81 with exactly two in (2.1, 2.2]
and none in (2.2, 2.35] mix on their own to 3.101190 (case27_k16f2E2.json),
3.102381 with their dense-search peaks, exact check running
(check_c27k16f2E2mt.log); the 144 with exactly one in (2.1, 2.2] and one to
three in (2.2, 2.35] mix on their own to 3.106891, and with that mix's
dense-search peaks the 108 with at least two in (2.05, 2.1] give 3.099618
(exact check running, check_c27k16f2E1cmt.log), the 18 with one there
3.103662, 2 above the level (re-mixed with the peaks), and the 18 with none
there up to 3.10741 (open). The mixes give 3.1197 at 16, 3.1312 at 18, 3.1446
at 21 and 3.1634 at 24 (floating point). The rounds for exactly 11 to 14 reach
3.08563, 3.09128, 3.09689 and 3.10236 on their own samples, about 0.0055 more
per close centre; a round aimed at exactly 15 (c27k15_d10.log) reaches 3.10766
on its samples, above the level, which suggests without proving that degree 10
in this form does not settle 15 or more close. This is stated in Remark 7.41.
