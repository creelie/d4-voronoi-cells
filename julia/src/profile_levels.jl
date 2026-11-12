# profile_levels.jl -- port of anc/profile_levels.py: the profile-level bound on U(Y) and the
# density bound it gives.  Every constant is an exact rational rounded UP to a multiple of
# 1e-13 and kept as its integer numerator (as in the Python); the recursion is integer
# arithmetic on these numerators (all below 1e16, far inside Int64; products and sums are
# overflow-checked).  The density comparison is in Arb balls at 200 bits.
#
# For a bin-count vector n = (n_0..n_6) (lower ends 2, 2.0161, 2.05, 2.1, 2.15, 2.2, 2.444;
# types AABBBBF):
#   G(n) = min( sum_b n_b S(e_b),  t_M/2 + m_0 n_0 + m_1 (n_1+..+n_6) for the five two-point
#               certificates,  F26(n) if n_6 <= 2,  min_b G(n - e_b) + S(e_b) ),
# and G(n) = sum_b n_b S(e_b) as soon as that sum is at most the cut 3.2.
module ProfileLevels
using ..Common
using ..Common: QQ
using ..Exact: fr, solve_exact, ldl_positive, arb_lo, arb_hi
using ..Radial: S_arb
using Arblib, JSON, NPZ, Printf
using Base.Checked: checked_add, checked_mul

const SC = big(10)^13
"numerator of the least multiple of 1e-13 that is >= x"
up(x) = Int64(cld(numerator(QQ(x))*SC, denominator(QQ(x))))
"an integer k with x < k/1e13 certified (same start as the Python: int(float(mid)*1e13) + 2)"
function arb_up2(x::Arb)
    k = Int64(trunc(Float64(Arblib.midref(x))*1.0e13)) + 2
    while !(x < Arb(QQ(k, SC); prec=200)); k += 1; end
    k
end
const EDGES = QQ[2, 20161//10000, 41//20, 21//10, 43//20, 11//5, 611//250]
const TYPES = "AABBBBF"
const POINT0 = ("0.10730785", "0.10730646", "-0.05347062", "-0.05363574", "-0.05425642", "-0.05347062", "-0.03946852")
const PAIR0 = Dict(zip(("AA", "AB", "AF", "BB", "BF", "FF"),
    fr.(("-.02327176", "-.01036330", "-.01503736", "-.001713338", "-.004512281", "-.003681549")) .+ QQ(1, 10^7)))
const TRIPLE0 = Dict(zip(("AAA", "AAB", "AAF", "ABB", "ABF", "AFF", "BBB", "BBF", "BFF"),
    fr.((".001497925", ".0009598990", ".001322105", ".0005966748", ".001043879", ".001145076", ".0004094797", ".0006564131", ".001043982")) .+ QQ(1, 10^7)))
const T26 = fr("8.282667610")

"(t/2, m_0, m_1) of a two-point certificate, rounded up; t = z^T A_0^{-1} z rounded up to 2^-48"
function two_point(path)
    leaf = JSON.parsefile(path)["leaves"][1]
    e = fr.(leaf["edges"])
    (e[1] == 2 && e[2] == 20161//10000 && length(leaf["m"]) == 2) || error("unexpected bins in $path")
    a0 = leaf["A"][1]; n = length(a0)
    A0 = QQ[fr(a0[i][j]) for i in 1:n, j in 1:n]; z = QQ[fr(x) for x in leaf["z"]]
    ldl_positive(A0) || error("A_0 is not positive definite in $path")
    tq = sum(z .* solve_exact(A0, z))
    t = QQ(cld(numerator(tq)*big(2)^48, denominator(tq)), big(2)^48)
    m = fr.(leaf["m"])
    (up(t/2), up(m[1]), up(m[2]))
end

"thresholds proved in a full PASS log of combo_case_check.py: printed value + 1e-6 of its size"
function thresholds_from_log(s)
    pr = Dict(String(m.captures[1]) => (v = fr(m.captures[2]); v + abs(v)/10^6) for m in eachmatch(r"pair (\w\w): float largest[^\n]*?threshold (\S+)", s))
    tr = Dict(String(m.captures[1]) => (v = fr(m.captures[2]); v + abs(v)/10^6) for m in eachmatch(r"triple (\w\w\w): float largest[^\n]*?threshold (\S+?);", s))
    pr, tr
end
full_pass(s) = length(collect(eachmatch(r"PASS", s))) == 1 && occursin("case26_all.json", s) && occursin("combo26_level_d8.npz", s)

struct Consts
    S::Vector{Int64}; TP::Vector{NTuple{3,Int64}}; POINT::Vector{Int64}
    PAIR::Dict{String,QQ}; TRIPLE::Dict{String,QQ}; ang::Dict{NTuple{3,Int},Int64}
end
function angular(C::Consts, a, b, f)
    get!(C.ang, (a, b, f)) do
        cnt = Dict('A' => a, 'B' => b, 'F' => f); tot = T26/2
        for (w, c) in Iterators.flatten((C.PAIR, C.TRIPLE))
            mult = big(1)
            for l in unique(w); mult *= binomial(big(cnt[l]), count(==(l), w)); end
            tot += mult*c
        end
        up(tot)
    end
end
dot7(n, v) = (s = Int64(0); for i in 1:7; s = checked_add(s, checked_mul(Int64(n[i]), v[i])); end; s)
function base(C::Consts, n)
    v = dot7(n, C.S)
    rest = Int64(n[2] + n[3] + n[4] + n[5] + n[6] + n[7])
    for (t2, m0, m1) in C.TP
        v = min(v, checked_add(t2, checked_add(checked_mul(m0, Int64(n[1])), checked_mul(m1, rest))))
    end
    if n[7] <= 2
        v = min(v, checked_add(angular(C, n[1] + n[2], n[3] + n[4] + n[5] + n[6], n[7]), dot7(n, C.POINT)))
    end
    v
end

"all n in N^7 with |n| = M"
function compositions(M)
    out = NTuple{7,Int8}[]
    for a in 0:M, b in 0:M-a, c in 0:M-a-b, d in 0:M-a-b-c, e in 0:M-a-b-c-d, f in 0:M-a-b-c-d-e
        push!(out, (a, b, c, d, e, f, M - a - b - c - d - e - f))
    end
    out
end
dec1(n, b) = ntuple(i -> i == b ? Int8(n[i] - 1) : n[i], 7)

"levels M = 0..Mmax by dynamic programming on |n| (G(n) needs only G(n - e_b)); returns per-M maxima with and without n_0 <= 24"
function levels(C::Consts, Mmax, cut)
    prev = Dict{NTuple{7,Int8},Int64}((Int8(0), Int8(0), Int8(0), Int8(0), Int8(0), Int8(0), Int8(0)) => 0)
    res = Dict{Int,Any}()
    for M in 1:Mmax
        cur = Dict{NTuple{7,Int8},Int64}(); sizehint!(cur, binomial(M + 6, 6))
        best = (Int64(-1), nothing); best24 = (Int64(-1), nothing); n24 = 0
        for n in compositions(M)
            cs = dot7(n, C.S)
            g = cs
            if cs > cut
                g = base(C, n)
                for b in 1:7
                    n[b] > 0 && (g = min(g, checked_add(prev[dec1(n, b)], C.S[b])))
                end
            end
            cur[n] = g
            g > best[1] && (best = (g, n))
            if n[1] <= 24
                n24 += 1
                g > best24[1] && (best24 = (g, n))
            end
        end
        res[M] = (best24, best, n24)
        prev = cur
    end
    res
end
f7(k) = @sprintf("%.7f", k/1.0e13)

# expected values of logs/density_levels.out and density_levels_tight.out (7 digits, as printed there)
const EXPECT = Dict(
  :base  => (["3.2395754", "3.2888648", "3.2876756", "3.2859021", "3.2586970", "3.2079814", "3.1999976"],
             ["3.2395754", "3.2888648", "3.2876756", "3.2859021", "3.2644542", "3.2102926", "3.1999976"],
             "3.2888648", (22, 0, 3, 0, 0, 0, 0), "0.631496", ["3.23958", "3.28887", "3.28768", "3.28591", "3.26446", "3.21030", "3.2"]),
  :tight => (["3.2395754", "3.2826739", "3.2809570", "3.2809595", "3.2586970", "3.2079814", "3.1999976"],
             ["3.2395754", "3.2826739", "3.2809570", "3.2809595", "3.2625398", "3.2102926", "3.1999976"],
             "3.2826739", (22, 0, 3, 0, 0, 0, 0), "0.630996", nothing))
const PROFILES = [593775, 736280, 906185, 1107540, 1344820, 1622950, 1947330]

"""
run(io; datadir, tightlog=nothing, log26=nothing)

`tightlog`: a full PASS log of combo_case_check.py on combo26_level_d8.npz / case26_all.json whose (smaller)
thresholds replace the constants PAIR/TRIPLE; the result then rests on that run.
`log26`: a full PASS log against which the built-in constants POINT/PAIR/TRIPLE/T26 are cross-checked.
"""
function run(io::IO=stdout; datadir=Common.datadir(), tightlog=nothing, log26=nothing)
    kind = tightlog === nothing ? :base : :tight
    R = Report("profile levels and density bound (anc/profile_levels.py)" * (kind == :tight ? ", TIGHT constants" : ""))
    S = [arb_up2(S_arb(e)) for e in EDGES]
    println(io, "S(e_b) upper bounds: ", [@sprintf("%.10f", s/1.0e13) for s in S])
    check!(R, io, "S(e_b) rounded up (certified S(e_b) < k/1e13) agree with the Python log",
           [@sprintf("%.10f", s/1.0e13) for s in S] == ["0.1389635009", "0.1273288867", "0.1046916513", "0.0758130327", "0.0521298076", "0.0333964703", "0.0000025325"] &&
           all(S_arb(e) < Arb(QQ(s, SC); prec=200) for (e, s) in zip(EDGES, S)))
    TP = [two_point(joinpath(datadir, "radial_certificates", "density_$M.json")) for M in (25, 27, 28, 29, 30)]
    tpexp = [("0.9095329909", "0.0970851006", "0.0878140183"), ("4.5238173034", "-0.0443547088", "-0.0501172844"), ("4.5338315966", "-0.0447180657", "-0.0504752608"),
             ("4.5482527029", "-0.0452186014", "-0.0510049757"), ("4.8498109771", "-0.0565351159", "-0.0565351159")]
    println(io, "two-point certificates (t/2, m0, m1): ", [map(x -> @sprintf("%.10f", x/1.0e13), tp) for tp in TP])
    check!(R, io, "two-point constants (t = z^T A_0^{-1} z exactly, rounded up; A_0 positive definite by exact LDL^T) agree with the Python log",
           [map(x -> @sprintf("%.10f", x/1.0e13), tp) for tp in TP] == tpexp)
    POINT = [up(fr(x)) for x in POINT0]; PAIR = copy(PAIR0); TRIPLE = copy(TRIPLE0)
    if tightlog !== nothing
        s = read(tightlog, String)
        pr, tr = thresholds_from_log(s)
        ok = full_pass(s) && Set(keys(pr)) == Set(keys(PAIR)) && Set(keys(tr)) == Set(keys(TRIPLE)) &&
             all(occursin("PAIR3_$k <= Pi + c2", s) for k in keys(PAIR)) && all(occursin("TRIPLE3_$k <= c3", s) for k in keys(TRIPLE)) &&
             all(pr[k] <= PAIR[k] for k in keys(PAIR)) && all(tr[k] <= TRIPLE[k] for k in keys(TRIPLE))
        check!(R, io, "the tight log is one full PASS on combo26_level_d8.npz / case26_all.json, names all 6 pair and 9 triple kinds, and its thresholds are at most the built-in ones", ok;
               got=join(["$k: " * @sprintf("%.9f", Float64(pr[k])) for k in sort(collect(keys(pr)))], ", "))
        ok || return R
        PAIR, TRIPLE = pr, tr
    end
    if log26 !== nothing
        s = read(log26, String)
        ok = full_pass(s)
        t = fr(match(r"\(t = ([0-9.]+)", s).captures[1]); ok &= t + QQ(1, 10^9) <= T26
        pt = Dict(String(m.captures[1]) => fr(m.captures[2]) for m in eachmatch(r"'(\w)': '([-0-9.]+)'", match(r"point terms: (\{.*?\})", s).captures[1]))
        m = npzread(joinpath(datadir, "combo", "combo26_level_d8.npz"), ["m"])["m"]
        for b in 1:7
            ok &= up(pt[string(TYPES[b])] + QQ(1, 10^6) + QQ(Float64(m[b])) + QQ(2, 10^6) + QQ(1, 10^9)) <= POINT[b] + up(QQ(1, 10^8))
        end
        pr, tr = thresholds_from_log(s)
        ok &= length(pr) == 6 && length(tr) == 9 && all(pr[k] <= PAIR[k] for k in keys(pr)) && all(tr[k] <= TRIPLE[k] for k in keys(tr))
        ok &= all(occursin("PAIR3_$k <= Pi + c2", s) for k in keys(PAIR)) && all(occursin("TRIPLE3_$k <= c3", s) for k in keys(TRIPLE))
        check!(R, io, "constants POINT/PAIR/TRIPLE/T26 of the typed M = 26 certificate are implied by the log $(basename(log26)) and the certificate file", ok; got="t = $(Float64(t)) + 1e-9 <= $(Float64(T26))")
    end
    C = Consts(S, TP, POINT, PAIR, TRIPLE, Dict{NTuple{3,Int},Int64}())
    cut = up(fr("3.2"))
    res = levels(C, 30, cut)
    e24, eall, etop, earg, edens, ecoarse = EXPECT[kind]
    worst = Int64(-1); warg = nothing
    for (i, M) in enumerate(24:30)
        b24, ball, n24 = res[M]
        worst < ball[1] && ((worst, warg) = ball)
        check!(R, io, "M = $M: max G over |n| = $M with n_0 <= 24 ($n24 profiles)", f7(b24[1]) == e24[i] && n24 == PROFILES[i];
               got="$(f7(b24[1])) at $(Int.(b24[2]))", expected=e24[i] * ", $(PROFILES[i]) profiles")
        okc = ecoarse === nothing || (M < 30 ? QQ(ball[1], SC) <= fr(ecoarse[i]) : QQ(ball[1], SC) < fr("3.2"))
        check!(R, io, "M = $M: max G over |n| = $M, no restriction on n_0", f7(ball[1]) == eall[i] && okc;
               got="$(f7(ball[1])) at $(Int.(ball[2]))", expected=eall[i] * (ecoarse === nothing ? "" : (M < 30 ? "; at most $(ecoarse[i])" : "; below 3.2")))
    end
    check!(R, io, "max G over 24 <= |n| <= 30", f7(worst) == etop && Tuple(Int.(warg)) == earg && QQ(worst, SC) <= fr(etop);
           got="$(worst)/1e13 at $(Int.(warg))", expected="$etop at $earg")
    L = QQ(worst + 1, SC); p = 200; PI = Arb(π; prec=p)
    v = 9PI^2/8 - Arb(L; prec=p) + Arb(QQ(1, 28600); prec=p)
    dens = (PI^2/2)/v
    println(io, "level L* = ", L, "; 9 pi^2/8 - L* + 1/28600 = ", v); println(io, "density bound = ", dens)
    check!(R, io, "9 pi^2/8 - $etop + 1/28600 > 0 and (pi^2/2)/(9 pi^2/8 - $etop + 1/28600) < $edens (Arb, 200 bits, with L = $etop exactly)",
           (vv = 9PI^2/8 - Arb(fr(etop); prec=p) + Arb(QQ(1, 28600); prec=p); Arblib.is_positive(vv) != 0 && (PI^2/2)/vv < Arb(fr(edens); prec=p));
           got=@sprintf("%.12f", arb_hi((PI^2/2)/(9PI^2/8 - Arb(fr(etop); prec=p) + Arb(QQ(1, 28600); prec=p)))))
    check!(R, io, "the same with L* = (max G + 1e-13) as in the Python", Arblib.is_positive(v) != 0 && dens < Arb(fr(edens); prec=p);
           got=@sprintf("density < %.15f", arb_hi(dens)), expected=kind == :base ? "0.631495029884106" : "0.630995132110131")
    kind == :base && check!(R, io, "density < 0.630996 is NOT given by the built-in constants (it needs the tight thresholds)", dens > Arb(fr("0.630996"); prec=p))
    R
end
end
