# gnbhd.jl -- the certified evaluations behind thm:G-nbhd (appendix app:G-nbhd, table
# tab:G-nbhd-bins) and the exact arithmetic of thm:hexagon-rigidity / cor:hexagon-close.
# Port of anc/gnbhd.py and anc/gnbhd_bins.py (python-flint, 200 bits) to Arblib at 200 bits: every
# constant is a ball containing its exact value, grid values are exact rationals, and a
# comparison counts only when Arb certifies it for the whole ball.
module GNbhd
using ..Common
using ..Common: QQ
using ..Exact: fr, arb_lo, arb_hi
using Arblib, Printf

const P = 200
Q(a, b=1) = Arb(QQ(a, b); prec=P)
Q(x::Rational) = Arb(QQ(x); prec=P)
mid(x::Arb) = Float64(Arblib.midref(x))
pos(x::Arb) = Arblib.is_positive(x) != 0
rt(n) = sqrt(Arb(n; prec=P))
PI() = Arb(π; prec=P)

"cap volume of B(R), R^2 = 3/2, beyond the hyperplane at height d/2 (closed form)"
function S(d::Arb)
    R2 = Q(3, 2); R = sqrt(R2); h = d/2
    4PI()/3*(3R2*R2/8*(PI()/2 - asin(h/R)) - h/8*(5R2 - 2h*h)*sqrt(R2 - h*h))
end
function G(h::Arb, x::Arb)
    R2 = Q(3, 2); t = tan(x)
    R2*R2*x - 2R2*h*h*t + h^4*(t + t^3/3)
end
"Pi(1,1,1/2) by lem:pair-closed; returns (value, side conditions hold)"
function omega_half()
    R = sqrt(Q(3, 2)); g = PI()/3; al = acos(1/R)
    L = g - al; U = al; pc = g/2
    ok = pos(U - L) && pos(pc - L) && pos(U - pc)
    one = Q(1)
    PI()/4*(G(one, pc - g) - G(one, L - g) + G(one, U) - G(one, pc)), ok
end
fshare(x::Arb) = (q = sqrt(2*(1 - x)/(1 + x)); (1 - q)^2*(2 + q)/4)
s_of(delta::Arb) = S(Q(2)) - S(2 + delta)

struct Base0; W12::Arb; Troot::Arb; Lip::Arb; ok::Bool; end
function base0()
    W12, ok = omega_half()
    Base0(W12, 9PI()^2/8 - 24S(Q(2)) + 96W12, PI()/(72rt(3)), ok)
end
"lem:T-near-roots: (pen1, bound1, bound2, side conditions hold)"
function lemmaT(B::Base0, d0::Arb)
    e0 = d0
    xbar = Q(1, 2) + 2*e0 + e0*e0
    ok = pos(Q(56, 100) - xbar) && pos(Q(1, 3) - 2*e0 - e0*e0) && pos(1 - 8fshare(Q(56, 100)))
    pen1 = B.Lip*(6d0*d0 + 6rt(2)*(1 + d0)*d0)
    pen2 = B.Lip*(6d0*d0 + 4rt(10)*(1 + d0)*d0)
    pen1, 8 - B.Troot + 8B.W12 + pen1, 8 - B.Troot + 16B.W12 + pen2, ok
end
"least n/step (n >= n0, n/step < 0.449) with f(n/step) certified positive"
function least(f, step, n0)
    for n in n0:(449*step ÷ 1000)
        pos(f(QQ(n, step))) && return QQ(n, step)
    end
    error("no grid point found")
end
pfun(x, y) = (x + y)^3*(Q(40, 3)*x + Q(8, 3)*y)

"thm:one-large-push: certified m with vol(V∩K) >= 8 + m S; returns (m, all side conditions hold, constants)"
function thmB(d0::Arb, eta1::Arb, eta2::Arb; Smin=nothing)
    s2 = rt(2); s3 = rt(3); s11 = rt(11); s24 = rt(24)
    e0 = d0
    r1 = s2/(1 - s2*e0); Rr = 1 + 3r1*e0
    c3 = (Q(1, 2) + e0)/(1 - e0*e0/2)
    Dmax = e0*sqrt(1 + Rr*Rr) + c3*(e0*e0/2 + Rr*e0)
    mu = 2Dmax/3
    k = 1/(s3/2 - (1 + c3)*e0)
    Lam = s3/24*(1 + 3mu)^4
    s1hat = (e0 + e0*e0/2)/(1 - e0*e0/2)
    kperp = sqrt((1 + s1hat)^2 + 1); Dperp = e0*kperp
    gbar = 1/(1 - e0*e0/2)^2
    sR = (Rr*e0 + e0*e0/2)/(1 - e0*e0/2)
    cond = Arb[1/s2 - e0, -(-(1 - sR)/2 + (1 + sR)*e0 + Rr/2 + Rr*e0 - 1), (1 - sR)*(1 - e0) - Rr*e0 + 1,
               s3/2 - (1 + c3)*e0, 1 - Dperp]
    Mbar = (k/2)*(16Lam*(1 + c3) + 2*(1 + 3mu)^3*(1 + c3)*e0) +
           (e0*k*k/2)*(16Lam*(1 + c3*c3) + 2s3*(1 + 3mu)^2*e0*e0*(1 + c3*c3)) +
           s24*kperp^3*e0*e0
    vbar = k*(4*(1 + 3mu)^3*e0*(1 + c3) + 2s3*(1 + 3mu)^2*(1 + c3)*e0*e0)
    CQ = -(Q(2, 3) - vbar/4 - Dperp^3/2 - gbar*Mbar)
    sl = (1 - 2*e0 - e0*e0)/(1 + 2*e0)
    push!(cond, sl - eta1)
    x = eta1/sl
    Psi = 1 - 3x/2 + x*x - x^3/4
    r = s2*(1 + eta1)/(1 - s2*e0)
    beta = s2*r + 2
    lam = (1 + eta1)^2/(1 - s2*e0)^2
    b = 1/(1 + eta1)
    push!(cond, (4 - 3lam)*s2*b*(1 - s2*e0) - 4*e0)
    push!(cond, eta1 - eta2 + Q(1, big(10)^30))
    al0 = s2*(eta1 + eta2); b0 = beta*e0
    alII = 2s2*eta2
    cII = s2*(20alII*alII + 28alII*b0 + Q(52, 3)*b0*b0)
    cI = (8al0^3 + 12al0*al0*b0)/eta1 + 16s2*b0*b0 + Q(4, 3)*s2*(1 - eta2/eta1)*(3al0*al0 + 3al0*b0 + b0*b0)
    cmix = Arb(prec=P); Arblib.max!(cmix, cI, cII, P)
    hullx = 6beta*cmix
    hullE = 12beta^4*e0*e0
    Gr = Q(4, 3)*Psi - 6pfun(eta1, eta2)/eta1 - 4Dperp^3
    Cx = k*2s2*(1 + 3mu)^3*(1 + c3) + hullx
    CE = CQ + k*2s3*(1 + 3mu)^2*(1 + c3)*eta1 + hullE
    c1 = (s11 + rt(27))/4
    c2 = 6*(s11 + rt(75))/8 + Q(1, 2)
    push!(cond, 1 - c2*d0)
    rho = c1*(1 + eta1)/(1 - c2*d0)
    # a certainly non-positive C_E is a favourable term and is dropped; otherwise the whole ball is used
    CEuse = Arblib.is_nonpositive(CE) != 0 ? Q(0) : CE
    # E^2 <= d0 rho S always; if S >= Smin, also E^2 <= d0^2 <= (d0^2/Smin) S
    fac = rho*d0
    if Smin !== nothing && Arblib.is_negative(d0*d0/Smin - fac) != 0
        fac = d0*d0/Smin
    end
    m = (Gr - d0*Cx)/2 - fac*CEuse
    m, all(pos, cond), (; r1, R1=Rr, c3, Dmax, mu, k, kperp, Dperp, Lam, CQ, sl, Gr, Cx, CE, rho, beta)
end

# tab:G-nbhd-bins of the paper: (x, x', y, m) in thousandths / ten-thousandths
const TABLE = [(97,101,97,1765),(101,105,93,1913),(105,109,90,2008),(109,113,86,2147),(113,117,82,2281),
    (117,121,79,2274),(121,125,75,2290),(125,129,72,2276),(129,133,69,2261),(133,137,66,2244),
    (137,141,63,2226),(141,145,60,2207),(145,149,57,2187),(149,153,54,2166),(153,157,52,2117),
    (157,161,49,2094),(161,165,46,2070),(165,169,44,2020),(169,173,42,1968),(173,177,39,1942),
    (177,181,37,1889),(181,185,35,1836),(185,189,33,1782),(189,193,30,1754),(193,197,28,1699),
    (197,201,26,1643),(201,205,25,1559),(205,209,23,1501),(209,212,21,1472)]

"quoted decimal `s` agrees with the ball x: |x - s| is certainly below one unit of the last digit"
function agrees(x::Arb, s::AbstractString)
    nd = occursin('.', s) ? length(split(s, '.')[2]) : 0
    pos(Q(1, big(10)^nd) - abs(x - Arb(fr(s); prec=P)))
end

function run(io::IO=stdout; N=48, step=250)
    R = Report("thm:G-nbhd: lem:T-near-roots and thm:one-large-push at e0 = 1/$N (Arb, $P bits)")
    B = base0(); d0 = Q(1, N)
    println(io, "S(2)       = ", S(Q(2))); println(io, "Pi(1,1,1/2)= ", B.W12)
    println(io, "T_root     = ", B.Troot); println(io, "L          = ", B.Lip)
    check!(R, io, "lem:pair-closed side conditions for Pi(1,1,1/2)", B.ok)
    for (nm, v, s) in [("S(2)", S(Q(2)), "0.1389635009"), ("Pi_0 = Pi(1,1,1/2)", B.W12, "0.0014454086"),
                       ("T_root", B.Troot, "7.9069401562"), ("L = pi/(72 sqrt 3)", B.Lip, "0.0251916578")]
        check!(R, io, "$nm agrees with the quoted value", agrees(v, s); got=@sprintf("%.12f", mid(v)), expected=s * "...")
    end
    pen, b1, b2, okT = lemmaT(B, d0)
    check!(R, io, "lem:T-near-roots side conditions (xbar < 0.56, 1/3 > 2 e0 + e0^2, 8 f(0.56) < 1)", okT)
    check!(R, io, "s(d_a) <= 0.1092348", b1 < Arb(fr("0.1092348"); prec=P) && b1 > Arb(fr("0.1092347"); prec=P); got=@sprintf("%.9f", mid(b1)), expected="0.1092348 (rounded up)")
    check!(R, io, "s(d_a)+s(d_b) <= 0.1230289", b2 < Arb(fr("0.1230289"); prec=P) && b2 > Arb(fr("0.1230288"); prec=P); got=@sprintf("%.9f", mid(b2)), expected="0.1230289 (rounded up)")
    check!(R, io, "s(2.212) > bound 1, so every delta_i < 0.212", pos(s_of(Q(212, 1000)) - b1); got=@sprintf("%.6f > %.6f", mid(s_of(Q(212, 1000))), mid(b1)), expected="s(2.212) = 0.109359")
    check!(R, io, "2 s(2.097) > bound 2, so at most one delta_i >= 0.097", pos(2s_of(Q(97, 1000)) - b2); got=@sprintf("%.6f > %.6f", mid(2s_of(Q(97, 1000))), mid(b2)), expected="2 s(2.097) = 0.123133")
    d1 = least(q -> s_of(Q(q)) - b1, 1000, 1)
    d2 = least(q -> 2s_of(Q(q)) - b2, 1000, 1)
    check!(R, io, "least grid thresholds (step 1/1000)", N != 48 || (d1 == 212//1000 && d2 == 97//1000); got="delta_max < $(Float64(d1)), second < $(Float64(d2))", expected="0.212, 0.097")
    # first case: every delta_i < d2
    m, ok, c = thmB(d0, Q(d2/2), Q(d2/2))
    println(io, "constants of prop:contact-cell-L2 / prop:hull-vertex and of the first case:")
    quoted = [(:r1, "1.45715"), (:R1, "1.09107"), (:c3, "0.52095"), (:Dmax, "0.04279"), (:mu, "0.02853"), (:k, "1.19855"),
              (:kperp, "1.42918"), (:Dperp, "0.02977"), (:Lam, "0.10023"), (:CQ, "0.93043"), (:sl, "0.91958"),
              (:Gr, "1.14379"), (:Cx, "36.2355"), (:CE, "2.85216"), (:rho, "2.78077"), (:beta, "4.16066")]
    if N == 48
        for (k, s) in quoted
            check!(R, io, "constant $k agrees with the quoted value", agrees(getfield(c, k), s); got=@sprintf("%.7f", mid(getfield(c, k))), expected=s)
        end
    end
    check!(R, io, "first case (all delta_i < $(Float64(d2)), eta1 = eta2 = $(Float64(d2/2))): side conditions hold and m > 0", ok && pos(m);
           got=@sprintf("m >= %.7f", arb_lo(m)), expected="m = 0.029208...")
    N == 48 && check!(R, io, "first case margin agrees with 0.029208", agrees(m, "0.029208") && pos(m - Q(292, 10000)); got=@sprintf("%.8f", mid(m)), expected="0.029208...")
    # bins on the largest push-out
    x = d2; nb = 0; worst = arb_lo(m); above = Int[]; worstbin = Inf
    while x < d1
        x2 = min(x + 1//step, d1)
        y = min(least(q -> s_of(Q(x)) + s_of(Q(q)) - b2, 1000, 0), x2)
        m, ok, _ = thmB(d0, Q(x2/2), Q(y/2); Smin=Q(x))
        nb += 1; worst = min(worst, arb_lo(m)); worstbin = min(worstbin, arb_lo(m))
        exp = ""; match = true
        if N == 48 && step == 250
            if nb <= length(TABLE)
                t = TABLE[nb]
                match = x == t[1]//1000 && x2 == t[2]//1000 && y == t[3]//1000 &&
                        pos(m - Q(t[4], 10000)) && pos(Q(1, 10000) - (m - Q(t[4], 10000)))   # the table prints m rounded down to 4 decimals
                pos(m - Q(t[4], 10000)) || push!(above, nb)
                exp = @sprintf("[%.3f,%.3f], y = %.3f, m >= %.4f", t[1]/1000, t[2]/1000, t[3]/1000, t[4]/10000)
            else
                match = false; exp = "no such row in the table"
            end
        end
        check!(R, io, @sprintf("bin %2d: delta_max in [%.3f,%.3f], others < %.3f: side conditions hold, m > 0%s", nb, Float64(x), Float64(x2), Float64(y), exp == "" ? "" : ", row matches the table"),
               ok && pos(m) && match; got=@sprintf("m >= %.6f", arb_lo(m)), expected=exp)
        x = x2
    end
    check!(R, io, "the bins cover [$(Float64(d2)), $(Float64(d1))]", x == d1 && (N != 48 || step != 250 || nb == 29); got="$nb bins", expected="29")
    check!(R, io, "m >= 0.147 on every bin (as used in the proof)", worstbin >= 0.147; got=@sprintf("%.6f", worstbin))
    check!(R, io, "every printed m is a certified lower bound", isempty(above); got="rows above the certified m: $(above)", expected="none")
    check!(R, io, "worst margin over all cases is positive", worst > 0; got=@sprintf("%.6f", worst), expected="0.0292 (first case); >= 0.147 on the bins")
    R
end

# ---------------------------------------------------------------------------------------------
# thm:hexagon-rigidity, exact rational arithmetic.
# F(eta,s)^2 = (sqrt(A) + 3 eta^2)^2 + eta^4/4 with A = 480 s^2 + 60 s eta^2 + 36 eta^4.
# F < eta  <=>  sqrt(A) + 3 eta^2 < sqrt(B), B = eta^2 - eta^4/4 >= 0
#          <=>  6 eta^2 sqrt(A) < C, C = B - A - 9 eta^4
#          <=>  C > 0 and 36 eta^4 A < C^2          (all quantities rational).
function F_below(eta::QQ, s::QQ)
    A = 480s^2 + 60s*eta^2 + 36eta^4
    Bq = eta^2 - eta^4/4
    C = Bq - A - 9eta^4
    (A >= 0 && Bq >= 0 && C > 0 && 36eta^4*A < C^2), C, C^2 - 36eta^4*A
end
function run_hexagon(io::IO=stdout)
    R = Report("thm:hexagon-rigidity and cor:hexagon-close: exact rational arithmetic")
    for (eta, s) in [(QQ(1, 48), QQ(17, 20000)), (QQ(7, 100), QQ(9, 5000))]
        ok, C, D = F_below(eta, s)
        check!(R, io, "F($eta, $s) < $eta (after squaring twice: C > 0 and 36 eta^4 A < C^2)", ok;
               got=@sprintf("C = %.6e, C^2 - 36 eta^4 A = %.6e", Float64(C), Float64(D)))
        Fv = sqrt((sqrt(Arb(480s^2 + 60s*eta^2 + 36eta^4; prec=P)) + 3Q(eta)^2)^2 + Q(eta)^4/4)
        check!(R, io, "the same in ball arithmetic", pos(Q(eta) - Fv); got=@sprintf("F = %.9f < eta = %.9f", mid(Fv), Float64(eta)))
    end
    # how sharp: the smallest failing multiples of 1/20000 and 1/5000 (information, a strict check each)
    okn, _, _ = F_below(QQ(1, 48), QQ(18, 20000))
    check!(R, io, "sharpness: F(1/48, 18/20000) < 1/48 fails (so 17/20000 is the last multiple of 1/20000)", !okn)
    okn, _, _ = F_below(QQ(7, 100), QQ(10, 5000))
    check!(R, io, "sharpness: F(7/100, 10/5000) < 7/100 fails", !okn)
    # s = 0: F(eta,0)^2 = (9 eta^2)^2 + eta^4/4 = (325/4) eta^4, so eta <= F forces eta = 0 or eta >= 2/sqrt(325)
    check!(R, io, "s = 0: 9^2 + 1/4 = 325/4, i.e. F(eta,0) = (sqrt(325)/2) eta^2", QQ(81) + 1//4 == 325//4)
    lo, hi = fr("0.110940"), fr("0.110941")
    check!(R, io, "2/sqrt(325) = 0.110940... (4/325 between 0.110940^2 and 0.110941^2)", lo^2 < 4//325 < hi^2;
           got=@sprintf("%.12f", mid(2/rt(325))), expected="0.110940...")
    check!(R, io, "2/sqrt(325) > 2/sqrt(397) (the radius of thm:local-uniqueness)", 325 < 397)
    for eta in (QQ(11, 100), QQ(1109, 10000))
        ok, _, _ = F_below(eta, QQ(0))
        check!(R, io, "s = 0: F($eta, 0) < $eta ($eta < 2/sqrt(325))", ok && eta^2 < 4//325)
    end
    # second inequality constants: (6+3)^2 + 1/4 = 325/4 again, and 2 sqrt(480)*sqrt(325)/2*... -> 156000 = 480*325
    check!(R, io, "480 * 325 = 156000 (discriminant of eta = sqrt(480) s + (sqrt(325)/2) eta^2)", 480*325 == 156000)
    check!(R, io, "eta_pm real for s <= 1/800: 4 * 156000 / 800^2 <= 1", 4*156000//800^2 <= 1; got=string(4*156000//800^2))
    # cor:hexagon-close: radius 2/sqrt(1-2s) at s = 17/20000
    s = QQ(17, 20000); r2 = 4/(1 - 2s)
    check!(R, io, "2/sqrt(1-2s) at s = 17/20000 lies in (2.0017021, 2.0017022), above 2.0017",
           fr("2.0017021")^2 < r2 < fr("2.0017022")^2 && r2 > fr("2.0017")^2; got=@sprintf("%.10f", mid(sqrt(Arb(r2; prec=P)))), expected="2.0017021...")
    check!(R, io, "1 - 2/rho^2 = 1/2 + s at rho = 2/sqrt(1-2s)", 1 - 2/r2 == 1//2 + s)
    R
end
end
