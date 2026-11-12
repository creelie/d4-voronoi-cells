# constants.jl -- ports of multi_cap/twentyfour_close_check.py (lem:twentyfour-close) and of
# the rigorous part of multi_cap/pushout_check.py (thm:G-pushout): constants in Arb ball
# arithmetic (200 and 128 bits) and exact rationals.  The qhull comparison at the end of
# pushout_check.py is floating point and not part of the proof; it is not ported.
module Constants
using ..Common
using ..Common: QQ
using ..Exact
using Arblib, Printf

mid(x::Arb) = Float64(Arblib.midref(x))
function Scap(d::Arb, prec)
    R2 = Arb(3; prec=prec)/2; R = sqrt(R2); PI = Arb(π; prec=prec)
    h = d/2
    anti = h/8*(5R2 - 2h*h)*sqrt(R2 - h*h) + 3R2*R2/8*asin(h/R)
    4PI/3*(3R2*R2/8*PI/2 - anti)
end
dec(s::AbstractString, prec) = Arb(fr(s); prec=prec)

function run_twentyfour_close(io::IO=stdout)
    R = Report("twentyfour_close_check (lem:twentyfour-close)")
    p = 200; PI = Arb(π; prec=p); S(s) = Scap(dec(s, p), p)
    tau = 9PI^2/8 - 8
    A = 9PI^2/8 - 23S("2")
    sD = S("2") - S("2.1648")
    B = dec("0.0929000002", p)
    s2444 = S("2.444")
    room = A + sD - 8
    @printf(io, "S(2) = %.10f  S(2.444) = %.10f  9pi^2/8 - 8 = %.10f\n", mid(S("2")), mid(s2444), mid(tau))
    check!(R, io, "(a) 22 S(2) + 11 S(2.444) < 9 pi^2/8 - 8", 22S("2") + 11s2444 < tau; got=@sprintf("%.10f < %.10f", mid(22S("2") + 11s2444), mid(tau)))
    check!(R, io, "s(D) < B = 0.0929000002, so min(s(D), B) = s(D)", sD < B; got=@sprintf("%.10f", mid(sD)))
    check!(R, io, "T(X) - 8 >= A* + s(D) - 8 > 0", room > 0; got=@sprintf("%.6e", mid(room)), expected="2.6492e-05 (room quoted in the paper)")
    for M in 25:33
        check!(R, io, "(b) M = $M: (M - 23) S(2.444) < A* + s(D) - 8", (M - 23)*s2444 < room; got=@sprintf("%.6e < %.6e", mid((M - 23)*s2444), mid(room)))
    end
    check!(R, io, "(b) at M = 33: T(Y) - 8 > 1.1e-6", room - 10s2444 > dec("1.1e-6", p); got=@sprintf("%.4e", mid(room - 10s2444)))
    check!(R, io, "(c) M = 34: 11 S(2.444) exceeds the room", 11s2444 > room; got=@sprintf("%.6e", mid(11s2444)))
    check!(R, io, "(d) 6 S(2.4) < 0.00368", 6S("2.4") < dec("0.00368", p); got=@sprintf("%.10f", mid(6S("2.4"))))
    R
end

function run_pushout(io::IO=stdout)
    R = Report("pushout_check (thm:G-pushout), rigorous part")
    p = 128; PI = Arb(π; prec=p); R2 = Arb(3; prec=p)/2; Rr = sqrt(R2)
    Aq(x) = Arb(QQ(x); prec=p)
    S(d) = Scap(Aq(d), p)
    G(h, x) = (tx = tan(x); R2^2*x - 2R2*h*h*tx + h^4*(tx + tx^3/3))
    function Pi(h1q, h2q, uq)        # lem:pair-closed, as a ball; the branch decisions must be definite
        h1, h2, ua = Aq(h1q), Aq(h2q), Aq(uq)
        g = acos(ua); a1, a2 = acos(h1/Rr), acos(h2/Rr)
        L = max(-a1, g - a2); U = min(a1, g + a2)
        L < U || error("branch not decided")
        pc = atan(h2 - h1*cos(g), h1*sin(g))
        (L < pc && pc < U) || error("branch not decided")
        PI/4*(G(h2, pc - g) - G(h2, L - g) + G(h1, U) - G(h1, pc))
    end
    DB = QQ(1971, 10000); EB = DB/2
    VB = 9PI^2/8
    Tb = VB - 23S(2) - S(2 + DB) + 88Pi(1, 1, 1//2) + 8Pi(1, 1 + DB/2, 1//2)
    check!(R, io, "T on the one-centre ray at delta = 0.1971 exceeds 8", Tb > 8; got=@sprintf("%.9f", mid(Tb)), expected="8.000015233")
    T0 = VB - 24S(2) + 96Pi(1, 1, 1//2)
    check!(R, io, "T at the root system is 7.906940...", abs(T0 - Arb(QQ(790694, 100000); prec=p)) < Arb(QQ(1, 100000); prec=p); got=@sprintf("%.9f", mid(T0)), expected="7.906940156")
    lam1 = (1 + EB)^2 - 1
    check!(R, io, "(1 + eta_B)^2 - 1 < 1/3", lam1 < 1//3; got=@sprintf("%.5f", Float64(lam1)), expected="0.20681")
    g_over(e) = QQ(1, 3)*(4 - 6e + 4e*e - e^3) - 768e^3
    check!(R, io, "g(eta)/eta decreases on [0, eta_B]: -2 + (8/3) eta_B < 0", -2 + QQ(8, 3)*EB < 0)
    check!(R, io, "g(eta)/eta >= 2/5 on [0, eta_B]", g_over(EB) >= 2//5; got=@sprintf("%.5f", Float64(g_over(EB))), expected="0.41379")
    check!(R, io, "pyramid: (2 + eta)/(1 + eta) <= 2", (2 + EB)/(1 + EB) <= 2)
    R
end
end
