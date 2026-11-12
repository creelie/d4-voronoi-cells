# covering.jl -- the rigorous part of multi_cap/covering_bound.py (sec:covering-bound):
# for m = 5..24, r_m is the root of 2r - sin 2r = 2 pi/m on (0, pi/2) and the bound is
# (pi m/3) tan^3 r_m (thm:covering-bound, tab:covering-bound); cor:m22 uses
# bound(22) >= 8 > bound(23).  The Python evaluates this in floating point (and adds
# Monte-Carlo and plotting cross-checks, which are not ported); here r_m is enclosed by a
# certified bisection (the left side is increasing) and the bound is an Arb ball.
module Covering
using ..Common
using ..Common: QQ
using ..Exact
using Arblib, Printf

function root_bracket(m::Int; prec=128, steps=100)
    PI = Arb(π; prec=prec)
    f(r::QQ) = (x = Arb(r; prec=prec); 2x - sin(2x) - 2PI/m)
    a, b = QQ(0), QQ(157, 100)                    # 1.57 < pi/2
    @assert f(a) < 0 && f(b) > 0
    for _ in 1:steps
        c = (a + b)/2; v = f(c)
        if v < 0; a = c; elseif v > 0; b = c; else; break; end   # undecided sign: stop (bracket stays valid)
    end
    (a, b)
end
"enclosure of (pi m/3) tan^3 r_m"
function bound(m::Int; prec=128)
    a, b = root_bracket(m; prec=prec)
    r = Arb((Arb(a; prec=prec), Arb(b; prec=prec)); prec=prec)
    Arb(π; prec=prec)*m/3*tan(r)^3
end
const EXPECTED = Dict(22 => 8.046376, 23 => 7.916728, 24 => 7.798989)
function run(io::IO=stdout)
    R = Report("covering_bound (rigorous table)")
    B = Dict(m => bound(m) for m in 5:24)
    for m in 5:24; @printf(io, "  m = %2d: bound in [%.9f, %.9f]\n", m, arb_lo(B[m]), arb_hi(B[m])); end
    check!(R, io, "the bound decreases in m on 5..24", all(B[m] > B[m+1] for m in 5:23))
    check!(R, io, "cor:m22: bound(22) >= 8 > bound(23)", B[22] > 8 && B[23] < 8; got=@sprintf("%.6f, %.6f", arb_lo(B[22]), arb_hi(B[23])), expected="8.046376, 7.916728")
    for (m, e) in sort(collect(EXPECTED))
        check!(R, io, "bound($m) agrees with the value of the Python run", abs(arb_lo(B[m]) - e) < 1e-6 && abs(arb_hi(B[m]) - e) < 1e-6; got=@sprintf("%.6f", arb_lo(B[m])), expected=@sprintf("%.6f", e))
    end
    dens = Arb(π; prec=128)^2/2/B[24]
    check!(R, io, "density pi^2/(2 bound(24)) agrees with the value of the Python run", abs(arb_lo(dens) - 0.632749) < 1e-6; got=@sprintf("%.6f", arb_lo(dens)), expected="0.632749")
    R
end
end
