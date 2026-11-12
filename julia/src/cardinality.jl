# cardinality.jl -- port of multi_cap/certify_cardinality.py: rigorous verification of a
# three-point certificate for A(4, t) <= N - 1 (thm:kissing-stable with N = 25,
# thm:twenty-six with N = 26).  Exact LDL^T of F_0..F_d, B = 1 + f(1) + F(1,1,1) exactly,
# P1 <= e1 on [-1, t] and P2 <= e2 on the ordered admissible domain by second-order Taylor
# forms in outward-rounded Float64 interval arithmetic (the same operations as the Python),
# and the final inequality in exact rational arithmetic.
module Cardinality
using ..Common
using ..Common: QQ
using ..Exact
using ..ThreePoint
using ..ThreePoint: P3, IPoly, Work, ieval, up, derivative, gegenbauer_S3, legendre_coeffs, chebyshev, phi_poly, univariate, S_entry, at_1uu, at_111
using ..Combo: triple_box_check
using NPZ, Printf

"the least double >= t"
top(t::QQ) = (x = Float64(t); QQ(x) >= t ? x : nextfloat(x))
"the largest double <= the decimal txt"
below(txt::AbstractString) = (q = fr(txt); x = Float64(q); QQ(x) <= q ? x : prevfloat(x))

function build(d, f, F)
    LC = legendre_coeffs(d); G = gegenbauer_S3(d); TC = chebyshev(d)
    P2 = P3(); P1 = P3(); F111 = QQ(0)
    for k in 1:d; P1 = ThreePoint.padd(P1, ThreePoint.pscale(univariate(G[k+1], 1), f[k])); end
    for k in 0:d
        phik = phi_poly(k, LC); n = d - k + 1
        for i in 0:n-1, j in 0:n-1
            F[k+1][i+1, j+1] == 0 && continue
            s = S_entry(i, j, phik, TC)
            P2 = ThreePoint.padd(P2, ThreePoint.pscale(s, F[k+1][i+1, j+1]))
            P1 = ThreePoint.padd(P1, ThreePoint.pscale(at_1uu(s), 3*F[k+1][i+1, j+1]))
            F111 += F[k+1][i+1, j+1]*at_111(s)
        end
    end
    P1 = ThreePoint.padd(P1, P3((0,0,0) => QQ(1)))
    (P1, P2, F111)
end

"P1(u) <= e1 on [-1, top(t)]; returns (ok, intervals, largest centre value)"
function verify_1d(P1::P3, t::QQ, e1::Float64; wmin=1e-7)
    Pp = IPoly(P1); D = IPoly(derivative(P1, 1)); H = IPoly(derivative(derivative(P1, 1), 1))
    W = Work()
    cur = [(-1.0, top(t))]; worst = -Inf; n_ok = 0
    while !isempty(cur)
        und = Tuple{Float64,Float64}[]
        for (lo, hi) in cur
            c = (lo + hi)/2; r = max(up(hi - c), up(c - lo))
            (_, p_c) = ieval(Pp, (c, 0.0, 0.0), (c, 0.0, 0.0), W)
            (gl, gh) = ieval(D, (c, 0.0, 0.0), (c, 0.0, 0.0), W); g = max(abs(gl), abs(gh))
            (hl, hh) = ieval(H, (lo, 0.0, 0.0), (hi, 0.0, 0.0), W); h = max(abs(hl), abs(hh))
            upper = up(up(p_c + up(g*r)) + up(h*up(r*r)/2))
            worst = max(worst, p_c)
            if upper <= e1; n_ok += 1; else; push!(und, (lo, hi)); end
        end
        (!isempty(und) && maximum(hi - lo for (lo, hi) in und) < wmin) && return (false, n_ok, worst)
        cur = vcat([(lo, (lo + hi)/2) for (lo, hi) in und], [((lo + hi)/2, hi) for (lo, hi) in und])
    end
    (true, n_ok, worst)
end

"""
    run(certpath; e1="1e-6", e2="1e-4", wmin=1e-5, N=25, expect=nothing, maxlevel)

`expect` = (B to 9 decimals as a string, intervals, boxes) from the log of the Python run.
"""
function run(certpath::AbstractString, io::IO=stdout; e1="1e-6", e2="1e-4", wmin=1e-5, N=25, expect=nothing, maxlevel=typemax(Int))
    R = Report("certify_cardinality: $(basename(certpath))")
    Z = npzread(certpath)
    d = length(Z["f"]); t = fr(string(Float64(Z["t"][])))
    e1f = below(e1); e2f = below(e2)
    f = [max(QQ(Float64(x)), QQ(0)) for x in Z["f"]]
    F = Matrix{QQ}[]
    for k in 0:d
        Mx = QQ.(Float64.(Z["F$k"])); n = size(Mx, 1)
        for i in 1:n, j in i+1:n
            m = (Mx[i,j] + Mx[j,i])/2; Mx[i,j] = m; Mx[j,i] = m
        end
        push!(F, Mx)
    end
    println(io, "certificate degree $d, t = $t = 1/2 + $(t - 1//2) exactly; the boxes run to $(top(t))")
    check!(R, io, "1. f_k >= 0 and F_0..F_$d positive definite (exact LDL^T)", all(ldl_positive(Mx) for Mx in F))
    P1, P2, F111 = build(d, f, F)
    B = 1 + sum(f) + F111
    Bs = @sprintf("%.9f", Float64(B))
    check!(R, io, "2. B = 1 + f(1) + F(1,1,1) (exact rational), value as in the Python log", expect === nothing || Bs == expect[1];
           got="$Bs; $(length(P2)) monomials in P2", expected=expect === nothing ? "" : expect[1])
    ok1, nint, w1 = verify_1d(P1, t, e1f)
    check!(R, io, "3. (i) P1 <= $e1 on [-1, t]", ok1 && (expect === nothing || nint == expect[2]);
           got=@sprintf("%d intervals, largest centre value %.3e", nint, w1), expected=expect === nothing ? "" : "$(expect[2]) intervals")
    tf = top(t)
    ok2, msg, nd, _ = triple_box_check(P2, (tf, tf, tf), 2, e2f; wmin=wmin, label="(ii)", maxlevel=maxlevel, io=io)
    smoke = maxlevel != typemax(Int)
    check!(R, io, "4. (ii) P2 <= $e2 on the ordered admissible domain", ok2 && (smoke || expect === nothing || nd == expect[3]);
           got=msg, expected=expect === nothing ? "" : "$(expect[3]) boxes")
    lhs = QQ(N - 1)*(1 - QQ(e1f)) - QQ((N - 1)*(N - 2))*QQ(e2f)
    check!(R, io, "5. at |C| = $N: (N-1)(1-e1) - (N-1)(N-2) e2 > B - 1 (exact)", lhs > B - 1;
           got=@sprintf("%.6f vs %.6f", Float64(lhs), Float64(B - 1)))
    println(io, Common.allpass(R) ? (smoke ? "PARTIAL (branch and bound cut short; not a proof)" : "PASS: A(4, $(Float64(t))) <= $(N-1)") : "FAILED")
    R
end
end
