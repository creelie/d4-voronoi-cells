# Ports of multi_cap/radial_count_check.py (thm:count31) and multi_cap/radial_case_check.py
# (thm:count30 cases 1-4, prop:C-radial, prop:G-radial, prop:levels).
# Exact Rational{BigInt}; Bernstein bounds in Float64 with an a-priori error allowance;
# S and Pi in Arb ball arithmetic (128 bits).
module Radial
using ..Common
using ..Common: QQ
using ..Exact
using ..Exact: issym
using Arblib, JSON

const Rsq = QQ(3,2)

# ---------- exact polynomials in s on [0,1] (index 1 = constant term)
function pmul(p::Vector{QQ}, q::Vector{QQ})
    out = zeros(QQ, length(p)+length(q)-1)
    for (i,a) in enumerate(p)
        a == 0 && continue
        for (j,b) in enumerate(q); out[i+j-1] += a*b; end
    end
    out
end
function padd(p::Vector{QQ}, q::Vector{QQ}, cp=1, cq=1)
    n = max(length(p), length(q)); out = zeros(QQ, n)
    for i in eachindex(p); out[i] += cp*p[i]; end
    for i in eachindex(q); out[i] += cq*q[i]; end
    out
end
function cheb_T(alpha::QQ, beta::QQ, n::Int)
    T = [QQ[1], QQ[alpha, beta]]
    for _ in 2:n; push!(T, padd(pmul(QQ[2alpha,2beta], T[end]), T[end-1], 1, -1)); end
    T[1:n+1]
end
function cheb_U(alpha::QQ, beta::QQ, n::Int)
    U = [QQ[1], QQ[2alpha, 2beta]]
    for _ in 2:n; push!(U, padd(pmul(QQ[2alpha,2beta], U[end]), U[end-1], 1, -1)); end
    U[1:n+1]
end
function bernstein(p::Vector{QQ}, n::Int)
    c = vcat(p, zeros(QQ, n+1-length(p)))
    [sum(QQ(binomial(big(i),big(j)), binomial(big(n),big(j)))*c[j+1] for j in 0:i) for i in 0:n]
end

"de Casteljau at the midpoint along `axis` of a Float64 array; the same operations as numpy's"
function split_mid(B::Array{Float64,N}, axis::Int) where N
    n = size(B, axis) - 1
    L = similar(B); Rt = similar(B)
    cur = copy(B)
    sel(A, i) = selectdim(A, axis, i)
    sel(L,1) .= sel(cur,1); sel(Rt,n+1) .= sel(cur,n+1)
    for k in 1:n
        for i in 1:n+1-k
            sel(cur,i) .= (sel(cur,i) .+ sel(cur,i+1)) .* 0.5
        end
        sel(L,k+1) .= sel(cur,1); sel(Rt,n+1-k) .= sel(cur,n+1-k)
    end
    (L, Rt)
end

# ---------- S and Pi in Arb
Rarb() = sqrt(A_(Rsq))
PIarb() = Arb(π; prec=PREC[])
function S_arb(d::QQ)
    d^2 >= 6 && return Arb(0; prec=PREC[])
    h = A_(d)/2; R2 = A_(Rsq)
    anti = h/8*(5R2 - 2h*h)*sqrt(R2 - h*h) + 3R2^2/8*asin(h/Rarb())
    4PIarb()/3*(3R2^2/8*PIarb()/2 - anti)
end
function Garb(h::Arb, x::Arb)
    tx = tan(x); R2 = A_(Rsq)
    R2^2*x - 2R2*h*h*tx + h^4*(tx + tx^3/3)
end
"a rigorous lower bound (Float64) for Pi(d1/2, d2/2, u); d1, d2, u rational"
function Pi_lower(d1::QQ, d2::QQ, u::QQ)
    (d1^2 >= 6 || d2^2 >= 6 || u <= -1) && return 0.0
    h1, h2, ua = A_(d1)/2, A_(d2)/2, A_(u)
    g = acos(ua); R = Rarb()
    a1, a2 = acos(h1/R), acos(h2/R)
    L = max(-a1, g - a2); U = min(a1, g + a2)
    (L < U) || return 0.0
    pc = atan(h2 - h1*cos(g), h1*sin(g))
    pc = min(max(pc, L), U)
    val = PIarb()/4*(Garb(h2, pc-g) - Garb(h2, L-g) + Garb(h1, U) - Garb(h1, pc))
    isfinite(val) || return 0.0
    max(arb_lo(val), 0.0)
end
# floating-point filter (never decides anything): truncated_search.pair
const R2f = 1.5; const Rf = sqrt(R2f); const R4f = R2f*R2f
Gf(h, x) = (tx = tan(x); R4f*x - 2*R2f*h*h*tx + h^4*(tx + tx^3/3))
function pair_float(h1::Float64, h2::Float64, u::Float64)
    g = acos(clamp(u,-1.0,1.0))
    a1 = acos(clamp(h1/Rf,-1.0,1.0)); a2 = acos(clamp(h2/Rf,-1.0,1.0))
    L = max(-a1, g-a2); U = min(a1, g+a2)
    pc = clamp(atan(h2 - h1*cos(g), h1*max(sin(g),1e-15)), L, U)
    val = Gf(h2, pc-g) - Gf(h2, L-g) + Gf(h1, U) - Gf(h1, pc)
    L < U ? pi/4*val : 0.0
end

# ---------- the pair branch and bound (step 2)
function root_coefficients(A::Vector{Matrix{QQ}}, D, r, c1, c2, dmax)
    alpha, beta = (4-c1)/c2, 2*(dmax-2)/c2
    umax = (2dmax^2-4)/(2dmax^2)
    Ts = [bernstein(p, r) for p in cheb_T(alpha, beta, r)]
    Us = [bernstein(p .// (k), D) for (k,p) in enumerate(cheb_U(QQ(-1), umax+1, D))]
    B0 = zeros(QQ, r+1, r+1, D+1)
    for k in 0:D, a in 0:r, b in 0:r
        A[k+1][a+1,b+1] == 0 && continue
        for i in 0:r
            ci = A[k+1][a+1,b+1]*Ts[a+1][i+1]
            ci == 0 && continue
            for j in 0:r
                cij = ci*Ts[b+1][j+1]
                cij == 0 && continue
                for l in 0:D; B0[i+1,j+1,l+1] += cij*Us[k+1][l+1]; end
            end
        end
    end
    (B0, umax)
end
"returns (ok, closed, dropped, arbcalls, message); `maxboxes` > 0 stops early (smoke tests)"
function pair_check(A, D, r, c1, c2, dmax; depth_max=80, maxboxes=0)
    B0, umax = root_coefficients(A, D, r, c1, c2, dmax)
    Bf = Float64.(B0)
    Cmax = maximum(abs, B0)
    err = Float64(Cmax)*2.0^-53*(1 + depth_max*max(r,D))*1.01
    dval(i,n) = 2 + (dmax-2)*QQ(i, big(2)^n)
    uval(i,n) = -1 + (umax+1)*QQ(i, big(2)^n)
    amaxq(d1,d2) = (d1*d1+d2*d2-4)/(2d1*d2)
    wu = (umax+1)/(dmax-2)/3
    stack = Tuple{NTuple{6,Int},Array{Float64,3},Int}[((0,0,0,0,0,0), Bf, 0)]
    closed = dropped = arbcalls = 0
    memo = Dict{NTuple{3,QQ},Float64}()
    while !isempty(stack)
        (i1,n1,i2,n2,i3,n3), B, depth = pop!(stack)
        dlo, dhi = dval(i1,n1), dval(i1+1,n1)
        elo, ehi = dval(i2,n2), dval(i2+1,n2)
        ulo = uval(i3,n3)
        if dlo > ehi || ulo > amaxq(dhi,ehi)
            dropped += 1; continue
        end
        kub = nextfloat(maximum(B) + err)
        if kub <= 0
            closed += 1; continue
        end
        pf = pair_float(Float64(dhi)/2, Float64(ehi)/2, Float64(ulo))
        if kub <= pf - 1e-9
            key = (dhi,ehi,ulo)
            v = get(memo, key, nothing)
            if v === nothing
                v = Pi_lower(dhi,ehi,ulo); memo[key] = v; arbcalls += 1
            end
            if kub <= v
                closed += 1; continue
            end
        end
        depth >= depth_max && return (false, closed, dropped, arbcalls, "depth limit at d in [$dlo, $dhi], d' in [$elo, $ehi], u from $ulo")
        (maxboxes > 0 && closed+dropped >= maxboxes) && return (true, closed, dropped, arbcalls, "stopped early (smoke test)")
        w = (QQ(1,big(2)^n1), QQ(1,big(2)^n2), QQ(1,big(2)^n3)*wu)
        ax = w[1] >= w[2] && w[1] >= w[3] ? 1 : (w[2] >= w[3] ? 2 : 3)   # first index of the maximum
        L, Rt = split_mid(B, ax)
        if ax == 1
            push!(stack, ((2i1,n1+1,i2,n2,i3,n3), L, depth+1)); push!(stack, ((2i1+1,n1+1,i2,n2,i3,n3), Rt, depth+1))
        elseif ax == 2
            push!(stack, ((i1,n1,2i2,n2+1,i3,n3), L, depth+1)); push!(stack, ((i1,n1,2i2+1,n2+1,i3,n3), Rt, depth+1))
        else
            push!(stack, ((i1,n1,i2,n2,2i3,n3+1), L, depth+1)); push!(stack, ((i1,n1,i2,n2,2i3+1,n3+1), Rt, depth+1))
        end
    end
    (true, closed, dropped, arbcalls, "$closed boxes closed, $dropped dropped, $arbcalls values of Pi in Arb")
end

function bracket_poly(A, z, D, r, c1, c2, dmax)
    alpha, beta = (4-c1)/c2, 2*(dmax-2)/c2
    Tq = cheb_T(alpha, beta, r)
    Asum = sum(A)
    q = QQ[0]
    for a in 0:r
        for b in 0:r
            Asum[a+1,b+1] != 0 && (q = padd(q, pmul(Tq[a+1],Tq[b+1]), 1, Asum[a+1,b+1]/2))
        end
        q = padd(q, Tq[a+1], 1, -z[a+1])
    end
    qe = bernstein(q, 2r)
    (Float64.(qe), Float64(maximum(abs, qe))*2.0^-53*(1+60*2r)*1.01)
end

parsemat(a) = QQ[fr(a[i][j]) for i in eachindex(a), j in eachindex(a[1])]
function positivity(A, z, D)
    good = all(ldl_psd(A[k+1])[1] for k in 1:D)
    g0, piv = ldl_psd(A[1])
    w = solve_exact(A[1], z)
    tq = sum(z .* w)
    t = QQ(cld(numerator(tq)*big(2)^48, denominator(tq)), big(2)^48)
    n = length(z)
    Zm = [A[1] z; permutedims(z) t]
    gz, pz = ldl_psd(Zm)
    (good, g0 && all(p > 0 for p in piv), gz, t, piv, pz)
end

# ---------- radial_count_check.py
function run_count(path::AbstractString, io::IO=stdout; expect=nothing)
    R = Common.Report("multi_cap/radial_count_check " * basename(path))
    c = JSON.parsefile(path)
    M, D, r = c["M"], c["D"], c["r"]
    c1, c2, dmax = fr(c["c1"]), fr(c["c2"]), fr(c["dmax"])
    A = [parsemat(a) for a in c["A"]]; z = QQ[fr(x) for x in c["z"]]; m = fr(c["m"])
    Common.check!(R, io, "dmax exceeds sqrt 6", dmax^2 > 6; got=Float64(dmax^2))
    Common.check!(R, io, "A_0, ..., A_$D symmetric", all(issym, A))
    good, pd0, gz, t, piv, pz = positivity(A, z, D)
    Common.check!(R, io, "A_1, ..., A_$D positive semidefinite (exact LDL^T)", good)
    Common.check!(R, io, "A_0 positive definite (exact LDL^T)", pd0; got="least pivot $(Float64(minimum(piv)))", expected=get(something(expect,Dict()), "pivot", ""))
    Common.check!(R, io, "[[A_0, z], [z^T, t]] positive semidefinite (exact LDL^T)", gz; got="t = $(Float64(t))", expected=get(something(expect,Dict()), "t", ""))
    ok, closed, dropped, calls, msg = pair_check(A, D, r, c1, c2, dmax)
    e = get(something(expect,Dict()), "pair", nothing)
    Common.check!(R, io, "pair inequality K <= Pi on the admissible domain", ok && (e === nothing || (closed,dropped)==e[1:2]); got="$closed closed, $dropped dropped, $calls Arb values", expected= e===nothing ? "" : "$(e[1]) closed, $(e[2]) dropped, $(e[3]) Arb values")
    qb, err1 = bracket_poly(A, z, D, r, c1, c2, dmax)
    stack = Tuple{Int,Int,Vector{Float64},Int}[(0,0,qb,0)]; nb = 0; gmax = -1e9; okb = true
    while !isempty(stack)
        i, n, B, depth = pop!(stack)
        dlo = 2 + (dmax-2)*QQ(i, big(2)^n)
        bound = S_arb(dlo) + Arb(nextfloat(maximum(B)+err1); prec=PREC[])
        if bound < A_(m)
            nb += 1; gmax = max(gmax, arb_hi(bound)); continue
        end
        if depth >= 60; okb = false; break; end
        L, Rt = split_mid(B, 1)
        push!(stack, (2i,n+1,L,depth+1)); push!(stack, (2i+1,n+1,Rt,depth+1))
    end
    e = get(something(expect,Dict()), "intervals", nothing)
    Common.check!(R, io, "S(d) + K(d, d, 1)/2 - z . p(d) <= m on [2, dmax]", okb && (e===nothing || nb==e); got="m = $(Float64(m)), largest bound $gmax, $nb intervals", expected= e===nothing ? "" : "$e intervals")
    Common.check!(R, io, "m < 0, so the bound falls with the count", m < 0)
    total = M*A_(m) + A_(t)/2; target = 9PIarb()^2/8 - 8
    Common.check!(R, io, "M m + t/2 < 9 pi^2/8 - 8", total < target; got="$(Float64(total)) < $(Float64(target))", expected=get(something(expect,Dict()), "total", ""))
    R
end

# ---------- radial_case_check.py
const RKS = QQ(20161,10000)
function arb_up(x::Arb)
    QQ(arb_hi(x))
end
parse_case(lst) = Dict{QQ,Tuple{Int,Int}}(fr(k)=>(lo,hi) for (k,lo,hi) in lst)
function implied(case, M, r)
    lo, hi = 0, M
    for (rj,(l,h)) in case
        rj <= r && (lo = max(lo,l))
        rj >= r && (hi = min(hi,h))
    end
    RKS >= r && (hi = min(hi,24))
    (lo,hi)
end
function check_tree(node, M, out)
    case = parse_case(node["case"])
    if !haskey(node, "children"); push!(out, case); return true; end
    r, c = fr(node["split"][1]), node["split"][2]
    lo, hi = implied(case, M, r)
    ch1, ch2 = parse_case(node["children"][1]["case"]), parse_case(node["children"][2]["case"])
    ok = true
    for (ch, want_lo, want_hi) in ((ch1, nothing, c-1), (ch2, c, nothing))
        rest = Dict(k=>v for (k,v) in ch if k != r); pr = Dict(k=>v for (k,v) in case if k != r)
        ok &= rest == pr && haskey(ch, r)
        haskey(ch, r) || return false
        l2, h2 = ch[r]
        if want_hi !== nothing
            ok &= h2 == want_hi && l2 <= lo
        else
            ok &= l2 == want_lo && h2 >= hi
        end
    end
    ok && check_tree(node["children"][1], M, out) && check_tree(node["children"][2], M, out)
end
function normalised(case, M)
    case = copy(case)
    lo, hi = get(case, RKS, (0,M))
    case[RKS] = (lo, min(hi,24))
    rs = sort(collect(keys(case))); out = Dict{QQ,Tuple{Int,Int}}()
    for (i,ri) in enumerate(rs)
        lo, hi = case[ri]
        for rj in rs[1:i-1]; lo = max(lo, case[rj][1]); end
        for rj in rs[i+1:end]; hi = min(hi, case[rj][2]); end
        out[ri] = (lo,hi)
    end
    out
end
"max of sum_b n_b w_b over integer n >= 0, sum n = M, lo_j <= n_0+...+n_j <= hi_j; exact. nothing if empty"
function dp_max(M, nrs, bounds, weights::Vector{QQ})
    best = Dict{Int,QQ}(0=>QQ(0))
    for b in 0:nrs
        new = Dict{Int,QQ}()
        for (cum,val) in best, n in 0:M-cum
            c2 = cum+n
            if b < nrs
                lo, hi = bounds[b+1]
                (lo <= c2 <= hi) || continue
            elseif c2 != M
                continue
            end
            v = val + n*weights[b+1]
            (!haskey(new,c2) || v > new[c2]) && (new[c2] = v)
        end
        best = new
    end
    get(best, M, nothing)
end
function bracket_check(A, z, D, r, c1, c2, dmax, edges, m)
    qb, err1 = bracket_poly(A, z, D, r, c1, c2, dmax)
    dval(i,n) = 2 + (dmax-2)*QQ(i, big(2)^n)
    stack = Tuple{Int,Int,Vector{Float64},Int}[(0,0,qb,0)]; nint = 0; worst = Inf
    while !isempty(stack)
        i, n, B, depth = pop!(stack)
        dlo, dhi = dval(i,n), dval(i+1,n)
        meet = [b for b in 1:length(edges)-1 if !(dhi < edges[b] || dlo > edges[b+1] || (b > 1 && dhi == edges[b]))]
        mb = minimum(m[b] for b in meet)
        bound = S_arb(dlo) + Arb(nextfloat(maximum(B)+err1); prec=PREC[])
        if bound < A_(mb)
            nint += 1; worst = min(worst, arb_lo(A_(mb)-bound)); continue
        end
        depth >= 60 && return (false, nint, worst)
        L, Rt = split_mid(B, 1)
        push!(stack, (2i,n+1,L,depth+1)); push!(stack, (2i+1,n+1,Rt,depth+1))
    end
    (true, nint, worst)
end

"expected values (closed, dropped, intervals, bound) per certified leaf, read from a Python log"
function expect_from_log(logpath)
    out = Tuple{Int,Int,Int,Float64}[]; cl = dr = ni = 0
    for line in eachline(logpath)
        m = match(r"(\d+) boxes closed, (\d+) dropped", line); m !== nothing && ((cl, dr) = parse.(Int, m.captures))
        m = match(r"(\d+) intervals", line); m !== nothing && (ni = parse(Int, m.captures[1]))
        m = match(r"below the target: \[([0-9.]+)", line); m !== nothing && push!(out, (cl, dr, ni, parse(Float64, m.captures[1])))
    end
    out
end

"""
run_case(path; expect, skip_pair) : radial_case_check.py on one certificate file.
`expect` is a vector with one entry per certified leaf: (closed, dropped, intervals, "bound-prefix").
"""
function run_case(path::AbstractString, io::IO=stdout; expect=nothing, maxboxes=0)
    R = Common.Report("multi_cap/radial_case_check " * basename(path))
    c = JSON.parsefile(path)
    M, D, r = c["M"], c["D"], c["r"]
    c1, c2, dmax = fr(c["c1"]), fr(c["c2"]), fr(c["dmax"])
    target = haskey(c, "level") ? A_(fr(c["level"])) : 9PIarb()^2/8 - 8
    Common.check!(R, io, "dmax exceeds sqrt 6", dmax^2 > 6)
    tl = Dict{QQ,Tuple{Int,Int}}[]
    ok = check_tree(c["tree"], M, tl)
    Common.check!(R, io, "the case tree covers every packing", ok && isempty(c["tree"]["case"]); got="$(length(tl)) leaves")
    lv = c["leaves"]
    Common.check!(R, io, "one certificate per leaf, in order", length(lv)==length(tl) && all(normalised(x,M)==parse_case(y["case"]) for (x,y) in zip(tl,lv)))
    ncert = 0; nres = 0
    for (n_, leaf) in enumerate(lv)
        case = parse_case(leaf["case"]); rs = sort(collect(keys(case)))
        edges = vcat(QQ(2), rs, dmax); bounds = [case[x] for x in rs]
        Sup = QQ[arb_up(S_arb(e)) for e in edges[1:end-1]]
        desc = join(["$lo <= N($(Float64(k))) <= $hi" for (k,(lo,hi)) in sort(collect(case); by=first)], ", ")
        println(io, "leaf $n_: $desc")
        if get(leaf, "residual", false)
            nres += 1; println(io, "  [----] residual case: no certificate"); continue
        end
        if get(leaf, "excluded", false)
            v = dp_max(M, length(rs), bounds, Sup)
            Common.check!(R, io, "leaf $n_ excluded: every count vector has sum of caps below the target", v === nothing || A_(v) < target; got= v===nothing ? "none" : Float64(v))
            continue
        end
        ncert += 1
        e = expect === nothing ? nothing : expect[ncert]
        A = [parsemat(a) for a in leaf["A"]]; z = QQ[fr(x) for x in leaf["z"]]; m = QQ[fr(x) for x in leaf["m"]]; s = fr(leaf["s"])
        Common.check!(R, io, "leaf $n_ edges as in the case", QQ[fr(x) for x in leaf["edges"]] == edges && length(m)==length(edges)-1)
        good, pd0, gz, t, _, _ = positivity(A, z, D)
        Common.check!(R, io, "leaf $n_ A_k and [[A_0, z], [z^T, t]] positive semidefinite (exact LDL^T)", good && pd0 && gz; got="t = $(round(Float64(t); digits=9))")
        okp, closed, dropped, calls, msg = pair_check(A, D, r, c1, c2, dmax; maxboxes=maxboxes)
        Common.check!(R, io, "leaf $n_ K <= Pi on the admissible domain", okp && (e===nothing || maxboxes>0 || (closed,dropped)==(e[1],e[2])); got="$closed closed, $dropped dropped, $calls Arb values"*(maxboxes>0 ? " (truncated smoke run)" : ""), expected= e===nothing ? "" : "$(e[1]) closed, $(e[2]) dropped")
        okb, nint, worst = bracket_check(A, z, D, r, c1, c2, dmax, edges, m)
        Common.check!(R, io, "leaf $n_ f <= m_b on every bin", okb && (e===nothing || nint==e[3]); got="$nint intervals, least room $worst", expected= e===nothing ? "" : "$(e[3]) intervals")
        Common.check!(R, io, "leaf $n_ s >= 0", s >= 0)
        wts = QQ[mb + s*sb for (mb,sb) in zip(m,Sup)]
        v = dp_max(M, length(rs), bounds, wts)
        total = A_(v) - A_(s)*target + A_(t)/2
        Common.check!(R, io, "leaf $n_ largest sum_b n_b m_b over the case, plus t/2, below the target", total < target && (e===nothing || abs(Float64(total)-e[4]) <= 1.01e-7); got="$(Float64(total)) < $(Float64(target))", expected= e===nothing ? "" : string(e[4]))
        if haskey(c, "level")
            un = A_(v) + A_(t)/2
            println(io, "  with no assumption on the count vector: U(Y) <= ", Float64(Arblib.ubound(un), RoundUp))
            push!(R.checks, Common.Check("leaf $n_ unconditional level", true, string(Float64(un)), ""))
        end
    end
    R
end
end
