# Ports of multi_cap/certificate_check.py (thm:certificate, at most 23 contacts) and
# multi_cap/labelled_certificate_check.py (thm:m23-labelled).
# Exact Rational{BigInt} for the certificate, its polynomial and its bound; Arb ball arithmetic
# (128 bits) for A_*, omega, omega', omega'' and the Gamma tables; outward-rounded Float64
# interval arithmetic (prevfloat/nextfloat after every operation) in the branch and bound.
#
# Difference from the Python original: there omega and its derivatives are evaluated by mpmath's
# 53-bit interval context from a sympy expression; here by Arb power series at 128 bits.  The
# tables therefore differ in their last bits (ours are tighter), and the number of boxes may
# differ from the Python log by a few units; both are valid certificates of the same inequality.
module ThreePoint
using ..Common
using ..Common: QQ
using ..Exact
using Arblib, NPZ

const N = 23
const SCALE = 1000
const THIRD = nextfloat(1.0/3.0)
down(x::Float64) = prevfloat(x)
up(x::Float64) = nextfloat(x)

# ---------------------------------------------------------------- exact polynomials in (u, v, t)
const P3 = Dict{NTuple{3,Int},QQ}
function padd(p::P3, q::P3)
    r = copy(p)
    for (k,c) in q; r[k] = get(r,k,QQ(0)) + c; end
    filter!(kv->kv[2]!=0, r)
end
pscale(p::P3, s) = s == 0 ? P3() : P3(k=>c*s for (k,c) in p)
function pmul(p::P3, q::P3)
    r = P3()
    for (k1,c1) in p, (k2,c2) in q
        k = (k1[1]+k2[1], k1[2]+k2[2], k1[3]+k2[3]); r[k] = get(r,k,QQ(0)) + c1*c2
    end
    filter!(kv->kv[2]!=0, r)
end
function ppow(p::P3, n)
    r = P3((0,0,0)=>QQ(1))
    for _ in 1:n; r = pmul(r,p); end
    r
end
const U_ = P3((1,0,0)=>QQ(1)); const V_ = P3((0,1,0)=>QQ(1)); const T_ = P3((0,0,1)=>QQ(1)); const ONE = P3((0,0,0)=>QQ(1))

function legendre_coeffs(d)
    P = [QQ[1], QQ[0,1]]
    for k in 1:d-1
        a = vcat(QQ(0), (2k+1) .* P[k+1]); b = vcat(P[k], zeros(QQ, length(a)-length(P[k])))
        push!(P, [(a[i]-k*b[i])/(k+1) for i in eachindex(a)])
    end
    P[1:d+1]
end
function gegenbauer_S3(d)
    Uc = [QQ[1], QQ[0,2]]
    for k in 1:d-1
        a = vcat(QQ(0), 2 .* Uc[k+1]); b = vcat(Uc[k], zeros(QQ, length(a)-length(Uc[k])))
        push!(Uc, a .- b)
    end
    [Uc[k+1] ./ (k+1) for k in 0:d]
end
function chebyshev(n)
    Tc = [QQ[1], QQ[0,1]]
    for k in 1:n-1
        a = vcat(QQ(0), 2 .* Tc[k+1]); b = vcat(Tc[k], zeros(QQ, length(a)-length(Tc[k])))
        push!(Tc, a .- b)
    end
    Tc[1:n+1]
end
function univariate(coeffs, var)
    p = P3()
    for (i,c) in enumerate(coeffs)
        c == 0 && continue
        k = [0,0,0]; k[var] = i-1; p[Tuple(k)] = c
    end
    p
end
function phi_poly(k, LC)
    x = padd(T_, pscale(pmul(U_,V_), -1))
    s2 = pmul(padd(ONE, pscale(pmul(U_,U_), -1)), padd(ONE, pscale(pmul(V_,V_), -1)))
    out = P3()
    for m in (k%2):2:k
        LC[k+1][m+1] != 0 && (out = padd(out, pscale(pmul(ppow(x,m), ppow(s2,(k-m)÷2)), LC[k+1][m+1])))
    end
    out
end
permute(p::P3, perm) = P3((e[perm[1]],e[perm[2]],e[perm[3]])=>c for (e,c) in p)
const PERMS = [(1,2,3),(1,3,2),(2,1,3),(2,3,1),(3,1,2),(3,2,1)]
function S_entry(i, j, phik, TC)
    y = pmul(pmul(univariate(TC[i+1],1), univariate(TC[j+1],2)), phik)
    s = P3()
    for perm in PERMS; s = padd(s, permute(y, perm)); end
    pscale(s, QQ(1,6))
end
function at_1uu(p::P3)
    r = P3()
    for ((a,b,c),co) in p; k = (b+c,0,0); r[k] = get(r,k,QQ(0)) + co; end
    filter!(kv->kv[2]!=0, r)
end
at_111(p::P3) = sum(values(p); init=QQ(0))

function build_P(d, f, F)
    LC = legendre_coeffs(d); G = gegenbauer_S3(d); TC = chebyshev(d)
    P = P3()
    for k in 0:d
        f[k+1] == 0 && continue
        for var in 1:3; P = padd(P, pscale(univariate(G[k+1], var), f[k+1])); end
    end
    F111 = QQ(0)
    for k in 0:d
        phik = phi_poly(k, LC); n = d-k+1
        for i in 0:n-1, j in 0:n-1
            F[k+1][i+1,j+1] == 0 && continue
            s = S_entry(i, j, phik, TC)
            P = padd(P, pscale(s, F[k+1][i+1,j+1]))
            s1 = at_1uu(s)
            for perm in ((1,2,3),(2,1,3),(3,2,1))
                P = padd(P, pscale(permute(s1, perm), F[k+1][i+1,j+1]/(N-2)))
            end
            F111 += F[k+1][i+1,j+1]*at_111(s)
        end
    end
    (P, F111)
end
function derivative(p::P3, var)
    r = P3()
    for (e,c) in p
        if e[var] > 0
            k = collect(e); k[var] -= 1; r[Tuple(k)] = get(r,Tuple(k),QQ(0)) + c*e[var]
        end
    end
    r
end
is_symmetric(P::P3) = all(get(P, (e[p[1]],e[p[2]],e[p[3]]), QQ(0)) == c for (e,c) in P for p in PERMS)

function load_certificate(path)
    Z = npzread(path)
    d = length(Z["f"]) - 1
    f = QQ[QQ(Float64(x)) for x in Z["f"]]
    F = [QQ[QQ(Float64(Z["F$k"][i,j])) for i in axes(Z["F$k"],1), j in axes(Z["F$k"],2)] for k in 0:d]
    (d, f, F)
end

# ---------------------------------------------------------------- closed forms in ball arithmetic
I_lin(r) = 1/cos(r)^4/4
I_r(r) = r/cos(r)^4/4 - (tan(r) + tan(r)^3/3)/4
I_sin2r(r) = 2*tan(r)^3/3
I_sin2(r) = tan(r)^4/4
piarb() = Arb(π; prec=PREC[])
r_star() = atan(1/sqrt(Arb(2; prec=PREC[])))
"omega as a power series in g about the ball g0 (degree 2): returns (omega, d omega/du, d^2 omega/du^2)"
function omega_jet(g0::Arb)
    rs = ArbSeries((r_star(),); degree=2, prec=PREC[])
    g = ArbSeries((g0, 1); degree=2, prec=PREC[])
    K = -g/4 + sin(g)/4 + tan(g/2)/2*sin(g/2)^2
    prim(r) = I_r(r)/2 + K*I_lin(r) - I_sin2r(r)/4 - tan(g/2)/2*I_sin2(r)
    om = 4*piarb()*(prim(rs) - prim(g/2))
    o0 = om[0]; og = om[1]; ogg = 2*om[2]
    sg = sin(g0); cg = cos(g0)
    (o0, -og/sg, ogg/sg^2 - og*cg/sg^3)
end
function A_star()
    pi_ = piarb(); rs = r_star()
    prim(r) = 4*(2*pi_^2*I_lin(r) - 23*pi_*(2*I_r(r) - I_sin2r(r)))
    (2*pi_^2 + prim(rs) - prim(Arb(0; prec=PREC[])))/4
end
linspace(a::Float64, b::Float64, n::Int) = (step = (b-a)/n; v = [i*step + a for i in 0:n]; v[end] = b; v)

struct Tables
    us::Vector{Float64}; olow::Vector{Float64}; ohigh::Vector{Float64}
    d1lo::Vector{Float64}; d1hi::Vector{Float64}; du::Float64; M2::Float64
end
"angles from 2 r_* down to g_lo; see certificate_check.omega_tables"
function omega_tables(g_lo::Float64; n=12000, n2=3000)
    g_hi = nextfloat(arb_hi(2*r_star()))
    gams = linspace(g_hi, g_lo, n)
    us = zeros(n+1); olow = zeros(n+1); ohigh = zeros(n+1); d1lo = zeros(n+1); d1hi = zeros(n+1)
    for (i,g) in enumerate(gams)
        gi = Arb(g; prec=PREC[])
        us[i] = up(arb_hi(cos(gi)))
        w, w1, _ = omega_jet(gi)
        olow[i] = max(down(arb_lo(w)), 0.0); ohigh[i] = max(up(arb_hi(w)), 0.0)
        d1lo[i] = down(arb_lo(w1)); d1hi[i] = up(arb_hi(w1))
    end
    @assert all(diff(us) .> 0) && all(diff(olow) .>= 0)
    du = up(maximum(diff(us)) + 1e-12)
    d1lo[1] = min(d1lo[1], -1e-300)
    M2 = 0.0
    edges = linspace(g_hi, g_lo, n2)
    for k in 1:n2
        a, b = edges[k+1], edges[k]
        _, _, w2 = omega_jet(Arb((a,b); prec=PREC[]))
        @assert isfinite(w2)
        M2 = max(M2, abs(arb_lo(w2)), abs(arb_hi(w2)))
    end
    Tables(us, olow, ohigh, d1lo, d1hi, du, up(M2))
end

# ---------------------------------------------------------------- outward-rounded Float64 intervals
@inline function imul(alo::Float64, ahi::Float64, blo::Float64, bhi::Float64)
    p1 = alo*blo; p2 = alo*bhi; p3 = ahi*blo; p4 = ahi*bhi
    (down(min(min(p1,p2),min(p3,p4))), up(max(max(p1,p2),max(p3,p4))))
end
struct IPoly
    exps::Vector{NTuple{3,Int}}
    clo::Vector{Float64}; chi::Vector{Float64}
    deg::NTuple{3,Int}
end
function IPoly(p::P3)
    items = sort(collect(p); by=first)
    exps = [e for (e,_) in items]; cf = [Float64(c) for (_,c) in items]
    IPoly(exps, down.(cf), up.(cf), ntuple(v->maximum(e[v] for e in exps), 3))
end
const PWMAX = 40
struct Work
    plo::Matrix{Float64}; phi::Matrix{Float64}
end
Work() = Work(zeros(PWMAX,3), zeros(PWMAX,3))
"interval evaluation of the polynomial on the box [lo, hi]; the operations and their order are those of Poly.eval"
function ieval(p::IPoly, lo::NTuple{3,Float64}, hi::NTuple{3,Float64}, w::Work)
    for v in 1:3
        w.plo[1,v] = 1.0; w.phi[1,v] = 1.0
        for k in 1:p.deg[v]
            (w.plo[k+1,v], w.phi[k+1,v]) = imul(w.plo[k,v], w.phi[k,v], lo[v], hi[v])
        end
    end
    slo = 0.0; shi = 0.0
    @inbounds for i in eachindex(p.exps)
        (a,b,c) = p.exps[i]
        (ml, mh) = imul(w.plo[a+1,1], w.phi[a+1,1], w.plo[b+1,2], w.phi[b+1,2])
        (ml, mh) = imul(ml, mh, w.plo[c+1,3], w.phi[c+1,3])
        (ml, mh) = imul(ml, mh, p.clo[i], p.chi[i])
        slo = down(slo + ml); shi = up(shi + mh)
    end
    (slo, shi)
end

# ---------------------------------------------------------------- the Gamma tables (labelled certificate)
const Dlab = QQ(21648,10000); const Hlab = Dlab/2; const TAU1 = QQ(1,100); const KC = 256; const MX = 1600
flo(x::QQ) = (f = Float64(x); while QQ(f) > x; f = prevfloat(f); end; f)
fhi(x::QQ) = (f = Float64(x); while QQ(f) < x; f = nextfloat(f); end; f)
flo(x::Arb) = down(arb_lo(x))
fhi(x::Arb) = up(arb_hi(x))
function Scap(d::QQ)
    R2 = A_(QQ(3,2)); h = A_(d)/2; s = R2 - h*h
    4*piarb()/3*(3*R2*R2/8*piarb()/2 - (h/8*(5*R2 - 2*h*h)*sqrt(s) + 3*R2*R2/8*asin(h/sqrt(R2))))
end
s_(d::QQ) = Scap(QQ(2)) - Scap(d)
function fr_hi(tau::QQ, x::QQ)
    R2 = A_(QQ(3,2)); t = A_(tau); xa = A_(x)
    q = (1 - t*xa)/(sqrt(1 - xa*xa)*sqrt(R2 - t*t))
    ql = Arb(Arblib.lbound(q); prec=PREC[])
    ql >= 1 && return 0.0
    fhi((1-ql)^2*(2+ql)/4)
end
amax(d1, d2) = (d1*d1 + d2*d2 - 4)/(2*d1*d2)

struct GammaTables
    ds_lo::Vector{Float64}; ds_hi::Vector{Float64}
    xs_lo::Vector{Float64}; third_lo::Float64
    FR::Matrix{Float64}; AMX::Matrix{Float64}
end
function GammaTables()
    taus = [1 + (Hlab-1)*QQ(k,KC) for k in 0:KC]
    Ss = [Scap(2t) for t in taus]
    ds = [Ss[k]-Ss[k+1] for k in 1:KC]
    ds_lo = flo.(ds); ds_hi = fhi.(ds)
    @assert all(ds_lo .> 0)
    xtop = QQ(fhi(amax(Dlab,Dlab)))
    xs = [QQ(1,3) + (xtop-QQ(1,3))*QQ(m,MX) for m in 0:MX]
    FR = zeros(KC, MX+2)
    for k in 1:KC, (m,x) in enumerate(xs)
        FR[k,m+1] = fr_hi(taus[k], x)
    end
    @assert all(diff(FR[:,2:end]; dims=2) .>= 0) && all(diff(FR; dims=1) .<= 0)
    AMX = [fhi(amax(2taus[k+1], 2taus[l+1])) for k in 1:KC, l in 1:KC]
    @assert all(diff(AMX; dims=2) .> 0)
    GammaTables(ds_lo, ds_hi, flo.(xs), flo(QQ(1,3)), FR, AMX)
end
function col(G::GammaTables, xhi::Float64)
    xhi <= G.third_lo && return 1
    m = searchsortedfirst(G.xs_lo, xhi) + 1           # 1-based column
    m > MX+2 && error("inner product beyond a_D")
    m
end
function gamma_env!(g::Vector{Float64}, G::GammaTables, ahi, bhi)
    ca = col(G, ahi); cb = col(G, bhi); before = 0.0
    @inbounds for k in 1:KC
        cut = up(G.ds_hi[k]*up(G.FR[k,ca] + G.FR[k,cb]))
        c_lo = down(down(G.ds_lo[k]/11) - cut)
        g[k] = before - cut - 1e-12
        before = k == 1 ? c_lo : before + c_lo          # np.cumsum
    end
    @inbounds for k in KC-1:-1:1; g[k] = min(g[k], g[k+1]); end
    g
end
kappa!(out::Vector{Int}, G::GammaTables, ylo::Float64) = (for k in 1:KC; out[k] = searchsortedfirst(view(G.AMX,k,:), ylo); end; out)   # 1-based cell, KC+1 if none
struct GainWork
    G1::Vector{Float64}; G2::Vector{Float64}; G3::Vector{Float64}; kt::Vector{Int}; kv::Vector{Int}; ku::Vector{Int}
end
GainWork() = GainWork(zeros(KC), zeros(KC+1), zeros(KC+1), zeros(Int,KC), zeros(Int,KC), zeros(Int,KC))
function gain(G::GammaTables, W::GainWork, L::NTuple{3,Float64}, Hh::NTuple{3,Float64})
    gamma_env!(W.G1, G, Hh[3], Hh[2]); gamma_env!(W.G2, G, Hh[3], Hh[1]); gamma_env!(W.G3, G, Hh[2], Hh[1])
    W.G2[KC+1] = Inf; W.G3[KC+1] = Inf
    kappa!(W.kt, G, L[3]); kappa!(W.kv, G, L[2])
    out = Inf
    if L[1] <= 0.5
        @inbounds for k in 1:KC
            out = min(out, W.G1[k] + W.G2[W.kt[k]] + W.G3[W.kv[k]])
        end
    else
        kappa!(W.ku, G, L[1])
        @inbounds for k1 in 1:KC, k2 in W.kt[k1]:KC
            out = min(out, W.G1[k1] + W.G2[k2] + W.G3[max(W.kv[k1], W.ku[k2])])
        end
    end
    down(SCALE*out - 1e-9)
end

# ---------------------------------------------------------------- branch and bound
"""
verify(P, tables, lo0, hi0, wmin; extra, counterexample, batch): the interval branch and bound of
certificate_check.verify_domain / labelled_certificate_check.verify on the ordered box lo0..hi0.
extra = (:none,), (:lin, c) or (:gain, G).  Returns (ok, verified, outside, levels, message).
"""
function verify(P::P3, tab::Tables, lo0, hi0, wmin; extra=(:none,), counterexample=false, batch=150000, io::IO=stdout, label="", maxlevel=typemax(Int))
    olow = down.(SCALE .* tab.olow); ohigh = up.(SCALE .* tab.ohigh); d1lo = down.(SCALE .* tab.d1lo); d1hi = up.(SCALE .* tab.d1hi)
    M2 = up(SCALE*tab.M2); dM = up(M2*tab.du); us = tab.us; nus = length(us)
    Pp = IPoly(P); Dp = [IPoly(derivative(P,v)) for v in 1:3]
    HK = [(1,1),(1,2),(1,3),(2,2),(2,3),(3,3)]
    Hs = [IPoly(derivative(derivative(P,v),w)) for (v,w) in HK]
    det = IPoly(P3((0,0,0)=>QQ(1), (1,1,1)=>QQ(2), (2,0,0)=>QQ(-1), (0,2,0)=>QQ(-1), (0,0,2)=>QQ(-1)))
    W = Work(); GW = GainWork()
    Box = Tuple{NTuple{3,Float64},NTuple{3,Float64}}
    cur = Box[(Tuple(Float64.(lo0)), Tuple(Float64.(hi0)))]
    level = 0; n_done = 0; n_out = 0; worst = Inf; t0 = time()
    while !isempty(cur)
        nxt = Box[]
        for s in 1:batch:length(cur)
            chunk = view(cur, s:min(s+batch-1, length(cur)))
            A = Box[]; B = Box[]; maxw = -Inf; anyund = false
            for (L,Hh) in chunk
                alive = (L[1] <= Hh[2]) && (L[2] <= Hh[3])
                (_, dhi) = ieval(det, L, Hh, W)
                alive &= dhi >= 0
                if !alive; n_out += 1; continue; end
                c = ntuple(v->(L[v]+Hh[v])/2, 3)
                r = ntuple(v->max(up(Hh[v]-c[v]), up(c[v]-L[v])), 3)
                q0 = 0.0; q0hi = 0.0
                glo = zeros(3); ghi = zeros(3); m2 = zeros(3)
                for v in 1:3
                    inside = L[v] >= THIRD
                    idx = clamp(searchsortedlast(us, c[v]), 1, nus)
                    q0 = down(q0 + (inside ? olow[idx] : 0.0))
                    q0hi = up(q0hi + (c[v] >= THIRD ? ohigh[min(idx+1, nus)] : 0.0))
                    glo[v] = inside ? down(d1lo[idx] - dM) : 0.0
                    ghi[v] = inside ? up(d1hi[idx] + dM) : 0.0
                    m2[v] = inside ? M2 : 0.0
                end
                (plo, phi) = ieval(Pp, c, c, W)
                q0 = down(q0 - phi); q0hi = up(q0hi - plo)
                if counterexample
                    (dclo, _) = ieval(det, c, c, W)
                    if dclo >= 0 && c[1] <= c[2] && c[2] <= c[3] && q0hi < 0
                        return (false, n_done, n_out, level, "(cert) is negative at the admissible triple $c: at most $(q0hi/SCALE)")
                    end
                end
                if extra[1] === :lin
                    cc = extra[2]::Float64
                    q0 = down(q0 + down(cc*(c[3]-0.5)))
                    glo[3] = down(glo[3] + cc); ghi[3] = up(ghi[3] + cc)
                end
                first = 0.0; contrib = zeros(3)
                for v in 1:3
                    (gl, gh) = ieval(Dp[v], c, c, W)
                    dq_lo = down(glo[v] - gh); dq_hi = up(ghi[v] - gl)
                    g_abs = max(abs(dq_lo), abs(dq_hi))
                    first = up(first + up(g_abs*r[v])); contrib[v] = up(g_abs*r[v])
                end
                second = 0.0
                for (k,(v,w)) in enumerate(HK)
                    (hl, hh) = ieval(Hs[k], L, Hh, W)
                    habs = max(abs(hl), abs(hh))
                    if v == w
                        term = up(up(habs + m2[v])*up(r[v]*r[v]))
                        second = up(second + term); contrib[v] = up(contrib[v] + term)
                    else
                        term = up(up(2*habs)*up(r[v]*r[w]))
                        second = up(second + term); contrib[v] = up(contrib[v] + term); contrib[w] = up(contrib[w] + term)
                    end
                end
                qlo = down(down(q0 - first) - up(second/2))
                if extra[1] === :gain
                    qlo < 0 && (qlo = down(qlo + gain(extra[2]::GammaTables, GW, L, Hh)))
                    for v in 1:3
                        contrib[v] = contrib[v] + (Hh[v] > 0.5 ? 130.0*(Hh[v]-L[v]) : 0.0)
                    end
                end
                worst = min(worst, qlo)
                if qlo >= 0; n_done += 1; continue; end
                anyund = true
                maxw = max(maxw, Hh[1]-L[1], Hh[2]-L[2], Hh[3]-L[3])
                ax = contrib[1] >= contrib[2] && contrib[1] >= contrib[3] ? 1 : (contrib[2] >= contrib[3] ? 2 : 3)
                mid = (L[ax]+Hh[ax])/2
                push!(A, (L, Base.setindex(Hh, mid, ax))); push!(B, (Base.setindex(L, mid, ax), Hh))
            end
            if anyund && maxw < wmin
                return (false, n_done, n_out, level, "undecided box of width below $wmin")
            end
            append!(nxt, A); append!(nxt, B)
        end
        level += 1; cur = nxt
        println(io, "  $label level $(lpad(level,2)): $n_done boxes verified, $n_out outside the domain, $(length(cur)) to bisect, least lower bound so far $(worst/SCALE)  [$(round(Int,time()-t0))s]"); flush(io)
        level >= maxlevel && !isempty(cur) && return (true, n_done, n_out, level, "stopped at level $level (smoke test)")
    end
    (true, n_done, n_out, level, "")
end

function certificate_and_bound(certpath, R, io; tag="step 1")
    d, f, F = load_certificate(certpath)
    f = vcat(f[1], [max(x, QQ(0)) for x in f[2:end]])
    Common.check!(R, io, "$tag: f_k >= 0 for k >= 1 and every F_k positive definite (exact LDL^T)", all(ldl_positive(F[k]) for k in 1:d+1); got="sizes $(size.(F,1))", expected="sizes [9, 8, 7, 6, 5, 4, 3, 2, 1]")
    P, F111 = build_P(d, f, F)
    B = N*(N*f[1] - sum(f))/2 - N*F111/(6*(N-2))
    Ast = A_star()
    gap_hi = QQ(up(arb_hi(8 - Ast)))       # a rational above 8 - A_*
    Bf = B/SCALE
    (P, Bf, Ast, gap_hi)
end

"certificate_check.py 8 1e-6"
function run_certificate(certpath, io::IO=stdout; wmin=1e-6, maxlevel=typemax(Int))
    R = Common.Report("multi_cap/certificate_check (thm:certificate)")
    P, Bf, Ast, gap_hi = certificate_and_bound(certpath, R, io)
    Common.check!(R, io, "step 2: exact bound B/1000 exceeds 8 - A_*", Bf > gap_hi;
        got="B/1000 = $(round(Float64(Bf); digits=9)); 8 - A_* in [$(round(arb_lo(8-Ast); digits=12)), $(round(arb_hi(8-Ast); digits=12))]; excess $(Float64(Bf-gap_hi))",
        expected="B/1000 = 0.092900000; 8 - A_* in [0.092855570294, 0.092855570294]; excess 4.443e-05")
    Common.check!(R, io, "step 3: P has 449 monomials, degree 16; P is symmetric", is_symmetric(P) && length(P)==449 && maximum(sum(e) for e in keys(P))==16;
        got="$(length(P)) monomials, degree $(maximum(sum(e) for e in keys(P)))", expected="449 monomials, degree 16")
    tab = omega_tables(Float64(pi/3))
    Common.check!(R, io, "omega tables", true;
        got="$(length(tab.us)) angles, du <= $(tab.du), omega(60deg) >= $(round(tab.olow[end]; digits=10)), omega' in [$(round(tab.d1lo[end]; digits=8)), $(round(tab.d1hi[end]; digits=8))], |omega''| <= $(round(tab.M2; digits=6))",
        expected="12001 angles, du <= 1.44e-05, omega >= 0.0014454086, omega' in [0.02519166, 0.02519166], |omega''| <= 0.289712")
    ok, nd, no, lev, msg = verify(P, tab, (-1.0,-1.0,-1.0), (0.5,0.5,0.5), wmin; counterexample=true, io=io, maxlevel=maxlevel)
    Common.check!(R, io, "step 3: (cert) holds on the admissible domain (interval branch and bound)", ok && msg=="";
        got="$nd boxes verified, $no outside the domain, $lev levels $msg", expected="313780 boxes verified, 13517 outside the domain, 44 levels (Python, mpmath 53-bit tables)")
    R
end

"labelled_certificate_check.py"
function run_labelled(certpath, io::IO=stdout; wmin=1e-7, skip_region_1=false)
    R = Common.Report("multi_cap/labelled_certificate_check (thm:m23-labelled)")
    P, Bf, Ast, gap_hi = certificate_and_bound(certpath, R, io)
    Common.check!(R, io, "step 1: B/1000 exceeds 8 - A_*; P is symmetric (449 monomials)", Bf > gap_hi && is_symmetric(P) && length(P)==449;
        got="B/1000 = $(round(Float64(Bf); digits=10))", expected="0.0929000002")
    gap_arb = A_(gap_hi)
    sD = s_(Dlab)
    Common.check!(R, io, "step 2: s(D) > 8 - A_* for D = 2.1648", Arb(Arblib.lbound(sD)) > gap_arb; got="s(D) = $(Float64(sD))", expected="0.0928820620010")
    aD = amax(Dlab,Dlab); aD_f = fhi(aD)
    mono = A_(aD)/(1 - A_(Hlab)*A_(aD))
    Common.check!(R, io, "step 2: the share fr decreases in both heights: a_D/(1 - (D/2) a_D) < 2", arb_hi(mono) < 2; got="a_D = $(Float64(aD)); ratio = $(Float64(mono))", expected="a_D = 0.5732294553; ratio = 1.5103410")
    f12 = fr_hi(QQ(1), QQ(1,2))
    Common.check!(R, io, "step 2: fr(1, 1/2) < 1/22", QQ(f12) < QQ(1,22); got=f12, expected="0.02371033")
    fs = fr_hi(QQ(1), QQ(1,2)+TAU1)
    rr = QQ(1,11) - 2*QQ(fs)
    kap = s_(2 + 4TAU1)/A_(4TAU1)
    c_slab = flo(4000*A_(rr)*kap)
    Common.check!(R, io, "step 2: slab constant c = 4000 r kappa > 0", rr > 0 && c_slab > 0;
        got="fr <= $(round(fs; digits=8)), r >= $(round(Float64(rr); digits=8)), kappa >= $(round(flo(kap); digits=8)), c >= $(round(c_slab; digits=6))",
        expected="fr <= 0.02650437, r >= 0.03790036, kappa >= 0.69636299, c >= 105.569628")
    g_lo = acos(aD_f) - 1e-6
    tab = omega_tables(g_lo; n=18000, n2=4000)
    Common.check!(R, io, "omega tables up to a_D", arb_lo(cos(Arb(g_lo; prec=PREC[]))) >= aD_f;
        got="$(length(tab.us)) angles up to u = $(round(tab.us[end]; digits=6)), du <= $(tab.du), |omega''| <= $(round(tab.M2; digits=6))",
        expected="18001 angles up to u = 0.573230, du <= 1.42e-05, |omega''| <= 0.401257")
    if !skip_region_1
        tab1 = omega_tables(Float64(pi/3))
        ok, nd, no, lev, msg = verify(P, tab1, (-1.0,-1.0,-1.0), (0.5,0.5,0.5), 1e-5; counterexample=true, io=io, label="I")
        Common.check!(R, io, "step 3: region I, u <= v <= t <= 1/2: the inequality (cert)", ok; got="$nd boxes verified, $no outside $msg", expected="313780 verified, 13517 outside (certificate_check log)")
    end
    t1 = Float64(QQ(1,2)+TAU1)
    ok, nd, no, lev, msg = verify(P, tab, (-1.0,-1.0,0.5), (t1,t1,t1), wmin; extra=(:lin, c_slab), batch=100000, io=io, label="II_s")
    Common.check!(R, io, "step 4: region II_s, 1/2 <= t <= 0.51: Q0 + c (t - 1/2) >= 0", ok; got="$nd boxes verified, $no outside the domain $msg", expected="30051 boxes verified, 3408 outside the domain")
    G = GammaTables()
    ok, nd, no, lev, msg = verify(P, tab, (-1.0,-1.0,t1), (aD_f,aD_f,aD_f), wmin; extra=(:gain, G), batch=100000, io=io, label="II_f")
    Common.check!(R, io, "step 5: region II_f, 0.51 <= t <= a_D: Q0 + 1000 sum Gamma_i >= 0", ok; got="$nd boxes verified, $no outside the domain $msg", expected="6482 boxes verified, 1242 outside the domain")
    R
end
end
