# combo.jl -- port of multi_cap/combo30_check.py, combo_case_check.py and
# combo_mixtype_check.py: the combined two-point (distance-labelled) and typed
# three-point certificates of statement (C) (thm:count30 residual case, prop:levels,
# the 24..29-centre cases of the gap closure).
#
# What is proved, and with which arithmetic (the same as the Python):
#   (a) exact LDL^T of the rounded blocks (Rational{BigInt});
#   (b) K + PAIR3 - c2 <= Pi on the admissible pairs: tensor Bernstein coefficients (exact,
#       then Float64 de Casteljau with an a-priori rounding allowance), Pi at box corners
#       in Arb;
#   (c) f <= m_b on every bin (Radial.bracket_check);
#   (d) TRIPLE3 <= c3 on the admissible triples: second-order Taylor forms in
#       outward-rounded Float64 interval arithmetic;
#   (e) the bound over all count vectors in exact rational arithmetic, compared with the
#       target in Arb.
#
# Two inputs of the Python are floating-point and cannot be reproduced bit for bit:
#   * psd_exact rounds numpy's eigh; here LinearAlgebra.eigen is used (the exact matrices
#     differ in the last bits; each is verified positive by exact LDL^T, which is all the
#     proof needs);
#   * the thresholds c2, c3 are "largest sampled value + margin" with numpy's generator
#     and scipy's SLSQP; here they are sampled on the same grids plus Julia's generator and
#     refined by a pattern search.  Any threshold that the branch and bound then proves is
#     valid, so this affects box counts and the last digits of the bound, not rigour.  If the
#     branch and bound of a triple inequality meets a point above its threshold (the sampling
#     missed a peak), the threshold is raised from that point and the inequality proved again
#     (`retries`); the Python stops there and is rerun by hand with TRIPLE_EXTRA.
# The Python options BORDER, LEVEL, TRIPLE_EXTRA, PAIR_MARGIN, SKIP_EMPTY, REFINE_TRIPLES,
# REFINE_PAIRS, SKIP_COVERED, REGION_TRIPLES, EXTRA_CASES and TRIPLE_TAYLOR=mono are
# keyword arguments of `run`.  PRECHECK, ONLY, CHECKPOINT, REFINE_DIVERSE and REFINE_GRID
# (bookkeeping and float search aids) are not ported.
module Combo
using ..Common
using ..Common: QQ
using ..Exact
using ..Radial
using ..Radial: cheb_T, cheb_U, bernstein, split_mid, Pi_lower, pair_float
using ..ThreePoint
using ..ThreePoint: P3, IPoly, Work, ieval, up, down, derivative, gegenbauer_S3, legendre_coeffs, chebyshev, phi_poly, univariate
using Arblib, NPZ, JSON, LinearAlgebra, Random, Printf

const D2 = 12
const R2 = 4
const DMAX = QQ(24494898, 10000000)
const C1 = QQ(1199745671238, 269636686585)     # radial_count_sdp.C1, C2
const C2 = QQ(121198924898, 269636686585)

amaxq(d1, d2) = (d1*d1 + d2*d2 - 4)/(2*d1*d2)
dyad(x::Float64) = QQ(round(BigInt, x*2.0^48), big(2)^48)
"the least double >= x, as a rational"
function above(x::QQ)
    y = Float64(x)
    QQ(y) < x && (y = nextfloat(y))
    QQ(y)
end
above(x::Float64) = QQ(x)
function below(x::QQ)
    y = Float64(x)
    QQ(y) > x && (y = prevfloat(y))
    QQ(y)
end

# ---------------------------------------------------------------- the case
struct Case
    name::String
    M::Int
    types::Vector{String}
    trange::Dict{String,Tuple{QQ,QQ}}
    edges::Vector{QQ}
    bintype::Vector{Union{String,Nothing}}
    counts::Vector{Vector{Int}}
    raw::Any
end
num(x) = x == "sqrt6" ? DMAX : fr(string(x))
function count_vectors(M, nb, constraints, exclude)
    out = Vector{Int}[]
    okc(v, cs) = all(lo <= sum(v[i+1] for i in bs) <= hi for (bs, lo, hi) in cs)
    function rec(prefix, left)
        if length(prefix) == nb - 1
            v = vcat(prefix, left)
            (okc(v, constraints) && !any(okc(v, X) for X in exclude)) && push!(out, v)
            return
        end
        for k in 0:left; rec(vcat(prefix, k), left - k); end
    end
    rec(Int[], M)
    out
end
function load_case(path)
    J = JSON.parsefile(path)
    types = String.(J["types"])
    trange = Dict(String(k) => (num(v[1]), num(v[2])) for (k, v) in J["trange"])
    bins = J["bins"]
    edges = vcat(num(bins[1][2]), [num(b[3]) for b in bins])
    bintype = Union{String,Nothing}[String(b[1]) for b in bins]
    for (i, b) in enumerate(bins)
        @assert num(b[2]) == edges[i] && trange[b[1]][1] <= num(b[2]) && num(b[3]) <= trange[b[1]][2]
    end
    cv = count_vectors(J["M"], length(bins), J["constraints"], get(J, "exclude", []))
    Case(basename(path), J["M"], types, trange, edges, bintype, cv, J)
end
"the fixed case of combo30_check.py"
function case30()
    tr = Dict("A" => (QQ(2), QQ(205,100)), "B" => (QQ(205,100), QQ(225,100)), "F" => (QQ(24,10), DMAX))
    Case("thm:count30 residual", 30, ["A","B","F"], tr, [QQ(2), QQ(205,100), QQ(215,100), QQ(225,100), QQ(24,10), DMAX],
         Union{String,Nothing}["A","B","B",nothing,"F"],
         [[22,1,1,0,6],[22,2,0,0,6],[23,0,1,0,6],[23,1,0,0,6],[24,0,0,0,6]], nothing)
end
function tcounts(C::Case, nb)
    out = Dict(s => 0 for s in C.types)
    for (ty, n) in zip(C.bintype, nb); ty === nothing || (out[ty] += n); end
    out
end
Npair(tc, s, t) = s == t ? binomial(tc[s], 2) : tc[s]*tc[t]
function Ntriple(tc, combo)
    out = 1
    for s in unique(combo); out *= binomial(tc[s], count(==(s), combo)); end
    out
end
function Nbins(nb, bs)
    out = 1
    for b in unique(bs); out *= binomial(nb[b], count(==(b), bs)); end
    out
end
tmax(C::Case, s, t) = amaxq(C.trange[s][2], C.trange[t][2])
combos2(T) = [(T[i], T[j]) for i in eachindex(T) for j in i:length(T)]
combos3(T) = [(T[i], T[j], T[k]) for i in eachindex(T) for j in i:length(T) for k in j:length(T)]

# ---------------------------------------------------------------- the certificate, made exact
function psd_exact(M::Matrix{Float64}, shift::QQ=QQ(1, big(2)^40))
    M = (M + M')/2
    E = eigen(Symmetric(M))
    M = (E.vectors .* max.(E.values, 0.0)') * E.vectors'; M = (M + M')/2
    n = size(M, 1)
    out = zeros(QQ, n, n)
    for i in 1:n, j in i:n
        out[i,j] = dyad(M[i,j]); out[j,i] = out[i,j]
    end
    for i in 1:n; out[i,i] += shift; end
    out
end
"the three-point blocks (typed3pt.Layout(m, d3)) from the stored vector: key (:f,k,0) or (:F,P,k), P zero-based"
function blocks_from_x3(x3::Vector{Float64}, d3::Int, m::Int)
    names = vcat([(:f, k, 0) for k in 1:d3], [(:F, P, k) for P in 0:m-1 for k in 0:d3])
    out = Dict{Tuple{Symbol,Int,Int},Matrix{Float64}}(); off = 0
    for name in names
        s = name[1] === :f ? m : m*(d3 - name[3] + 1)
        n = s*(s+1) ÷ 2
        v = x3[off+1:off+n]; off += n
        Mx = zeros(s, s)
        for i in 1:s; Mx[i,i] = v[i]; end
        k = s
        for i in 1:s, j in i+1:s
            k += 1; Mx[i,j] = v[k]; Mx[j,i] = v[k]
        end
        out[name] = Mx
    end
    @assert off == length(x3) (off, length(x3))
    out
end
const Blocks = Dict{Tuple{Symbol,Int,Int},Matrix{QQ}}
"one three-point certificate and the certificate type (zero-based) it gives to every case type"
struct Group
    Bk::Blocks
    tmap::Vector{Int}
end

# ---------------------------------------------------------------- the three-point polynomials
function point_term(Bk::Blocks, d3, s)
    v = sum(Bk[(:f,k,0)][s+1,s+1] for k in 1:d3)
    mk = d3 + 1; blk = Bk[(:F,s,0)]
    v + sum(blk[s*mk+i, s*mk+j] for i in 1:mk, j in 1:mk)
end
ladd(p, q, cp=1, cq=1) = Radial.padd(p, q, cp, cq)
lmul(p, q) = Radial.pmul(p, q)
function pair3_poly(Bk::Blocks, d3, s, t)
    G = gegenbauer_S3(d3); TC = chebyshev(d3)
    p = QQ[0]
    for k in 1:d3; p = ladd(p, G[k+1], 1, 2*Bk[(:f,k,0)][s+1,t+1]); end
    for (P, Q) in ((s, t), (t, s))
        mk = d3 + 1; blk = Bk[(:F,P,0)]
        for i in 0:mk-1, j in 0:mk-1
            c = blk[P*mk+i+1, Q*mk+j+1]
            c != 0 && (p = ladd(p, TC[j+1], 1, 2c))
        end
        for k in 0:d3
            mk = d3 - k + 1; blk = Bk[(:F,P,k)]
            w = QQ[1]
            for _ in 1:k; w = lmul(w, QQ[1, 0, -1]); end
            q = QQ[0]
            for i in 0:mk-1, j in 0:mk-1
                c = blk[Q*mk+i+1, Q*mk+j+1]
                c != 0 && (q = ladd(q, lmul(TC[i+1], TC[j+1]), 1, c))
            end
            p = ladd(p, lmul(w, q))
        end
    end
    while length(p) > 1 && p[end] == 0; pop!(p); end
    p
end
"TRIPLE3 in (u12, u13, u23) for the certificate types s1, s2, s3 (zero-based)"
function triple3_poly(Bk::Blocks, d3, s1, s2, s3)
    LC = legendre_coeffs(d3); TC = chebyshev(d3)
    out = P3()
    for (P, (Qa, va), (Qb, vb), vt) in ((s1, (s2, 1), (s3, 2), 3), (s2, (s1, 1), (s3, 3), 2), (s3, (s1, 2), (s2, 3), 1))
        for k in 0:d3
            mk = d3 - k + 1; blk = Bk[(:F,P,k)]
            core = P3()
            for i in 0:mk-1, j in 0:mk-1
                c = blk[Qa*mk+i+1, Qb*mk+j+1]
                c != 0 && (core = ThreePoint.padd(core, ThreePoint.pscale(ThreePoint.pmul(univariate(TC[i+1], 1), univariate(TC[j+1], 2)), c)))
            end
            isempty(core) && continue
            y = ThreePoint.pmul(core, phi_poly(k, LC))
            perm = (va, vb, vt)
            ren = P3()
            for (e, c) in y
                ee = [0, 0, 0]
                for src in 1:3; ee[perm[src]] += e[src]; end
                kk = (ee[1], ee[2], ee[3]); ren[kk] = get(ren, kk, QQ(0)) + c
            end
            out = ThreePoint.padd(out, ThreePoint.pscale(ren, 2))
        end
    end
    out
end
g_point(G::Vector{Group}, d3, i) = sum(point_term(g.Bk, d3, g.tmap[i]) for g in G)
function g_pair(G::Vector{Group}, d3, i, j)
    p = QQ[0]
    for g in G; p = ladd(p, pair3_poly(g.Bk, d3, g.tmap[i], g.tmap[j])); end
    while length(p) > 1 && p[end] == 0; pop!(p); end
    p
end
function g_triple(G::Vector{Group}, d3, i, j, k)
    out = P3()
    for g in G; out = ThreePoint.padd(out, triple3_poly(g.Bk, d3, g.tmap[i], g.tmap[j], g.tmap[k])); end
    out
end

# ---------------------------------------------------------------- (b) the pairs
"K(d, d', u) + Pu(u) - c2 <= Pi(d/2, d'/2, u) on the box of distances, -1 <= u <= a(d, d') (d <= d' if symmetric)"
function pair_box_check(A::Vector{Matrix{QQ}}, Pu::Vector{QQ}, c2::QQ, box, symmetric::Bool; depth_max=70, maxboxes=0)
    (alo, ahi), (blo, bhi) = box
    umax = amaxq(ahi, bhi)
    Du = max(D2, length(Pu) - 1)
    Ta = [bernstein(p, R2) for p in cheb_T((2alo - C1)/C2, 2*(ahi - alo)/C2, R2)]
    Tb = [bernstein(p, R2) for p in cheb_T((2blo - C1)/C2, 2*(bhi - blo)/C2, R2)]
    Us = [bernstein(p .// k, Du) for (k, p) in enumerate(cheb_U(QQ(-1), umax + 1, D2))]
    q = QQ[0]; pw = QQ[1]
    for c in Pu
        q = ladd(q, pw, 1, c); pw = lmul(pw, QQ[-1, umax + 1])
    end
    q[1] -= c2
    Bq = bernstein(q, Du)
    B0 = zeros(QQ, R2+1, R2+1, Du+1)
    for i in 1:R2+1, j in 1:R2+1, l in 1:Du+1; B0[i,j,l] = Bq[l]; end
    for k in 0:D2, a in 0:R2, b in 0:R2
        A[k+1][a+1,b+1] == 0 && continue
        for i in 0:R2
            ci = A[k+1][a+1,b+1]*Ta[a+1][i+1]
            ci == 0 && continue
            for j in 0:R2
                cij = ci*Tb[b+1][j+1]
                cij == 0 && continue
                for l in 0:Du; B0[i+1,j+1,l+1] += cij*Us[k+1][l+1]; end
            end
        end
    end
    Cmax = maximum(abs, B0)
    Bf = Float64.(B0)
    err = Float64(Cmax)*2.0^-53*(1 + depth_max*max(R2, Du))*1.01
    fa, wa, fb, wb, fu = Float64(alo), Float64(ahi - alo), Float64(blo), Float64(bhi - blo), Float64(umax + 1)
    corner(I, L) = (alo + (ahi - alo)*QQ(I[1] + 1, big(2)^L[1]), blo + (bhi - blo)*QQ(I[2] + 1, big(2)^L[2]),
                    -1 + (umax + 1)*QQ(I[3], big(2)^L[3]))
    stack = Tuple{NTuple{3,Int},NTuple{3,Int},Array{Float64,3}}[((0,0,0), (0,0,0), Bf)]
    closed = dropped = arbcalls = 0
    memo = Dict{Tuple{NTuple{3,Int},NTuple{3,Int}},Float64}()
    while !isempty(stack)
        I, L, B = pop!(stack)
        maximum(L) > depth_max && return (false, closed, dropped, arbcalls, "depth limit reached")
        sc = (2.0^(-L[1]), 2.0^(-L[2]), 2.0^(-L[3]))
        dhi = fa + wa*(I[1] + 1)*sc[1]; ehi = fb + wb*(I[2] + 1)*sc[2]; ulo = -1.0 + fu*I[3]*sc[3]
        drop = false
        if symmetric
            m = max(L[1], L[2])
            drop = (big(I[1]) << (m - L[1])) > (big(I[2] + 1) << (m - L[2]))
        end
        if !drop
            am = (dhi*dhi + ehi*ehi - 4)/(2*dhi*ehi)
            if ulo > am + 1e-12
                drop = true
            elseif ulo > am - 1e-12
                dh, eh, ul = corner(I, L)
                drop = ul > amaxq(dh, eh)
            end
        end
        if drop; dropped += 1; continue; end
        kub = nextfloat(maximum(B) + err)
        ok = kub <= 0
        if !ok && kub <= pair_float(dhi/2, ehi/2, ulo) - 1e-9
            v = get(memo, (I, L), nothing)
            if v === nothing
                v = Pi_lower(corner(I, L)...); memo[(I, L)] = v; arbcalls += 1
            end
            ok = kub <= v
        end
        if ok; closed += 1; continue; end
        (maxboxes > 0 && closed + dropped >= maxboxes) && return (true, closed, dropped, arbcalls, "stopped early (smoke test)")
        w = (wa*sc[1], wb*sc[2], fu*sc[3]/3)
        ax = w[1] >= w[2] && w[1] >= w[3] ? 1 : (w[2] >= w[3] ? 2 : 3)
        Lt, Rt = split_mid(B, ax)
        Ln = Base.setindex(L, L[ax] + 1, ax)
        push!(stack, (Base.setindex(I, 2I[ax], ax), Ln, Lt)); push!(stack, (Base.setindex(I, 2I[ax] + 1, ax), Ln, Rt))
    end
    (true, closed, dropped, arbcalls, "$closed boxes closed, $dropped set aside, Pi in Arb at $arbcalls corners")
end

# ---------------------------------------------------------------- (d) the triples
"P <= cmax on {bottoms_i <= x_i <= tops_i, Gram det >= 0}, x1 <= x2 if sym >= 1, x2 <= x3 if sym == 2"
function triple_box_check(P::P3, tops, sym::Int, cmax::Float64; label="", wmin=1e-6, batch=150000, bottoms=(-1.0,-1.0,-1.0),
                          mono=false, maxlevel=typemax(Int), io::IO=stdout, progress=true)
    Pp = IPoly(P); Dv = [IPoly(derivative(P, v)) for v in 1:3]
    HK = [(1,1),(1,2),(1,3),(2,2),(2,3),(3,3)]
    Hs = [IPoly(derivative(derivative(P, v), w)) for (v, w) in HK]
    det = IPoly(P3((0,0,0)=>QQ(1), (1,1,1)=>QQ(2), (2,0,0)=>QQ(-1), (0,2,0)=>QQ(-1), (0,0,2)=>QQ(-1)))
    W = Work()
    Box = Tuple{NTuple{3,Float64},NTuple{3,Float64}}
    cur = Box[(Tuple(Float64.(bottoms)), Tuple(Float64.(tops)))]
    n_done = 0; n_out = 0; worst = -Inf; level = 0; t0 = time()
    while !isempty(cur)
        nxt = Box[]
        for s in 1:batch:length(cur)
            chunk = view(cur, s:min(s + batch - 1, length(cur)))
            Abox = Box[]; Bbox = Box[]; maxw = -Inf; anyund = false
            for (L, Hh) in chunk
                alive = true
                sym >= 1 && (alive &= L[1] <= Hh[2])
                sym == 2 && (alive &= L[2] <= Hh[3])
                (_, dhi) = ieval(det, L, Hh, W)
                alive &= dhi >= 0
                if !alive; n_out += 1; continue; end
                Lt, Ht = L, Hh
                if mono
                    for v in 1:3
                        (gblo, gbhi) = ieval(Dv[v], L, Hh, W)
                        gblo > 0 && (Lt = Base.setindex(Lt, Hh[v], v))
                        gbhi < 0 && (Ht = Base.setindex(Ht, L[v], v))
                    end
                end
                c = ntuple(v -> (Lt[v] + Ht[v])/2, 3)
                r = ntuple(v -> max(up(Ht[v] - c[v]), up(c[v] - Lt[v])), 3)
                (_, p_c) = ieval(Pp, c, c, W)
                (dclo, _) = ieval(det, c, c, W)
                inside = dclo >= 0
                sym >= 1 && (inside &= c[1] <= c[2])
                sym == 2 && (inside &= c[2] <= c[3])
                if inside
                    worst = max(worst, p_c)
                    worst > cmax && return (false, @sprintf("value %.6e above %.6e", worst, cmax), n_done, worst, c)
                end
                first = 0.0; contrib = zeros(3)
                for v in 1:3
                    (gl, gh) = ieval(Dv[v], c, c, W)
                    g = max(abs(gl), abs(gh))
                    term = up(g*r[v]); first = up(first + term); contrib[v] = term
                end
                second = 0.0
                for (k, (v, w)) in enumerate(HK)
                    (hl, hh) = ieval(Hs[k], Lt, Ht, W)
                    habs = max(abs(hl), abs(hh))
                    if v == w
                        hdiag = mono ? max(hh, 0.0) : habs
                        term = up(hdiag*up(r[v]*r[v])); second = up(second + term); contrib[v] = up(contrib[v] + term)
                    else
                        term = up(up(2*habs)*up(r[v]*r[w])); second = up(second + term)
                        contrib[v] = up(contrib[v] + term); contrib[w] = up(contrib[w] + term)
                    end
                end
                upper = up(up(p_c + first) + up(second/2))
                if upper <= cmax; n_done += 1; continue; end
                anyund = true
                wd = (Hh[1] - L[1], Hh[2] - L[2], Hh[3] - L[3])
                maxw = max(maxw, wd[1], wd[2], wd[3])
                ax = contrib[1] >= contrib[2] && contrib[1] >= contrib[3] ? 1 : (contrib[2] >= contrib[3] ? 2 : 3)
                if mono && p_c > cmax
                    ax = wd[1] >= wd[2] && wd[1] >= wd[3] ? 1 : (wd[2] >= wd[3] ? 2 : 3)
                end
                mid = (L[ax] + Hh[ax])/2
                push!(Abox, (L, Base.setindex(Hh, mid, ax))); push!(Bbox, (Base.setindex(L, mid, ax), Hh))
            end
            (anyund && maxw < wmin) && return (false, "undecided box below width $wmin", n_done, worst, nothing)
            append!(nxt, Abox); append!(nxt, Bbox)
        end
        level += 1; cur = nxt
        if progress && level % 10 == 0
            println(io, "      $label level $level: $n_done boxes verified, $(length(cur)) to bisect [$(round(Int, time() - t0)) s]"); flush(io)
        end
        (level >= maxlevel && !isempty(cur)) && return (true, "stopped at level $level with $(length(cur)) boxes left (smoke test), $n_done boxes verified", n_done, worst, nothing)
    end
    (true, @sprintf("%d boxes, largest centre value %.4e", n_done, worst), n_done, worst, nothing)
end
"the variables of the branch and bound, as positions of (u12, u13, u23), and the symmetry used"
function triple_layout(labels)
    s1, s2, s3 = labels
    s1 == s2 == s3 && return ((1,2,3), 2)
    s1 == s2 && return ((2,3,1), 1)
    s2 == s3 && return ((1,2,3), 1)
    ((1,2,3), 0)
end

# ---------------------------------------------------------------- thresholds (floating point; proved afterwards)
linspace(a, b, n) = n == 1 ? [a] : (v = [a + (b - a)*i/(n - 1) for i in 0:n-1]; v[end] = b; v)
function kfloat(Af, d, e, u)
    x = (2d - Float64(C1))/Float64(C2); y = (2e - Float64(C1))/Float64(C2)
    Tx = zeros(R2+1); Ty = zeros(R2+1); Tx[1] = 1; Ty[1] = 1; Tx[2] = x; Ty[2] = y
    for k in 3:R2+1; Tx[k] = 2x*Tx[k-1] - Tx[k-2]; Ty[k] = 2y*Ty[k-1] - Ty[k-2]; end
    Um2 = 1.0; Um1 = 2u; s = dot(Tx, Af[1], Ty) + Um1/2*dot(Tx, Af[2], Ty)
    for k in 2:D2
        Uk = 2u*Um1 - Um2; s += Uk/(k+1)*dot(Tx, Af[k+1], Ty); Um2 = Um1; Um1 = Uk
    end
    s
end
polyval(p, u) = (s = 0.0; for c in Iterators.reverse(p); s = s*u + c; end; s)
"pattern search (all 3^n - 1 directions) for a local maximum of f on the set where feas holds;
`repair` maps a trial point to a list of candidates on the feasible set (so that the search can
follow a curved boundary and its intersections with the faces of the box)"
function ascend(f, feas, x0; h0=0.02, hmin=1e-11, repair=y -> (y,))
    x = copy(x0); fx = f(x); h = h0; n = length(x)
    dirs = [collect(d) for d in Iterators.product(ntuple(_ -> (-1, 0, 1), n)...) if any(!=(0), d)]
    while h > hmin
        moved = false
        for d in dirs, y in repair(x .+ h .* d)
            feas(y) || continue
            fy = f(y)
            if fy > fx; x = y; fx = fy; moved = true; end
        end
        moved || (h /= 2)
    end
    fx
end
function pair_setup(Af, Pu, box, margin, rng)
    (a0, a1), (b0, b1) = Float64.(box[1]), Float64.(box[2])
    pc = Float64.(Pu)
    f(x) = kfloat(Af, x[1], x[2], x[3]) + polyval(pc, x[3]) - pair_float(x[1]/2, x[2]/2, x[3])
    top(p, q) = (p*p + q*q - 4)/(2p*q)
    feas(x) = a0 <= x[1] <= a1 && b0 <= x[2] <= b1 && -1 <= x[3] <= top(x[1], x[2])
    pts = Vector{Float64}[]
    n = 200000; k6 = n ÷ 6
    for i in 1:n
        p = a0 + (a1 - a0)*rand(rng); q = b0 + (b1 - b0)*rand(rng)
        i <= k6 ? (p = a0) : i <= 2k6 ? (q = b0) : i <= 3k6 ? (p = a1) : i <= 4k6 ? (q = b1) : nothing
        t = top(p, q)
        u = (4k6 < i <= 5k6) ? t - 3e-3*rand(rng) : -1 + (t + 1)*sqrt(rand(rng))
        push!(pts, [p, q, u])
    end
    sv = vcat(0.0, linspace(0.0, 1.0, 80) .^ 2, 1.0)
    for gd in linspace(a0, a1, 25), ge in linspace(b0, b1, 25), s in sv
        push!(pts, [gd, ge, -1 + (top(gd, ge) + 1)*(1 - s)])
    end
    vals = [f(x) for x in pts]
    vmax = maximum(vals)
    best = -Inf
    for i in partialsortperm(vals, 1:40; rev=true)
        feas(pts[i]) && (best = max(best, ascend(f, feas, pts[i]; repair = y -> ([clamp(y[1], a0, a1), clamp(y[2], b0, b1), clamp(y[3], -1.0, top(clamp(y[1], a0, a1), clamp(y[2], b0, b1)))],))))
    end
    (above(max(vmax, best) + margin), vmax)
end
struct FPoly
    exps::Vector{NTuple{3,Int}}; cf::Vector{Float64}; deg::NTuple{3,Int}
end
FPoly(P::P3) = (it = sort(collect(P); by=first); FPoly([e for (e, _) in it], [Float64(c) for (_, c) in it],
                                                     ntuple(v -> maximum(e[v] for (e, _) in it; init=0), 3)))
function (p::FPoly)(x)
    pw = ntuple(v -> (a = ones(p.deg[v] + 1); for k in 1:p.deg[v]; a[k+1] = a[k]*x[v]; end; a), 3)
    s = 0.0
    @inbounds for i in eachindex(p.exps)
        (a, b, c) = p.exps[i]; s += p.cf[i]*pw[1][a+1]*pw[2][b+1]*pw[3][c+1]
    end
    s
end
gramdet(x) = 1 + 2x[1]*x[2]*x[3] - x[1]^2 - x[2]^2 - x[3]^2
"move a trial point into the box [Bf, Tf], then one coordinate into the interval that the Gram
determinant allows given the other two (each of the three choices is a candidate)"
function gram_repair(Bf, Tf)
    function repair(y)
        y = [clamp(y[i], Bf[i], Tf[i]) for i in 1:3]
        gramdet(y) >= 0 && return (y,)
        out = Vector{Float64}[]
        for (i, j, k) in ((1, 2, 3), (1, 3, 2), (2, 3, 1))
            sroot = sqrt(max((1 - y[i]^2)*(1 - y[j]^2), 0.0))*(1 - 1e-13)
            lo = max(y[i]*y[j] - sroot, Bf[k]); hi = min(y[i]*y[j] + sroot, Tf[k])
            lo <= hi || continue
            z = copy(y); z[k] = clamp(y[k], lo, hi)
            push!(out, z)
        end
        out
    end
end
"local ascent of P (in the variables of the branch and bound) from the point x0"
function triple_refine_from(P::P3, tops, bottoms, x0)
    Tf = Float64.(collect(tops)); Bf = Float64.(collect(bottoms)); f = FPoly(P)
    feas(x) = all(Bf[i] <= x[i] <= Tf[i] for i in 1:3) && gramdet(x) >= 0
    ascend(f, feas, collect(Float64, x0); repair=gram_repair(Bf, Tf))
end
"sample and refine P on {B <= x <= T, Gram det >= 0}; returns (threshold, sampled max, #samples, refined max)"
function triple_sample(P::P3, T, B, margin, rng; nrand=200000)
    Tf = Float64.(collect(T)); Bf = Float64.(collect(B))
    f = FPoly(P)
    feas(x) = all(Bf[i] <= x[i] <= Tf[i] for i in 1:3) && gramdet(x) >= 0
    repair = gram_repair(Bf, Tf)
    pts = Vector{Float64}[]
    for dim in (4, 3, 2)
        got = 0
        for _ in 1:30nrand
            X = randn(rng, 3, dim)
            for i in 1:3; X[i, :] ./= norm(X[i, :]); end
            g = [dot(X[1, :], X[2, :]), dot(X[1, :], X[3, :]), dot(X[2, :], X[3, :])]
            if all(Bf[i] <= g[i] <= Tf[i] for i in 1:3)
                push!(pts, g); got += 1
                got >= nrand && break
            end
        end
    end
    ax = [vcat(Bf[i], linspace(Bf[i], Tf[i], 40), Tf[i]) for i in 1:3]
    for x in ax[1], y in ax[2], z in ax[3]
        g = [x, y, z]
        gramdet(g) >= 0 && push!(pts, g)
    end
    vals = [f(x) for x in pts]
    vmax = maximum(vals)
    best = -Inf
    for i in partialsortperm(vals, 1:min(40, length(vals)); rev=true)
        best = max(best, ascend(f, feas, pts[i]; repair=repair))
    end
    (above(max(vmax, best) + margin), vmax, length(pts), best)
end
"lay P(u12, u13, u23) out for the branch and bound: returns (Pr, tops, sym)"
function layout(P::P3, T, labels)
    perm, sym = triple_layout(labels)
    pos(src) = findfirst(==(src), perm)
    tops = zeros(3)
    for src in 1:3; tops[pos(src)] = Float64(above(T[src])); end
    Pr = P3()
    for (e, co) in P
        ee = [0, 0, 0]
        for src in 1:3; ee[pos(src)] += e[src]; end
        k = (ee[1], ee[2], ee[3]); Pr[k] = get(Pr, k, QQ(0)) + co
    end
    (Pr, tops, sym)
end

# ---------------------------------------------------------------- main
"""
    run(certpath, d3; case=case30() | path, mg2, mg3, mgm, options...) -> Report

`case` is a Case (case30() for combo30_check.py) or the path of a case JSON file
(combo_case_check.py).  With `mixtype=true` the certificate holds x3g/gtyp
(combo_mixtype_check.py; skip_empty is then set).  Bins in refine_* / region are numbered
from 1 as in the Python options.  `maxlevel` / `maxboxes` cut the branch and bounds short
for smoke tests, `bb=false` skips them altogether (the PRECHECK of the Python) and
`only=["FFF", ...]` runs the branch and bound of the named triple kinds alone (ONLY); the
result is then not a proof of the case and the last line says PARTIAL.  `triples=false` runs
everything but the triple branch and bounds.  The thresholds do not depend on these switches
(the sampling sequence is fixed by `seed`), so a run with `triples=false` and runs with `only`
that together name every triple kind prove the case, as with ONLY in the Python.
"""
function run(certpath::AbstractString, d3::Int, io::IO=stdout; case=case30(), mg2=2e-5, mg3=2e-6, mgm=2e-6,
             border=false, level=nothing, triple_extra=Dict{String,Float64}(), pair_margin=Dict{String,Float64}(),
             skip_empty=false, refine_triples=Vector{Vector{Int}}(), refine_pairs=Vector{Vector{Int}}(),
             skip_covered=false, region=Tuple{Vector{Int},QQ,QQ}[], extra_cases=String[], mono=false, mixtype=false,
             maxlevel=typemax(Int), maxboxes=0, seed=7, expect_bound=nothing, expect_tol=2e-5, bb=true, only=String[], retries=3, triples=true)
    t0 = time(); R = Report("combo: $(basename(certpath))")
    C = case isa Case ? case : load_case(case)
    smoke = maxlevel != typemax(Int) || maxboxes > 0 || !bb || !isempty(only) || !triples
    mixtype && (skip_empty = true)
    Z = npzread(certpath, mixtype ? ["A", "z", "t", "m", "x3g", "gtyp"] : ["A", "z", "t", "m", "x3"])   # the other arrays (samples, labels) are not used
    rng = Xoshiro(seed)
    target = level === nothing ? 9*Arb(π; prec=PREC[])^2/8 - 8 : A_(fr(level))
    tname = level === nothing ? "9 pi^2/8 - 8" : level
    T = C.types; nT = length(T); tidx = Dict(s => i for (i, s) in enumerate(T))
    println(io, "statement (C) at $(C.M) centres, case $(C.name): certificate $(basename(certpath)), two-point degree $D2/$R2, three-point degree $d3; $(length(C.counts)) count vectors"); flush(io)
    check!(R, io, "dmax exceeds sqrt 6", DMAX^2 > 6)
    # (a) positivity
    Af = Float64.(Z["A"]); nA = size(Af, 2)
    sh30 = QQ(1, big(2)^30)
    local A0, z, t, Zm
    if border
        Bd = zeros(nA+1, nA+1)
        Bd[1:nA, 1:nA] = Af[1, :, :]; Bd[1:nA, nA+1] = Z["z"]; Bd[nA+1, 1:nA] = Z["z"]; Bd[nA+1, nA+1] = Z["t"][]
        Zm = psd_exact(Bd, sh30)
        A0 = Zm[1:nA, 1:nA]; z = Zm[1:nA, nA+1]; t = Zm[nA+1, nA+1]
    else
        A0 = psd_exact(Af[1, :, :], sh30); z = dyad.(Float64.(Z["z"]))
    end
    A = vcat([A0], [psd_exact(Af[k+1, :, :], sh30) for k in 1:D2])
    good = all(ldl_psd(A[k+1])[1] for k in 1:D2)
    g0, piv = ldl_psd(A0)
    if !border
        wv = solve_exact(A0, z); tq = sum(z .* wv)
        t = QQ(cld(numerator(tq)*big(2)^48, denominator(tq)), big(2)^48)
        Zm = [A0 z; permutedims(z) t]
    end
    gz, _ = ldl_psd(Zm)
    check!(R, io, "A_1..A_D and [[A_0, z], [z^T, t]] positive semidefinite (exact LDL^T)", good && g0 && all(p > 0 for p in piv) && gz;
           got=@sprintf("t = %.9f, float t %.9f", Float64(t), Z["t"][]))
    local groups::Vector{Group}
    if mixtype
        GT = Int.(Z["gtyp"]); NG = size(GT, 1)
        @assert all(b !== nothing for b in C.bintype) && [tidx[b] for b in C.bintype] == collect(1:nT) && size(GT, 2) == nT "the case must give each bin a type of its own, in bin order"
        groups = [Group(Blocks(k => psd_exact(Mx) for (k, Mx) in blocks_from_x3(Float64.(Z["x3g"][g, :]), d3, 3)), [GT[g, i] for i in 1:nT]) for g in 1:NG]
        println(io, "  mixed typing: $NG groups, typings ", [join("ABF"[GT[g, i]+1] for i in 1:nT) for g in 1:NG])
    else
        groups = [Group(Blocks(k => psd_exact(Mx) for (k, Mx) in blocks_from_x3(Float64.(vec(Z["x3"])), d3, nT)), collect(0:nT-1))]
    end
    nblocks = sum(length(g.Bk) for g in groups)
    check!(R, io, "three-point blocks positive definite (exact LDL^T)", all(ldl_positive(Mx) for g in groups for Mx in values(g.Bk)); got="$nblocks blocks")
    pt = Dict(s => g_point(groups, d3, i) for (i, s) in enumerate(T))
    println(io, "  point terms: ", join(["$s: " * @sprintf("%.6f", Float64(pt[s])) for s in T], ", ")); flush(io)
    # (b) pairs
    Af2 = [Float64.(a) for a in A]
    c2 = Dict{Tuple{String,String},QQ}(); c2order = Tuple{String,String}[]
    pjobs = Any[]
    for (s, tt) in combos2(T)
        (skip_empty && all(Npair(tcounts(C, nb), s, tt) == 0 for nb in C.counts)) && continue
        Pu = g_pair(groups, d3, tidx[s], tidx[tt])
        box = (C.trange[s], C.trange[tt])
        c2[(s, tt)], vmax = pair_setup(Af2, Pu, box, get(pair_margin, s*tt, mg2), rng)
        push!(c2order, (s, tt))
        @printf(io, "  pair %s%s: float largest K + PAIR3 - Pi %.6e; threshold %.6e\n", s, tt, vmax, Float64(c2[(s, tt)])); flush(io)
        push!(pjobs, ((s, tt), "K + PAIR3_$s$tt <= Pi + c2 on the admissible pairs", Pu, c2[(s, tt)], box, s == tt))
    end
    c2r = Dict{Vector{Int},QQ}(); c2rorder = Vector{Int}[]
    prng = Xoshiro(13)
    for bs0 in refine_pairs
        bs = sort(bs0); s, tt = C.bintype[bs[1]], C.bintype[bs[2]]
        box = ((C.edges[bs[1]], C.edges[bs[1]+1]), (C.edges[bs[2]], C.edges[bs[2]+1]))
        Pu = g_pair(groups, d3, tidx[s], tidx[tt])
        c2r[bs], vmax = pair_setup(Af2, Pu, box, get(pair_margin, s*tt, mg2), prng); push!(c2rorder, bs)
        name = "bins " * join(bs, "-")
        @printf(io, "  pair %s%s on %s: float largest %.6e; threshold %.6e against %.6e for the kind\n", s, tt, name, vmax, Float64(c2r[bs]), Float64(c2[(s, tt)]))
        push!(pjobs, (nothing, "K + PAIR3_$s$tt <= Pi + its threshold on the pairs of $name", Pu, c2r[bs], box, bs[1] == bs[2]))
    end
    covered2(st, vecs) = (rs = [bs for bs in c2rorder if (C.bintype[bs[1]], C.bintype[bs[2]]) == st];
                          (!isempty(rs) && all(Npair(tcounts(C, nb), st...) == sum(Nbins(nb, bs) for bs in rs) for nb in vecs), rs))
    skip2 = Set{Tuple{String,String}}()
    if skip_covered
        for st in c2order
            cov, rs = covered2(st, C.counts)
            if cov
                push!(skip2, st); c2[st] = max(c2[st], maximum(c2r[bs] for bs in rs))
                println(io, "  pair $(st[1])$(st[2]): covered by the listed pairs of bins in every count vector; its branch and bound is skipped")
            end
        end
    end
    for (st, label, Pu, thr, box, symm) in pjobs
        (st !== nothing && st in skip2) && continue
        (bb && isempty(only)) || continue
        ok, cl, dr, ac, msg = pair_box_check(A, Pu, thr, box, symm; maxboxes=maxboxes)
        check!(R, io, label, ok; got=msg)
    end
    # (c) bins
    mf = Float64.(vec(Z["m"]))
    m = QQ[]; im = 0
    for b in C.bintype
        if b === nothing; push!(m, QQ(10^6)); else; im += 1; push!(m, above(mf[im] + mgm)); end
    end
    @assert im == length(mf)
    if bb && isempty(only)
        okb, nint, room = Radial.bracket_check(A, z, D2, R2, C1, C2, DMAX, C.edges, m)
        check!(R, io, "f <= m_b on every bin", okb; got=@sprintf("%d intervals, least room %.2e", nint, room))
    end
    # (d) triples
    c3 = Dict{NTuple{3,String},QQ}(); c3order = NTuple{3,String}[]
    tjobs = Any[]     # (kind or nothing, label, P, tops, sym, cmax, bottoms)
    for combo in combos3(T)
        all(Ntriple(tcounts(C, nb), collect(combo)) == 0 for nb in C.counts) && continue
        P = g_triple(groups, d3, tidx[combo[1]], tidx[combo[2]], tidx[combo[3]])
        Tt = (tmax(C, combo[1], combo[2]), tmax(C, combo[1], combo[3]), tmax(C, combo[2], combo[3]))
        kind = join(combo)
        thr, vmax, ns, best = triple_sample(P, Tt, (-1, -1, -1), mg3 + get(triple_extra, kind, 0.0), rng)
        Pr, tops, sym = layout(P, Tt, combo)
        c3[combo] = thr; push!(c3order, combo)
        @printf(io, "  triple %s: float largest %.6e on %d samples (refined %.6e); threshold %.6e; %d monomials, symmetry %d\n", kind, vmax, ns, best, Float64(thr), length(P), sym); flush(io)
        push!(tjobs, (kind, kind, Pr, tops, sym, thr, (-1.0, -1.0, -1.0), v -> (c3[combo] = v), mg3 + get(triple_extra, kind, 0.0)))
    end
    c3r = Dict{Vector{Int},QQ}(); c3rorder = Vector{Int}[]
    rrng = Xoshiro(11)
    for bs0 in refine_triples
        bs = sort(bs0); combo = Tuple(String(C.bintype[b]) for b in bs)
        @assert haskey(c3, combo)
        hi = [C.edges[b+1] for b in bs]
        Tt = (amaxq(hi[1], hi[2]), amaxq(hi[1], hi[3]), amaxq(hi[2], hi[3]))
        P = g_triple(groups, d3, tidx[combo[1]], tidx[combo[2]], tidx[combo[3]])
        ex = join(bs, "-"); name = "bins " * ex
        thr, vmax, ns, best = triple_sample(P, Tt, (-1, -1, -1), mg3 + get(triple_extra, ex, 0.0), rrng)
        Pr, tops, sym = layout(P, Tt, bs)
        c3r[bs] = thr; push!(c3rorder, bs)
        @printf(io, "  triple %s on %s: float largest %.6e on %d samples; threshold %.6e against %.6e for the kind; symmetry %d\n", join(combo), name, vmax, ns, Float64(thr), Float64(c3[combo]), sym)
        push!(tjobs, (nothing, join(combo) * " on " * name, Pr, tops, sym, thr, (-1.0, -1.0, -1.0), v -> (c3r[bs] = v), mg3 + get(triple_extra, ex, 0.0)))
    end
    c3g = Dict{Vector{Int},Tuple{Int,Int,QQ,QQ,Int,Int}}()
    grng = Xoshiro(13)
    PAIRIDX = Dict((1,2) => 1, (1,3) => 2, (2,3) => 3)
    for (bs0, s1, s2) in region
        bs = sort(bs0)
        @assert haskey(c3r, bs) && !haskey(c3g, bs)
        i, j, k = (bs[1] == bs[2] != bs[3]) ? (1, 2, 3) : (2, 3, 1)
        @assert bs[i] == bs[j] != bs[k]
        xi, yi = PAIRIDX[minmax(i, k)], PAIRIDX[minmax(j, k)]
        b, p = bs[i], bs[k]
        combo = Tuple(String(C.bintype[x]) for x in bs); name = "bins " * join(bs, "-")
        P = g_triple(groups, d3, tidx[combo[1]], tidx[combo[2]], tidx[combo[3]])
        swap(e) = begin v = collect(e); v[xi], v[yi] = v[yi], v[xi]; (v[1], v[2], v[3]) end
        check!(R, io, "TRIPLE3 on $name symmetric under exchange of the two centres of bin $b (exact coefficients)", all(get(P, swap(e), QQ(0)) == co for (e, co) in P))
        hi = [C.edges[x+1] for x in bs]
        Tt = [amaxq(hi[1], hi[2]), amaxq(hi[1], hi[3]), amaxq(hi[2], hi[3])]
        c0 = amaxq(C.edges[b+1], C.edges[b+1]); Tx = Tt[xi]
        @assert s1 < 0 < s2 < Tx < 1 && s1*s1 > c0 && s2*s2 < c0
        Ka = Int(floor((1 - c0)/(s1*s1 - c0)))
        g = (c0 - s2*s2)/(1 - Tx*Tx)
        Kt = 2
        while (1 + g)/2 >= (1 - QQ(2, Kt + 1))^2; Kt += 1; end
        Bm = QQ[-1, -1, -1]; Tm = copy(Tt); Bm[xi] = below(s1); Tm[xi] = above(s2)
        Bx = QQ[-1, -1, -1]; Txx = copy(Tt); Txx[xi] = above(s1); Bx[yi] = below(s2)
        ex = join(bs, "-")
        cm, vm, nm, bm = triple_sample(P, Tm, Bm, mg3 + get(triple_extra, ex * "/mid", 0.0), grng)
        cx, vx, nx, bx = triple_sample(P, Txx, Bx, mg3 + get(triple_extra, ex * "/mix", 0.0), grng)
        c3g[bs] = (b, p, cm, cx, Ka, Kt)
        @printf(io, "  triple %s on %s split at [%s, %s]: middle threshold %.6e, mixed threshold %.6e, against %.6e for the bins; Ka = %d, Kt = %d\n",
                join(combo), name, string(s1), string(s2), Float64(cm), Float64(cx), Float64(c3r[bs]), Ka, Kt)
        push!(tjobs, (nothing, join(combo) * " on " * name * " middle", P, Float64.(Tm), 0, cm, Tuple(Float64.(Bm)),
                      v -> (c3g[bs] = Base.setindex(c3g[bs], v, 3)), mg3 + get(triple_extra, ex * "/mid", 0.0)))
        push!(tjobs, (nothing, join(combo) * " on " * name * " mixed", P, Float64.(Txx), 0, cx, Tuple(Float64.(Bx)),
                      v -> (c3g[bs] = Base.setindex(c3g[bs], v, 4)), mg3 + get(triple_extra, ex * "/mix", 0.0)))
    end
    covered3(combo, vecs) = (rs = [bs for bs in c3rorder if Tuple(String(C.bintype[b]) for b in bs) == combo];
                             (!isempty(rs) && all(Ntriple(tcounts(C, nb), collect(combo)) == sum(Nbins(nb, bs) for bs in rs) for nb in vecs), rs))
    skip3 = Set{String}()
    if skip_covered
        for combo in c3order
            cov, rs = covered3(combo, C.counts)
            if cov
                push!(skip3, join(combo)); c3[combo] = max(c3[combo], maximum(c3r[bs] for bs in rs))
                println(io, "  triple $(join(combo)): covered by the listed combinations of bins in every count vector; its branch and bound is skipped")
            end
        end
    end
    nbox = 0
    for (kind, label, Pr, tops, sym, thr0, bottoms, set!, margin) in tjobs
        (kind !== nothing && kind in skip3) && continue
        (bb && triples && (isempty(only) || kind in only)) || continue
        thr = thr0; local ok, msg, nd
        for attempt in 1:(retries + 1)
            ok, msg, nd, wv, wc = triple_box_check(Pr, tops, sym, Float64(thr); label=label, bottoms=bottoms, mono=mono, maxlevel=maxlevel, io=io)
            (ok || wc === nothing || attempt > retries) && break
            # the sampling missed a peak: the branch and bound met a point of the domain above the
            # threshold.  Climb from it, raise the threshold of this inequality and prove it again;
            # the bound below uses the raised threshold.
            best = max(wv, triple_refine_from(Pr, tops, bottoms, wc))
            thr = above(best + max(margin, mg3)); set!(thr)
            @printf(io, "  %s: %s; the sampled threshold was too low, raised to %.9e (largest value found %.9e) and checked again\n", label, msg, Float64(thr), best)
        end
        check!(R, io, occursin(" on ", label) ? "TRIPLE3_$label <= its threshold on the admissible triples" : "TRIPLE3_$label <= c3 on the admissible triples", ok; got=msg)
        nbox += nd
    end
    for combo in c3order          # thresholds raised by a retry: keep the covered kinds at the largest of their bins
        join(combo) in skip3 && (c3[combo] = max(c3[combo], maximum(c3r[bs] for bs in covered3(combo, C.counts)[2])))
    end
    println(io, "  triple boxes in all: $nbox"); flush(io)
    # (e) the bound, over every count vector (exact)
    kindof2(bs) = (String(C.bintype[bs[1]]), String(C.bintype[bs[2]]))
    kindof3(bs) = Tuple(String(C.bintype[b]) for b in bs)
    refined_gain(nb) = sum(Nbins(nb, bs)*max(c3[kindof3(bs)] - c3r[bs], 0) for bs in c3rorder; init=QQ(0)) +
                       sum(Nbins(nb, bs)*max(c2[kindof2(bs)] - c2r[bs], 0) for bs in c2rorder; init=QQ(0))
    function region_gain(nb)
        out = QQ(0)
        for (bs, (b, p, cm, cx, Ka, Kt)) in c3g
            n = nb[b]; ch = c3r[bs]
            (nb[p] == 0 || n < 2) && continue
            cm = min(cm, ch); cx = min(cx, ch); tot = binomial(n, 2); best = nothing
            for ka in 0:min(Ka, n), kt in 0:min(Kt, n - ka)
                hot = binomial(ka, 2) + binomial(kt, 2); mix = ka*kt
                v = hot*ch + mix*cx + (tot - hot - mix)*cm
                (best === nothing || v > best) && (best = v)
            end
            out += nb[p]*(tot*ch - best)
        end
        out
    end
    cache = Dict{Vector{Int},QQ}()
    function bound_over(vecs)
        worst = nothing; wvec = nothing
        for nb in vecs
            tc = tcounts(C, nb); key = [tc[s] for s in T]
            base = get!(cache, key) do
                t/2 + sum(Npair(tc, s, tt)*c2[(s, tt)] for (s, tt) in c2order; init=QQ(0)) +
                    sum(Ntriple(tc, collect(cb))*c3[cb] for cb in c3order; init=QQ(0))
            end
            val = base + sum(nb[b]*(m[b] + pt[C.bintype[b]]) for b in eachindex(nb) if nb[b] != 0; init=QQ(0))
            val -= refined_gain(nb) + region_gain(nb)
            (worst === nothing || val > worst) && ((worst, wvec) = (val, nb))
        end
        (worst, wvec)
    end
    worst, wvec = bound_over(C.counts)
    @printf(io, "  %d count vectors, %d type-count vectors; largest bound %.6f at %s\n", length(C.counts), length(cache), Float64(worst), string(Tuple(wvec)))
    for extra in extra_cases
        E = JSON.parsefile(extra)
        @assert (E["M"], E["types"], E["trange"], E["bins"]) == (C.raw["M"], C.raw["types"], C.raw["trange"], C.raw["bins"]) extra
        vecs = count_vectors(E["M"], length(E["bins"]), E["constraints"], get(C.raw, "exclude", []))
        if !isempty(skip2) || !isempty(skip3)
            check!(R, io, "the skipped kinds are covered by the listed bins in the count vectors of $(basename(extra))",
                   all(covered2(st, vecs)[1] for st in skip2) && all(covered3(cb, vecs)[1] for cb in c3order if join(cb) in skip3))
        end
        ew, ev = bound_over(vecs)
        check!(R, io, "with these thresholds, largest bound over the $(length(vecs)) count vectors of $(basename(extra)) below $tname", A_(ew) < target;
               got=@sprintf("%.6f at %s", Float64(ew), string(Tuple(ev))))
    end
    check!(R, io, "largest bound over the count vectors below $tname", A_(worst) < target;
           got=@sprintf("%.6f < %s", Float64(worst), string(target)))
    if expect_bound !== nothing
        check!(R, io, "bound agrees with the Python log (thresholds are sampled independently, tolerance $expect_tol)", abs(Float64(worst) - expect_bound) <= expect_tol;
               got=@sprintf("%.6f", Float64(worst)), expected=@sprintf("%.6f", expect_bound))
    end
    println(io, (Common.allpass(R) ? (smoke ? "PARTIAL (branch and bounds skipped, restricted or cut short; not a proof of the case)" : "PASS") : "FAILED"),
            @sprintf(": case %s, U(Y) <= %.6f against %s [%.0f s]", C.name, Float64(worst), tname, time() - t0))
    R
end

end # module
