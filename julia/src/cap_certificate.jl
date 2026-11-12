# Port of cap_certificate/cap_inequality_certificate.py (sec:cap-theorem), exact steps
# (A), (B), (D), (E), (F), (F').  The floating-point cross-checks (C), (G) of the Python
# script are not relied on by the proof and are not ported; in exchange vol(Q)=25/3 and
# vol(24-cell)=8, which the Python script obtains from a floating-point convex hull, are
# computed here exactly (barycentric triangulation over the face lattice, rational determinants).
module CapCertificate
using ..Common
using ..Common: QQ
using LinearAlgebra: det

dotq(a,b) = sum(a[i]*b[i] for i in eachindex(a))

function integer_roots()
    out = NTuple{4,Int}[]
    for i in 1:4, j in i+1:4, si in (1,-1), sj in (1,-1)
        v = [0,0,0,0]; v[i]=si; v[j]=sj; push!(out, Tuple(v))
    end
    out
end

"exact rank of a rational matrix"
function qrank(M::Matrix{QQ})
    A = copy(M); r = 0; (m,n) = size(A)
    for c in 1:n
        p = findfirst(i->A[i,c]!=0, r+1:m)
        p === nothing && continue
        p += r; r += 1
        A[r,:], A[p,:] = A[p,:], A[r,:]
        for i in r+1:m
            if A[i,c] != 0
                f = A[i,c]/A[r,c]; A[i,:] .-= f .* A[r,:]
            end
        end
        r == m && break
    end
    r
end
"solve a square rational system exactly; nothing if singular"
function qsolve(M::Matrix{QQ}, b::Vector{QQ})
    n = size(M,1); A = hcat(copy(M), b)
    for c in 1:n
        p = findfirst(i->A[i,c]!=0, c:n)
        p === nothing && return nothing
        p += c-1
        A[c,:], A[p,:] = A[p,:], A[c,:]
        A[c,:] ./= A[c,c]
        for i in 1:n
            i==c && continue
            A[i,c] != 0 && (A[i,:] .-= A[i,c] .* A[c,:])
        end
    end
    A[:,n+1]
end
function qdet(M::Matrix{QQ})
    A = copy(M); n = size(A,1); d = QQ(1)
    for c in 1:n
        p = findfirst(i->A[i,c]!=0, c:n)
        p === nothing && return QQ(0)
        p += c-1
        if p != c; A[c,:], A[p,:] = A[p,:], A[c,:]; d = -d; end
        d *= A[c,c]
        for i in c+1:n
            A[i,c] != 0 && (A[i,:] .-= (A[i,c]/A[c,c]) .* A[c,:])
        end
    end
    d
end

"vertices of {z : <z,n_i> <= 1} by exact enumeration"
function vertices_of(normals)
    V = Vector{Vector{QQ}}()
    m = length(normals)
    for a in 1:m, b in a+1:m, c in b+1:m, d in c+1:m
        M = QQ[normals[i][k] for i in (a,b,c,d), k in 1:4]
        z = qsolve(M, QQ[1,1,1,1])
        z === nothing && continue
        all(dotq(o,z) <= 1 for o in normals) || continue
        z in V || push!(V, z)
    end
    V
end
affdim(pts) = length(pts) <= 1 ? length(pts)-1 : qrank(QQ[pts[i][k]-pts[1][k] for i in 2:length(pts), k in 1:4])
centroid(pts) = sum(pts) .// length(pts)

"exact volume of {z : <z,n_i> <= 1} (origin interior) from its vertices"
function exact_volume(normals, V)
    facets = Vector{Vector{Int}}()
    for n in normals
        idx = [i for i in eachindex(V) if dotq(n,V[i])==1]
        (length(idx) >= 4 && affdim(V[idx])==3) && !(idx in facets) && push!(facets, idx)
    end
    vol = QQ(0)
    for F in facets
        cF = centroid(V[F])
        ridges = Set{Vector{Int}}()
        for G0 in facets
            G0 === F && continue
            G = intersect(F,G0)
            length(G) >= 3 && affdim(V[G])==2 && push!(ridges, G)
        end
        for G in ridges
            cG = centroid(V[G])
            edges = Set{Vector{Int}}()
            for H0 in facets
                E = intersect(G,H0)
                length(E) == 2 && push!(edges, E)
            end
            for E in edges
                cE = centroid(V[E])
                for v in E
                    M = QQ[x[k] for x in (cF,cG,cE,V[v]), k in 1:4]
                    vol += abs(qdet(M))//24
                end
            end
        end
    end
    vol
end

# Laurent polynomials in t = sqrt(a), exact
const LP = Dict{Int,QQ}
dda(p::LP) = LP(k-2 => c*k//2 for (k,c) in p if k != 0)      # d/da = (1/(2t)) d/dt
evalp(p::LP, t::QQ) = sum(c*t^k for (k,c) in p; init=QQ(0))
clean(p::LP) = LP(k=>c for (k,c) in p if c != 0)

# dense polynomials with rational coefficients (index 1 = constant term)
const P = Vector{QQ}
ptrim(p::P) = (while length(p) > 1 && p[end]==0; p = p[1:end-1]; end; p)
function pmul(a::P,b::P); c = zeros(QQ, length(a)+length(b)-1); for i in eachindex(a), j in eachindex(b); c[i+j-1] += a[i]*b[j]; end; ptrim(c); end
function padd(a::P,b::P); n=max(length(a),length(b)); c=zeros(QQ,n); c[1:length(a)].+=a; c[1:length(b)].+=b; ptrim(c); end
pscale(a::P,s) = ptrim(a .* QQ(s))
ppow(a::P,n) = n==0 ? QQ[1] : pmul(a, ppow(a,n-1))
peval(a::P,x::QQ) = foldr((c,acc)->c+x*acc, a; init=QQ(0))
pder(a::P) = length(a)==1 ? QQ[0] : QQ[(i-1)*a[i] for i in 2:length(a)]
function pdivrem(a::P,b::P)
    a = copy(a); q = zeros(QQ, max(length(a)-length(b)+1,1))
    while length(a) >= length(b) && !(length(a)==1 && a[1]==0)
        k = length(a)-length(b); f = a[end]/b[end]; q[k+1] = f
        for i in eachindex(b); a[i+k] -= f*b[i]; end
        a = a[1:end-1]; isempty(a) && (a = QQ[0])
    end
    (ptrim(q), ptrim(a))
end
function sturm_chain(p::P)
    ch = [p, pder(p)]
    while !(length(ch[end])==1 && ch[end][1]==0)
        (_,r) = pdivrem(ch[end-1], ch[end]); push!(ch, pscale(r,-1))
    end
    pop!(ch); ch
end
function sign_changes(ch, x::QQ)
    s = [sign(peval(p,x)) for p in ch]; s = filter(!=(0), s)
    count(i->s[i]!=s[i+1], 1:length(s)-1)
end
"number of distinct real roots in (a,b], a<b rational (p(a) != 0 required)"
function count_roots(p::P, a::QQ, b::QQ)
    ch = sturm_chain(p)
    sign_changes(ch,a) - sign_changes(ch,b)
end

# confluent divided differences of h(x) = (x-1/4)_+^4  (C^3)
function hderiv(x::QQ, k::Int)   # h^(k)(x)/k!
    x <= 1//4 && return QQ(0)
    binomial(4,k) * (x-1//4)^(4-k)
end
function h_divided(v)
    xs = sort(collect(QQ.(v))); n = length(xs)
    T = Matrix{QQ}(undef, n, n)
    for i in 1:n; T[i,1] = hderiv(xs[i],0); end
    for k in 2:n, i in 1:n-k+1
        T[i,k] = xs[i]==xs[i+k-1] ? hderiv(xs[i],k-1) : (T[i+1,k-1]-T[i,k-1])/(xs[i+k-1]-xs[i])
    end
    T[1,n]
end

# (F') box certificate
g_at(t::QQ) = t <= 1//2 ? QQ(0) : t*t*(2t-1)^4
function divided_difference(ts)
    a = [x*x for x in ts]; v = [g_at(x) for x in ts]
    for k in 1:3
        v = [(v[i+1]-v[i])/(a[i+k]-a[i]) for i in 1:4-k]
    end
    v[1]
end
const EPS = QQ(1, 10^6); const DEN = big(10)^6
function corner_bound(lo, hi)
    a3 = 1 - lo[1]-lo[2]-lo[3]
    a3 < 0 && (a3 = QQ(0))
    highs = sort([hi[1],hi[2],hi[3],a3]; rev=true)
    highs = [x + (i-1)*EPS for (i,x) in enumerate(highs)]
    ts = [QQ(isqrt(floor(BigInt, x*DEN*DEN))+1, DEN) for x in highs]
    sort!(ts; rev=true)
    for i in 2:4
        ts[i] >= ts[i-1] && (ts[i] = ts[i-1] - QQ(1,10^7))
    end
    divided_difference(ts)
end
function feasible(lo,hi)
    lo[1]+lo[2]+lo[3] > 1 && return false
    hi[1]+hi[2]+hi[3]+hi[3] < 1 && return false
    hi[2] < 1//4 && return false
    (hi[1] < lo[2] || hi[2] < lo[3]) && return false
    true
end
function box_certificate()
    N0 = 6; STEP = QQ(3, 4N0)
    stack = Tuple{Vector{QQ},Vector{QQ},Int}[]
    for i in 0:N0-1, j in 0:N0-1, k in 0:N0-1
        push!(stack, (QQ[i*STEP,j*STEP,k*STEP], QQ[(i+1)*STEP,(j+1)*STEP,(k+1)*STEP], 0))
    end
    n = 0; largest = QQ(0); unresolved = false
    while !isempty(stack)
        lo,hi,depth = pop!(stack)
        feasible(lo,hi) || continue
        v = corner_bound(lo,hi); n += 1
        if v <= 1
            v > largest && (largest = v); continue
        end
        if depth >= 25; unresolved = true; break; end
        w = hi .- lo; m = argmax(w)       # first index of the maximum, as Python's max(range(3), key=...)
        mid = (lo[m]+hi[m])/2
        lo2,hi2 = copy(lo),copy(hi); hi2[m] = mid
        lo3,hi3 = copy(lo),copy(hi); lo3[m] = mid
        push!(stack,(lo2,hi2,depth+1)); push!(stack,(lo3,hi3,depth+1))
    end
    (n, largest, unresolved)
end

function run(io::IO=stdout)
    R = Common.Report("cap_certificate/cap_inequality_certificate")
    ROOTS = integer_roots()
    Common.check!(R, io, "(A) root system has 24 elements of squared length 2", length(ROOTS)==24 && all(dotq(r,r)==2 for r in ROOTS))
    alpha0 = ROOTS[1]; others = ROOTS[2:end]
    Qv = vertices_of(others)
    Common.check!(R, io, "(A) Q has 25 vertices", length(Qv)==25; got=length(Qv), expected=25)
    volQ = 4*exact_volume(others, Qv)      # x = sqrt(2) z
    Common.check!(R, io, "(A) vol(Q) = 25/3 (exact)", volQ==25//3; got=string(volQ), expected="25//3")
    Cv = vertices_of(ROOTS)
    volC = 4*exact_volume(ROOTS, Cv)
    Common.check!(R, io, "(A) vol(24-cell) = 8 and vol(Q) - vol(24-cell) = 1/3 (exact)", volC==8 && volQ-volC==1//3; got=string(volC), expected="8//1")
    orth = [r for r in ROOTS if dotq(r,alpha0)==0]; frame = NTuple{4,Int}[]
    for r in orth
        all(dotq(r,f)==0 for f in frame) && push!(frame, r)
        length(frame)==3 && break
    end
    E = vcat([alpha0], frame)
    Common.check!(R, io, "(B) frame is pairwise orthogonal", length(E)==4 && all(dotq(E[i],E[j])==0 for i in 1:4 for j in i+1:4))
    # l^1 norm = (sum_k |<E_k,v>|)/sqrt2 <= sqrt2  <=>  sum_k |<E_k,v>| <= 2
    l1 = [sum(abs(dotq(E[k],v)) for k in 1:4) for v in Qv]
    Common.check!(R, io, "(B) every vertex has l^1 norm at most sqrt(2), attained", all(x<=2 for x in l1) && maximum(l1)==2; got="max sqrt2*l1 = $(maximum(l1))", expected="2")
    # (D)
    g = LP(); for (k,c) in enumerate(Common.QQ[1,-8,24,-32,16]); g[k+1] = c; end   # t^2 (2t-1)^4
    g = LP(2=>1, 3=>-8, 4=>24, 5=>-32, 6=>16)
    d1 = dda(g); d2 = dda(d1); d3 = dda(d2); d4 = dda(d3); d5 = dda(d4)
    Common.check!(R, io, "(D) g''''(t^2) = 3(20t^2 - 3)/(2 t^5)", clean(d4)==LP(-3=>30, -5=>-9//2))
    Common.check!(R, io, "(D) nonnegative at t = 1/2 (and 20t^2-3 increasing)", evalp(d4,QQ(1,2)) >= 0; got=string(evalp(d4,QQ(1,2))), expected="96")
    Common.check!(R, io, "(D) g'''(1/4) = 0", evalp(d3,QQ(1,2))==0)
    # (E)
    c = QQ[0,1]; one = QQ[1]
    e1 = padd(pmul(ppow(c,2), padd(pscale(ppow(c,2),2), QQ[-1])), pscale(ppow(QQ[-1,2],4),-1))
    e2 = padd(pmul(padd(ppow(c,2),QQ[-1//4]), padd(pscale(ppow(c,2),2),QQ[-3//4])), pscale(ppow(QQ[-1,2],4),-1))
    e3 = padd(ppow(QQ[1//2,1],3), pscale(pmul(QQ[-1//2,1], ppow(c,2)), -8))
    # rational enclosures of sqrt3/2 and 1/sqrt2 widened by 1e-9 as in the Python tolerance
    s3lo, s3hi = QQ(866025403,10^9), QQ(866025404,10^9); @assert 4s3lo^2 < 3 < 4s3hi^2
    s2lo, s2hi = QQ(707106781,10^9), QQ(707106782,10^9); @assert 2s2lo^2 < 1 < 2s2hi^2
    tol = QQ(1,10^9)
    cases = [("c >= sqrt(3)/2", e1, s3lo-tol, QQ(1)+tol), ("1/sqrt(2) <= c <= sqrt(3)/2", e2, s2lo-tol, s3hi+tol), ("1/2 <= c <= 1/sqrt(2)", e3, QQ(1,2)-tol, s2hi+tol)]
    allE = true
    for (name, ex, lo, hi) in cases
        q = ex; k1 = 0
        while peval(q, QQ(1)) == 0; (q, r) = pdivrem(q, QQ[-1,1]); @assert r==QQ[0]; k1 += 1; end
        nroots = count_roots(q, lo, hi)
        mid = (lo+hi)/2
        ok = nroots==0 && peval(ex, mid) > 0 && peval(q,lo) != 0
        allE &= ok
        Common.check!(R, io, "(E) no sign change on $name", ok; got="roots other than c=1 in widened interval: $nroots; (c-1)-multiplicity $k1")
    end
    Common.check!(R, io, "(E) case complete", allE)
    # (F)
    Common.check!(R, io, "(F) g^(5)(a) = (45/4) a^(-7/2) (1 - 4a)", clean(d5)==LP(-7=>45//4, -5=>-45))
    Common.check!(R, io, "(F) g''''(1/4) = 96", evalp(d4,QQ(1,2))==96)
    ineq = [QQ[1,-1,0,0], QQ[0,1,-1,0], QQ[0,0,1,-1], QQ[0,0,0,1], QQ[0,1,0,0]]; rhs = QQ[0,0,0,0,1//4]
    verts = Set{Vector{QQ}}()
    for a in 1:5, b in a+1:5, d in b+1:5
        M = QQ[(i<=3 ? ineq[(a,b,d)[i]][k] : QQ(1)) for i in 1:4, k in 1:4]
        z = qsolve(M, QQ[rhs[a],rhs[b],rhs[d],1])
        z === nothing && continue
        all(dotq(ineq[i],z) >= rhs[i] for i in 1:5) && push!(verts, z)
    end
    expectV = Set([QQ[3//4,1//4,0,0], QQ[1//2,1//2,0,0], QQ[1//2,1//4,1//4,0], QQ[1//3,1//3,1//3,0], QQ[1//4,1//4,1//4,1//4]])
    Common.check!(R, io, "(F) R has the five vertices of the paper", verts==expectV; got=length(verts), expected=5)
    vals = sort([h_divided(v) for v in verts]; rev=true)
    Common.check!(R, io, "(F) vertex values of h are 2/9, 3/16, 1/8, 27/256, 0", vals==QQ[2//9,3//16,1//8,27//256,0]; got=join(string.(vals),", "), expected="2//9, 3//16, 1//8, 27//256, 0")
    Common.check!(R, io, "(F) so F <= 4 * 2/9 = 8/9 < 1 on R", 4*maximum(vals) < 1; got=string(4*maximum(vals)), expected="8//9")
    n, largest, unres = box_certificate()
    lg = string(round(Float64(largest); digits=8))
    Common.check!(R, io, "(F') every corner bound is at most 1: 303 boxes, largest 0.99755050", !unres && n==303 && lg=="0.9975505"; got="boxes=$n largest=$lg", expected="boxes=303 largest=0.99755050")
    R
end
end
