# Ports of multi_cap/llm24_certificate_check.py (prop:verified, steps 1, 2, 4, 6, 7 of the
# verification of the certificate of de Laat, Leijenhorst and de Muinck Keizer) and of
# multi_cap/explicit_eps0.py (thm:m24-robust, thm:near-contact: epsilon_0 = 4e-26).
# Exact Rational{BigInt}; Arb ball arithmetic at 256 bits where the Python uses python-flint's arb.
module LLM24
using ..Common
using ..Common: QQ
using ..Exact
using ..CapCertificate: count_roots, pdivrem, peval
import ..CapCertificate
const CPmul = CapCertificate.pmul
using Arblib, SHA, Downloads, JSON

const DOI = "10.4121/74ce1c25-6fca-4680-8a36-e9c18e7e9594"
const ZIPNAME = "LasserreSphericalCodes.zip"
const SHA256 = "2f806fac9f5f8644c53f9ad1590d920d346b61549d5d8860729ae6125a7568cc"
const ZIPBYTES = 152065368

sha256_of(path) = open(io->bytes2hex(sha256(io)), path)

"""
fetch(dest): download LasserreSphericalCodes.zip from 4TU.ResearchData (DOI $DOI), check its
SHA-256 and unpack it under `dest`.  Returns the folder proofs/4_24.  An already present and
correct archive is not downloaded again.
"""
function fetch(dest::AbstractString; io::IO=stdout)
    mkpath(dest)
    zip = joinpath(dest, ZIPNAME)
    if !(isfile(zip) && sha256_of(zip) == SHA256)
        # the record's identifier at 4TU.ResearchData is the suffix of the DOI
        api = "https://data.4tu.nl/v2/articles/" * last(split(DOI, "/")) * "/files"
        files = JSON.parse(String(take!(Downloads.download(api, IOBuffer()))))
        isempty(files) && error("the 4TU API lists no file for $DOI")
        f = only(filter(x->x["name"] == ZIPNAME, files))
        println(io, "downloading ", f["download_url"], " (", f["size"], " bytes)")
        Downloads.download(f["download_url"], zip)
    end
    got = sha256_of(zip)
    got == SHA256 || error("SHA-256 mismatch for $zip: got $got, expected $SHA256")
    println(io, "SHA-256 ", got, " matches (", ZIPNAME, ")")
    out = joinpath(dest, "LasserreSphericalCodes", "proofs", "4_24")
    isdir(out) || run(`unzip -q -o $zip -d $dest`)
    out
end

# ---------------------------------------------------------------- reading
function parse_q(s::AbstractString)
    i = findfirst('/', s)
    i === nothing && return QQ(parse(BigInt, s))
    j = findlast('/', s)
    QQ(parse(BigInt, SubString(s, 1, i-1)), parse(BigInt, SubString(s, j+1)))
end
function load_dense(path)
    open(path) do f
        r, c = parse.(Int, split(readline(f)))
        @assert r > 0 && c > 0
        M = Matrix{QQ}(undef, r, c); i = 0
        for line in eachline(f)
            isempty(strip(line)) && continue
            i += 1; t = split(line); @assert length(t) == c
            for j in 1:c; M[i,j] = parse_q(t[j]); end
        end
        @assert i == r
        M
    end
end
const MP = Dict{Vector{Int},QQ}       # multivariate polynomial
function load_low_rank(path)
    open(path) do f
        nvecs, veclength = parse.(Int, split(readline(f)))
        @assert nvecs > 0 && veclength > 0
        function read_poly()
            d = MP()
            while !eof(f)
                line = readline(f)
                isempty(strip(line)) && break
                t = split(line); e = parse.(Int, t[2:end])
                d[e] = get(d, e, QQ(0)) + parse_q(t[1])
            end
            filter!(kv->kv[2] != 0, d)
        end
        pref = [read_poly() for _ in 1:nvecs]
        vecs = [[read_poly() for _ in 1:veclength] for _ in 1:nvecs]
        (pref, vecs)
    end
end
blockname(fname, prefix) = Tuple(startswith(p, "sos") ? p : parse(Int, p) for p in split(fname[length(prefix)+1:end-4], "_"))

# ---------------------------------------------------------------- polynomials
mp_const(nv, c=QQ(1)) = MP(zeros(Int, nv) => QQ(c))
mp_gen(nv, i) = (e = zeros(Int, nv); e[i] = 1; MP(e => QQ(1)))
function mp_add(p::MP, q::MP, s=1)
    r = copy(p)
    for (e,c) in q; r[e] = get(r, e, QQ(0)) + s*c; end
    filter!(kv->kv[2] != 0, r)
end
function mp_mul(p::MP, q::MP)
    r = MP()
    for (e1,c1) in p, (e2,c2) in q
        e = e1 .+ e2; r[e] = get(r, e, QQ(0)) + c1*c2
    end
    filter!(kv->kv[2] != 0, r)
end
mp_scale(p::MP, s) = MP(e => c*s for (e,c) in p)
mp_deg(p::MP) = maximum(sum(e) for e in keys(p))
function mp_det(M)
    n = length(M)
    n == 1 && return M[1][1]
    res = MP()
    for j in 1:n
        minor = [vcat(row[1:j-1], row[j+1:end]) for row in M[2:end]]
        res = mp_add(res, mp_mul(M[1][j], mp_det(minor)), isodd(j) ? 1 : -1)
    end
    res
end
function gram_vars(k, nv)
    X = [Vector{MP}(undef, k) for _ in 1:k]; idx = 0
    for i in 1:k
        X[i][i] = mp_const(nv)
        for j in i+1:k
            idx += 1; X[i][j] = X[j][i] = mp_gen(nv, idx)
        end
    end
    X
end
function subsets(n, r)
    out = Vector{Vector{Int}}()
    rec(start, cur) = length(cur) == r ? push!(out, copy(cur)) : (for i in start:n; push!(cur, i); rec(i+1, cur); pop!(cur); end)
    rec(1, Int[]); out
end
function elementary(polys, nv)
    [reduce(mp_add, [reduce(mp_mul, polys[S]) for S in subsets(length(polys), r)]) for r in 1:length(polys)]
end
const NV = [1, 1, 3, 6]
"the weights of the domain of k points, with names: e<r>_pair, det, e<r>_minor"
function domain_weights(k, costheta::QQ)
    nv = NV[k]; ws = Tuple{String,MP}[]
    pairw(i) = mp_mul(mp_add(mp_gen(nv,i), mp_const(nv)), mp_add(mp_const(nv, costheta), mp_gen(nv,i), -1))
    if k == 2
        push!(ws, ("e1_pair", pairw(1)))
    elseif k >= 3
        X = gram_vars(k, nv)
        push!(ws, ("det", mp_det(X)))
        for (r,w) in enumerate(elementary([pairw(i) for i in 1:nv], nv)); push!(ws, ("e$(r)_pair", w)); end
        if k == 4
            mins = [mp_det([[X[a][b] for b in S] for a in S]) for S in subsets(4, 3)]
            for (r,w) in enumerate(elementary(mins, nv)); push!(ws, ("e$(r)_minor", w)); end
        end
    end
    ws
end
"(name, c) with e = c * weight, or (\"const\", c), (\"zero\", 0); nothing if e is not a multiple of a weight"
function classify(e::MP, weights)
    isempty(e) && return ("zero", QQ(0))
    mp_deg(e) == 0 && return ("const", first(values(e)))
    for (name, w) in weights
        (mp_deg(w) == mp_deg(e) && keys(w) == keys(e)) || continue
        ratios = Set(e[m]/w[m] for m in keys(e))
        length(ratios) == 1 && return (name, first(ratios))
    end
    nothing
end

# ---------------------------------------------------------------- admissible tuples (Section 2.3 of the cited paper)
function lambdas(t, d, lmax=d)
    d < 0 && return Vector{Int}[]
    (d == 0 || t == 0) && return [zeros(Int, t)]
    out = Vector{Int}[]
    for l1 in 0:min(d, lmax), aux in lambdas(t-1, d-l1, l1)
        push!(out, vcat(l1, aux))
    end
    out
end
countels(lam, j, k) = j == 1 ? lam[1]-k : lam[2]+k
function isadmissible(lam, i, j, k)
    i == 0 && return sum(lam) == 0 && j == 0 && k == 0
    i == 1 && return lam[2] == 0 && j == 0 && k == 0
    countels(lam, 2, k) % 2 == 0
end
function admissible_count(lam, d1, d2)
    ws = [k for k in 0:lam[1]-lam[2] if countels(lam, 2, k) % 2 == 0]
    dt = (d2 - sum(lam)) ÷ 2
    count(isadmissible(lam, i, j, k) for i in 0:2 for j in 0:dt for k in ws)
end

# ---------------------------------------------------------------- positivity
function ldlt_positive(M::Matrix{QQ})
    n = size(M,1); A = copy(M)
    for k in 1:n
        p = A[k,k]
        p <= 0 && return false
        for i in k+1:n
            f = A[k,i]/p
            f == 0 && continue
            for j in i:n; A[i,j] -= f*A[k,j]; end
        end
    end
    true
end
"rigorous Cholesky in ball arithmetic: every pivot a ball inside (0, inf). Returns (ok, precision, lower bound of the least pivot)"
function cholesky_positive(M::Matrix{QQ}; prec=256, maxprec=2^14)
    n = size(M,1)
    while prec <= maxprec
        A = ArbMatrix(n, n; prec=prec)
        for i in 1:n, j in 1:n; A[i,j] = Arb(M[i,j]; prec=prec); end
        L = ArbMatrix(n, n; prec=prec)
        if Arblib.cho!(L, A) != 0
            least = nothing; ok = true
            for k in 1:n
                d = L[k,k]^2
                ok &= Arblib.ispositive(L[k,k]) && Arblib.ispositive(d)
                lb = Arb(Arblib.lbound(d); prec=prec)
                least = least === nothing || lb < least ? lb : least
            end
            ok && return (true, prec, least)
        end
        prec *= 2
    end
    (false, prec, nothing)
end

"p_2 as a coefficient vector (index 1 = constant term), built exactly from the data"
function two_point_polynomial(folder, haveblock)
    coeffs = Dict{Int,QQ}()
    for fn in sort(readdir(folder))
        startswith(fn, "poly_2_") || continue
        name = fn[8:end-4]
        haveblock(name) || continue
        X = load_dense(joinpath(folder, "dense_$name.txt")); B = load_dense(joinpath(folder, "transform_$name.txt"))
        M = B*X*permutedims(B); m = size(M,1)
        pref, vecs = load_low_rank(joinpath(folder, fn))
        for (e, v) in zip(pref, vecs)
            @assert length(v) == m
            q = MP()
            for i in 1:m, j in 1:m
                M[i,j] == 0 && continue
                q = mp_add(q, mp_scale(mp_mul(v[i], v[j]), M[i,j]))
            end
            for (mon, c) in mp_mul(e, q); coeffs[mon[1]] = get(coeffs, mon[1], QQ(0)) + c; end
        end
    end
    deg = maximum(keys(coeffs))
    QQ[get(coeffs, i, QQ(0)) for i in 0:deg]
end

function run_check(folder::AbstractString, io::IO=stdout; expected_p2=nothing)
    R = Common.Report("multi_cap/llm24_certificate_check (prop:verified, steps 1, 2, 4, 6, 7)")
    lines = readlines(joinpath(folder, "metadata.txt"))
    t = split(lines[2]); n, d1, d2, delta = parse(Int,t[1]), parse(Int,t[4]), parse(Int,t[5]), parse(Int,t[6])
    Ncode = parse_q(t[2]); costheta = parse_q(t[3]); ips = parse_q.(split(lines[3]))
    Common.check!(R, io, "metadata: n = 4, N = 24, cos theta = 1/2, d1 = 14, d2 = 16, inner products -1, -1/2, 0, 1/2",
        lines[1] == "codebound" && n == 4 && Ncode == 24 && costheta == 1//2 && d1 == 14 && d2 == 16 && delta == 16 && ips == QQ[-1, -1//2, 0, 1//2]; got=join(lines, " | "), expected="codebound | 4 24 1//2 14 16 16 | -1 -1//2 0 1//2")
    files = sort(readdir(folder))
    dense = [fn for fn in files if startswith(fn, "dense")]
    names = Dict(fn => blockname(fn, "dense_") for fn in dense)
    # step 1 and step 2 together, one block at a time (the blocks are large)
    okfmt = true; allpos = true; precs = Dict{Int,Int}(); least_all = nothing; exact_blocks = 0
    sizes = Dict{Any,Int}(); trows = Dict{Any,Int}(); obj = nothing; t1 = time()
    for fn in dense
        b = names[fn]
        tf = "transform_" * fn[7:end]
        okfmt &= tf in files
        X = load_dense(joinpath(folder, fn)); B = load_dense(joinpath(folder, tf))
        okfmt &= size(X,1) == size(X,2) == size(B,2)
        sizes[b] = size(X,1); trows[b] = size(B,1)
        ok, prec, least = cholesky_positive(X)
        if ok
            precs[prec] = get(precs, prec, 0) + 1
            least_all = least_all === nothing || least < least_all ? least : least_all
        else
            allpos = false; println(io, "   block $b: Cholesky failed up to precision $prec bits")
        end
        if size(X,1) <= 16
            exact_blocks += 1
            ldlt_positive(X) || (allpos = false; println(io, "   block $b: exact LDL^T finds a nonpositive pivot"))
        end
        b == (0,0) && (obj = (B*X*permutedims(B))[1,1])
    end
    nirr = count(b->b[1] isa Int, keys(sizes)); cnt(s) = count(b->b[1] == s, keys(sizes))
    Common.check!(R, io, "step 1: 127 dense blocks X_i read, each square and matching its transformation B_i", length(sizes) == 127 && okfmt;
        got="irreducible representations: $nirr blocks; sos2: $(cnt("sos2")); sos3: $(cnt("sos3")); sos4: $(cnt("sos4")); total dimension $(sum(values(sizes))), largest $(maximum(values(sizes)))",
        expected="irreducible representations: 60 blocks; sos2: 2; sos3: 15; sos4: 50; total dimension 3726, largest 350")
    irreps = [lam for lam in lambdas(2, d1) if lam[1] != lam[2] || lam[1] % 2 == 0]
    expd = Dict(Tuple(lam) => admissible_count(lam, d1, delta) for lam in irreps)
    irb = [b for b in keys(sizes) if b[1] isa Int]
    Common.check!(R, io, "the block of every lambda has as many rows as there are admissible tuples (i, j, k)", all(haskey(expd, b) && trows[b] == expd[b] for b in irb);
        got="$(length(irb)) of the $(length(expd)) signatures", expected="60 of the 60 signatures")
    polyfiles = [fn for fn in files if startswith(fn, "poly")]
    pk(fn) = parse(Int, split(fn[1:end-4], "_")[2])
    pb(fn) = Tuple(startswith(p, "sos") ? p : parse(Int, p) for p in split(fn[1:end-4], "_")[3:end])
    counts = [count(fn->pk(fn) == k, polyfiles) for k in 1:4]
    Common.check!(R, io, "step 1: 76 sum-of-squares files (1 + 2 + 15 + 58 for the 1-, 2-, 3-, 4-point constraints)", counts == [1,2,15,58]; got=counts, expected="[1, 2, 15, 58]")
    Common.check!(R, io, "step 2: every one of the 127 blocks X_i is positive definite (Cholesky in ball arithmetic; exact LDL^T for size <= 16)", allpos;
        got="precision (bits => blocks) $(sort(collect(precs))); least pivot >= $(Float64(least_all)); $exact_blocks blocks also by exact LDL^T; $(round(Int, time()-t1)) s",
        expected="[(256, 127)]; least pivot >= 1.379798072887995561e-15; 81 blocks by exact LDL^T")
    Common.check!(R, io, "step 6: the objective (B X B^T)_{11} for lambda = (0,0) equals 24 exactly", obj == Ncode; got=string(obj), expected="24")
    # step 4
    ok = true; pc = Dict{Int,Int}()
    for k in 1:4
        W = domain_weights(k, costheta)
        for fn in polyfiles
            (pk(fn) == k && haskey(sizes, pb(fn))) || continue
            pref, _ = load_low_rank(joinpath(folder, fn))
            for e in pref
                pc[k] = get(pc, k, 0) + 1
                c = classify(e, W)
                (c !== nothing && c[2] >= 0) || (ok = false; println(io, "   constraint $k block $(pb(fn)): bad prefactor"))
            end
        end
    end
    Common.check!(R, io, "step 4: every sum-of-squares prefactor is a nonnegative constant or a nonnegative multiple of a domain weight", ok; got=sort(collect(pc)), expected="{2: 2, 3: 20, 4: 103}")
    # step 7
    P2 = two_point_polynomial(folder, name->isfile(joinpath(folder, "dense_$name.txt")))
    deg = length(P2)-1
    Common.check!(R, io, "step 7: the two-point polynomial p_2(u) vanishes at u = -1, -1/2, 0, 1/2", all(peval(P2, ip) == 0 for ip in ips); got="degree $deg", expected="degree 16")
    if expected_p2 !== nothing
        ref = QQ[]
        for line in eachline(expected_p2)
            (startswith(line, "#") || isempty(strip(line))) && continue
            push!(ref, fr(split(line, limit=2)[2]))
        end
        Common.check!(R, io, "p_2 equals the coefficients written by the Python run (llm24_out/llm24_p2.txt), exactly", ref == P2; got="$(length(P2)) coefficients", expected="$(length(ref)) coefficients")
    end
    Q = P2; mult = Int[]
    for ip in ips
        m = 0
        while peval(Q, ip) == 0; (Q, r) = pdivrem(Q, QQ[-ip, 1]); m += 1; end
        push!(mult, m)
    end
    nz = count_roots(Q, QQ(-1), costheta)
    Common.check!(R, io, "step 7: after removing those four zeros, p_2 has no further zero on [-1, 1/2] (Sturm sequence over Q)", nz == 0 && peval(Q, QQ(-1)) != 0 && peval(Q, costheta) != 0;
        got="multiplicities $mult; zeros of the quotient: $nz", expected="multiplicities [1, 2, 2, 1]; zeros of the quotient: 0")
    Common.check!(R, io, "p_2 is nonnegative at 401 rational points of [-1, 1/2]", all(peval(P2, QQ(-1) + QQ(3,2)*QQ(i,400)) >= 0 for i in 0:400))
    R
end

# ================================================================ explicit_eps0.py
const P256 = 256
aq(q) = Arb(QQ(q); prec=P256)
upb(x::Arb) = Arb(Arblib.ubound(x); prec=P256)
const _MC = Dict{Int,Dict{Int,QQ}}()
function mono_cheb(n)
    get!(_MC, n) do
        out = Dict{Int,QQ}()
        for j in 0:n
            k = abs(n-2j); out[k] = get(out, k, QQ(0)) + QQ(binomial(big(n), big(j)), big(2)^n)
        end
        out
    end
end
"l1 norm of the tensor Chebyshev coefficients: a bound for sup |poly| on [-1,1]^n"
function cheb_l1(poly::MP)
    out = MP()
    for (e, c) in poly
        terms = MP(Int[] => c)
        for n in e
            new = MP()
            for (key, v) in terms, (k, a) in mono_cheb(n)
                kk = vcat(key, k); new[kk] = get(new, kk, QQ(0)) + v*a
            end
            terms = new
        end
        for (key, v) in terms; out[key] = get(out, key, QQ(0)) + v; end
    end
    s = Arb(0; prec=P256)
    for v in values(out); v != 0 && (s += abs(aq(v))); end
    s
end
function lambda_max_bound(X::Matrix{QQ})
    n = size(X,1); tr = Arb(0; prec=P256); rowmax = Arb(0; prec=P256)
    for i in 1:n
        rs = Arb(0; prec=P256)
        for j in 1:n
            a = aq(X[i,j]); rs += abs(a)
            j == i && (tr += a)
        end
        u = upb(rs)
        (u - rowmax) > 0 && (rowmax = u)
    end
    (tr - rowmax) < 0 ? upb(tr) : rowmax
end
function part_A(folder)
    S = Dict{Tuple{Int,Int},Arb}()
    for k in (3,4)
        W = domain_weights(k, QQ(1,2))
        for fn in sort(readdir(folder))
            startswith(fn, "poly_$(k)_") || continue
            name = fn[length("poly_$(k)_")+1:end-4]
            isfile(joinpath(folder, "dense_$name.txt")) || continue
            pref, vecs = load_low_rank(joinpath(folder, fn))
            cls = [classify(pf, W) for pf in pref]
            @assert all(c !== nothing && c[2] >= 0 for c in cls)
            any(endswith(c[1], "_pair") for c in cls) || continue
            lam = lambda_max_bound(load_dense(joinpath(folder, "dense_$name.txt")))
            T = load_dense(joinpath(folder, "transform_$name.txt")); m = size(T,2)
            for ((wname, c), v) in zip(cls, vecs)
                @assert endswith(wname, "_pair")
                s2 = Arb(0; prec=P256)
                for b in 1:m
                    acc = MP()
                    for a in 1:size(T,1)
                        tab = T[a,b]
                        (tab == 0 || isempty(v[a])) && continue
                        for (e, q) in v[a]; acc[e] = get(acc, e, QQ(0)) + tab*q; end
                    end
                    l = cheb_l1(acc); s2 += l*l
                end
                key = (k, parse(Int, wname[2:2]))
                S[key] = get(S, key, Arb(0; prec=P256)) + aq(c)*lam*s2
            end
        end
    end
    Dict(k => upb(v) for (k,v) in S)
end
horner(c::Vector{Arb}, x::Arb) = foldr((a, acc)->a + x*acc, c; init=Arb(0; prec=P256))

const S11 = QQ(33166248, 10^7)
const RIG_L = S11/4 + QQ(12991,10000)
const RIG_M = 6*S11/8 + 6*QQ(10826,10000) + QQ(1,2)
function rigidity(S::QQ)
    x1 = RIG_L*(1+S)*S/(1 - RIG_M/48)
    @assert x1 <= QQ(1,48)
    RIG_L*(1+S)*S/(1 - RIG_M*x1)
end
function local_bracket(S::QQ)
    s2, s3, s3l, s5, s24 = QQ(14143,10000), QQ(17321,10000), QQ(1732,1000), QQ(22361,10000), QQ(4899,1000)
    d = S; E = rigidity(S); e = E
    r1 = s2*(1+d/2)/(1-s2*e)
    eta1 = d/2 + r1*e; eta2 = eta1 + d/2 + r1*e
    R = 1 + 2eta2; Rp = R + e/2
    c3 = (QQ(1,2)+e)/(1-e*e/2)
    k1 = 2/(s3l*(1-e))
    Dmax = Rp*e + d/2 + c3*(d/2 + Rp*e)
    thmax = k1*Dmax
    c0 = 1/s3 - thmax
    a_in = 3*s3/4*(1-c0*c0)
    lam = 1 + s3*thmax
    a_out = 3*s3/4*(lam*lam - QQ(1,3))
    gbar = 1/(1-e*e/2)^2
    A1 = k1*8*(Rp + c3*Rp); B1 = k1*4*(1+c3)
    nd = S
    s_e = A1*E*E + B1*E*nd; s_d = A1*nd*E + B1*nd*nd
    loss_v = QQ(1,2)*a_in*(QQ(1,2)*s_d + QQ(1,2)*e*s_e)
    loss_m = QQ(1,2)*gbar*(a_in + R*a_out)*s_e
    Dperp = s5*e + e*e/2 + d/2 + e*(d/2 + e*e/2 + e)/(1-e*e/2)
    loss_c = QQ(1,4)*(QQ(1,2)*S + QQ(1,2)*E*E + gbar*s24*E)*4*Dperp^3
    Theta = e + d
    hull = 24*2*4*(QQ(3,2)+QQ(142,100))^4*Theta^4
    alpha = QQ(1,2)*gbar*(a_in + R*a_out)*A1 + QQ(1,4)*a_in*e*A1
    @assert alpha > QQ(2,3)
    total = QQ(2,3)*S + QQ(2,3)*E*E - loss_v - loss_m - loss_c - hull
    (total/S, (R=R, theta=thmax, a_in=a_in, a_out=a_out, v=loss_v/S, m=loss_m/S, c=loss_c/S, hull=hull/S))
end
rel(a, b; tol=2e-5) = abs(a-b) <= tol*abs(b)

function run_eps0(folder::AbstractString, io::IO=stdout)
    R = Common.Report("multi_cap/explicit_eps0 (thm:near-contact, epsilon_0 = 4e-26)")
    ok0 = S11^2 > 11 && QQ(12991,10000)^2 > QQ(27,16) && QQ(10826,10000)^2 > QQ(75,64)
    S = part_A(folder)
    expS = Dict((3,1)=>1.87162e7, (3,2)=>1.80813e7, (3,3)=>2.65100e6, (4,1)=>1.21897e9, (4,2)=>2.03959e6, (4,3)=>88995.9, (4,4)=>956.266, (4,5)=>50.5478, (4,6)=>2.34368)
    Common.check!(R, io, "A. the sizes S_{k,r} of the sum-of-squares terms whose weight can change sign", ok0 && keys(S) == keys(expS) && all(rel(Float64(S[k]), expS[k]) for k in keys(S));
        got=join(["S$k<=$(Float64(S[k]))" for k in sort(collect(keys(S)))], " "), expected=join(["S$k<=$(expS[k])" for k in sort(collect(keys(expS)))], " "))
    p = two_point_polynomial(folder, name->isfile(joinpath(folder, "dense_$name.txt")))
    x = QQ[0,1]; h = QQ(1,2)
    mul = CPmul
    den = mul(mul(mul(QQ[1,1], mul(QQ[h,1],QQ[h,1])), mul(x,x)), QQ[-h,1])
    q, r = pdivrem(p, den)
    Common.check!(R, io, "B. p_2 = (u+1)(u+1/2)^2 u^2 (u-1/2) q(u) exactly, deg q = 10", r == QQ[0] && length(p) == 17 && length(q) == 11; got="deg q = $(length(q)-1)")
    A = [aq(-c) for c in q]; dA = [k*A[k+1] for k in 1:length(A)-1]
    Nn = 30000; lo = nothing
    for k in 0:Nn-1
        a = QQ(-1) + QQ(3,2)*QQ(k,Nn); b = QQ(-1) + QQ(3,2)*QQ(k+1,Nn)
        m = aq((a+b)/2); rad = (b-a)/2
        ball = Arb((a, b); prec=P256); pm = Arb((-rad, rad); prec=P256)
        val = horner(A, m) + horner(dA, ball)*pm
        lb = Arb(Arblib.lbound(val); prec=P256)
        lo = lo === nothing || lb < lo ? lb : lo
    end
    qmin = lo
    Common.check!(R, io, "B. -q >= q_min > 0 on [-1, 1/2]", qmin > 0 && rel(Float64(qmin), 0.000262008); got=Float64(qmin), expected="0.000262008")
    delta = QQ(3,10000)
    d = aq(delta)
    row = sqrt(Arb(8; prec=P256) + 3*(1+2d)^2); okC = true
    for k in 1:5
        err = (row + 2d*sqrt(Arb(k; prec=P256)))^k - row^k
        okC &= err < 1
    end
    Common.check!(R, io, "C. every k x k minor (k <= 5) of 2g is within less than 1 of that of 2G", okC)
    Dn = 23d
    dW = upb(2Dn/sqrt(6-Dn) + sqrt(Arb(24; prec=P256))*Dn/6)
    Common.check!(R, io, "C. then d(W) < 1/48", dW < Arb(1; prec=P256)/48 && rel(Float64(dW), 0.0112709); got=Float64(dW), expected="0.0112709")
    f(u) = (u+1)*(u+QQ(1,2))^2*u^2*(QQ(1,2)-u)
    fmin = min(f(QQ(-1)+delta), f(QQ(-1,2)-delta), f(QQ(-1,2)+delta), f(-delta), f(delta), f(QQ(1,2)-delta))
    need = qmin*aq(fmin)
    b = Arb(9; prec=P256)/16
    B3 = sum(binomial(3,r)*b^(r-1)*S[(3,r)] for r in 1:3 if haskey(S,(3,r)))
    B4 = sum(binomial(6,r)*b^(r-1)*S[(4,r)] for r in 1:6 if haskey(S,(4,r)))
    Pa = [aq(c) for c in p]; dP = [k*Pa[k+1] for k in 1:16]
    L2 = upb(abs(horner(dP, Arb((QQ(1,2), QQ(1,2)+QQ(1,1000)); prec=P256))))
    Efun(m, kap) = (a = kap*(Arb(3; prec=P256)/2 + kap); binomial(m,2)*kap*L2 + binomial(m,3)*a*B3 + binomial(m,4)*a*B4)
    kap = nothing
    for k in 20:39
        for j in 9:-1:1
            cand = Arb(j; prec=P256)/Arb(10; prec=P256)^k
            if Efun(24, cand) < need; kap = (j,k); break; end
        end
        kap === nothing || break
    end
    kstar = Arb(kap[1]; prec=P256)/Arb(10; prec=P256)^kap[2]
    Common.check!(R, io, "E. constants: need, B_3, B_4, |p_2'|, E(24, kappa*)", rel(Float64(need), 2.94494e-12) && rel(Float64(B3), 8.7499e7; tol=2e-4) && rel(Float64(B4), 7.3316e9; tol=2e-4) && rel(Float64(L2), 0.22572; tol=2e-4) && rel(Float64(upb(Efun(24,kstar))), 2.34247e-12);
        got="need $(Float64(need)), B_3 $(Float64(B3)), B_4 $(Float64(B4)), L2 $(Float64(L2)), E $(Float64(upb(Efun(24,kstar))))", expected="need 2.94494e-12, B_3 8.7499e+7, B_4 7.3316e+9, L2 0.22572, E 2.34247e-12")
    Common.check!(R, io, "E. kappa* = 2e-26 and E(24, kappa*) < need", kap == (2,26) && Efun(24,kstar) < need; got="kappa* = $(kap[1])e-$(kap[2])", expected="kappa* = 2e-26")
    Common.check!(R, io, "E. E(25, kappa*) < 1: no 25 directions with inner products <= 1/2 + kappa*", Efun(25,kstar) < 1)
    Common.check!(R, io, "E. kappa* < delta", kstar < d)
    eps0 = 2kstar
    Common.check!(R, io, "E. epsilon_0 = 2 kappa* = 4e-26: 1/2 - 2/(2 + epsilon_0)^2 <= kappa*", Arb(1; prec=P256)/2 - 2/(2+eps0)^2 <= kstar)
    Smax = 24*QQ(2kap[1], big(10)^kap[2])
    S4 = QQ(4,1000)
    Common.check!(R, io, "D. estimate (a): ||eps|| <= 2.38 (1 + delta) S for S <= 4e-3", RIG_L < QQ(21283,10000) && RIG_M < QQ(949,100) && rigidity(S4) <= QQ(238,100)*S4*(1+S4);
        got="$(round(Float64(RIG_L); digits=4)), $(round(Float64(RIG_M); digits=3)), $(round(Float64(rigidity(S4)/(S4*(1+S4))); digits=4))", expected="2.1283, 9.483, 2.3674")
    br, _ = local_bracket(Smax)
    Common.check!(R, io, "D. at every S in (0, 24 epsilon_0] the bracket is positive", br > 0; got=round(Float64(br); digits=6), expected="0.666667")
    br, parts = local_bracket(S4)
    Common.check!(R, io, "D. the same whenever d(W) <= 1/48 and 0 < S <= 4e-3 (bracket >= 0.1)", rigidity(S4)+S4 <= QQ(14,1000) && br >= QQ(1,10);
        got="bracket $(round(Float64(br); digits=4)), Theta <= $(round(Float64(rigidity(S4)+S4); digits=5)), R $(round(Float64(parts.R); digits=5)), theta $(round(Float64(parts.theta); digits=5)), a_in $(round(Float64(parts.a_in); digits=5)), a_out $(round(Float64(parts.a_out); digits=5))",
        expected="bracket 0.1355, Theta <= 0.01351, R 1.06263, theta 0.02138, a_in 0.89755, a_out 0.96404")
    vals = [local_bracket(QQ(k, 10^5))[1] for k in 1:400]
    Common.check!(R, io, "D. the bracket decreases in S on a grid of step 1e-5 up to 4e-3", all(vals[i] >= vals[i+1] for i in 1:399))
    got = [round(Float64(local_bracket(fr(s))[1]); digits=4) for s in ("1e-3","2e-3","3e-3","4.5e-3","4.6e-3")]
    Common.check!(R, io, "D. bracket at S = 1e-3, 2e-3, 3e-3, 4.5e-3, 4.6e-3", got == [0.5868, 0.4829, 0.3392, 0.0030, -0.0264]; got=got, expected="[0.5868, 0.4829, 0.3392, 0.003, -0.0264]")
    R
end
end
