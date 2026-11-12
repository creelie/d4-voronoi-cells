# Exact rational linear algebra, polynomial helpers, parsing; and Arb helpers.
module Exact
using ..Common: QQ
using Arblib
export fr, ldl_positive, ldl_psd, solve_exact, A_, arb_lo, arb_hi, PREC

"parse a rational: integer, float (exact binary value), 'a/b', or decimal string"
fr(x::Integer) = QQ(x)
fr(x::Rational) = QQ(x)
fr(x::AbstractFloat) = QQ(x)            # exact
function fr(s::AbstractString)
    s = strip(s)
    if occursin('/', s)
        a, b = split(s, '/'); return QQ(parse(BigInt, a), parse(BigInt, b))
    end
    e = 0
    m = match(r"^([^eE]+)[eE]([+-]?\d+)$", s)
    m !== nothing && ((s, e) = (m.captures[1], parse(Int, m.captures[2])))
    neg = startswith(s, "-"); s = lstrip(s, ['+','-'])
    if occursin('.', s)
        ip, fp = split(s, '.'); v = QQ(parse(BigInt, ip*fp), big(10)^length(fp))
    else
        v = QQ(parse(BigInt, s))
    end
    v *= QQ(10)^e
    neg ? -v : v
end

"exact LDL^T; true iff symmetric and every pivot is positive"
function ldl_positive(M::Matrix{QQ})
    n = size(M,1); A = copy(M)
    for i in 1:n, j in 1:n; A[i,j] == A[j,i] || return false; end
    for k in 1:n
        piv = A[k,k]
        piv <= 0 && return false
        for i in k+1:n
            m = A[i,k]/piv
            m == 0 && continue
            for j in k+1:n; A[i,j] -= m*A[k,j]; end
        end
    end
    true
end
"exact LDL^T; (true, pivots) if positive semidefinite (zero pivots need zero columns)"
function ldl_psd(M::Matrix{QQ})
    n = size(M,1); A = copy(M); piv = QQ[]
    for k in 1:n
        p = A[k,k]
        p < 0 && return (false, piv)
        if p == 0
            any(A[i,k] != 0 for i in k+1:n) && return (false, piv)
            push!(piv, p); continue
        end
        push!(piv, p)
        for i in k+1:n
            f = A[i,k]/p
            f == 0 && continue
            for j in k+1:n; A[i,j] -= f*A[k,j]; end
        end
    end
    (true, piv)
end
function solve_exact(M::Matrix{QQ}, b::Vector{QQ})
    n = size(M,1); A = hcat(copy(M), b)
    for k in 1:n
        p = findfirst(i->A[i,k]!=0, k:n) + k - 1
        if p != k; A[k,:], A[p,:] = A[p,:], A[k,:]; end
        for i in 1:n
            if i != k && A[i,k] != 0
                f = A[i,k]/A[k,k]; A[i,:] .-= f .* A[k,:]
            end
        end
    end
    [A[i,n+1]/A[i,i] for i in 1:n]
end
issym(M) = all(M[i,j]==M[j,i] for i in axes(M,1), j in axes(M,2))

const PREC = Ref(128)
A_(x::QQ; prec=PREC[]) = Arb(x; prec=prec)
A_(x; prec=PREC[]) = Arb(QQ(x); prec=prec)
"a double that is <= every point of the ball"
arb_lo(x::Arb) = Float64(Arblib.lbound(x), RoundDown)
"a double that is >= every point of the ball"
arb_hi(x::Arb) = Float64(Arblib.ubound(x), RoundUp)
"python '%s' % arb.str(8)-like display"
arbstr(x::Arb, d=8) = Arblib.string_nice(x, d)
end
