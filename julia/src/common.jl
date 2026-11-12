# Shared reporting helpers and exact rational interval arithmetic.
module Common
export Report, check!, QQ, RI
const QQ = Rational{BigInt}

struct Check
    name::String
    pass::Bool
    got::String
    expected::String
end
mutable struct Report
    title::String
    checks::Vector{Check}
    seconds::Float64
    Report(t) = new(t, Check[], 0.0)
end
function check!(R::Report, io::IO, name, pass::Bool; got="", expected="")
    c = Check(String(name), pass, string(got), string(expected))
    push!(R.checks, c)
    print(io, pass ? "PASS " : "FAIL ", name)
    c.got != "" && print(io, "  [got ", c.got, c.expected != "" ? "; expected " * c.expected : "", "]")
    println(io); flush(io)
    pass
end
allpass(R::Report) = all(c.pass for c in R.checks)

datadir() = normpath(joinpath(@__DIR__, "..", "data"))

# ---------------------------------------------------------------------------
# Exact interval arithmetic with rational endpoints (no rounding at all).
struct RI
    lo::QQ
    hi::QQ
end
RI(x::Union{Integer,Rational}) = RI(QQ(x), QQ(x))
enclose(x::RI) = x
enclose(x) = RI(x)
Base.:+(a::RI, b::RI) = RI(a.lo+b.lo, a.hi+b.hi)
Base.:+(a::RI, b) = a + enclose(b)
Base.:+(a, b::RI) = enclose(a) + b
Base.:-(a::RI) = RI(-a.hi, -a.lo)
Base.:-(a::RI, b::RI) = a + (-b)
Base.:-(a::RI, b) = a + (-enclose(b))
Base.:-(a, b::RI) = enclose(a) + (-b)
function Base.:*(a::RI, b::RI)
    p = (a.lo*b.lo, a.lo*b.hi, a.hi*b.lo, a.hi*b.hi)
    RI(minimum(p), maximum(p))
end
Base.:*(a::RI, b) = a * enclose(b)
Base.:*(a, b::RI) = enclose(a) * b
function Base.:/(a::RI, b::RI)
    @assert b.lo > 0
    a * RI(1//b.hi, 1//b.lo)
end
Base.:/(a::RI, b) = a / enclose(b)
Base.:/(a, b::RI) = enclose(a) / b

"floor-sqrt enclosure on the grid 10^-24 (as in anc/cap_comparison.py)"
function sqrt_enclosure(v::QQ)
    @assert v >= 0
    grid = big(10)^24
    n = isqrt(fld(numerator(v)*grid^2, denominator(v)))
    RI(n//grid, (n+1)//grid)
end
function atan_endpoint(v::QQ; terms=48)
    @assert 0 <= v < 1
    lower = sum((-1)^k * v^(2k+1)//(2k+1) for k in 0:terms-1)
    RI(lower, lower + v^(2terms+1)//(2terms+1))
end
atan_enclosure(v::RI) = RI(atan_endpoint(v.lo).lo, atan_endpoint(v.hi).hi)
atan_enclosure(v) = atan_enclosure(enclose(v))
const PI_RI = Ref{RI}()
function pi_ri()
    isassigned(PI_RI) || (PI_RI[] = 16*atan_enclosure(QQ(1,5)) - 4*atan_enclosure(QQ(1,239)))
    PI_RI[]
end
function decimal_interval(x::RI, digits=8)
    s = big(10)^digits
    lo = fld(numerator(x.lo*s), denominator(x.lo*s)); hi = cld(numerator(x.hi*s), denominator(x.hi*s))
    f(n) = string(fld(n,s), ".", lpad(string(mod(n,s)), digits, '0'))
    "[" * f(lo) * ", " * f(hi) * "]"
end
end
