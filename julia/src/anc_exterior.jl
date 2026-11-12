# Port of anc/exterior_test.py : exact exterior-cut test in R^4 over Q and Q(sqrt d).
module AncExterior
using ..Common
using ..Common: QQ
const V = NTuple{4,QQ}
const ZERO = (QQ(0),QQ(0),QQ(0),QQ(0))
const AXES = ntuple(j->ntuple(i->QQ(i==j),4),4)
dot(u,v) = u[1]*v[1]+u[2]*v[2]+u[3]*v[3]+u[4]*v[4]
shift(u,v,a) = ntuple(i->u[i]+a*v[i],4)
scale(v,a) = ntuple(i->a*v[i],4)
sgn(x) = (x>0) - (x<0)
function radical_sign(a,b,d)
    (b==0 || d==0) && return sgn(a)
    (a==0 || sgn(a)==sgn(b)) && return sgn(b)
    sgn(a)*sgn(a*a-b*b*d)
end
function project(row, normals)
    for (n,len) in normals
        c = dot(row,n)
        c != 0 && (row = shift(row,n,-c/len))
    end
    row
end
function affine_space(eqs)
    point = ZERO; normals = Tuple{V,QQ}[]
    for (row,rhs) in eqs
        res = rhs - dot(row,point)
        n = project(row, normals)
        if n == ZERO
            res != 0 && return nothing
            continue
        end
        len = dot(n,n)
        point = shift(point,n,res/len)
        push!(normals,(n,len))
    end
    (point, Tuple(normals))
end
function stratum(eqs)
    planes = Tuple{V,QQ}[(b,c) for (a,b,c) in eqs if a==0]
    spheres = Tuple{V,QQ}[(scale(b,1/a), c/a) for (a,b,c) in eqs if a!=0]
    local row, rhs
    if !isempty(spheres)
        row, rhs = spheres[1]
        for (b,c) in spheres[2:end]; push!(planes,(shift(b,row,-1), c-rhs)); end
    end
    aff = affine_space(planes)
    aff === nothing && return nothing
    point, normals = aff
    if isempty(spheres)
        return length(normals)==4 ? (point, (), QQ(0)) : nothing
    end
    centre = scale(row, QQ(-1,2))
    diff = shift(centre, point, -1)
    tangent = project(diff, normals)
    perp = shift(diff, tangent, -1)
    r2 = rhs + dot(centre,centre) - dot(perp,perp)
    (r2 < 0 || (r2 > 0 && length(normals)==4)) && return nothing
    (shift(point,tangent,1), r2 != 0 ? normals : (), r2)
end
struct ExteriorTest
    constraints::Vector{Tuple{QQ,V,QQ}}
    strata::Set{Any}
end
function ExteriorTest(shell)
    cons = Tuple{QQ,V,QQ}[(QQ(1), ZERO, QQ(2,3))]
    for y in shell; push!(cons, (4-dot(y,y), scale(y,QQ(4)), QQ(4))); end
    strata = Set{Any}()
    n = length(cons)
    for i in 1:n
        d = stratum([cons[i]]); d===nothing || push!(strata,d)
        for j in i+1:n
            d = stratum([cons[i],cons[j]]); d===nothing || push!(strata,d)
            for k in j+1:n
                d = stratum([cons[i],cons[j],cons[k]]); d===nothing || push!(strata,d)
                for l in k+1:n
                    d = stratum([cons[i],cons[j],cons[k],cons[l]]); d===nothing || push!(strata,d)
                end
            end
        end
    end
    ExteriorTest(cons, strata)
end
function feasible(t::ExteriorTest, p, v, d)
    sq = dot(p,p)+dot(v,v)*d; cr = 2dot(p,v)
    for (a,b,c) in t.constraints
        radical_sign(c-a*sq-dot(b,p), -a*cr-dot(b,v), d) < 0 && return false
    end
    true
end
function cutting_witness(t::ExteriorTest, x)
    for (p,normals,r2) in t.strata
        local dirs, d
        if r2 == 0
            dirs = (ZERO,); d = QQ(0)
        else
            v = ZERO
            for f in (x, AXES...)
                v = project(f, normals)
                v != ZERO && break
            end
            d = r2/dot(v,v)
            dirs = length(normals)==3 ? (v, scale(v,QQ(-1))) : (v,)
        end
        excess = dot(x,p)-1
        for v in dirs
            radical_sign(excess, dot(x,v), d) > 0 && feasible(t,p,v,d) && return (p,v,d)
        end
    end
    nothing
end
function verify_packing(shell)
    all(4 <= dot(y,y) < 6 for y in shell) || return false
    for i in 1:length(shell), j in i+1:length(shell)
        d = shift(shell[i], shell[j], -1)
        dot(d,d) >= 4 || return false
    end
    true
end
function rational_roots()
    r = V[]
    for ax in AXES, s in (-1,1); push!(r, scale(ax, QQ(s))); end
    for s in Iterators.product((-1,1),(-1,1),(-1,1),(-1,1)); push!(r, ntuple(i->QQ(s[i],2),4)); end
    r
end
function layer_shell(n)
    eq = V[]
    for i in 1:3, j in i+1:3, s in (-1,1), t in (-1,1)
        push!(eq, shift(scale(AXES[i], QQ(17s,10)), AXES[j], QQ(17t,10)))
    end
    # itertools.product order: last coordinate varies fastest
    sg = [(a,b,c) for a in (-1,1) for b in (-1,1) for c in (-1,1)]
    upper = V[(QQ(s[1]),QQ(s[2]),QQ(s[3]),QQ(17,10)) for s in sg]
    lower = V[(QQ(s[1]),QQ(s[2]),QQ(s[3]),QQ(-17,10)) for s in sg]
    vcat(eq, upper, lower[1:n-20])
end
function run(io::IO=stdout)
    R = Common.Report("anc/exterior_test")
    a = QQ(1001,1000); c = QQ(99999,100000)
    shell = [scale(u,2a) for u in rational_roots()]
    ok = verify_packing(shell)
    t = ExteriorTest(shell)
    vertex = (a, 2c-a, 2(a-c), QQ(0))
    ok &= all(2dot(vertex,y) <= dot(y,y) for y in shell)
    ok &= cutting_witness(t, vertex) === nothing
    ok &= cutting_witness(t, (a,a,QQ(0),QQ(0))) !== nothing
    Common.check!(R, io, "24 neighbours: clipped orbit protected; untrimmed vertex cut", ok; got="strata=$(length(t.strata))")
    lower = 8a^4 - 192(a-c)^4
    delta = QQ(1,2000)
    Common.check!(R, io, "24 neighbours: exact clipped volume > 8.032; perturbations of 1/2000 give > 8.007",
        lower > QQ(8032,1000) && (1-QQ(3,2)*delta)^4*lower > QQ(8007,1000); got=string(Float64(lower)), expected="> 8.032")
    for n in (25,26,27)
        shell = layer_shell(n)
        ok = verify_packing(shell)
        t = ExteriorTest(shell)
        near = scale(AXES[4], QQ(5,4)); far = scale(AXES[4], QQ(3,2))
        ok &= all(2dot(far,y) <= dot(y,y) for y in shell)
        ok &= cutting_witness(t, near) === nothing
        ok &= cutting_witness(t, far) !== nothing
        ok &= all(dot(y,y) >= (QQ(12,5)+QQ(1,500))^2 for y in shell)
        ok &= QQ(1,2)*QQ(314,100)^2*QQ(6,5)^4 > QQ(51,5)
        Common.check!(R, io, "$n neighbours: x=(0,0,0,5/4) protected; x=(0,0,0,3/2) cut; perturbations of 1/500 give volume > 10.2", ok; got="strata=$(length(t.strata))")
    end
    R
end
end
