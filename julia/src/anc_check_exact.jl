# Port of anc/check_exact.py : exact rational/integer identities.
module AncCheckExact
using ..Common
const F = Rational{BigInt}
function run(io::IO=stdout)
    R = Common.Report("anc/check_exact")
    roots = Set{NTuple{4,Int}}()
    for i in 1:4, j in i+1:4, a in (-1,1), b in (-1,1)
        r = [0,0,0,0]; r[i]=a; r[j]=b; push!(roots, Tuple(r))
    end
    vertices = Set{NTuple{4,F}}()
    for i in 1:4, s in (-1,1)
        v = F[0,0,0,0]; v[i]=s; push!(vertices, Tuple(v))
    end
    for s in Iterators.product((-1,1),(-1,1),(-1,1),(-1,1))
        push!(vertices, Tuple(F(x,2) for x in s))
    end
    dot4(a,b) = sum(a[i]*b[i] for i in 1:4)
    ok = length(roots)==24 && length(vertices)==24
    for v in vertices
        ok &= maximum(dot4(v,r) for r in roots)==1
        ok &= count(r->dot4(v,r)==1, roots)==6
    end
    for r in roots
        ok &= count(v->dot4(v,r)==1, vertices)==6
    end
    Common.check!(R, io, "24 root directions, 24 dual facets, six incidences at each end", ok)
    had = ((1,1,1,1),(1,1,-1,-1),(1,-1,1,-1),(1,-1,-1,1))
    img = Set(Tuple(sum(F(h,2)*rj for (h,rj) in zip(row,r)) for row in had) for r in roots)
    q(t) = F(4,3)-2t+F(4,3)*t*t-F(2305,3)*t*t*t
    ok = img == Set(Tuple(F.(r)) for r in roots) && q(F(1,10))==F(227,600) && F(227,600)>F(3,8)
    Common.check!(R, io, "Hadamard root symmetry and exact radial coefficient 227/600", ok; got=string(q(F(1,10))), expected="227//600")
    coefficient = F(3,8)-F(1424,12500)
    Common.check!(R, io, "angular comparison coefficient > 1/4", coefficient>F(1,4); got=string(coefficient), expected="6527//25000")
    ok = true; n=0
    for k in 0:11
        a = F(1, big(2)^k); an = a/2
        for j in 1:100
            ok &= an*max(an,F(1,j)) <= a*max(a,F(1,j)); n+=1
        end
    end
    Common.check!(R, io, "transfer comparison monotonicity at 1200 rational entries", ok && n==1200; got=n, expected=1200)
    r = F(1,2); v=(1,-1)
    Common.check!(R, io, "cross-cell mask has a negative rational quadratic form", 2r*v[1]*v[2] < 0)
    orbs = Set{Any}()
    for a in 0:10, b in 0:10, c in 0:10
        t=(a,b,c)
        orbit = Set([t,(b,c,a),(c,a,b)])
        d = sum((x+1)*(y+1)*(z+1) for (x,y,z) in orbit)
        d<=10 && push!(orbs, (Tuple(sort(collect(orbit))), d))
    end
    expected = Set{Any}([(((0,0,0),),1), (((1,1,1),),8), (((0,0,1),(0,1,0),(1,0,0)),6), (((0,0,2),(0,2,0),(2,0,0)),9)])
    Common.check!(R, io, "complete small highest-weight orbit enumeration", orbs==expected; got=length(orbs), expected=4)
    ok = true; cnt = 0
    for n in 1:10
        for m in 0:(2^n-1)
            signs = [((m>>(j-1))&1)==1 ? 1 : -1 for j in 1:n]
            pre = cumsum(vcat(0,signs))
            lhs = 0
            for t in -n:n-1
                s = 0
                for j in 1:n
                    (t<j<=t+n) && (s+=signs[j])
                end
                lhs += s^2
            end
            rhs = sum(pre[j+1]^2+(pre[n+1]-pre[j+1])^2 for j in 0:n-1)
            ok &= lhs==rhs && 2lhs>=n*pre[n+1]^2; cnt+=1
        end
    end
    Common.check!(R, io, "signed-window energy and prefix identity for all 2046 sign strings of lengths 1-10", ok && cnt==2046; got=cnt, expected=2046)
    R
end
end
