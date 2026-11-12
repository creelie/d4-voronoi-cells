# Ports of anc/annular_constants.py, anc/cap_comparison.py, anc/cap_profiles_26.py
# Exact rational arithmetic; transcendental quantities are enclosed in rational intervals.
module AncAnnular
using ..Common
using ..Common: RI, QQ, sqrt_enclosure, atan_enclosure, pi_ri, decimal_interval
const F = QQ

function run_annular_constants(io::IO=stdout)
    R = Common.Report("anc/annular_constants")
    atanb(x) = (l = sum((-1)^j * x^(2j+1)//(2j+1) for j in 0:7); (l, l + x^17//17))
    a, a_upper = F(12247,10000), F(122475,100000)
    q_lower, q_upper = F(12649,10000), F(1265,1000)
    sl, su = F(77,100), F(31,40)
    ok = a^2 < F(3,2) < a_upper^2 && q_lower^2 < F(8,5) < q_upper^2 && sl^2 < F(3,5) < su^2
    rho, b = F(5,4), F(1,4)
    loss2 = (F(2,5)*F(96,95)^2)^2 * (2rho*(rho-a))^3 / (a^6*b^6)
    exp1 = big"22535460391546257408000000000"//big"439736697911706067259903971009"
    Common.check!(R, io, "exterior removal fraction < 23/100", ok && loss2==exp1 && loss2 < F(23,100)^2; got=string(loss2), expected=string(exp1))
    tail = 2*F(314,100)*(q_lower-a_upper)^4/su^3
    exp2 = big"65277247497517"//big"1861937500000000000"
    Common.check!(R, io, "vertex-cone slope < 1/20 and volume correction > 1/28600",
        (q_upper-a)/(a*sl) < F(1,20) && tail==exp2 && tail > F(1,28600); got=string(tail), expected=string(exp2))
    l5,u5 = atanb(F(1,5)); l239,u239 = atanb(F(1,239))
    pl, pu = F(31415926,10000000), F(31415927,10000000)
    den = 9*pl^2/8 - F(33274,10000) + F(1,28600)
    dens = pl^2/(2den)
    Common.check!(R, io, "3.1415926 < pi < 3.1415927 and density bound < 0.634625",
        pl < 16l5-4u239 && 16u5-4l239 < pu && 0 < den < 8 && dens < F(634625,1000000); got=string(Float64(dens)), expected="< 0.634625")
    R
end

function single_cap(radius2::F, h::F)
    r4 = radius2^2
    perp = sqrt_enclosure(radius2 - h^2)
    angle = atan_enclosure(perp / h)
    (4*pi_ri()/3) * (F(3,8)*r4*angle - h*(5radius2-2h^2)*perp/8)
end
function cap_data(radius2::F, h::F)
    r4 = radius2^2
    tangent = sqrt_enclosure(radius2-h^2)/h
    angle = atan_enclosure(tangent)
    prim(arg, theta) = r4*theta - 2radius2*h^2*arg + h^4*(arg + arg*arg*arg/3)
    pair = pi_ri()/2 * (prim(tangent, angle) - prim(sqrt_enclosure(F(1,3)), pi_ri()/6))
    (pi_ri()*pi_ri()*r4/2, single_cap(radius2,h), pair)
end

function run_cap_comparison(io::IO=stdout)
    R = Common.Report("anc/cap_comparison")
    dil = F(53,40)
    ok = (dil^2+8dil)/9 < F(11,8) && F(32*27*6^3, 105*4^7)/F(314,100) < F(7,200) && 5*F(7,200)/3==F(7,120) && F(29,9)<F(15,4)
    Common.check!(R, io, "pair dilation, triple-to-pair ratio and five-completion constants", ok)
    expected = Dict(F(1)=>"T=[7.90694015, 7.90694016] one-radius lower expression=[7.91466037, 7.91466038] two-radius lower expression=[7.93603476, 7.93603477]",
                    F(10031,10000)=>"T=[7.99621535, 7.99621536] one-radius lower expression=[8.00086993, 8.00086994] two-radius lower expression=[8.02919724, 8.02919725]")
    for h in (F(1), F(10031,10000))
        b0,s0,p0 = cap_data(F(3,2),h); b1,s1,p1 = cap_data(F(25,16),h)
        trunc_ = b0 - 24*s0 + 96*p0
        common = b1 - b0 - 24*(s1-s0)
        beta = dil^3 - 1 - F(7,120)
        coarse = trunc_ + F(77,100)*(common + 96*beta*p0)
        sharp = trunc_ + F(77,100)*(common + 96*(p1 - F(127,120)*p0))
        ok = h==1 ? sharp.hi < F(79361,10000) : (trunc_.hi < F(79963,10000) && coarse.lo > F(80008,10000) && sharp.lo > F(80291,10000))
        got = "T=$(decimal_interval(trunc_)) one-radius lower expression=$(decimal_interval(coarse)) two-radius lower expression=$(decimal_interval(sharp))"
        Common.check!(R, io, "h=$(h) thresholds and logged intervals", ok && got==expected[h]; got=got, expected=expected[h])
    end
    R
end

function compositions(total::Int, bins::Int)
    out = Vector{Vector{Int}}()
    function rec(pre, rem, b)
        if b==1; push!(out, vcat(pre, rem)); return; end
        for f in 0:rem; rec(vcat(pre,f), rem-f, b-1); end
    end
    rec(Int[], total, bins); out
end

function run_cap_profiles_26(io::IO=stdout)
    R = Common.Report("anc/cap_profiles_26")
    dec(s) = parse(BigInt, replace(s, "."=>"", "-"=>"")) // big(10)^(length(s)-findfirst('.',s)) * (startswith(s,"-") ? -1 : 1)
    edges = dec.(["2.0","2.0161","2.05","2.1","2.15","2.2","2.444"])
    POINT = dec.(["0.10730785","0.10730646","-0.05347062","-0.05363574","-0.05425642","-0.05347062","-0.03946852"])
    PAIR = Dict(k=>dec(v)+F(1,10^7) for (k,v) in zip(("AA","AB","AF","BB","BF","FF"),("-0.02327176","-0.01036330","-0.01503736","-0.001713338","-0.004512281","-0.003681549")))
    TRIPLE = Dict(k=>dec(v)+F(1,10^7) for (k,v) in zip(("AAA","AAB","AAF","ABB","ABF","AFF","BBB","BBF","BFF"),
        ("0.001497925","0.0009598990","0.001322105","0.0005966748","0.001043879","0.001145076","0.0004094797","0.0006564131","0.001043982")))
    T = dec("8.282667610")
    cache = Dict{NTuple{3,Int},F}()
    function angular(a,b,f)
        get!(cache,(a,b,f)) do
            @assert f <= 2
            cnt = Dict('A'=>a,'B'=>b,'F'=>f)
            tot = T/2
            for (w,c) in Iterators.flatten((PAIR,TRIPLE))
                m = big(1)
                for l in Set(w); m *= binomial(big(cnt[l]), count(==(l), w)); end
                tot += m*c
            end
            tot
        end
    end
    grid = big(10)^10
    roundup(v) = cld(numerator(v)*grid, denominator(v))//grid
    caps = [roundup(single_cap(F(3,2), e/2).hi) for e in edges]
    target = (F(9,8)*pi_ri()*pi_ri() - 8).lo
    best = (F(0), Int[]); np = 0; excl = 0
    for n in compositions(26,7)
        cum = cumsum(n)
        (n[1] > 24 || any(cum[k+1] < m for (k,m) in enumerate((15,20,22,23,24)))) && continue
        np += 1
        a,b,f = cum[2], cum[6]-cum[2], n[7]
        cert = angular(a,b,f) + sum(x*m for (x,m) in zip(n,POINT))
        co = sum(x*c for (x,c) in zip(n,caps))
        j = min(cert,co)
        j < target && (excl += 1)
        j > best[1] && (best = (j, n))
    end
    Common.check!(R, io, "U < 3.288 on all 17845 count profiles; largest at (20,0,5,0,0,0,1)",
        np==17845 && best[1] < dec("3.288") && best[2]==[20,0,5,0,0,0,1]; got="profiles=$np max=$(Float64(best[1])) at $(Tuple(best[2]))", expected="17845 profiles, max at (20, 0, 5, 0, 0, 0, 1)")
    Common.check!(R, io, "profiles clearing T > 8 by these two bounds alone", excl==12378; got=excl, expected=12378)
    pm = (15,0,9,0,0,0,2)
    cs = sum(x*c for (x,c) in zip(pm,caps))
    Common.check!(R, io, "former maximising profile (15,0,9,0,0,0,2) has cap sum < 3.027", cs < dec("3.027"); got=Float64(cs), expected="< 3.027")
    l25 = dec("1.819065983")/2 + 23*dec("0.097085101") + 2*dec("0.087814019")
    Common.check!(R, io, "25-neighbour level < 3.3182 if at most 23 lie within 2.0161", l25 < dec("3.3182"); got=Float64(l25), expected="< 3.3182")
    R
end
end
