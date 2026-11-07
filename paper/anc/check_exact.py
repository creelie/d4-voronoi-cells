from fractions import Fraction as F
from itertools import product, combinations, permutations

roots=set()
for i,j in combinations(range(4),2):
    for a,b in product((-1,1),repeat=2):
        r=[0]*4
        r[i],r[j]=a,b
        roots.add(tuple(r))
vertices=set()
for i in range(4):
    for s in (-1,1):
        v=[F(0)]*4
        v[i]=F(s)
        vertices.add(tuple(v))
vertices.update(tuple(F(s,2) for s in signs) for signs in product((-1,1),repeat=4))
assert len(roots)==len(vertices)==24
for v in vertices:
    assert max(sum(a*b for a,b in zip(v,r)) for r in roots)==1
    assert sum(sum(a*b for a,b in zip(v,r))==1 for r in roots)==6
for r in roots:
    assert sum(sum(a*b for a,b in zip(v,r))==1 for v in vertices)==6
hadamard=((1,1,1,1),(1,1,-1,-1),(1,-1,1,-1),(1,-1,-1,1))
assert {tuple(sum(F(h,2)*rj for h,rj in zip(row,r)) for row in hadamard) for r in roots}==roots
q=lambda t:F(4,3)-2*t+F(4,3)*t*t-F(2305,3)*t*t*t
assert q(F(1,10))==F(227,600)>F(3,8)
coefficient=F(3,8)-F(1424,12500)
assert coefficient>F(1,4)
print('PASS 24 root directions, 24 dual facets, six incidences at each end')
print('PASS Hadamard root symmetry and exact radial coefficient 227/600')
print('PASS angular comparison coefficient',coefficient,'> 1/4')

for k in range(12):
    a=F(1,2**k)
    a_next=a/2
    for j in range(1,101):
        assert a_next*max(a_next,F(1,j)) <= a*max(a,F(1,j))
r=F(1,2)
v=(1,-1)
masked_q=2*r*v[0]*v[1]
assert masked_q<0
print('PASS transfer comparison monotonicity at 1200 rational entries')
print('PASS cross-cell mask has a negative rational quadratic form')

orbit_dimensions=set()
for t in product(range(11),repeat=3):
    orbit={t,t[1:]+t[:1],t[2:]+t[:2]}
    dimension=sum((a+1)*(b+1)*(c+1) for a,b,c in orbit)
    if dimension<=10:
        orbit_dimensions.add((tuple(sorted(orbit)),dimension))
expected={(((0,0,0),),1),(((1,1,1),),8),(((0,0,1),(0,1,0),(1,0,0)),6),(((0,0,2),(0,2,0),(2,0,0)),9)}
assert orbit_dimensions==expected
print('PASS complete small highest-weight orbit enumeration')

for n in range(1,11):
    for signs in product((-1,1),repeat=n):
        prefixes=[0]
        for sign in signs:
            prefixes.append(prefixes[-1]+sign)
        lhs=sum(sum(signs[j-1] for j in range(1,n+1) if t<j<=t+n)**2 for t in range(-n,n))
        rhs=sum(prefixes[j]**2+(prefixes[n]-prefixes[j])**2 for j in range(n))
        assert lhs==rhs
        assert 2*lhs>=n*prefixes[n]**2
print('PASS signed-window energy and prefix identity for all 2046 sign strings of lengths 1–10')
