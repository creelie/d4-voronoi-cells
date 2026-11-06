"""c27_config_search.py -- local search for packing sets with prescribed bin counts: centres in
R^4 with |y| in given ranges and pairwise distance >= 2, by penalised L-BFGS from starts near a
root system (even tries) and random starts (odd tries).  Floating point; a positive best
penalty is evidence, not a proof, that no such set exists.
usage: python3 c27_config_search.py COUNT:LO:HI[,COUNT:LO:HI...] TRIES [SEED]
   e.g. 10:2:2.0161,11:2.0161:2.05,3:2.05:2.1,3:2.35:9 200 1   (HI is capped at sqrt6)"""
import numpy as np, sys
from scipy.optimize import minimize
spec=[tuple(map(float,g.split(':'))) for g in sys.argv[1].split(',')]  # count:lo:hi
tries=int(sys.argv[2]); seed=int(sys.argv[3]) if len(sys.argv)>3 else 0
R6=np.sqrt(6.0)
lo=np.concatenate([[a]*int(c) for c,a,b in spec]); hi=np.concatenate([[min(b,R6)]*int(c) for c,a,b in spec]); N=len(lo)
n1=int(sum(c for c,a,b in spec if b<2.3)); n2=N-n1
iu=np.triu_indices(N,1)
def f(x):
    X=x.reshape(N,4)
    D=X[:,None,:]-X[None,:,:]; d2=(D**2).sum(-1)
    v=np.maximum(0,4-d2)[iu]
    r2=(X**2).sum(1)
    a=np.maximum(0,lo**2-r2); b=np.maximum(0,r2-hi**2)
    val=(v**2).sum()+(a**2).sum()+(b**2).sum()
    # grad
    G=np.zeros_like(X)
    W=np.zeros((N,N)); W[iu]=v; W=W+W.T
    G+= -4*(W[:,:,None]*D).sum(1)
    G+= (-4*a+4*b)[:,None]*X
    return val,G.ravel()
rng=np.random.default_rng(seed)
# 24-cell roots
import itertools
roots=[]
for i,j in itertools.combinations(range(4),2):
    for s in (1,-1):
        for t in (1,-1):
            v=np.zeros(4); v[i]=s; v[j]=t; roots.append(v*np.sqrt(2))
roots=np.array(roots)
best=1e9
for t in range(tries):
    if t%2==0 and n1<=24:
        X1=roots[rng.permutation(24)[:n1]]*(1+0.05*rng.random())+0.3*rng.standard_normal((n1,4))
    else:
        X1=rng.standard_normal((n1,4)); X1*=2.05/np.linalg.norm(X1,axis=1)[:,None]
    X2=rng.standard_normal((n2,4)); X2*=2.4/np.linalg.norm(X2,axis=1)[:,None]
    x0=np.vstack([X1,X2]); x0=x0/np.linalg.norm(x0,axis=1)[:,None]*((lo+hi)/2)[:,None]; X1=x0[:n1]; X2=x0[n1:]
    x0=np.vstack([X1,X2]).ravel()
    r=minimize(f,x0,jac=True,method='L-BFGS-B',options={'maxiter':20000,'ftol':1e-16,'gtol':1e-14})
    if r.fun<best:
        best=r.fun
        X=r.x.reshape(N,4); d=np.sqrt(((X[:,None]-X[None])**2).sum(-1))[iu].min(); rr=np.linalg.norm(X,axis=1)
        print(t,'penalty %.3e mindist %.6f'%(r.fun,d),np.round(np.sort(rr),4).tolist(),flush=True)
        if r.fun<1e-20: np.save('found_%s.npy'%sys.argv[1].replace(':','_').replace(',','-'),X); break
print("best",best)
