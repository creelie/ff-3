import sys, itertools
from fractions import Fraction as Fr
import numpy as np

import sympy as sp
from mf import E
from Gmu import forms, mul
def theta3(d,Q):
    v=[0]*Q
    for n in range(-40,41):
        if d*n*n<Q: v[d*n*n]+=1
    return [Fr(x) for x in v]
def theta2_4(d,Q):  # theta2(4 d tau)^2/4 = q^d (sum q^{4d n(n+1)})^2 ; generalised: sum over n of q^{d(2n+1)^2/ ...}
    # theta2(4d tau) = sum_n q^{d(2n+1)^2/2}? theta2(t)=sum q^{(n+1/2)^2/2}; at 4d tau: q^{2d(n+1/2)^2}=q^{d(2n+1)^2/2}
    # so theta2(4d tau)^2 has powers d((2a+1)^2+(2b+1)^2)/2 integral; /4 normalises leading
    v=[0]*Q
    for a in range(-40,41):
        for b in range(-40,41):
            e=d*((2*a+1)**2+(2*b+1)**2)//2
            if e<Q: v[e]+=1
    return [Fr(x,4) for x in v]
def Eodd(d,Q):
    v=[Fr(0)]*Q; v[0]=Fr(1)
    for n in range(1,Q):
        if n%d==0:
            m=n//d; v[n]=Fr(24*sum(e for e in range(1,m+1) if m%e==0 and e%2))
    return v
def blocks(N,Q):
    dN=[d for d in range(1,N+1) if N%d==0]
    w1=[]  # weight-1 theta-type with chi_{-4}
    for d in dN: w1.append(("t3(%d)^2"%d,mul(theta3(d,Q),theta3(d,Q),Q)))
    for d in dN: w1.append(("t2(4*%d)^2/4"%d,theta2_4(d,Q)))
    w2=[]
    for i in range(len(w1)):
        for j in range(i,len(w1)):
            w2.append((w1[i][0]+"*"+w1[j][0],mul(w1[i][1],w1[j][1],Q)))
    for d in [d for d in range(1,2*N+1) if (2*N)%d==0]: w2.append(("Ee(%d)"%d,Eodd(d,Q)))
    w4=[]
    for i in range(len(w2)):
        for j in range(i,len(w2)):
            w4.append((w2[i][0]+" | "+w2[j][0],mul(w2[i][1],w2[j][1],Q)))
    for d in [d for d in range(1,4*N+1) if (4*N)%d==0]: w4.append(("E4(%d)"%d,E(4,d,Q)))
    return w4
def solve(N,target,Q,K):
    B=blocks(N,Q)
    A=np.array([[float(b[1][n]) for b in B] for n in range(K)])
    t=np.array([float(x) for x in target[:K]])
    # scale rows
    s=np.maximum(1,np.abs(A).max(axis=1)); A2=A/s[:,None]; t2=t/s
    from scipy.optimize import linprog
    r=linprog(np.ones(len(B)),A_eq=A2,b_eq=t2,bounds=(0,None),method="highs")
    if r.status!=0: return None,r.message
    supp=[i for i,x in enumerate(r.x) if x>1e-9]
    M=sp.Matrix([[sp.Rational(B[i][1][n].numerator,B[i][1][n].denominator) for i in supp] for n in range(Q)])
    tt=sp.Matrix([sp.Rational(x.numerator,x.denominator) for x in target[:Q]])
    try: sol,par=M.gauss_jordan_solve(tt)
    except ValueError: return None,"exact inconsistent on support"
    if par.shape[0]: sol=sol.subs({p:0 for p in par})
    return [(B[i][0],sol[k]) for k,i in enumerate(supp)],"ok"
if __name__=="__main__":
    N=int(sys.argv[1]); which=sys.argv[2]; Q=60
    b,c,G0,G1=forms(N,Q)
    tgt={"G0":G0,"G1":G1,"c":c}[which]
    res,msg=solve(N,tgt,Q,40)
    print(msg)
    if res:
        for nm,v in res: print("  ",v,nm)
        print("all nonneg:",all(v>=0 for _,v in res))
