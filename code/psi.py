import sys, math
from fractions import Fraction as Fr
from collections import defaultdict
from chl import CHL
# series: dict n -> dict j -> int ; truncated at q-order Q
def mul(A,B,Q):
    C=defaultdict(lambda: defaultdict(int))
    for n1,p1 in A.items():
        for n2,p2 in B.items():
            n=n1+n2
            if n>Q: continue
            c=C[n]
            for j1,v1 in p1.items():
                for j2,v2 in p2.items(): c[j1+j2]+=v1*v2
    return {n:{j:v for j,v in p.items() if v} for n,p in C.items()}
def powfac(l,j,a,Q):
    # (1 - q^l y^j)^(-a), l>=1
    out={0:{0:1}}
    if a==0: return out
    m=1; 
    # binomial series: sum_k C(a+k-1,k) x^k  (for a>0); general integer a: use (1-x)^(-a)
    k=1; coef=1
    while l*k<=Q:
        coef=coef*(a+k-1)//k if a>0 else None
        if a<0:
            coef=math.comb(-a,k)*(-1)**k
        out.setdefault(l*k,{})[j*k]=out.get(l*k,{}).get(j*k,0)+coef
        k+=1
    return out
def build(N,Q):
    C=CHL(N,Fr(Q+4))
    A0={}; A1={}; A2={}
    for l in range(0,Q+2):
        J=int(math.isqrt(4*l+4))+1
        for j in range(-J,J+1):
            if l>=1: 
                v=C.a(Fr(0),l,j)
                if v: A0[(l,j)]=v
            for (dic,kp) in ((A1,Fr(1,N)),(A2,Fr(2,N))):
                if 4*kp*l-j*j>=-1:
                    v=C.a(kp,l,j)
                    if v: dic[(l,j)]=v
    # P0 = q^{-1} prod_{l>=1}(1-q^l y^j)^{-a0}  (shifted: store q^{n} for n>=0 meaning q^{n-1})
    P={0:{0:1}}
    for (l,j),a in A0.items():
        if l<=Q+1: P=mul(P,powfac(l,j,a,Q+1),Q+1)
    X1=defaultdict(dict); X2=defaultdict(dict)
    for (l,j),a in A1.items():
        X1[l][j]=X1[l].get(j,0)+a
        if 2*l<=Q+1: X2[2*l][2*j]=X2[2*l].get(2*j,0)+Fr(a,2)
    for (l,j),a in A2.items(): X2[l][j]=X2[l].get(j,0)+a
    X1=dict(X1); X2=dict(X2)
    X11=mul(X1,X1,Q+1)
    S={}
    for n in set(X2)|set(X11):
        d=defaultdict(Fr)
        for j,v in X2.get(n,{}).items(): d[j]+=v
        for j,v in X11.get(n,{}).items(): d[j]+=Fr(v,2)
        S[n]={j:v for j,v in d.items() if v}
    PsiT=mul(P,S,Q+1)   # = psi_1*(y-2+1/y), index n means q^{n-1}
    # D(tau): set y=1
    Dser={n:sum(p.values()) for n,p in PsiT.items()}
    # polar*(y-2+1/y) = D * sum_s q^{N(s^2+s)} y^{2s} (1-y)^2/(1-q^{Ns}y)^2
    Pol=defaultdict(lambda: defaultdict(Fr))
    s_rng=range(-Q-2,Q+3)
    for s in s_rng:
        e0=N*(s*s+s)
        if e0>Q+1: continue
        # expand y^{2s}(1-y)^2/(1-q^{Ns}y)^2
        terms=defaultdict(lambda: defaultdict(Fr))
        if s==0:
            terms[0][0]+=1
        elif s>0:
            # sum_k (k+1) q^{Nsk} y^k
            k=0
            while e0+N*s*k<=Q+1:
                for jj,cc in ((0,1),(1,-2),(2,1)):
                    terms[N*s*k][2*s+k+jj]+=(k+1)*cc
                k+=1
        else:
            # 1/(1-q^{Ns}y)^2 = q^{-2Ns}y^{-2}/(1-q^{-Ns}y^{-1})^2 = sum_k (k+1) q^{-Ns(k+2)} y^{-(k+2)}
            k=0
            while e0-N*s*(k+2)<=Q+1:
                for jj,cc in ((0,1),(1,-2),(2,1)):
                    terms[-N*s*(k+2)][2*s+jj-(k+2)]+=(k+1)*cc
                k+=1
        for n1,p in terms.items():
            for n2,dv in Dser.items():
                n=e0+n1+n2
                if n<=Q+1:
                    for j,v in p.items(): Pol[n][j]+=v*dv
    return PsiT,Pol,Dser,A1
def finite(N,Q):
    PsiT,Pol,Dser,A1=build(N,Q)
    F={}
    for n in range(0,Q+1):
        d=defaultdict(Fr)
        for j,v in PsiT.get(n,{}).items(): d[j]+=v
        for j,v in Pol.get(n,{}).items(): d[j]-=v
        d={j:v for j,v in d.items() if v}
        # divide by (y-2+1/y)
        if not d: F[n-1]={}; continue
        lo=min(d); hi=max(d); q={}; r=dict(d)
        for e in range(hi,lo+1,-1):
            c=r.get(e,0)
            if c:
                q[e-1]=c; r[e]=r.get(e,0)-c; r[e-1]=r.get(e-1,0)+2*c; r[e-2]=r.get(e-2,0)-c
        if any(v for v in r.values()): 
            F[n-1]=None; continue
        F[n-1]=q
    return F,Dser
if __name__=="__main__":
    N=int(sys.argv[1]); Q=int(sys.argv[2])
    F,D=finite(N,Q)
    print("D:",[D.get(n,0) for n in range(6)])
    for n in range(-1,Q):
        if F[n] is None: print(n,"NOT DIVISIBLE (truncation?)"); continue
        print(n,[ (l,(-1)**(l+1)*F[n].get(l,0)) for l in range(0,6) if 4*n-N*l*l>0 or F[n].get(l,0)])
