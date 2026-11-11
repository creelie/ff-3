#!/usr/bin/env python3
"""Exact checks for "Single-Centred Black Holes in the Composite CHL Models Have Positive Index".

V1  product formula reproduces Sen's tables for N=2,3 at P^2=2
V2  polar coefficient ghat(1)=16,7,6 and f = eta product of the cycle shape
V3  closed form f psi^F = b phi01(N tau) + c' phim2(N tau) - 2 ghat H(N tau) to q^Q
V4  phi01 = 12 theta2(z)^2/theta2^2 - 2 Eodd phim2 to q^Q
V5  G_mu equal the nonnegative theta/Eisenstein combinations of reps.json to q^100
V6  (-1)^(l+1) c^F(n,l) > 0 for every Delta>0 with n <= Q (product formula)
V7  domination inequality for mu=1, n <= 3000, and exact P_1 for small n
"""
import sys, json, math
from fractions import Fraction as Fr
from collections import defaultdict
from psi import finite
from closed import run, H
from chl import weak_tables
from Gmu import forms, mul
from lp import blocks
from mf import eta_prod
fails=[]
def check(name,ok):
    print(("ok   " if ok else "FAIL ")+name)
    if not ok: fails.append(name)
SHAPE={2:{1:8,2:8},3:{1:6,3:6},4:{1:4,2:2,4:4},6:{1:2,2:2,3:2,6:2},8:{1:2,2:1,4:1,8:2}}
GHAT={2:52,3:27,4:16,6:7,8:6}
Q=int(sys.argv[1]) if len(sys.argv)>1 else 30
# V1
F2,_=finite(2,4); F3,_=finite(3,4)
check("V1 Z2 (Q^2,P^2)=(1,2): -c^F(1,0)=2164, c^F(1,1)=360", -F2[1].get(0)==2164 and F2[1].get(1)==360)
check("V1 Z3 (Q^2,P^2)=(2/3,2): -c^F(1,0)=540, c^F(1,1)=27", -F3[1].get(0)==540 and F3[1].get(1)==27)

# V4: phi01 = 12 theta2(tau,z)^2/theta2(tau)^2 - 2 Eodd(tau) phim2, on theta components:
# component mu (y^l, l = mu mod 2) of the right side: sum_{u != l mod 2} q^{(l^2+u^2)/4} [12/theta2^2 + (-1)^(l+1)(-2 Eodd)/eta^6]
def v4(Qv=40):
    c01,cm2=weak_tables(Qv+2)
    from Gmu import mul as m_
    K=4*Qv
    # quarter-power series: eta^6 = q^{1/4} prod(1-q^n)^6 ; theta2^2 = 4 q^{1/4} prod (1-q^{2n})^4/(1-q^n)^2 ... use 3 theta4^4 /4 = eta^6/theta2^2
    th4=[Fr(0)]*(Qv+1)
    for n in range(-20,21):
        if n*n<=Qv: th4[n*n]+=(1 if n%2==0 else -1)
    t44=m_(m_(th4,th4,Qv+1),m_(th4,th4,Qv+1),Qv+1)
    Eo=[Fr(1)]+[Fr(24*sum(d for d in range(1,n+1) if n%d==0 and d%2)) for n in range(1,Qv+1)]
    # A_mu * eta^6 = theta-half * (3 theta4^4 -/+ ...): mu=1: 3t44-2Eo ; mu=0: 3t44+2Eo  (times sum over u)
    eta6=[Fr(0)]*(Qv+1); eta6[0]=Fr(1)
    for k in range(1,Qv+1):
        for _ in range(6):
            for m in range(Qv,k-1,-1): eta6[m]-=eta6[m-k]
    inv=[Fr(0)]*(Qv+1); inv[0]=Fr(1)
    for n in range(1,Qv+1): inv[n]=-sum(eta6[k]*inv[n-k] for k in range(1,n+1))
    ok=True
    for mu,sg in ((1,-2),(0,2)):
        W=[3*t44[n]+sg*Eo[n] for n in range(Qv+1)]
        W=m_(W,inv,Qv+1)
        # theta half: sum_{u != mu} q^{(u^2 - (1-mu))/4 ...}; compare c01(D) for D = 4n - mu
        half=[Fr(0)]*(Qv+1)
        for u in range(-20,21):
            if (u-mu)%2:
                e=(u*u-(1-mu))//4   # mu=1: u even -> u^2/4 ; mu=0: u odd -> (u^2-1)/4
                if 0<=e<=Qv: half[e]+=1
        S=m_(half,W,Qv+1)
        for n in range(Qv-2):
            D=4*n-1 if mu==1 else 4*n
            if S[n]!=c01[D]: ok=False
    return ok
check("V4 phi01 = 12 theta2(z)^2/theta2^2 - 2 Eodd phim2 (theta components, D<160)", v4())

for N in (4,6,8):
    print("--- N =",N)
    F,D=finite(N,Q)
    b,c,f=run(N,Q)   # run asserts exact divisibility (V3 structure) and returns b, c', f
    ep=eta_prod(SHAPE[N],Q+2)
    check("V2 ghat(1)=%d and f = eta product to q^%d"%(GHAT[N],Q), D[0]==GHAT[N] and all(f[i]==ep[i+1] for i in range(Q)))
    bb,cc,G0,G1=forms(N,max(Q,101))
    check("V3 b, c' equal the stated Eisenstein/eta expressions to q^%d"%(Q-1), all(b[n]==bb[n] for n in range(Q)) and all(c[n]==cc[n] for n in range(Q)))
    reps=json.load(open("reps.json"))
    B=dict(blocks(N,101))
    for mu,G in ((0,G0),(1,G1)):
        rep=reps["%d_G%d"%(N,mu)]
        s=[Fr(0)]*101
        for nm,v in rep:
            v=Fr(v); assert v>=0
            for n in range(101): s[n]+=v*B[nm][n]
        check("V5 G_%d = nonnegative combination (%d terms) to q^100"%(mu,len(rep)), all(s[n]==G[n] for n in range(101)))
    # V6
    ok=True; cnt=0
    for n in range(0,Q-1):
        for l,v in F[n].items():
            if 4*n-N*l*l>0:
                cnt+=1
                if (-1)**(l+1)*v<=0: ok=False
    check("V6 (-1)^(l+1) c^F(n,l) > 0 for all %d coefficients with Delta>0, n<=%d"%(cnt,Q-2), ok)
    # V7
    mmin={4:2,6:2,8:3}[N]
    g=GHAT[N]
    gmin=min(G1[m] for m in range(mmin,mmin+N))
    check("V7 G_1(m)>0 for %d<=m<=100"%mmin, all(G1[m]>0 for m in range(mmin,101)) and all(G1[m]==0 for m in range(mmin)))
    ok=True
    for n in range(mmin,3001):
        j=(n-mmin)//N
        if math.comb(j+5,5)*gmin <= 2*g*Fr(4*n,N):
            # exact: E1 = sum_j a_j G1(n-Nj), a_j coefficients of theta3 * prod(1-q^k)^-6
            pass
    nb=64*g*N//gmin
    check("V7 binomial bound beats 2 ghat D for %d<n<=3000 (analytic for n>%d)"%(max(n for n in range(mmin,3001) if math.comb((n-mmin)//N+5,5)*gmin<=2*g*Fr(4*n,N)),nb),
          all(math.comb((n-mmin)//N+5,5)*gmin>2*g*Fr(4*n,N) for n in range(60,3001)))
    # exact P_1 for n<60
    K=60
    a=[Fr(0)]*K; a[0]=Fr(1)
    for k in range(1,K):
        for _ in range(6):
            for m in range(k,K): a[m]+=a[m-k]
    th=[Fr(0)]*K
    for u in range(-10,11):
        if u*u<K: th[u*u]+=1
    a=mul(th,a,K)
    ok=True
    for n in range(mmin,K):
        E1=sum(a[j]*G1[n-N*j] for j in range(0,n//N+1))
        Dl=4*n-N
        Ht=2*g*H(Dl//N) if n%N==0 else 0   # class-number term only when N | m
        if E1-Ht<=0: ok=False
    check("V7 P_1(Delta)=E_1-2 ghat H(Delta/N) > 0 exactly for n<%d"%K, ok)
print("ALL CHECKS PASSED" if not fails else "FAILED: %s"%fails)
