import sys, math
from fractions import Fraction as Fr
from collections import defaultdict
from psi import finite
sys.path.insert(0,'.')
from chl import weak_tables
def H(n):
    # Hurwitz class number, H(0)=-1/12
    if n==0: return Fr(-1,12)
    if n<0 or n%4 in (1,2): return Fr(0)
    tot=Fr(0)
    b=n%2
    while b*b<=n//3:
        ac=(b*b+n)//4
        a=max(b,1)
        while a*a<=ac:
            if ac%a==0:
                c=ac//a
                if a==b or a==c or b==0: w=Fr(1)
                else: w=Fr(2)
                if a==b==c: w=Fr(1,3)
                elif b==0 and a==c: w=Fr(1,2)
                tot+=w
            a+=1
        b+=2
    return tot
def run(N,Q):
    F,D=finite(N,Q)
    g=D[0]  # D = ghat/f ; D[0]=ghat (index shifted: D[n] coef of q^{n-1})
    # f = ghat / D(tau) ; D(tau)=q^{-1} sum D[n] q^n -> f = q * ghat/ sum D[n]q^n
    Dl=[D.get(n,0) for n in range(Q+2)]
    inv=[Fr(1,1)/Dl[0]]+[0]*(Q+1)
    for n in range(1,Q+2): inv[n]=-sum(Dl[k]*inv[n-k] for k in range(1,n+1))/Dl[0]
    f=[g*v for v in inv]  # f = q * sum f[n] q^n
    # f*psiF : coefficient of q^m = sum f[i]*F[m-1-i]
    L=int(math.isqrt(4*(Q+1)))+3
    G=defaultdict(lambda: defaultdict(Fr))
    for m in range(0,Q):
        for i in range(0,m+1):
            Fm=F.get(m-1-i)
            if Fm is None: continue
            for l,v in Fm.items(): G[m][l]+=f[i]*v
    # add 2 ghat H(N tau,z): sum H(4n-r^2) q^{Nn} y^r
    for n in range(0,Q//N+1):
        if N*n<Q:
            r=0
            for r in range(-L,L+1):
                h=H(4*n-r*r)
                if h: G[N*n][r]+=2*g*h
    # Phi0, Phi2 at N tau
    c01,cm2=weak_tables(Q//N+3)
    def ph(tab,n,r):
        D_=4*n-r*r
        return tab.get(D_,0) if D_>=-1 else 0
    # b(tau): z=0 value /12
    b=[sum(G[m].values())/12 for m in range(Q)]
    # remainder R = G - b*Phi0 ; c = R/Phi2
    R=defaultdict(lambda: defaultdict(Fr))
    for m in range(Q):
        for l,v in G[m].items(): R[m][l]+=v
        for i in range(0,m+1):
            if b[i]==0: continue
            nn=m-i
            if nn%N: continue
            n=nn//N
            for r in range(-L,L+1):
                v=ph(c01,n,r)
                if v: R[m][r]-=b[i]*v
    # divide R by Phi2(N tau,z) as series: Phi2 = sum_n P_n(y) q^{Nn}, P_0 = y-2+1/y
    c=[]
    Rem={m:dict(R[m]) for m in range(Q)}
    for m in range(Q):
        p={k:v for k,v in Rem[m].items() if v}
        # p should equal c_m*(y-2+1/y) + ...: leading
        cm=-p.get(0,0)/2 if p else Fr(0)
        c.append(cm)
        for nn in range(0,Q-m):
            if nn%N: continue
            n=nn//N
            for r in range(-L,L+1):
                v=ph(cm2,n,r)
                if v and m+nn<Q: Rem[m+nn][r]=Rem[m+nn].get(r,0)-cm*v
        assert all(v==0 for v in Rem[m].values()), (m,{k:v for k,v in Rem[m].items() if v})
    return b,c,f
if __name__=="__main__":
    N=int(sys.argv[1]);Q=int(sys.argv[2])
    b,c,f=run(N,Q)
    print("b",[str(x) for x in b]);print("c",[str(x) for x in c])
