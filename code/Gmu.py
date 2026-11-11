import sys
from fractions import Fraction as Fr
from mf import E, eta_prod
def sig(n,k): return sum(d**k for d in range(1,n+1) if n%d==0)
def scale(v,N,Q):
    w=[Fr(0)]*Q
    for n in range(Q):
        if n*N<Q: w[n*N]=v[n]
    return w
def mul(a,b,Q):
    c=[Fr(0)]*Q
    for i in range(Q):
        if a[i]:
            for j in range(Q-i): c[i+j]+=a[i]*b[j]
    return c
def lin(d,Q):
    out=[Fr(0)]*Q
    for v,co in d:
        for n in range(Q): out[n]+=co*v[n]
    return out
COEF={4:({2:Fr(1,72),4:Fr(5,72)},{1:Fr(1,120),2:Fr(-7,360),4:Fr(7,30)},None),
      6:({2:Fr(1,288),3:Fr(1,108),6:Fr(13,864)},{1:Fr(7,3600),2:Fr(-29,7200),3:Fr(-7,3600),6:Fr(81,800)},({1:2,2:2,3:2,6:2},Fr(1,5))),
      8:({2:Fr(0),4:Fr(1,288),8:Fr(1,96)},{1:Fr(7,11520),2:Fr(-1,768),4:Fr(1,1440),8:Fr(1,12)},({2:4,4:4},Fr(3,16)))}
def forms(N,Q):
    bd,cd,cu=COEF[N]
    E21=E(2,1,Q)
    b=lin([([E21[n]-d*E(2,d,Q)[n] for n in range(Q)],co) for d,co in bd.items()],Q)
    c=lin([(E(4,d,Q),co) for d,co in cd.items()]+([(eta_prod(cu[0],Q),cu[1])] if cu else []),Q)
    th4=[Fr(0)]*Q
    for n in range(-20,21):
        if n*n<Q: th4[n*n]+=(1 if n%2==0 else -1)
    th4N=scale(th4,N,Q); th44=mul(mul(th4N,th4N,Q),mul(th4N,th4N,Q),Q)
    Eodd=[Fr(1)]+[Fr(24*sum(d for d in range(1,n+1) if n%d==0 and d%2)) for n in range(1,Q)]
    EN=scale(Eodd,N,Q)
    bth=mul(b,th44,Q); bE=mul(b,EN,Q)
    G1=[3*bth[n]+c[n]-2*bE[n] for n in range(Q)]
    G0=[-3*bth[n]+c[n]-2*bE[n] for n in range(Q)]
    return b,c,G0,G1
if __name__=="__main__":
    N=int(sys.argv[1]);Q=int(sys.argv[2])
    b,c,G0,G1=forms(N,Q)
    print("G0",[str(x) for x in G0]);print("G1",[str(x) for x in G1])
