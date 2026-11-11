import sys, math
from fractions import Fraction as Fr
from mpmath import mp, mpc, exp, pi, nint
mp.dps=60
# ---------- weak Jacobi discriminant tables ----------
def weak_tables(NQ2):
    """Discriminant tables c01[D], cm2[D] of phi_{0,1}, phi_{-2,1} for D = 4n - j^2, n <= NQ2."""
    from collections import defaultdict
    Q=NQ2
    def mul(A,B):
        C=defaultdict(int)
        for (n1,j1),v1 in A.items():
            for (n2,j2),v2 in B.items():
                if n1+n2<=Q: C[(n1+n2,j1+j2)]+=v1*v2
        return {k:v for k,v in C.items() if v}
    def powser(n,j,e):
        out={(0,0):1}; k=1
        while n*k<=Q:
            out[(n*k,j*k)]=math.comb(e,k)*(-1)**k if e>0 else math.comb(-e+k-1,k)
            k+=1
        return out
    R={(0,0):1}
    for n in range(1,Q+1):
        for t in (powser(n,1,2),powser(n,-1,2),powser(n,0,-4)): R=mul(R,t)
    pm=mul({(0,1):1,(0,0):-2,(0,-1):1},R)
    S=defaultdict(int)
    for n in range(1,Q+1):
        for d in range(1,n+1):
            if n%d==0: S[(n,d)]+=d; S[(n,-d)]+=d; S[(n,0)]-=2*d
    phi=defaultdict(int)
    for k,v in pm.items(): phi[k]+=v
    for k,v in R.items(): phi[k]+=12*v
    for k,v in mul(pm,dict(S)).items(): phi[k]+=12*v
    c01={};cm2={}
    for n in range(Q+1):
        for j in (0,1):
            c01[4*n-j*j]=phi.get((n,j),0); cm2[4*n-j*j]=pm.get((n,j),0)
    return c01,cm2
def sigma(n): return sum(d for d in range(1,n+1) if n%d==0)
CLASS={'1A':(24,[]),'2A':(8,[(2,Fr(4,3))]),'3A':(6,[(3,Fr(3,2))]),
       '4B':(4,[(2,Fr(-1,3)),(4,Fr(2))]),'6A':(2,[(2,Fr(-1,6)),(3,Fr(-1,2)),(6,Fr(5,2))]),
       '8A':(2,[(4,Fr(-1,2)),(8,Fr(7,3))])}
POW={2:{0:'1A',1:'2A'},3:{0:'1A',1:'3A',2:'3A'},
     4:{0:'1A',1:'4B',2:'2A',3:'4B'},
     6:{0:'1A',1:'6A',2:'3A',3:'2A',4:'3A',5:'6A'},
     8:{0:'1A',1:'8A',2:'4B',3:'8A',4:'2A',5:'8A',6:'4B',7:'8A'}}
SHAPE={2:{1:8,2:8},3:{1:6,3:6},4:{1:4,2:2,4:4},6:{1:2,2:2,3:2,6:2},8:{1:2,2:1,4:1,8:2}}
def egcd(a,b):
    if b==0: return (a,1,0)
    g,x,y=egcd(b,a%b); return (g,y,x-(a//b)*y)
def E2series(alpha,beta,delta,N,Kmax):
    # E2((alpha tau+beta)/delta) as dict t->coef with t in (1/N)Z, t = alpha*n/delta ; up to t<=Kmax
    out={Fr(0):mpc(1)}
    n=1
    while Fr(alpha*n,delta)<=Kmax:
        t=Fr(alpha*n,delta)
        out[t]=out.get(t,mpc(0))-24*sigma(n)*exp(2j*pi*n*beta/delta)
        n+=1
    return out
def phiM_slash(M,a,b,c,d,N,Kmax):
    # (1/(M-1)) [ (M/delta^2) E2((alpha tau+beta)/delta) - E2(tau) ]
    A=(M*a,M*b,c,d)
    g,x,y=egcd(A[0],A[2])  # x*Ma + y*c = g
    alpha=abs(g); s=1 if g>0 else -1
    p,r=A[0]//g, A[2]//g   # gamma'' first column (p,r)
    # need second column (u,v) with p*v - r*u = 1 : from x*p*g + y*r*g = g -> x*p + y*r = 1 -> v=x, u=-y
    u,v=-y,x
    beta = v*A[1]-u*A[3]; delta=-r*A[1]+p*A[3]
    if s<0: alpha=-alpha  # not expected
    assert alpha>0 and delta>0 and alpha*delta==M, (alpha,delta,M)
    beta%=delta
    E=E2series(alpha,beta,delta,N,Kmax); E1=E2series(1,0,1,N,Kmax)
    out={}
    for t,v_ in E.items(): out[t]=out.get(t,0)+Fr(M,delta*delta)*v_/(M-1)
    for t,v_ in E1.items(): out[t]=out.get(t,0)-v_/(M-1)
    return out
def find_gamma(r,s,dd,N):
    # (0,dd)*(a b;c d) = (dd*c, dd*d) == (r,s) mod N
    for c in range(0,N*4):
        for d_ in range(-N*4,N*4):
            if (dd*c-r)%N==0 and (dd*d_-s)%N==0 and math.gcd(c,d_)==1:
                g,x,y=egcd(d_,c)  # x d - ? ; need a d - b c =1
                # solve a*d_ - b*c = 1
                g2,X,Y=egcd(d_,-c) if c!=0 else (d_,1 if d_==1 else -1,0)
                # X*d_ + Y*(-c) = g2 = ±1
                if abs(g2)!=1: continue
                a,b=X*g2, Y*g2
                assert a*d_-b*c==1
                return a,b,c,d_
    raise Exception("no gamma")
class CHL:
    def __init__(s,N,Kmax):
        s.N=N; s.Kmax=Kmax
        s.c01,s.cm2=weak_tables(int(Kmax)+3)
        s.T={}  # (r,s)-> dict of T|gamma ; chi
        for r in range(N):
            for t in range(N):
                dd=math.gcd(math.gcd(r,t),N)
                cls=POW[N][dd%N]
                chi,beta=CLASS[cls]
                if r==0:
                    a,b,c,d=1,0,0,1
                    # need (0,dd)->(0,t): F^{(0,t)} is twining by g^t of same class as g^dd: use identity gamma
                else:
                    a,b,c,d=find_gamma(r,t,dd,N)
                Tg={}
                for (M,bm) in beta:
                    for tt,v in phiM_slash(M,a,b,c,d,N,Kmax).items():
                        Tg[tt]=Tg.get(tt,0)+bm*v
                s.T[(r,t)]=(chi,Tg)
    def coef(s,r,t,n,j):
        # coefficient of q^n y^j in F^{(r,t)}, n in (1/N)Z
        chi,Tg=s.T[(r,t)]
        D=4*n-j*j
        val=mpc(0)
        if D.denominator==1 and int(D) in s.c01:
            val+=Fr(chi,12)*s.c01[int(D)]
        for tt,v in Tg.items():
            m=n-tt
            if m<0 or m.denominator!=1: continue
            D2=4*int(m)-j*j
            if D2 in s.cm2: val+=v*s.cm2[D2]
            elif D2>=-1: raise Exception("cm2 range %d"%D2)
        return val/s.N
    def a(s,kp,l,j):
        # exponent: sum_s e^{2 pi i s l/N} c^{(r,s)}(4k'l-j^2), k' in (1/N)Z
        N=s.N; r=int(kp*N)%N
        D=4*kp*l-j*j
        j0=j%2; n0=(D+j0*j0)/4
        tot=mpc(0)
        for t in range(N):
            tot+=exp(-2j*pi*t*l/N)*s.coef(r,t,n0,j0)
        v=complex(tot); iv=int(round(v.real))
        assert abs(v.real-iv)<1e-8 and abs(v.imag)<1e-8, (kp,l,j,v)
        return iv
if __name__=="__main__":
    N=int(sys.argv[1]); C=CHL(N,Fr(6))
    for r in range(N):
        print(r,[C.a(Fr(r,N) if r else Fr(0),l,1) for l in range(0,4)],[C.a(Fr(r,N) if r else Fr(0),l,0) for l in range(0,4)])
