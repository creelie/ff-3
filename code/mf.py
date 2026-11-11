import sys
from fractions import Fraction as Fr
import sympy as sp
from closed import run
def sig(n,k): return sum(d**k for d in range(1,n+1) if n%d==0)
def E(k,d,Q):
    c={2:-24,4:240,6:-504}[k]
    v=[Fr(0)]*Q; v[0]=Fr(1)
    for n in range(1,Q):
        if n%d==0: v[n]=Fr(c*sig(n//d,k-1))
    return v
def eta_prod(shape,Q):
    # prod eta(d tau)^e as q-series (assume total q-power integer, returns shifted series starting q^{sum d e/24})
    s=sum(d*e for d,e in shape.items()); assert s%24==0
    v=[Fr(0)]*Q; v[s//24 if s//24<Q else 0]=Fr(1) if s//24<Q else 0
    for d,e in shape.items():
        for n in range(1,Q):
            if d*n>=Q: break
            # multiply by (1-q^{dn})^e
            for _ in range(abs(e)):
                if e>0:
                    for m in range(Q-1,d*n-1,-1): v[m]-=v[m-d*n]
                else:
                    for m in range(d*n,Q): v[m]+=v[m-d*n]
    return v
