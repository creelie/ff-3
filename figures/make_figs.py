"""Figures for the paper. Run from the repository root: python3 figures/make_figs.py"""
import sys, os, math, json
here=os.path.dirname(os.path.abspath(__file__)); root=os.path.dirname(here)
sys.path.insert(0,os.path.join(root,'code'))
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from fractions import Fraction as Fr
from psi import finite
from Gmu import forms
from closed import H
out=os.path.join(root,'paper')
GH={4:16,6:7,8:6}
cols={4:'#1f6fb4',6:'#d1495b',8:'#2e933c'}
# Figure 1: 3D view of log10 d* for N=4,6,8
Q=36
fig=plt.figure(figsize=(12,4.2))
data={}
for k,N in enumerate((4,6,8)):
    F,_=finite(N,Q)
    xs=[];ys=[];zs=[]
    for n in range(0,Q-1):
        for l,v in F[n].items():
            if l>=0 and 4*n-N*l*l>0:
                d=Fr((-1)**(l+1))*v
                xs.append(n); ys.append(l); zs.append(math.log10(float(d)))
                data.setdefault(str(N),[]).append([n,l,str(d)])
    ax=fig.add_subplot(1,3,k+1,projection='3d')
    ax.scatter(xs,ys,zs,c=zs,cmap='viridis',s=14,depthshade=True)
    ax.set_xlabel(r'$m=NQ^2/2$'); ax.set_ylabel(r'$\ell=Q\cdot P$'); ax.set_zlabel(r'$\log_{10} d^*$')
    ax.set_title(r'$\mathbb{Z}_%d$ CHL model, $P^2=2$'%N)
    ax.view_init(24,-58); ax.set_box_aspect(None,zoom=0.82)
plt.tight_layout(); plt.savefig(os.path.join(out,'fig_index3d.png'),dpi=220); plt.close()
json.dump(data,open(os.path.join(root,'code','data','dstar_P2_2.json'),'w'))
# Figure 2: G_mu(n)/n^3
fig,axs=plt.subplots(1,2,figsize=(10,3.8))
for N in (4,6,8):
    b,c,G0,G1=forms(N,101)
    ns=range(1,101)
    axs[0].plot(ns,[float(G0[n])/n**3 for n in ns],'.-',ms=3,lw=0.7,color=cols[N],label=r'$N=%d$'%N)
    ns1=range(2,101)
    axs[1].plot(ns1,[float(G1[n])/n**3 for n in ns1],'.-',ms=3,lw=0.7,color=cols[N],label=r'$N=%d$'%N)
for a,t in zip(axs,(r'$G_0(n)/n^3$',r'$G_1(n)/n^3$')):
    a.set_xlabel(r'$n$'); a.set_title(t); a.axhline(0,color='k',lw=0.5); a.legend(frameon=False)
plt.tight_layout(); plt.savefig(os.path.join(out,'fig_Gmu.png'),dpi=220); plt.close()
# Figure 3: margin of mu=1 over the class-number term
fig,ax=plt.subplots(figsize=(6,3.8))
from Gmu import mul
for N in (4,6,8):
    b,c,G0,G1=forms(N,121)
    K=121
    a=[Fr(0)]*K; a[0]=Fr(1)
    for k in range(1,K):
        for _ in range(6):
            for m in range(k,K): a[m]+=a[m-k]
    th=[Fr(0)]*K
    for u in range(-11,12):
        if u*u<K: th[u*u]+=1
    a=mul(th,a,K)
    xs=[];ys=[]
    for n in range(2,K):
        Dl=4*n-N
        if Dl<=0 or n%N: continue
        D=Dl//N
        E1=sum(a[j]*G1[n-N*j] for j in range(0,n//N+1))
        xs.append(Dl); ys.append(math.log10(float(E1/(2*GH[N]*H(D)))))
    ax.plot(xs,ys,'o-',ms=3,lw=0.8,color=cols[N],label=r'$N=%d$'%N)
ax.set_xlabel(r'$\Delta=4m-N\ell^2$'); ax.set_ylabel(r'$\log_{10}\,E_1(\Delta)/(2\hat g\,H(\Delta/N))$')
ax.axhline(0,color='k',lw=0.5); ax.legend(frameon=False)
plt.tight_layout(); plt.savefig(os.path.join(out,'fig_margin.png'),dpi=220); plt.close()
print("figures written")
