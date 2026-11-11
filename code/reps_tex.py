import json, re
from fractions import Fraction as Fr
reps=json.load(open('reps.json'))
def tok(t):
    t=t.strip()
    m=re.fullmatch(r't3\((\d+)\)\^2',t)
    if m: return ('t',int(m.group(1)))
    m=re.fullmatch(r't2\(4\*(\d+)\)\^2/4',t)
    if m: return ('s',int(m.group(1)))
    m=re.fullmatch(r'Ee\((\d+)\)',t)
    if m: return ('e',int(m.group(1)))
    m=re.fullmatch(r'E4\((\d+)\)',t)
    if m: return ('E',int(m.group(1)))
    raise ValueError(t)
def mono(name):
    parts=[tok(x) for x in re.split(r'\*(?=[tE])|\|',name)]
    cnt={}
    for p in parts: cnt[p]=cnt.get(p,0)+1
    order={'t':0,'s':1,'e':2,'E':3}
    out=''
    for (k,a),c in sorted(cnt.items(),key=lambda x:(order[x[0][0]],x[0][1])):
        base={'t':'t_{%d}'%a,'s':'s_{%d}'%a,'e':'e_{%d}'%a,'E':('E_4(%d\\tau)'%a if a>1 else 'E_4(\\tau)')}[k]
        out+=base+('^{%d}'%c if c>1 else '')
    return out
def frac(v):
    v=Fr(v)
    return '%d'%v.numerator if v.denominator==1 else '\\tfrac{%d}{%d}'%(v.numerator,v.denominator)
lines=[]
for N in (4,6,8):
    for mu in (0,1):
        terms=reps['%d_G%d'%(N,mu)]
        s=' + '.join(('' if Fr(v)==1 else frac(v)+'\\,')+mono(k) for k,v in terms)
        lines.append((N,mu,s))
for N,mu,s in lines:
    print('%%N=%d mu=%d'%(N,mu)); print(s)
