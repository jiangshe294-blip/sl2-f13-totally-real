"""Exact outward-rounded integer complex boxes and moving Rouche tubes.

Mathematical functions retained verbatim from the original exact checker.
All inputs are rational strings; no floating-point arithmetic is used.
"""
import math
from fractions import Fraction
S=10**100
class Box:
    __slots__=('rl','rh','il','ih')
    def __init__(self,rl,rh,il=0,ih=0):self.rl,self.rh,self.il,self.ih=rl,rh,il,ih
    @staticmethod
    def rational(q):
        f=Fraction(str(q));n,d=f.numerator*S,f.denominator
        return Box(n//d,-((-n)//d))
    @staticmethod
    def pair(pair):
        r,i=map(Box.rational,pair)
        return Box(r.rl,r.rh,i.rl,i.rh)
    def __add__(a,b):
        if not isinstance(b,Box):b=Box.rational(b)
        return Box(a.rl+b.rl,a.rh+b.rh,a.il+b.il,a.ih+b.ih)
    __radd__=__add__
    def __neg__(a):return Box(-a.rh,-a.rl,-a.ih,-a.il)
    def __sub__(a,b):return a+-b
    def __rsub__(a,b):return -a+b
    @staticmethod
    def rm(a,b,c,d):
        terms=(a*c,a*d,b*c,b*d)
        return min(terms)//S,-((-max(terms))//S)
    def __mul__(a,b):
        if not isinstance(b,Box):
            if isinstance(b,int):
                if b>=0:return Box(a.rl*b,a.rh*b,a.il*b,a.ih*b)
                return -a*(-b)
            b=Box.rational(b)
        r0,r1=Box.rm(a.rl,a.rh,b.rl,b.rh)
        r2,r3=Box.rm(a.il,a.ih,b.il,b.ih)
        i0,i1=Box.rm(a.rl,a.rh,b.il,b.ih)
        i2,i3=Box.rm(a.il,a.ih,b.rl,b.rh)
        return Box(r0-r3,r1-r2,i0+i2,i1+i3)
    __rmul__=__mul__
    def upper(a):return max(abs(a.rl),abs(a.rh))+max(abs(a.il),abs(a.ih))
    def lower(a):
        def lb(lo,hi):return min(abs(lo),abs(hi)) if lo*hi>0 else 0
        return max(lb(a.rl,a.rh),lb(a.il,a.ih))
    def widen(a,r):return Box(a.rl-r,a.rh+r,a.il-r,a.ih+r)
ZERO=Box(0,0);ONE=Box(S,S)
def add(*terms):
    out=[ZERO]*max(map(len,terms))
    for a in terms:
        for i,v in enumerate(a):out[i]=out[i]+v
    return out

def scale(a,c):return [v*c for v in a]
def midpoint(a,b):return (a+b)*Box.rational('1/2')

def ceildiv(a,b):return -((-a)//b)

def taylor(coeff,c):
    out=[coeff[-1]]
    for v in coeff[-2::-1]:
        out=[out[0]*c+v]+[out[j]*c+out[j-1] for j in range(1,len(out))]+[out[-1]]
    return out

def g_at(data,D,t):return add(data['R'],scale(data['P'],-2*t),scale(D,t*t))

def rouche_at(coeff,c,rho):
    tay=taylor(coeff,c);lhs=tay[0].upper();power=rho
    for k in range(2,len(tay)):
        power=ceildiv(power*rho,S)
        lhs+=ceildiv(tay[k].upper()*power,S)
    rhs=tay[1].lower()*rho//S
    assert lhs<rhs,('base Rouche',lhs,rhs)
    return lhs,rhs

def tube_bounds(data,D,t0,t1,c0,c1):
    tt=midpoint(t0,t1);dt=t1-t0;cc=midpoint(c0,c1);dc=c1-c0
    aa=taylor(g_at(data,D,tt),cc)
    bb=taylor(add(scale(data['P'],-2*dt),scale(D,2*tt*dt)),cc)
    dd=taylor(scale(D,dt*dt),cc)
    powers=[ONE]
    for j in range(14):powers.append(powers[-1]*dc)
    upper=[];lower1=0
    for k in range(15):
        row=[]
        for j in range(15-k):
            value=aa[k+j]*powers[j]*math.comb(k+j,k)
            if j>=1 and k+j-1<len(bb):value=value+bb[k+j-1]*powers[j-1]*math.comb(k+j-1,k)
            if j>=2 and k+j-2<len(dd):value=value+dd[k+j-2]*powers[j-2]*math.comb(k+j-2,k)
            row.append(value)
        upper.append(sum(ceildiv(v.upper(),2**j) for j,v in enumerate(row)))
        if k==1:lower1=row[0].lower()-sum(ceildiv(v.upper(),2**j) for j,v in enumerate(row) if j)
    return upper,lower1

def one_segment(data,D,start,end):
    t0,t1=Box.pair(start['t']),Box.pair(end['t'])
    c0,c1=[list(map(Box.pair,frame['centers'])) for frame in [start,end]]
    assert len(c0)==len(c1)==14
    mid=[midpoint(x,y) for x,y in zip(c0,c1)];delta=[y-x for x,y in zip(c0,c1)]
    distances=[[0]*14 for _ in range(14)]
    for i in range(14):
        for j in range(i):
            value=(mid[i]-mid[j]).lower()-ceildiv((delta[i]-delta[j]).upper(),2)
            assert value>0,('tube centers not separated',i,j)
            distances[i][j]=distances[j][i]=value
    radii=[min(distances[i][j] for j in range(14) if j!=i)//16 for i in range(14)]
    ratios=[]
    for i,rho in enumerate(radii):
        assert rho>0
        upper,lower=tube_bounds(data,D,t0,t1,c0[i],c1[i])
        assert lower>0,('derivative lower bound',i)
        lhs=upper[0];power=rho
        for k in range(2,15):
            power=ceildiv(power*rho,S)
            lhs+=ceildiv(upper[k]*power,S)
        rhs=lower*rho//S
        assert lhs<rhs,('uniform Rouche',i,Fraction(lhs,max(1,rhs)))
        ratios.append(Fraction(lhs,rhs))
    assert all(radii[i]+radii[j]<distances[i][j] for i in range(14) for j in range(i))
    return radii,max(ratios)

def point(pair):return tuple(map(Fraction,pair))

def on_segment(q,a,b):
    return (b[0]-a[0])*(q[1]-a[1])==(b[1]-a[1])*(q[0]-a[0]) and all(min(a[k],b[k])<=q[k]<=max(a[k],b[k]) for k in range(2))
