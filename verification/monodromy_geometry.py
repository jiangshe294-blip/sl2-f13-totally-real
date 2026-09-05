"""Exact branch geometry and finite group enumeration.

Functions retained verbatim from the original rational-cover checker.
"""
import itertools
from fractions import Fraction as F
import sympy as S
from monodromy_intervals import point,on_segment
def expected_vertices(name,manifest):
    base=point(manifest['basepoint'])
    b=next(v for v in manifest['branches'] if v['name']==name)
    c,r=F(b['center']),F(b['radius'])
    stem=[base,(c,base[1]),(c,r)]
    raw=stem+[(c-r,r),(c-r,-r),(c+r,-r),(c+r,r),(c,r)]+stem[-2::-1]
    answer=[]
    for v in raw:
        if not answer or v!=answer[-1]:answer.append(v)
    return answer

def check_manifest(manifest,k):
    t=S.symbols('t');B=S.Poly(t**3-3*k*t-2*k,t,domain=S.QQ)
    assert B.degree()==3 and S.gcd(B,B.diff()).degree()==0
    branches=manifest['branches'];assert [v['name'] for v in branches]==['b1','b2','b3']
    base=point(manifest['basepoint']);assert base[1]>0
    spans=[];reports=[]
    for b in branches:
        a,z,c,r=[F(b[key]) for key in ['left','right','center','radius']]
        assert r>0 and base[1]>r and c-r<a<z<c+r
        assert B.eval(S.Rational(a))!=0 and B.eval(S.Rational(z))!=0
        assert B.count_roots(S.Rational(a),S.Rational(z))==1
        spans.append((c-r,c+r))
        reports.append({'name':b['name'],'left':str(a),'right':str(z),'exact_root_count':1})
    assert all(spans[i][1]<spans[i+1][0] for i in range(2))
    return reports

def verify_path(name,frames,manifest):
    vertices=expected_vertices(name,manifest)
    assert len(frames)>=2
    assert point(frames[0]['t'])==vertices[0] and point(frames[-1]['t'])==vertices[-1]
    edge=0;previous=vertices[0]
    for frame in frames[1:]:
        q=point(frame['t'])
        while edge<len(vertices)-2 and previous==vertices[edge+1]:edge+=1
        a,b=vertices[edge:edge+2]
        assert on_segment(q,a,b),('polygon mismatch',name,edge,q)
        axis=0 if a[0]!=b[0] else 1
        assert (q[axis]-previous[axis])*(b[axis]-a[axis])>0
        previous=q
    assert edge==len(vertices)-2

def cycle_type(p):
    unseen=set(range(len(p)));sizes=[]
    while unseen:
        a=next(iter(unseen));v=a;n=0
        while v in unseen:unseen.remove(v);n+=1;v=p[v]
        assert v==a;sizes.append(n)
    return sorted(sizes)

def identify_psl(permutations):
    """Construct and check new labels, with a hard finite enumeration cap."""
    assert len(permutations)==3
    assert all(sorted(p)==list(range(14)) for p in permutations)
    assert all(cycle_type(p)==[1,1,2,2,2,2,2,2] for p in permutations)
    def mul(a,b):return tuple(a[b[i]] for i in range(14))
    identity=tuple(range(14));known={identity};queue=[identity]
    for g in queue:
        for h in permutations:
            q=mul(g,h)
            if q not in known:
                known.add(q);queue.append(q)
                assert len(known)<=1092,'Generated group exceeds PSL order.'
    assert len(known)==1092
    unipotent=next(g for g in queue if cycle_type(g)==[1,13])
    fixed=next(i for i in range(14) if unipotent[i]==i)
    phi=[None]*14;phi[fixed]=13
    v=next(i for i in range(14) if i!=fixed)
    for j in range(13):phi[v]=j;v=unipotent[v]
    assert sorted(phi)==list(range(14))
    standard={}
    for a,b,c,d in itertools.product(range(13),repeat=4):
        if (a*d-b*c)%13!=1:continue
        p=[]
        for z in range(14):
            num,den=(a,c) if z==13 else ((a*z+b)%13,(c*z+d)%13)
            p.append(13 if den==0 else num*pow(den,-1,13)%13)
        standard.setdefault(tuple(p),(a,b,c,d))
    assert len(standard)==1092
    matrices=[]
    for p in permutations:
        conjugate=[None]*14
        for i in range(14):conjugate[phi[i]]=phi[p[i]]
        assert tuple(conjugate) in standard,'New generators are not conjugate to the standard PSL action.'
        matrices.append(list(standard[tuple(conjugate)]))
    return {'order':1092,'root_label_to_P1_F13':phi,'generator_SL2_matrices':matrices,
            'independently_enumerated_standard_projective_matrix_actions':1092}
