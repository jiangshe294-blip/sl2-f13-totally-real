"""Small standalone exact g/U and finite-field Spin checks; standard library only.

No archive assertions or saved PASS statuses are imported. The published input
coefficients, factors and local matrices are treated as untrusted proposals.
"""
from fractions import Fraction as Q
from math import gcd,lcm,isqrt
def trim(a):
    while len(a)>1 and a[-1]==0:a.pop()
    return a
def remainder(a,b,p=None):
    a=trim(a[:]);b=trim(b[:]);assert b!=[0]
    while len(a)>=len(b) and a!=[0]:
        c=a[-1]*pow(b[-1],-1,p)%p if p else Q(a[-1])/b[-1];s=len(a)-len(b)
        for i,v in enumerate(b):a[s+i]=(a[s+i]-c*v)%p if p else a[s+i]-c*v
        trim(a)
    return a
def times(a,b,p):
    out=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]=(out[i+j]+x*y)%p
    return trim(out)
def pgcd(a,b,p):
    while b!=[0]:a,b=b,remainder(a,b,p)
    return [(x*pow(a[-1],-1,p))%p for x in a]
def power(a,n,m,p):
    out=[1]
    while n:
        if n&1:out=remainder(times(out,a,p),m,p)
        n//=2
        if n:a=remainder(times(a,a,p),m,p)
    return out
def subtract_x(a,p):
    out=a[:]+[0]*max(0,2-len(a));out[1]=(out[1]-1)%p;return trim(out)
def primitive(a):
    d=1
    for q in a:d=lcm(d,Q(q).denominator)
    z=[int(Q(q)*d) for q in a];c=gcd(*z)
    return [x//c for x in z] if c else [0]
def signs_at(seq,minus):
    return [(-1 if f[-1]<0 else 1)*(-1 if minus and (len(f)-1)%2 else 1) for f in seq]
def changes(s):return sum(a!=b for a,b in zip(s,s[1:]))
def check_g(core,local):
    g=list(map(Q,core['g_Q_monic']));Z=list(map(int,core['g_Z']));assert len(g)==15 and g[-1]==1 and all(Q(z,Z[-1])==x for z,x in zip(Z,g)) and gcd(*Z)==1
    c=core['cover'];t=Q(c['t_star']);P=list(map(Q,c['P0']));R=list(map(Q,c['R0']));D=list(map(Q,c['D']));at=lambda a,i:a[i] if i<len(a) else Q(0)
    G=[at(R,i)-2*t*at(P,i)+t*t*at(D,i) for i in range(15)];assert [x/G[-1] for x in G]==g
    seq=[Z,primitive([i*Z[i] for i in range(1,15)])]
    while len(seq[-1])>1:
        r=remainder(seq[-2],seq[-1]);assert r!=[0];seq.append(primitive([-x for x in r]))
    minus=changes(signs_at(seq,True));plus=changes(signs_at(seq,False));assert minus-plus==14
    patterns=[]
    for cert in local['modular_certificates']:
        p=cert['prime'];assert p in [11,17];assert all(p%d for d in range(2,isqrt(p)+1));red=[z%p for z in Z];assert red[-1]!=0
        monic=[x*pow(red[-1],-1,p)%p for x in red];factors=cert['monic_factors_ascending'];product=[1]
        assert pgcd(monic,[(i*monic[i])%p for i in range(1,15)],p)==[1]
        for f in factors:
            n=len(f)-1;assert f[-1]==1 and n in [1,7,13]
            if n>1:
                xp=power([0,1],p,f,p);assert pgcd(f,subtract_x(xp,p),p)==[1]
                x=[0,1]
                for _ in range(n):x=power(x,p,f,p)
                assert x==[0,1]
            product=times(product,f,p)
        assert product==monic;patterns.append(sorted(len(f)-1 for f in factors))
    assert sorted(patterns)==[[1,13],[7,7]]
    # Reducible factors would have degrees which are subset sums of both
    # patterns. No proper positive degree can be in both sets.
    possible={0}
    for n in patterns[0]:possible|={s+n for s in list(possible)}
    other={0}
    for n in patterns[1]:other|={s+n for s in list(other)}
    assert possible&other=={0,14}
    traces=[Q(14)]
    for k in range(1,27):
        if k<=14:v=-k*g[14-k]-sum(g[14-i]*traces[k-i] for i in range(1,k))
        else:v=-sum(g[14-i]*traces[k-i] for i in range(1,15))
        traces.append(v)
    U=[list(map(Q,row)) for row in core['trace_isometry_U']];assert len(U)==14 and all(len(row)==14 for row in U)
    TU=[[sum(traces[i+k]*U[k][j] for k in range(14)) for j in range(14)] for i in range(14)]
    assert all(sum(U[k][i]*TU[k][j] for k in range(14))==int(i==j) for i in range(14) for j in range(14))
    return {'status':'PASS_EXACT_COVER_SPECIALIZATION_REAL_IRREDUCIBLE_G_AND_TRACE_U','degree':14,'real_roots':minus-plus,'sturm_variations':[minus,plus],'factor_patterns':patterns,'rational_isometry_coordinates':196,'PSL_containment_is_separate_written_theorem':True}
def action(m):
    a,b,c,d=m;out=[]
    for x in [None]+list(range(13)):
        num,den=(a,c) if x is None else ((a*x+b)%13,(c*x+d)%13)
        out.append(0 if den==0 else (num*pow(den,-1,13)%13)+1)
    return tuple(out)
def group():
    ans=set();lifts=0
    for a in range(13):
        for b in range(13):
            for c in range(13):
                ds=[(1+b*c)*pow(a,-1,13)%13] if a else range(13) if (-b*c)%13==1 else []
                for d in ds:lifts+=1;ans.add(action((a,b,c,d)))
    assert lifts==2184 and len(ans)==1092;return ans
def compose(a,b):return tuple(a[b[i]] for i in range(14))
def bareiss(A):
    a=[r[:] for r in A];n=len(a);sign=1;previous=1
    for k in range(n-1):
        if a[k][k]==0:
            pivot=next((i for i in range(k+1,n) if a[i][k]),None)
            if pivot is None:return 0
            a[k],a[pivot]=a[pivot],a[k];sign=-sign
        p=a[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                q,r=divmod(p*a[i][j]-a[i][k]*a[k][j],previous);assert r==0;a[i][j]=q
            a[i][k]=0
        previous=p
    return sign*a[-1][-1]
def spin_setup(core,local):
    s=local['spin'];f=local['frame'];p=f['prime'];assert p==21767 and all(p%d for d in range(2,isqrt(p)+1))
    perms=list(map(tuple,s['group_permutations']));assert len(perms)==1092 and set(perms)==group()
    roots=f['projective_root_images'];assert roots==core['local']['root_images'] and len(set(roots))==14
    def modq(x):q=Q(x);return q.numerator*pow(q.denominator,-1,p)%p
    gm=list(map(modq,core['g_Q_monic']));assert all(sum(c*pow(r,i,p) for i,c in enumerate(gm))%p==0 for r in roots)
    U=[[modq(q) for q in row] for row in core['trace_isometry_U']];assert U==f['U_mod_p']
    A=[[sum(pow(roots[i],k,p)*U[k][j] for k in range(14))%p for j in range(14)] for i in range(14)];assert A==f['A_mod_p']
    assert all(sum(A[k][i]*A[k][j] for k in range(14))%p==int(i==j) for i in range(14) for j in range(14));assert f['D_signs']==[-1]+[1]*13
    return p,perms,A,f['D_signs']
def spin_slice(core,local,k):
    p,perms,A,signs=spin_setup(core,local);start=156*k;end=start+156;values=[]
    for g in perms[start:end]:values.append(bareiss([[(int(i==j)+A[g[i]][j]*signs[j])%p for j in range(14)] for i in range(14)])%p)
    assert values==local['spin']['deltas_mod_p'][start:end]
    return {'status':'PASS_RECOMPUTED_INTEGER_BAREISS_SPIN_SLICE','start':start,'end':end,'deltas_mod_p':values,'embedding_semantics_requires_characteristic_zero_root_theorem':True}
def spin_conclusion(core,local,slices):
    p,perms,A,signs=spin_setup(core,local);deltas=[]
    for k,s in enumerate(slices):assert s['start']==156*k and s['end']==156*(k+1);deltas+=s['deltas_mod_p']
    idx={g:i for i,g in enumerate(perms)};H={tuple([0]+[(s*x+t)%13+1 for x in range(13)]) for s in [1,3,9] for t in range(13)};assert len(H)==39 and H<=set(perms)
    assert all(compose(h,k) in H for h in H for k in H);covered=set();betas=[]
    for record in local['spin']['right_cosets']:
        g=perms[record['representative_index']];members={idx[compose(g,h)] for h in H};assert len(members)==39 and members==set(record['member_indices']) and not(covered&members);covered|=members
        v=1
        for i in members:v=v*deltas[i]%p
        assert v==record['beta_mod_p'];betas.append(v)
    assert len(betas)==28 and covered==set(range(1092)) and betas==local['spin']['beta_values_mod_p'] and len(set(betas))==28
    coeff=[1]
    for b in betas:coeff=times(coeff,[-b%p,0,1],p)
    assert coeff==local['spin']['P56_coefficients_ascending'];zero=[i for i,b in enumerate(betas) if b==0];assert zero==[8]
    assert all(pow(b,(p-1)//2,p)==p-1 for b in betas if b)
    derivative=[i*coeff[i]%p for i in range(1,57)];assert pgcd(coeff,trim(derivative),p)==[0,1]
    return {'status':'PASS_RECOMPUTED_1092_SPIN_DETERMINANTS_AND_28_NORMS','p':p,'distinct_beta_residues':28,'one_zero_beta_index':8,'nonzero_beta_nonsquares':27,'mod56_derivative_gcd':'X','mod56_is_not_squarefree':True,'PSL_field_and_Spin_theorem_are_separate':True}
