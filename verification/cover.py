"""Exact cover identities and open conditions from public rational coefficients.

SymPy is used only for exact rational polynomial operations and modular
resultants. The declared kappa is input; no Belyi-family history is loaded.
"""
import time
import sympy as S
x,T=S.symbols("x T")
def scalar_poly(seq):
    return S.Poly.from_list([S.Rational(c) for c in seq[::-1]],x,domain=S.QQ)

def fiber_poly(rows):
    assert all(len(row)==3 for row in rows)
    return S.Poly(sum(S.Rational(c)*x**i*T**j for i,row in enumerate(rows) for j,c in enumerate(row)),x,T,domain=S.QQ)

def t_remainder(f,k):
    """Reduce with T^3=3kT+2k, using coefficient triples of Q[x]."""
    terms=[S.Poly(0,x,domain=S.QQ) for _ in range(max(3,f.degree(T)+1))]
    for (i,j),c in f.terms():terms[j]+=S.Poly(c*x**i,x,domain=S.QQ)
    for j in range(len(terms)-1,2,-1):
        terms[j-2]+=terms[j].mul_ground(3*k)
        terms[j-3]+=terms[j].mul_ground(2*k)
        terms[j]=S.Poly(0,x,domain=S.QQ)
    return terms[:3]

def cubic_norm(f,k):
    cs=t_remainder(f,k)
    cols=[t_remainder(f*S.Poly(T**j,x,T,domain=S.QQ),k) for j in range(3)]
    return (cols[0][0]*(cols[1][1]*cols[2][2]-cols[1][2]*cols[2][1])
            -cols[1][0]*(cols[0][1]*cols[2][2]-cols[0][2]*cols[2][1])
            +cols[2][0]*(cols[0][1]*cols[1][2]-cols[0][2]*cols[1][1]))

def reduction(f,p):
    terms={m:int(c.p)*pow(int(c.q),-1,p)%p for m,c in f.terms()}
    return S.Poly.from_dict(terms,f.gens,modulus=p)
def check(core):
    started=time.monotonic();c=core["cover"]
    d={"P":c["P0"],"Q":c["Q0"],"F":c["F"],"R":c["R0"],"D":c["D"],"U":c["U0_by_X_then_T"],"V":c["V0_by_X_then_T"]}
    P,Q,F,R,D=[scalar_poly(d[n]) for n in ['P','Q','F','R','D']]
    U,V=fiber_poly(d['U']),fiber_poly(d['V'])
    k=S.Rational(c["kappa"])
    assert k not in [0,1]
    B=S.Poly(T**3-3*k*T-2*k,T,domain=S.QQ)
    assert [f.degree() for f in [P,Q,F,R,D]]==[8,6,4,14,2]
    assert D==S.Poly(x*x+3,x,domain=S.QQ) and F.LC()==1
    assert U.degree(x)==6 and U.coeff_monomial(x**6)==1
    assert all(U.coeff_monomial(x**6*T**j)==0 for j in (1,2))
    assert V.degree(x)==2
    norm=P*P-Q*Q*F-D*R
    assert norm.is_zero,'Norm identity failed'
    fiber=S.Poly(R.as_expr()-2*P.as_expr()*T+D.as_expr()*T*T,x,T,domain=S.QQ)-U*U*V
    remainders=t_remainder(fiber,k)
    assert all(f.is_zero for f in remainders),'Symmetric fibers failed'
    normU=cubic_norm(U,k)
    AA=D*P.diff()-D.diff()*P
    CC=2*F*(D*Q.diff()-D.diff()*Q)+D*Q*F.diff()
    critical=CC*CC-4*AA*AA*F+144*F.LC()*R.LC()*D*D*normU
    assert critical.is_zero,'Critical divisor identity failed'
    # Scalar open conditions over Q.
    assert all(S.gcd(f,f.diff()).degree()==0 for f in [F,Q,D])
    assert all(S.gcd(f,g).degree()==0 for f,g in [(F,Q),(F,D),(Q,D)])
    # A single good reduction proves the following resultants have no common
    # zero with the branch cubic, simultaneously for every characteristic-zero
    # embedding of its etale algebra. All relevant leading degrees are retained.
    certificate=None
    for p in [41]:
        try:
            up,vp=reduction(U,p),reduction(V,p)
            bp=reduction(B,p)
            scalars={n:reduction(f,p) for n,f in [('Q',Q),('F',F),('D',D)]}
        except ValueError:continue
        if S.gcd(bp,bp.diff()).degree()>0:continue
        if up.degree(x)!=6 or vp.degree(x)!=2:continue
        # Nonvanishing V leading coefficient in each cubic embedding.
        vlead=S.Poly(sum(c*T**j for (i,j),c in vp.terms() if i==2),T,modulus=p)
        if S.gcd(vlead,bp).degree()>0:continue
        relations=[('U_Uprime',up.as_expr(),S.diff(up.as_expr(),x)),
                   ('V_Vprime',vp.as_expr(),S.diff(vp.as_expr(),x)),
                   ('U_V',up.as_expr(),vp.as_expr())]
        relations += [('U_'+n,up.as_expr(),f.as_expr()) for n,f in scalars.items()]
        residues={};passed=True
        for name,f,g in relations:
            z=S.Poly(S.resultant(f,g,x),T,modulus=p).rem(bp)
            if S.gcd(z,bp).degree()>0:passed=False;break
            residues[name]=[int(z.nth(i))%p for i in range(3)]
        if passed:
            certificate={'prime':p,'B_ascending':[int(bp.nth(i))%p for i in range(4)],
                         'nonzero_norm_resultant_residues_mod_B':residues,
                         'V_leading_coefficient_ascending':[int(vlead.nth(i))%p for i in range(3)]}
            break
    assert certificate is not None,'No good-reduction open-condition certificate found'
    discR=R.discriminant()
    rn=S.integer_nthroot(abs(discR.p),2);rd=S.integer_nthroot(discR.q,2)
    assert discR>0 and rn[1] and rd[1]
    return {"status":"PASS_EXACT_COVER_IDENTITIES_AND_OPEN_CONDITIONS","norm_identity":True,"symmetric_branch_fibers":True,"critical_divisor_identity":True,"squarefree_disjoint_F_Q_D":True,"open_conditions":certificate,"R_discriminant_positive_square":True,"seconds":time.monotonic()-started}
