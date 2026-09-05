#!/usr/bin/env python3
"""Direct characteristic-zero reference checks for the full splitting field.

Run with SageMath: sage -python reference_bridge.py --data DIR --phase roots
This intentionally simple reference uses generic exact number-field arithmetic.
It is NOT the optimized p-adic/norm-separation verifier used in the research
calculation. Its large-field phases may be very expensive. Shipping this source
does not mean that those phases have been timed or completed on a new platform.
The public reproduction report distinguishes actually executed checks from this
reference implementation. No historical receipt, binary, path or network access
is needed. Inputs are ordinary exact rational coefficient tables.
"""
import argparse
import gzip
import json
import sys
import time
from pathlib import Path


def load(path):
    raw = path.read_bytes()
    if path.suffix == '.gz':
        raw = gzip.decompress(raw)
    return json.loads(raw)


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--data', type=Path, default=Path(__file__).resolve().parent.parent)
    ap.add_argument('--phase', choices=['tower', 'roots', 'delta', 'cubic',
                                       'embedding', 'descent', 'all'], required=True)
    args = ap.parse_args()
    # Imported only after argument parsing so --help needs no Sage installation.
    from sage.all import QQ, PolynomialRing, NumberField, matrix, identity_matrix, diagonal_matrix
    started = time.monotonic()
    core = load(args.data / 'math_core.json')
    tower = load(args.data / 'tower.json.gz')
    root_data = load(args.data / 'root_frame.json.gz')['roots']
    bridge = load(args.data / 'delta_and_embedding.json.gz')
    norm13 = load(args.data / 'norm13_input.json.gz')
    gq = list(map(QQ, core['g_Q_monic']))
    require(len(gq) == 15 and gq[-1] == 1, 'degree-fourteen monic input')

    def ev(coefficients, point):
        out = point.parent()(0)
        for coefficient in reversed(coefficients):
            out = out * point + coefficient
        return out

    QX = PolynomialRing(QQ, 'X')
    g = QX(gq)
    # check=True deliberately asks Sage to verify the field presentations.
    # The paper also proves their irreducibility from the projective action.
    E = NumberField(g, 'a', check=True)
    a = E.gen()
    EY = PolynomialRing(E, 'Y')
    q13c = [sum((gq[i] * a**(i - 1 - j) for i in range(j + 1, 15)), E(0))
            for j in range(14)]
    q13 = EY(q13c)
    require(q13.degree() == 13 and q13.is_monic(), 'q13 degree')
    V = E.extension(q13, 'b', check=True)
    b = V.gen()
    av = V(a)

    def in_V(block, aa, bb):
        rows = block['numerator_rows_x_then_y']
        require(len(rows) == 14 and all(len(row) == 13 for row in rows),
                '14 by 13 coefficient block')
        denominator = QQ(block['denominator'])
        require(denominator > 0, 'positive V denominator')
        # Nested Horner evaluation fixes the convention a^i*b^j.
        by_a = [ev(list(map(QQ, row)), bb) for row in rows]
        return ev(by_a, aa) / denominator

    VZ = PolynomialRing(V, 'Z')
    qplus = VZ([in_V(v, av, b) for v in tower['q6_plus_ascending']])
    qminus = VZ([in_V(v, av, b) for v in tower['q6_minus_ascending']])
    Z = VZ.gen()
    require(qplus.degree() == qminus.degree() == 6, 'two sextic degrees')
    require(qplus.is_monic() and qminus.is_monic(), 'monic sextics')
    require(qplus * qminus * (Z - av) * (Z - b) == VZ(gq),
            'exact two-sextic factorization')
    K = V.extension(qplus, 'z', check=True)
    ak, bk, zk = K(a), K(b), K.gen()
    require(K.absolute_degree() == 1092, 'splitting tower absolute degree')

    def in_K(blocks, aa=ak, bb=bk, zz=zk):
        require(len(blocks) == 6, 'six z coefficient blocks')
        return ev([in_V(v, aa, bb) for v in blocks], zz)

    def root_at(record, aa=ak, bb=bk, zz=zk):
        if 'monomial_exponents_z_a_b' in record:
            k, i, j = record['monomial_exponents_z_a_b']
            return zz**k * aa**i * bb**j
        return in_K(record['coefficients_by_z_degree_ascending'], aa, bb, zz)

    labels = ['infinity'] + [str(i) for i in range(13)]
    require(set(root_data) == set(labels), 'all fourteen labels supplied')
    roots = [root_at(root_data[label]) for label in labels]
    require(roots[:3] == [ak, bk, zk], 'normalized initial triple')

    def delta_at(aa=ak, bb=bk, zz=zk):
        pairs = bridge['delta_numerator_denominator_pairs']
        require(len(pairs) == 1092, '1092 delta coordinates')
        blocks = []
        for k in range(6):
            by_a = []
            for i in range(14):
                coefficients = []
                for j in range(13):
                    numerator, denominator = pairs[182*k + 13*i + j]
                    require(QQ(denominator) > 0, 'positive delta denominator')
                    coefficients.append(QQ(numerator) / QQ(denominator))
                by_a.append(ev(coefficients, bb))
            blocks.append(ev(by_a, aa))
        return ev(blocks, zz)

    def permutation(function):
        return tuple(0 if function(t) is None else function(t) + 1
                     for t in [None] + list(range(13)))

    translation = permutation(lambda t: None if t is None else (t + 1) % 13)
    inversion = permutation(lambda t: 0 if t is None else
                            (None if t == 0 else (-pow(t, -1, 13)) % 13))
    scaling = permutation(lambda t: None if t is None else 3*t % 13)

    def check_action(p):
        aa, bb, zz = [roots[p[i]] for i in range(3)]
        require(ev(gq, aa) == 0, 'image of first relation')
        q13_at = [sum((gq[i]*aa**(i-1-j) for i in range(j+1, 15)), K(0))
                  for j in range(14)]
        require(ev(q13_at, bb) == 0, 'image of q13 relation')
        q6_at = [in_V(c, aa, bb) for c in tower['q6_plus_ascending']]
        require(ev(q6_at, zz) == 0, 'image of selected sextic relation')
        # These relation tests define a Q-embedding K -> K, hence an
        # automorphism. Test its action on EVERY supplied root, not just a,b,z.
        for index, label in enumerate(labels):
            require(root_at(root_data[label], aa, bb, zz) == roots[p[index]],
                    'actual root action at label ' + label)

    def check_roots():
        for i, root in enumerate(roots):
            require(ev(gq, root) == 0, 'root relation ' + str(i))
            for earlier in roots[:i]:
                require(root != earlier, 'distinct roots')
        for p in (translation, inversion, scaling):
            check_action(p)
        identity = tuple(range(14))
        group, queue = {identity}, [identity]
        for p in queue:
            for h in (translation, inversion):
                q = tuple(p[h[i]] for i in range(14))
                if q not in group:
                    group.add(q)
                    queue.append(q)
                    require(len(group) <= 1092, 'projective group order bound')
        require(len(group) == 1092, 'PSL2 projective generators')

    def check_delta():
        W = matrix(K, [[root**j for j in range(14)] for root in roots])
        U = matrix(QQ, [list(map(QQ, row)) for row in core['trace_isometry_U']])
        I = identity_matrix(K, 14)
        D0 = diagonal_matrix(K, [-1] + [1]*13)
        A = W*U*D0
        require(A.transpose()*A == I, 'full trace isometry at actual roots')
        require(A.det() == 1, 'orientation of actual orthogonal matrix')
        delta = delta_at()
        require(delta != 0, 'nonzero delta')
        require((I + A).det() == delta, 'direct exact Spin determinant')

    def check_cubic():
        delta = delta_at()
        # tau fixes infinity and 0 and sends label 1 to label 3.
        n3 = delta * delta_at(ak, bk, roots[4]) * delta_at(ak, bk, roots[10])
        require(n3 == in_K(bridge['norm3_delta_by_z_degree']),
                'delta*tau(delta)*tau^2(delta)')

    def embedding():
        C = [1, 3, 9]
        star = K(0)
        for t in range(13):
            theta = sum((roots[1+(t+c) % 13] - roots[1+(t-c) % 13]
                         for c in C), K(0))
            star += roots[t+1]**2 * theta
        require(star == in_K(bridge['s_star_by_z_degree']), 'exact s-star expression')
        t_of_a = ev(list(map(QQ, bridge['t_coefficients_ascending'])), ak)
        require(star != 0 and t_of_a != 0, 'nonzero quadratic embedding')
        w = star/t_of_a
        require(w*w == ev(list(map(QQ, core['r'])), ak), 'w^2=r(a)')
        return w

    def check_descent():
        w = embedding()
        n3 = in_K(bridge['norm3_delta_by_z_degree'])
        value = (in_V(norm13['A_in_V'], ak, bk)
                 + w*in_V(norm13['B_in_canonical_w_basis'], ak, bk))
        require(n3 == value, 'same cubic norm as degree-thirteen input')

    phases = {'tower': lambda: None, 'roots': check_roots, 'delta': check_delta,
              'cubic': check_cubic, 'embedding': embedding, 'descent': check_descent}
    selected = list(phases) if args.phase == 'all' else [args.phase]
    for name in selected:
        phase_start = time.monotonic()
        phases[name]()
        print(json.dumps({'phase': name, 'status': 'PASS_EXACT_EQUALITIES',
                          'seconds': time.monotonic()-phase_start}), flush=True)
    unexecuted = [name for name in phases if name not in selected and name != 'tower']
    print(json.dumps({'status': 'PASS_SELECTED_REFERENCE_PHASES', 'phases': selected,
                      'unexecuted_reference_phases': unexecuted,
                      'complete_reference_bridge_checked_in_this_invocation': args.phase == 'all',
                      'seconds': time.monotonic()-started,
                      'remaining': 'Unexecuted reference phases are not implied. Cover/monodromy, degree13 norm, Horner and group-theoretic proofs are also separate.'}))


if __name__ == '__main__':
    if not __debug__:
        raise SystemExit('Run without Python optimization.')
    if hasattr(sys, 'set_int_max_str_digits'):
        sys.set_int_max_str_digits(0)
    main()
