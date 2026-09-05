# Concrete mathematical data

These files give exact coefficients, not a list of references to undisclosed
objects. No network access, historical state engine, or machine-specific path
is needed to read them. Fractions are decimal strings `n` or `n/d` with positive
denominators. Unless otherwise specified polynomial coefficients are ascending.
The paper and the separate verification programs provide the mathematical
arguments; file hashes alone are not proofs of irreducibility or a Galois group.

## The actual polynomial

`f56_Z.txt.gz` decompresses to exactly 57 ASCII lines. Line i, with zero-based
indexing, is the complete signed decimal integer z_i. Put

    Z(X) = sum(i=0,...,56) z_i X^i,
    D = z_56,
    f(X) = Z(X)/D,
    P28(T) = sum(i=0,...,28) (z_(2i)/D) T^i.

Thus f(X)=P28(X^2), f is monic, and D is positive. All 28 odd coefficients
are zero, all 29 even coefficients are nonzero. Z has content one; that is
an arithmetic assertion checked by the independent encoding verifier, not
something inferred from gzip or its checksum. D has 1,088,453 decimal digits.
Keeping Z alone avoids storing a second 30 MB decimal expansion of the same
rational coefficients.

    python3 read_coefficients.py
    python3 read_coefficients.py --form Q --coefficient 56
    python3 read_coefficients.py --form P28 --all

The reader uses only the standard library and never converts the huge decimal
strings to Python integers. Its Q output `n/D` need not be in lowest terms;
it is exactly the desired rational number. The reader checks file hashes and
encoding, but does not pretend that these tests establish the paper's theorem.

## Finite-algebra construction inputs

`math_core.json` contains the following conventional objects:

| Key | Mathematical meaning |
| --- | --- |
| g_Z | primitive degree-14 g_0, all 15 integer coefficients |
| g_Q_monic | g=g_0/lc(g_0); E=Q[a]/(g(a)) |
| r | r(a) in E; M=E[w]/(w^2-r(a)) |
| beta_x, beta_y | beta=x(a)+w y(a), 14 rational coefficients each |
| cover.P0, Q0, F, R0, D | actual cover polynomials; D=X^2+3 |
| cover.kappa | B(T)=T^3-3 kappa T-2 kappa |
| cover.t_star | the actual rational specialization, not an unspecified parameter |
| cover.U0_by_X_then_T, V0_by_X_then_T | coefficient of X^i is a degree<=2 polynomial in T, modulo B(T) |
| trace_isometry_U | all 196 rational entries of U |
| local | prime 21767 and projectively labelled root residues |
| quadratic_layer_formula | exact two-torsion formula data for the chosen r |

The cover map on y^2=F(X) is (P0(X)+Q0(X)y)/D(X). It satisfies
P0^2-Q0^2 F=D R0. Define G(X,T)=R0(X)-2 T P0(X)+T^2 D(X).
The supplied scalar in `cover.G_to_g_Z_scalar` gives
g_0(X)=scalar*G(X,t_star) exactly.

`tower.json.gz` gives both sextics q6plus and q6minus in V[Z],
V=E[b]/(q13(a,b)), with

    q13(a,Y) = sum(i=1,...,14) g[i] sum(j=0,...,i-1) a^(i-1-j) Y^j.

Each V coefficient is encoded by a positive common `denominator` and a
14-by-13 `numerator_rows_x_then_y` array: entry (i,j) means a^i b^j.
Then K=V[z]/(q6plus(a,b,z)), a 1092-dimensional Q-algebra in the proof.

`root_frame.json.gz` gives all 14 roots r_infinity,r_0,...,r_12 of g in K.
Here r_infinity=a, r_0=b, r_1=z; the remaining roots have explicit K coordinates.
A K vector is six V blocks, ascending in the power of z. The exceptional
`monomial_exponents_z_a_b` array means precisely z^k a^i b^j.
This data determines the matrix W with row labelled t equal to
(1,r_t,...,r_t^13). With D0=diag(-1,1,...,1), define

    delta = det(I14 + W U D0).

`delta_and_embedding.json.gz` supplies the exact 1092 rational coordinates
of delta, the exact K vector s_star, the 14 coefficients of t(a), and the
cubic relative norm of delta. The explicit embedding of M into K is
w=s_star/t(a). Flat delta entry 182*k+13*i+j is the coefficient of z^k a^i b^j.

`norm13_input.json.gz` supplies A,B in V such that the cubic norm of delta is
A(a,b)+w B(a,b). This provides a compact fresh-computation input to the
13-step relative norm verifier. The computed result is the beta in math_core.
`local_spin.json.gz` provides modular factorization certificates for g,
the finite projective action, modular spin determinants, cosets and norm data.
These are mathematical values only; portable verifiers recompute their claims.

`monodromy_tubes.json.gz` adds the exact rational complex centres and polygon
frames for all 506 certified continuation segments (183,179,144 on the three
loops), the common base centres, branch isolating intervals, and loop geometry.
Each stored decimal means an exact rational, not a floating-point assertion.
The integer-interval verifier uses scale 10^100 and base radius 10^80 in scaled
coordinates. This file is 409,491 bytes and is part of the public manifest.

The final construction can therefore be stated without any JSON-pointer
definition: beta=Norm_(K/M)(delta), and f(X)=Norm_(M/Q)(X^2-beta).
All objects in these formulas are specified above by ordinary exact algebraic
coefficient data. The file formats are simply serializations of those objects.

## Size and reproducibility boundary

The compressed 57-line polynomial is 14,800,139 bytes. The complete explicit
root frame is an additional 10,788,760 bytes. All eight numerical-data files
total 38,123,949 bytes before outer compression of math_core.json; the source
code and documentation add a small amount. A smaller core delivery may omit
the explicit root-frame file only if it clearly states that omission and
supplies its constructive regeneration procedure; it is not then the full
explicit-coordinate collection. Do not equate this compact collection with
the separately preserved approximately 1 GB historical scientific archive.

`PUBLIC_MANIFEST.json` binds each distributed numerical file, source and
document, including uncompressed numerical bytes where applicable. Historical
exporters, author-review reports, compiled binaries, TeX files and research
logs are not part of this ancillary archive. No public coefficient file carries
a historical pending flag. The separate publisher archive retains the old
records unchanged.

Run `python3 verify_distribution.py` after unpacking. It independently reads
every public file, decompresses the supplied numerical gzip files, checks
stored and unpacked SHA-256 values, and checks the 57-line polynomial encoding.
Its scope is explicitly integrity, not Galois-theoretic proof. The source-only
mathematical commands are described in `verification/README.md` and
`verification/MONODROMY.md`.

The paper's separately supplied TeX source includes an explicit table of all
15 coefficients of the primitive degree-fourteen input polynomial. Line
wrapping inside each decimal integer is purely typographic; every digit is
present. The ancillary computation archive itself contains no TeX files.
