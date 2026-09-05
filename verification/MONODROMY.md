# Exact cover and monodromy verification

Run from this directory, with Python 3.9 or later and SymPy installed:

    python3 monodromy.py --data .. --loop all

The only mathematical inputs are `math_core.json` and
`monodromy_tubes.json.gz` in the data directory. The input directory can be
anywhere; no historical archive, receipt, executable, or network access is used.
The result is printed as JSON. The supplied `monodromy_example_run.json` is
an example of an actual fresh execution, not an input trusted by the checker.

`--loop b1`, `--loop b2`, and `--loop b3` are available for individual loop
diagnostics. They certify only the selected loop. There is deliberately no
receipt-merging mode: use `--loop all` to obtain the complete geometric group
verification. On the release machine that command took 18.98 seconds. Each
invocation has a 35-second alarm; slower machines may require an explicitly
reviewed resource adjustment before an all-loop run. No arithmetic conclusion
is inferred from a timeout.

## Mathematical checks

`cover.py` uses exact SymPy polynomials over Q and modulo 41 to verify:

1. Degrees, monicity and D=X^2+3.
2. P0^2-Q0^2 F=D R0.
3. R0-2 T P0+T^2 D=U0^2 V0 in Q[X,T]/(T^3-3 kappa T-2 kappa).
4. The critical-divisor identity using the cubic norm of U0.
5. Squarefreeness and pairwise disjointness of F,Q0,D over Q.
6. Separability of the branch cubic, preserved leading degrees, nonzero
   V0 leading coefficient, and all six requisite resultant gcd tests modulo 41.
7. The discriminant of R0 is a positive rational square.

Kappa is the explicit rational in the public coefficient data. No upstream
family parametrization is required to check these identities for this fixed
cover. The separate paper explains the degree, genus and Riemann--Hurwitz
argument turning these verified identities into an exhaustive branch inventory.

`monodromy_intervals.py` uses integer rectangles with common denominator 10^100.
Every multiplication rounds lower endpoints down and upper endpoints up.
The upper complex norm bound is l1; the lower bound is the maximum of the
individual real and imaginary distance-to-zero bounds. These bounds enclose
the Euclidean modulus and do not rely on floating-point arithmetic.

For every moving centre and every segment, Taylor coefficients are enclosed
uniformly on [-1/2,1/2]. The checker proves the strict Rouché inequality and
pairwise separation of the fourteen moving discs. At a shared endpoint the
centre is identical and the radii are nested, so the uniquely enclosed roots
agree. Strict containment identifies the initial and final base discs.
No proposed permutation from the input is read.

`monodromy_geometry.py` verifies three disjoint branch isolating intervals,
the prescribed positively oriented squares and upper-half-plane stems, and
the monotone subdivision of every polygon edge. It closes the group generated
by the three resulting permutations, obtaining 1092 elements. It independently
enumerates all determinant-one projective matrix actions over F13 and builds
a common labelling carrying the three generators to those actions.

The output therefore certifies the geometric PSL2(F13) group of this fixed
cover together with its required exact cover preconditions. It does not on
its own establish the arithmetic specialization group or the SL2 central lift;
those are separate arguments and verification components in the paper.

## Porting and validation

The 22 mathematical functions handling boxes, tubes, path geometry, group
enumeration, and cover polynomial reductions were copied without changing
their abstract syntax from the previously used exact checker. The public
driver changes only input loading, result serialization and resource bounds.
The cover driver uses the explicitly supplied kappa, fixes the successful
prime 41 rather than searching primes, and asserts the discriminant-square
condition. Unrelated source-family and source-j exploration is omitted.

Fresh validation checked all 506 segments (183,179,144) and all cover conditions.
Perturbing the constant coefficient of P0 was rejected, as was forcing two
moving tube centres to coincide. A separate reviewer compared the public
implementation with the mathematical preconditions and checked that the
reported scope does not include the arithmetic specialization or central lift.
