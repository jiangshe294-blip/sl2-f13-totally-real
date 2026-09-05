# Public source-only arithmetic verification

This directory uses only its C++/Python sources and the mathematical data files
in its parent directory. It does **not** read the historical research archive,
run an archived executable, import a producer's saved powers, or require an old
binary hash. The containing directory may have any name: there is no hard-coded
`anc/`, project path, or original computer path.

## Requirements

- Python 3.11 or later, with assertions enabled (never `-O` or `PYTHONOPTIMIZE`).
- A C++17 compiler and the GMP C and C++ development libraries (`gmp`, `gmpxx`).
- A POSIX system. The tested platform is macOS arm64; Linux uses the same source
  interfaces, but a separate Linux test is not claimed by the included receipt.
- No Sage, SymPy, NumPy, internet connection, or Python package installation is
  needed for this particular runner. Separate reference verifiers may have
  additional requirements.

The recommended build interface is `run.py build`. It uses `pkg-config gmpxx`
when available and otherwise the compiler's library search path. If necessary,
pass `--gmp-prefix PREFIX` to supply `PREFIX/include` and `PREFIX/lib`. Optional
`CMakeLists.txt` gives the same two native targets for readers who prefer CMake;
the Python runner's build receipt is still produced by `run.py build`.

## Execute in phases

From the directory containing `math_core.json`:

```sh
python3 verification/run.py build --work ./my-verification-work
python3 verification/run.py g --work ./my-verification-work
python3 verification/run.py spin0 --work ./my-verification-work
```

Continue with `spin1` through `spin6`, then `spin-conclusion`. Next execute
`norm01` through `norm13` **in ascending order**. Next execute `h27` through
`h00` **in descending order**. Finish with `encoding` and `conclusion`.
`python3 verification/run.py list` displays the phase families. The phases are
separate invocations so no native arithmetic call can silently exceed its
35-second resource cap. Each Python phase has a 42-second wall-clock cap.

For convenience, a POSIX shell can run the already enumerated phases:

```sh
for i in 0 1 2 3 4 5 6; do
  python3 verification/run.py spin$i --work ./my-verification-work || exit 1
done
python3 verification/run.py spin-conclusion --work ./my-verification-work
for i in 01 02 03 04 05 06 07 08 09 10 11 12 13; do
  python3 verification/run.py norm$i --work ./my-verification-work || exit 1
done
for i in 27 26 25 24 23 22 21 20 19 18 17 16 15 14 13 12 11 10 09 08 07 06 05 04 03 02 01 00; do
  python3 verification/run.py h$i --work ./my-verification-work || exit 1
done
python3 verification/run.py encoding --work ./my-verification-work
python3 verification/run.py conclusion --work ./my-verification-work
```

The loops do not change the per-phase caps. A nonzero exit, timeout, missing
checkpoint, failed exact comparison, or altered dependency is not a PASS.
Work is append-only: already completed phase files are never overwritten.
After an interrupted sequence, execute the next missing phase; for an entirely
new run choose a new work directory. The runtime residuals are much larger than
the compressed input distribution; allow approximately 1 GB of working space.
They are generated locally rather than shipped as opaque verification inputs.

## What is recomputed

`g` reconstructs the degree-fourteen specialization from the published rational
cover coefficients and target. It recomputes the exact rational Sturm sequence,
obtaining fourteen real roots. At 11 and 17 it recomposes the supplied proposed
factors and independently applies the prime-degree Rabin criterion. The factor
degree patterns are 7+7 and 1+13; their proper subset-degree sets are disjoint,
proving rational irreducibility. It also recomputes all trace moments needed for
the full exact rational identity U-transpose T U = I.

The seven `spin` slices independently enumerate SL2(F13) and its 1,092 projective
permutations, check the published root residues against g, recompute WU and its
orthogonality, and calculate all 1,092 determinants by exact integer Bareiss
elimination (156 per slice). `spin-conclusion` checks the subgroup H of order 39,
all 28 cosets, their norms and the resulting degree-56 polynomial modulo 21767.
It correctly records the **one zero beta residue**: the modular polynomial is
not squarefree, its derivative gcd is X, and the 27 nonzero beta residues are
distinct nonsquares. No squarefreeness claim is made at this prime.

The `norm` phases start only from the published u=A+B*w in

    F = (Q[a]/g(a))[b,w] / (q13(a,b), w^2-r(a)),
    q13(a,Y) = g(Y)/(Y-a).

The code computes u^k in an integer pseudo-remainder representation. It derives
the q13 power traces directly from g, calculates Tr(u^k), and applies Newton's
identities to the characteristic coefficients. Stage 13 compares the newly
computed norm against the published beta in all 28 rational coordinates.
No historical Newton-stage coefficients or expected intermediate powers are
loaded. Only the final beta is treated as the proposed output to be checked.

The `h` phases execute a new exact Horner calculation P28(beta), starting from
1, where P28(T) is obtained from the **actual complete** published coefficient
file by P28[T^i] = Z[X^(2i)] / Z[X^56]. Each new 28-coordinate residual is saved
locally. The last phase requires all 28 coordinates to be exactly zero. A
separate integer gcd check verifies that the 57-line Z file has content one;
positive leading coefficient and zero odd coefficients are also checked. Thus
the publication file is a concrete polynomial, not an instruction to search
for one or an unexpanded norm placeholder.

The reduction algorithms in `exact_V.hpp` and `norm28.cpp` are source-only
copies of separately audited exact arithmetic, with the historical main entry
and archive coupling removed where appropriate. `norm13.cpp` uses that
independent arithmetic to **recompute**, rather than compare supplied, powers
and Newton coefficients. These public wrappers are new packaging code, not a
claim of another independent mathematical author.

## Explicit proof boundary

The final runner PASS is a precisely named arithmetic-spine result, not a
standalone assertion that it has reproved every inverse-Galois theorem.
Identification of F with the intended field, the canonical M embedding, and
u with Norm(K/F)(delta) require the separate exact field-bridge proof/checks.
The geometric/arithmetic PSL cover theorem, characteristic-zero root labels,
and the Spin central-extension and stabilizer arguments remain part of the
paper and its explicitly separated verification layers. In particular,
P28(beta)=0 identifies P28 with Norm(M/Q)(T-beta) once the proved degree-28
stabilizer statement is supplied. A local residue frame alone is never used as
a proof that all fourteen characteristic-zero coordinates are actual roots.

The complete generic field-bridge reference supplied alongside this directory
is intended to expose a direct exact alternative to the historical integral
overorder/p-adic separation certificates. Its expensive phases must be reported
separately from the phases actually executed in a fresh directory here.

Its precise command, from the parent data directory, is:

```sh
sage -python verification/reference_bridge.py --data . --phase roots
```

The available reference phases are `tower`, `roots`, `delta`, `cubic`,
`embedding`, `descent`, and `all`. Each invocation constructs the exact field
tower with field-presentation checks enabled. For the complete direct bridge,
use `--phase all`. This generic implementation may require substantial time and
memory and is not represented as meeting the optimized runner's 35-second
arithmetic slices. **Only static source review and Python syntax compilation
were performed for this generic reference in the publication preparation;
the full generic Sage computation was not run.** This is distinct from the
actually executed source-only arithmetic-spine and monodromy checks.

The separate `monodromy.py` verifier, supplied in this same directory, requires
SymPy and `monodromy_tubes.json.gz`. From the parent data directory, run:

```sh
python3 verification/monodromy.py --loop all
```

It checks the complete fixed-cover norm/fiber/critical identities and open
conditions, the positive-square discriminant of R, all 506 exact continuation
segments, and the geometric PSL2(F13) permutation-group identification. Its
fresh-run receipt and separate-author source review are independent of the
53-phase arithmetic-spine summary. Individual `--loop b1`, `b2`, `b3` runs
check only the named loop; they do not by themselves trigger the aggregate PSL
identification. The geometric result still needs the paper's arithmetic
normalizer/discriminant argument and its characteristic-zero field bridge.

## Reproducibility receipts

Each new checkpoint binds the input bytes, executable built from the current
source bytes, and predecessor checkpoints. A freshly compiled binary may have
a different hash on another machine: that new hash is recorded rather than
compared against an old machine's executable. No normalization or assertion
suppression is necessary. The public clean-directory test summary records the
exact tested source/input hashes and every completed phase. A same-machine
empty-directory test is not described as an independent second-machine test.
