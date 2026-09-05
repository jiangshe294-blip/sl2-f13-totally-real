# Certified totally real `SL2(F13)` realization

This repository contains the exact computational data and source-only
verification programs accompanying Jiang Yu's preprint

> *A certified degree-56 totally real realization of SL2(F13)*.

The paper constructs an irreducible polynomial `f in Q[x]` of degree 56 whose
roots are all real and whose splitting field has Galois group isomorphic to
`SL2(F13)`.  The complete primitive integral coefficient list is stored in
[`f56_Z.txt.gz`](f56_Z.txt.gz).  The paper PDF and its arXiv-style LaTeX source
archive are deposited with this same data release on Zenodo:

- Version DOI: <https://doi.org/10.5281/zenodo.22317659>
- All-versions DOI: <https://doi.org/10.5281/zenodo.22317658>

## Read the exact polynomial online

The complete degree-56 coefficient list can be inspected without downloading
or opening a 31.6 MB line of digits:

- Public HTML reader: <https://sl2-f13-polynomial.jiangshe294.chatgpt.site>
- Editable LaTeX companion: <https://sl2-f13-polynomial.jiangshe294.chatgpt.site/f56_reader_source.tex>
- Reproducible website source: [`readable-web/`](readable-web/)

The reader loads one coefficient at a time, displays numbered 100-digit blocks,
and exposes the per-coefficient digest and source-line metadata.  Its generated
data reconstruct the canonical `f56_Z.txt` payload byte for byte.

The construction reported in the paper relies on substantive contributions
from GPT-5.6 Sol, as stated in the abstract and Zenodo description.

## Repository contents

- [`MATHEMATICAL_DATA.md`](MATHEMATICAL_DATA.md) specifies the polynomial and
  every exact algebraic data format in mathematical terms.
- [`PUBLIC_MANIFEST.json`](PUBLIC_MANIFEST.json) records the SHA-256 digest and
  byte size of every distribution file.
- [`verify_distribution.py`](verify_distribution.py) checks byte integrity,
  gzip payloads, and the 57-line polynomial encoding.
- [`verification/`](verification/) contains the exact C++/Python arithmetic
  verifier, the monodromy verifier, documentation, and fresh-run receipts.
- The compressed JSON and text files contain the full explicit coefficient,
  field-tower, root-frame, Spin-determinant, and monodromy data.

No historical research logs, machine-specific paths, compiled executables, or
opaque cached powers are included.

## Quick integrity check

Python 3.11 or later is sufficient for the distribution check:

```sh
python3 verify_distribution.py
```

The exact arithmetic verification additionally requires a C++17 compiler and
GMP/GMPXX.  From the repository root:

```sh
python3 verification/run.py build --work ./my-verification-work
python3 verification/run.py g --work ./my-verification-work
python3 verification/run.py spin0 --work ./my-verification-work
```

Continue with the phase order documented in
[`verification/README.md`](verification/README.md).  The monodromy certificate
has its own instructions in
[`verification/MONODROMY.md`](verification/MONODROMY.md).

## Proof and verification boundary

The source-only programs certify the exact arithmetic and geometric layers
described in their documentation.  They do not replace the mathematical
arguments in the paper: in particular, the characteristic-zero field bridge,
the central extension from `PSL2(F13)` to `SL2(F13)`, the stabilizer argument,
irreducibility, and total reality must be read together with the paper.

The generic Sage field-bridge verifier is supplied as source, but the full
generic Sage run was not executed during publication preparation.  This is
explicitly disclosed in [`verification/README.md`](verification/README.md).

## Citation and license

Citation metadata is provided in [`CITATION.cff`](CITATION.cff).  The deposited
materials are released under the Creative Commons Attribution 4.0 International
license; see [`LICENSE.md`](LICENSE.md).
