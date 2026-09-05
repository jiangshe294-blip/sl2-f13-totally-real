# Public exact-coefficient reader

This directory is the reproducible source for the public HTML reader:

<https://sl2-f13-polynomial.jiangshe294.chatgpt.site>

The interface renders its mathematical notation from LaTeX using KaTeX.  It
loads one exact coefficient at a time and paginates each million-digit integer
into numbered 100-digit blocks.  An editable standalone LaTeX companion is in
[`public/f56_reader_source.tex`](public/f56_reader_source.tex).

The generated files under `public/data/` are intentionally not committed.
They are deterministically reconstructed from the repository's canonical
[`../f56_Z.txt.gz`](../f56_Z.txt.gz), whose uncompressed SHA-256 is

```text
5058792bf79dd594034393954aac53bf801a57f7a36e989c55efefc9c270fd50
```

## Reproduce locally

Install Node.js 22.13 or later and Python 3.11 or later, then run:

```sh
cd readable-web
npm ci
npm run prepare:data
npm run dev
```

The preparation script refuses to proceed unless the decompressed canonical
file has the expected digest and exactly 57 valid decimal-integer lines.  It
also reconstructs the original byte stream from the generated files and checks
exact equality before writing the index.

For a production build:

```sh
npm run lint
npm run build
```

The live deployment is public and the accompanying mathematical materials are
released under CC BY 4.0.
