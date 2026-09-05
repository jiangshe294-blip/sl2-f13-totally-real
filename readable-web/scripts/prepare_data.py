#!/usr/bin/env python3
"""Generate browser-friendly coefficient files from the canonical gzip payload."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
from pathlib import Path


EXPECTED_SOURCE_SHA256 = (
    "5058792bf79dd594034393954aac53bf801a57f7a36e989c55efefc9c270fd50"
)
INTEGER_PATTERN = re.compile(r"-?(?:0|[1-9][0-9]*)\Z")


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def main() -> None:
    repo_root = Path(__file__).resolve().parents[2]
    default_input = repo_root / "f56_Z.txt.gz"
    default_output = Path(__file__).resolve().parents[1] / "public" / "data"

    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=default_input)
    parser.add_argument("--output", type=Path, default=default_output)
    args = parser.parse_args()

    with gzip.open(args.input, "rb") as source:
        canonical = source.read()

    source_hash = sha256(canonical)
    if source_hash != EXPECTED_SOURCE_SHA256:
        raise SystemExit(
            "canonical coefficient digest mismatch: "
            f"expected {EXPECTED_SOURCE_SHA256}, got {source_hash}"
        )
    if not canonical.endswith(b"\n"):
        raise SystemExit("canonical coefficient file must end with a newline")

    lines = canonical[:-1].decode("ascii").split("\n")
    if len(lines) != 57:
        raise SystemExit(f"expected 57 coefficient lines, got {len(lines)}")

    coefficient_dir = args.output / "coefficients"
    coefficient_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []

    for index, value in enumerate(lines):
        if INTEGER_PATTERN.fullmatch(value) is None:
            raise SystemExit(f"line {index + 1} is not a canonical decimal integer")
        is_zero = value == "0"
        sign = "0" if is_zero else ("-" if value.startswith("-") else "+")
        unsigned = value[1:] if value.startswith("-") else value
        encoded = value.encode("ascii")
        (coefficient_dir / f"z{index}.txt").write_bytes(encoded + b"\n")
        rows.append(
            {
                "index": index,
                "line": index + 1,
                "sign": sign,
                "digits": len(unsigned),
                "prefix": value[:24],
                "suffix": value[-24:],
                "sha256": sha256(encoded),
                "is_zero": is_zero,
            }
        )

    rebuilt = b"".join(
        (coefficient_dir / f"z{index}.txt").read_bytes() for index in range(57)
    )
    if rebuilt != canonical:
        raise SystemExit("generated coefficient files do not reconstruct the source")

    index_document = {
        "schema": "f56-Z-readable-index-v1",
        "source": "f56_Z.txt",
        "source_bytes": len(canonical),
        "source_sha256": source_hash,
        "coefficient_order": "ascending; line i+1 is z_i",
        "rows": rows,
    }
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "f56_Z_index.json").write_text(
        json.dumps(index_document, indent=2) + "\n", encoding="utf-8"
    )
    print(f"generated 57 exact coefficients; source SHA-256 {source_hash}")


if __name__ == "__main__":
    main()
