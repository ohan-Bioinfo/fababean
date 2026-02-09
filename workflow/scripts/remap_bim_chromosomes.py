#!/usr/bin/env python3
"""Remap PLINK BIM chromosome names to numeric codes."""

from __future__ import annotations

import argparse
from pathlib import Path


def load_mapping(path: Path) -> dict[str, str]:
    mapping: dict[str, str] = {}
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            mapping[parts[0]] = parts[1]
    if not mapping:
        raise ValueError(f"No mapping entries found in {path}")
    return mapping


def remap_bim(input_bim: Path, output_bim: Path, mapping: dict[str, str]) -> int:
    converted = 0
    with input_bim.open("r", encoding="utf-8") as src, output_bim.open(
        "w", encoding="utf-8"
    ) as dst:
        for raw in src:
            fields = raw.rstrip("\n").split()
            if not fields:
                continue
            chrom = fields[0]
            if chrom in mapping:
                fields[0] = mapping[chrom]
                converted += 1
            dst.write("\t".join(fields) + "\n")
    return converted


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mapping", required=True, type=Path)
    parser.add_argument("--input-bim", required=True, type=Path)
    parser.add_argument("--output-bim", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    mapping = load_mapping(args.mapping)
    converted = remap_bim(args.input_bim, args.output_bim, mapping)
    print(f"Remapped {converted} BIM rows: {args.input_bim} -> {args.output_bim}")


if __name__ == "__main__":
    main()
