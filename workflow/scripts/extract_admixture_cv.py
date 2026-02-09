#!/usr/bin/env python3
"""Extract ADMIXTURE CV errors from log_K*.out files."""

from __future__ import annotations

import argparse
import re
from pathlib import Path


CV_RE = re.compile(r"CV error \(K=(\d+)\):\s*([0-9]*\.?[0-9]+)")
LOG_K_RE = re.compile(r"log_K(\d+)\.out$")


def parse_log(log_path: Path) -> tuple[int, float | None, str]:
    match = LOG_K_RE.search(log_path.name)
    if not match:
        raise ValueError(f"Unexpected log filename: {log_path.name}")
    k = int(match.group(1))
    text = log_path.read_text(encoding="utf-8", errors="ignore")
    cv_match = CV_RE.search(text)
    if cv_match:
        return k, float(cv_match.group(2)), "ok"
    if "Error opening .bed file" in text:
        return k, None, "bed_open_error"
    return k, None, "cv_not_found"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--detailed-output", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    logs = sorted(args.input_dir.glob("log_K*.out"), key=lambda p: int(LOG_K_RE.search(p.name).group(1)))  # type: ignore[arg-type]
    if not logs:
        raise FileNotFoundError(f"No log_K*.out files found in {args.input_dir}")

    rows = [parse_log(log) for log in logs]
    ok_rows = [(k, cv) for k, cv, status in rows if status == "ok" and cv is not None]

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as out:
        out.write("K CV_Error\n")
        for k, cv in sorted(ok_rows):
            out.write(f"{k} {cv:.5f}\n")

    args.detailed_output.parent.mkdir(parents=True, exist_ok=True)
    with args.detailed_output.open("w", encoding="utf-8") as out:
        out.write("K\tCV_Error\tStatus\n")
        for k, cv, status in sorted(rows):
            cv_str = "" if cv is None else f"{cv:.5f}"
            out.write(f"{k}\t{cv_str}\t{status}\n")

    if ok_rows:
        best_k, best_cv = min(ok_rows, key=lambda item: item[1])
        print(f"Best K: {best_k} (CV error: {best_cv:.5f})")
    else:
        print("No valid CV values were found.")


if __name__ == "__main__":
    main()
