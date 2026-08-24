#!/usr/bin/env python3
"""Compare two pipeline OUT_ROOT trees for output parity.

Walks every CSV under the baseline root, finds its counterpart under the
other root, and verifies identical shape, columns, row order, exact object
values, and numeric values within tolerance. Path-valued cells are compared
after normalizing each run's own root prefix, and timing/log artifacts are
excluded, so two runs of the same inputs must compare clean regardless of
worker counts or cache state.

Usage:
    python tools/compare_run_outputs.py BASELINE_ROOT CANDIDATE_ROOT \
        [--rtol 1e-9] [--atol 1e-12]

Exit code 0 = parity; 1 = any mismatch or missing counterpart.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

EXCLUDED_NAMES = {"stage_timings.csv"}


def _normalize_paths(series: pd.Series, root: Path) -> pd.Series:
    """Replace this run's root prefix in string cells so both runs align."""
    text = series.astype(str)
    return text.str.replace(str(root.resolve()), "<OUT_ROOT>", regex=False).str.replace(
        str(root), "<OUT_ROOT>", regex=False
    )


def compare_csv(baseline_path: Path, candidate_path: Path, baseline_root: Path, candidate_root: Path,
                rtol: float, atol: float) -> list[str]:
    """Return a list of mismatch descriptions for one CSV pair (empty = parity)."""
    problems: list[str] = []
    try:
        baseline = pd.read_csv(baseline_path, low_memory=False)
    except pd.errors.EmptyDataError:
        baseline = pd.DataFrame()
    try:
        candidate = pd.read_csv(candidate_path, low_memory=False)
    except pd.errors.EmptyDataError:
        candidate = pd.DataFrame()
    if list(baseline.columns) != list(candidate.columns):
        problems.append(f"column mismatch: {list(baseline.columns)} vs {list(candidate.columns)}")
        return problems
    if baseline.shape != candidate.shape:
        problems.append(f"shape mismatch: {baseline.shape} vs {candidate.shape}")
        return problems
    for col in baseline.columns:
        b = baseline[col]
        c = candidate[col]
        b_num = pd.to_numeric(b, errors="coerce")
        c_num = pd.to_numeric(c, errors="coerce")
        both_numeric = b_num.notna().equals(c_num.notna()) and (
            b_num.notna() | b.isna()
        ).all() and (c_num.notna() | c.isna()).all()
        if both_numeric:
            b_vals = b_num.to_numpy(dtype=float)
            c_vals = c_num.to_numpy(dtype=float)
            close = np.isclose(b_vals, c_vals, rtol=rtol, atol=atol, equal_nan=True)
            if not close.all():
                worst = int(np.argmax(~close))
                problems.append(
                    f"numeric column '{col}' differs at row {worst}: {b_vals[worst]!r} vs {c_vals[worst]!r}"
                )
        else:
            b_txt = _normalize_paths(b, baseline_root)
            c_txt = _normalize_paths(c, candidate_root)
            equal = (b_txt == c_txt) | (b.isna() & c.isna())
            if not equal.all():
                worst = int(np.argmax(~equal.to_numpy()))
                problems.append(
                    f"text column '{col}' differs at row {worst}: {b_txt.iloc[worst]!r} vs {c_txt.iloc[worst]!r}"
                )
    return problems


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("baseline_root", type=Path)
    parser.add_argument("candidate_root", type=Path)
    parser.add_argument("--rtol", type=float, default=1e-9)
    parser.add_argument("--atol", type=float, default=1e-12)
    args = parser.parse_args()

    baseline_root = args.baseline_root
    candidate_root = args.candidate_root
    if not baseline_root.is_dir() or not candidate_root.is_dir():
        print(f"Both roots must exist: {baseline_root} {candidate_root}", file=sys.stderr)
        return 2

    n_compared = 0
    n_failed = 0
    for baseline_path in sorted(baseline_root.rglob("*.csv")):
        if baseline_path.name in EXCLUDED_NAMES:
            continue
        rel = baseline_path.relative_to(baseline_root)
        candidate_path = candidate_root / rel
        if not candidate_path.exists():
            print(f"[MISSING] {rel}")
            n_failed += 1
            continue
        problems = compare_csv(
            baseline_path, candidate_path, baseline_root, candidate_root, args.rtol, args.atol
        )
        n_compared += 1
        if problems:
            n_failed += 1
            print(f"[DIFF   ] {rel}")
            for problem in problems[:5]:
                print(f"          - {problem}")
        else:
            print(f"[OK     ] {rel}")
    print(f"\ncompared={n_compared} failed={n_failed}")
    return 1 if n_failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
