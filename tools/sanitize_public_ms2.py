#!/usr/bin/env python3
"""Strip raw observed-MS2-ion data from shared report sidecars for delivery.

Public reports must not contain observed MS2 ions. The shared per-report sidecars
(`*_top_candidates.csv`, `*_candidate_score_summary.csv`) ship in every report's
`data/` dir and carry observed ion/loss bins + matched-peak m/z. This removes any
column whose name indicates raw observed MS2 ion data, in place. The derived MS2
fit score (used for confidence, internal-only display) is not an ion list and is
kept in `top_candidates.csv` pipeline score columns; only raw-ion columns go.

Usage:
    python3 tools/sanitize_public_ms2.py <reports_root>
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

csv.field_size_limit(10 ** 7)

# Substrings identifying raw observed-MS2-ion columns (not derived scores).
RAW_ION_PATTERNS = (
    "observed_ion_bins",
    "observed_loss_bins",
    "observed_ion_intensities",
    "observed_loss_intensities",
    "ms2_pairs",
    "observed_ms2",
    "spectrum_matched_peaks",
)
TARGET_SUFFIXES = ("_top_candidates.csv", "_candidate_score_summary.csv")


def _raw_cols(header: list[str]) -> list[str]:
    return [c for c in header if any(p in c.lower() for p in RAW_ION_PATTERNS)]


def sanitize_csv(path: Path) -> int:
    with path.open(newline="") as h:
        reader = csv.reader(h)
        try:
            header = next(reader)
        except StopIteration:
            return 0
        drop = set(_raw_cols(header))
        if not drop:
            return 0
        keep_idx = [i for i, c in enumerate(header) if c not in drop]
        rows = [[row[i] for i in keep_idx if i < len(row)] for row in reader]
    with path.open("w", newline="") as h:
        writer = csv.writer(h)
        writer.writerow([header[i] for i in keep_idx])
        writer.writerows(rows)
    return len(drop)


def main() -> None:
    root = Path(sys.argv[1] if len(sys.argv) > 1 else "_audit_pull/reports")
    files = 0
    cols = 0
    for d in sorted(p for p in root.iterdir() if p.is_dir()):
        data = d / "data"
        if not data.is_dir():
            continue
        for csv_path in data.glob("*.csv"):
            if not any(csv_path.name.endswith(s) for s in TARGET_SUFFIXES):
                continue
            dropped = sanitize_csv(csv_path)
            if dropped:
                files += 1
                cols += dropped
    print(f"[sanitize] stripped raw-ion columns from {files} sidecar file(s) ({cols} column-instances)")


if __name__ == "__main__":
    main()
