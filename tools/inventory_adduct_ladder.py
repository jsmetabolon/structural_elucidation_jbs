#!/usr/bin/env python3
"""Compare CHRO_LIB hypothesis labels to the Excel adduct ladder (read-only).

Does not rerun inference. Prints how often stored masses would move if
``--derive-neutral-mass-from-adduct`` were turned on, and the analog regen
command if those pools are still used.

Usage:
    python3 tools/inventory_adduct_ladder.py
    python3 tools/inventory_adduct_ladder.py --limit-per-file 200
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from structural_elucidation.chemistry.adduct_ladder import (  # noqa: E402
    derive_target_neutral_mass,
    row_for_observed_label,
)

csv.field_size_limit(10**7)
MOVE_DA = 0.005


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--feature-rows-dir",
        type=Path,
        default=ROOT / "data" / "feature_rows",
    )
    parser.add_argument("--limit-per-file", type=int, default=500)
    args = parser.parse_args()

    paths = sorted(args.feature_rows_dir.glob("*_neutral_mass_hypotheses.csv"))
    if not paths:
        print(f"No hypothesis CSVs under {args.feature_rows_dir}")
        return

    mapped = 0
    derived = 0
    moved = 0
    unmapped = Counter()
    n = 0
    for path in paths:
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            for i, row in enumerate(reader):
                if i >= int(args.limit_per_file):
                    break
                n += 1
                ion_mode = str(row.get("ion_mode") or "negative")
                label = str(row.get("adduct_label") or "")
                try:
                    mz = float(row.get("parent_mz") or "nan")
                    stored = float(row.get("neutral_mass") or "nan")
                except ValueError:
                    continue
                hit = row_for_observed_label(label, ion_mode)
                if hit is not None:
                    mapped += 1
                mass, source, shift = derive_target_neutral_mass(
                    mz, ion_mode, label, stored
                )
                if source == "label":
                    derived += 1
                    if shift is not None and abs(shift) > MOVE_DA:
                        moved += 1
                elif hit is None:
                    unmapped[label.split(":")[-1].strip() or label] += 1

    print("Excel ladder vs CHRO_LIB hypothesis sample")
    print(f"  rows_read={n} files={len(paths)} limit_per_file={args.limit_per_file}")
    print(f"  mapped_to_ladder={mapped} derived_non_default={derived}")
    print(f"  would_move_gt_{MOVE_DA}_Da={moved}")
    print("  top unmapped labels:")
    for label, count in unmapped.most_common(15):
        print(f"    {count:5d}  {label}")
    print()
    print("Do not rerun run_inference.sh for this inventory.")
    print("Analog regen only if those pools are still consumed:")
    print("  python -m structural_elucidation.evaluation.analog")
    print("Optional NEG-only polar filter: add --drop-polar-adducts")
    print("Candidate derivation (then full candidates-to-rerank + McNemar):")
    print("  python -m structural_elucidation.pipeline.candidates --derive-neutral-mass-from-adduct")


if __name__ == "__main__":
    main()
