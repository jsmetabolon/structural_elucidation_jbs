#!/usr/bin/env python3
"""Sort the run-library compounds still missing from the 10 ppm summary.

Groups, each written as an ID list under ``--out-dir``:

- ``library10``: scored in the main pool, which searched the library at 5 ppm.
  These get a 10 ppm library re-search. Inputs are sliced into
  ``<out-dir>/library10_inputs`` for ``scripts/run_ppm10_extra_job.sh``.
- ``searched10_over``: scored by the 10 ppm recovery run and still over 10 ppm.
- ``nohit10``: sent to the 10 ppm recovery run with inputs and given no library
  candidate. These are the de novo set.
- ``no_inputs``: no feature row and hypothesis in the job or recovered-MS2 files.
- ``excluded``: Exclude or Artifact in the catalog.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from prepare_ppm10_recovery_inputs import (
    assign_inputs,
    catalog_targets,
    load_input_ids,
    scored_ids,
    write_inputs,
)


def classify(
    kept: list[str],
    excluded: set[str],
    published: set[str],
    main_scored: set[str],
    recovery_selected: set[str],
    recovery_scored: set[str],
) -> dict[str, list[str]]:
    """Split unpublished run-library IDs by why they are still missing."""

    groups: dict[str, list[str]] = {
        "library10": [],
        "searched10_over": [],
        "nohit10": [],
        "no_inputs": [],
    }
    for cid in kept:
        if cid in published:
            continue
        if cid in recovery_selected:
            groups["searched10_over" if cid in recovery_scored else "nohit10"].append(cid)
        elif cid in main_scored:
            groups["library10"].append(cid)
        else:
            groups["no_inputs"].append(cid)
    groups["excluded"] = sorted(excluded - published)
    return groups


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog", type=Path, required=True)
    parser.add_argument("--main-root", type=Path, required=True)
    parser.add_argument("--recovery-root", type=Path, required=True)
    parser.add_argument("--recovery-selection", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True, help="exec_summary_by_compound.csv")
    parser.add_argument("--feature-dir", type=Path, required=True)
    parser.add_argument("--recovered-dir", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args(argv)

    kept, excluded = catalog_targets(pd.read_csv(args.catalog, low_memory=False))
    summary = pd.read_csv(args.summary, usecols=["x_compound_id"], low_memory=False)
    published = set(summary["x_compound_id"].astype(str))
    selection = json.loads(args.recovery_selection.read_text(encoding="utf-8"))
    recovery_selected = {row["control_id"] for row in selection.get("assignments", [])}
    groups = classify(
        kept,
        excluded,
        published,
        scored_ids(args.main_root),
        recovery_selected,
        scored_ids(args.recovery_root),
    )

    args.out_dir.mkdir(parents=True, exist_ok=True)
    for name, ids in groups.items():
        (args.out_dir / f"{name}_ids.txt").write_text("\n".join(ids) + ("\n" if ids else ""), encoding="utf-8")

    inputs = load_input_ids(args.feature_dir, args.recovered_dir)
    assignments = assign_inputs(groups["library10"], *inputs)
    inputs_dir = args.out_dir / "library10_inputs"
    used = write_inputs(assignments, args.feature_dir, args.recovered_dir, inputs_dir)
    library10_counts = {
        "candidates": len(groups["library10"]),
        "selected": len(assignments),
        "no_inputs": len(groups["library10"]) - len(assignments),
    }
    (inputs_dir / "selection.json").write_text(
        json.dumps({"counts": library10_counts, "run_libs": ",".join(used), "assignments": assignments}, indent=2),
        encoding="utf-8",
    )

    counts = {name: len(ids) for name, ids in groups.items()}
    counts["catalog_run_libs"] = len(kept) + len(excluded)
    counts["published"] = len(published & (set(kept) | excluded))
    (args.out_dir / "remaining.json").write_text(
        json.dumps({"counts": counts, "library10": library10_counts}, indent=2), encoding="utf-8"
    )
    print("[remaining]", json.dumps(counts))
    print("[remaining] library10", json.dumps(library10_counts), f"libs={','.join(used) or 'none'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
