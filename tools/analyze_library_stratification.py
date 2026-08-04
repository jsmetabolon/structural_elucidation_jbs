#!/usr/bin/env python3
"""Test whether pathway agreement varies across the chromatographic libraries.

The manuscript states that the delivered set is imbalanced across libraries and flags
library bias as a limitation, but never measures it. If agreement differed sharply
between libraries, the pooled figures would be an average over dissimilar populations
and the pooled interval would understate the uncertainty.

Only the inferred assignments are stratified, since direct lookups copy a reference
label and would measure the catalogue rather than the predictor.

The structural cohort is reported as counts only. It cannot be stratified: two of the
libraries contribute fewer than ten truth-verified controls, so a per-library recovery
rate would be noise. That thinness is itself the evidence for the imbalance caveat.

Inputs are the committed per-control audit outputs, so this runs in any checkout.

Usage:
    python3 tools/analyze_library_stratification.py
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter
from pathlib import Path
from typing import Any

# Below this many controls a per-library rate is not worth reporting as a comparison.
MIN_COHORT = 30


def wilson(correct: int, n: int, z: float = 1.959963984540054) -> dict[str, Any]:
    if n <= 0:
        return {"n": 0, "correct": 0, "accuracy": None, "ci95": [None, None]}
    p = correct / n
    denom = 1.0 + z * z / n
    center = (p + z * z / (2.0 * n)) / denom
    half = z * math.sqrt((p * (1.0 - p) + z * z / (4.0 * n)) / n) / denom
    return {
        "n": n,
        "correct": correct,
        "accuracy": round(p, 4),
        "ci95": [round(max(0.0, center - half), 4), round(min(1.0, center + half), 4)],
    }


def _truthy(value: str) -> bool:
    return value.strip().lower() == "true"


def intervals_overlap(blocks: list[dict[str, Any]]) -> bool:
    """True when every pair of confidence intervals overlaps."""

    usable = [block for block in blocks if block["accuracy"] is not None]
    if len(usable) < 2:
        return True
    return max(block["ci95"][0] for block in usable) <= min(block["ci95"][1] for block in usable)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--pathway-accuracy",
        type=Path,
        default=Path("docs/paper/data/pathway/audit_pathway_accuracy.csv"),
    )
    parser.add_argument(
        "--structure-controls",
        type=Path,
        default=Path("docs/paper/data/structure/audit_structure_controls.csv"),
    )
    parser.add_argument("--out", type=Path, default=Path("docs/paper/data/library_stratification.json"))
    args = parser.parse_args()

    with args.pathway_accuracy.open(newline="") as handle:
        pathway_rows = list(csv.DictReader(handle))
    with args.structure_controls.open(newline="") as handle:
        structure_rows = list(csv.DictReader(handle))

    inferred = [row for row in pathway_rows if row.get("source") == "neighbor_consensus"]
    by_library: dict[str, dict[str, Any]] = {}
    underpowered: dict[str, int] = {}
    for library in sorted({row["lib"].strip() for row in inferred if row["lib"].strip()}):
        cohort = [row for row in inferred if row["lib"].strip() == library]
        if len(cohort) < MIN_COHORT:
            underpowered[library] = len(cohort)
            continue
        by_library[library] = {
            "super": wilson(sum(1 for row in cohort if _truthy(row["correct_super"])), len(cohort)),
            "sub": wilson(sum(1 for row in cohort if _truthy(row["correct_sub"])), len(cohort)),
        }
    unlabelled = sum(1 for row in inferred if not row["lib"].strip())
    if unlabelled:
        underpowered["unlabelled"] = unlabelled

    structure_counts = Counter(row["lib"].strip() or "unlabelled" for row in structure_rows)

    summary = {
        "method": (
            "Inferred pathway assignments stratified by chromatographic library. Direct "
            "lookups are excluded because they copy a reference label rather than predict "
            f"one. Libraries contributing fewer than {MIN_COHORT} controls are listed as "
            "underpowered rather than scored."
        ),
        "min_cohort": MIN_COHORT,
        "pathway_inferred_by_library": by_library,
        "underpowered_cohorts": underpowered,
        "super_intervals_overlap": intervals_overlap([block["super"] for block in by_library.values()]),
        "sub_intervals_overlap": intervals_overlap([block["sub"] for block in by_library.values()]),
        "structure_truth_controls_by_library": dict(sorted(structure_counts.items())),
        "structure_stratification_supported": all(
            count >= MIN_COHORT for count in structure_counts.values()
        ),
        "scope_note": (
            "Consistency across libraries speaks to chromatographic method, not to "
            "chemical class. The libraries differ in ion mode and size, and all of them "
            "are scored against the same independent label source, so they share its "
            "annotation ceiling."
        ),
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2) + "\n")
    for library, block in by_library.items():
        print(
            f"[library-stratification] library {library}: n={block['super']['n']}, "
            f"super {block['super']['accuracy']:.1%}, sub {block['sub']['accuracy']:.1%}"
        )
    print(
        f"[library-stratification] super intervals overlap: {summary['super_intervals_overlap']}; "
        f"underpowered: {summary['underpowered_cohorts']}"
    )
    print(f"[library-stratification] wrote {args.out}")


if __name__ == "__main__":
    main()
