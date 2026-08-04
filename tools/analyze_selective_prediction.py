#!/usr/bin/env python3
"""Measure whether reporting only the most confident calls improves their accuracy.

Both documents describe the workflow as an instrument for prioritization, which is a
claim about ranking: if a user works down the list from most to least confident, the
calls they act on first should be the ones most likely to be right. Aggregate accuracy
cannot show that. This computes the operating curve directly, by sorting each cohort
from most to least confident and scoring successively larger prefixes of it.

The two cohorts are ranked by different quantities because they are graded by different
mechanisms: structural calls carry a leave-one-out identification-confidence probability,
whereas inferred pathway assignments carry the similarity of the closest reference
compound. Both are the value a user would actually sort on.

Inputs are the committed per-control audit outputs, so this runs in any checkout and
needs no controlled report bundle.

Usage:
    python3 tools/analyze_selective_prediction.py
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

# Cumulative prefixes of each cohort, most confident first.
COVERAGE_GRID = (0.2, 0.4, 0.6, 0.8, 1.0)


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
    return value.strip().lower() in {"true", "1", "yes"}


def _float(value: str) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(number) else number


def curve(rows: list[dict[str, str]], rank_by: str, correct_column: str) -> list[dict[str, Any]]:
    """Score cumulative prefixes of ``rows`` sorted by descending ``rank_by``."""

    scored = [row for row in rows if _float(row.get(rank_by, "")) is not None]
    scored.sort(key=lambda row: -float(row[rank_by]))
    points = []
    for coverage in COVERAGE_GRID:
        kept = scored[: max(1, round(len(scored) * coverage))]
        correct = sum(1 for row in kept if _truthy(row.get(correct_column, "")))
        point = wilson(correct, len(kept))
        point["coverage"] = coverage
        points.append(point)
    return points


def discrimination(points: list[dict[str, Any]]) -> float | None:
    """Accuracy at the most confident prefix minus accuracy over the whole cohort.

    Positive means confidence ranks usefully; near zero means it does not rank at all.
    """

    if not points or points[0]["accuracy"] is None or points[-1]["accuracy"] is None:
        return None
    return round(points[0]["accuracy"] - points[-1]["accuracy"], 4)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--structure-controls", type=Path, default=Path("docs/paper/data/structure/audit_structure_controls.csv"))
    parser.add_argument("--pathway-accuracy", type=Path, default=Path("docs/paper/data/pathway/audit_pathway_accuracy.csv"))
    parser.add_argument("--out", type=Path, default=Path("docs/paper/data/selective_prediction.json"))
    args = parser.parse_args()

    with args.structure_controls.open(newline="") as handle:
        structure_rows = list(csv.DictReader(handle))
    with args.pathway_accuracy.open(newline="") as handle:
        pathway_rows = list(csv.DictReader(handle))

    # Only the inferred path is a prediction; direct lookups copy a reference label and
    # would flatter the curve without testing anything.
    inferred = [row for row in pathway_rows if row.get("source") == "neighbor_consensus"]

    structure_points = curve(structure_rows, "loo_probability_full", "top1_correct_0p70")
    super_points = curve(inferred, "best_similarity", "correct_super")
    sub_points = curve(inferred, "best_similarity", "correct_sub")

    summary = {
        "method": (
            "Each cohort is sorted from most to least confident and scored over cumulative "
            "prefixes. A rising curve means the confidence signal ranks usefully; a flat "
            "curve means it does not, whatever its absolute accuracy."
        ),
        "coverage_grid": list(COVERAGE_GRID),
        "structure_rank1": {
            "cohort": "truth-verified controls",
            "ranked_by": "leave-one-out identification-confidence probability",
            "correctness": "rank-1 truth Tanimoto >= 0.70",
            "n": len(structure_rows),
            "points": structure_points,
            "discrimination": discrimination(structure_points),
        },
        "pathway_super": {
            "cohort": "neighbor-consensus controls with independent known_X labels",
            "ranked_by": "best-reference structural similarity",
            "correctness": "predicted super-pathway equals known_X label",
            "n": len(inferred),
            "points": super_points,
            "discrimination": discrimination(super_points),
        },
        "pathway_sub": {
            "cohort": "neighbor-consensus controls with independent known_X labels",
            "ranked_by": "best-reference structural similarity",
            "correctness": "predicted subpathway equals known_X label",
            "n": len(inferred),
            "points": sub_points,
            "discrimination": discrimination(sub_points),
        },
        "scope_note": (
            "Retrospective selective prediction on fixed cohorts. It shows how accuracy "
            "would change if low-confidence calls were withheld, not the accuracy of a "
            "deployed abstention rule, which would need a prospectively chosen threshold."
        ),
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2) + "\n")
    for key in ("structure_rank1", "pathway_super", "pathway_sub"):
        block = summary[key]
        first, last = block["points"][0], block["points"][-1]
        print(
            f"[selective-prediction] {key}: top {COVERAGE_GRID[0]:.0%} "
            f"{first['accuracy']:.1%} (n={first['n']}) vs all {last['accuracy']:.1%} "
            f"(n={last['n']}); discrimination {block['discrimination']:+.4f}"
        )
    print(f"[selective-prediction] wrote {args.out}")


if __name__ == "__main__":
    main()
