#!/usr/bin/env python3
"""Separate taxonomy-vocabulary disagreement from prediction error in pathway results.

Pathway assignments are scored against an independent label source that uses a
different vocabulary, so the raw agreement figures understate the predictor by an
unknown amount. This measures that amount two ways and, importantly, measures the
same correction against the annotation ceiling, because a correction that lifts the
achievable rate as much as it lifts the observed rate has explained nothing.

Two corrections are applied, both to the super-pathway and subpathway levels
respectively:

* Super-pathway: amino acid and peptide are treated as one class. They are a
  structural continuum rather than distinct categories, and dipeptides filed under
  one heading in one scheme and the other heading in the other are the single
  largest source of disagreement between the two vocabularies.
* Subpathway: rows whose truth label is not a taxonomy label at all are excluded.
  The independent label column carries internal review notes ("Ask ...") and explicit
  non-answers ("Unknown ..."), against which no prediction can ever score.

A maximal collapse of every commonly confused class is also computed. It is reported
only as an upper bound and is deliberately not adopted, because merging four of the
eight super-pathways would destroy the distinction the label exists to make.

Inputs are the committed per-control audit output, so this runs in any checkout.

Usage:
    python3 tools/analyze_taxonomy_harmonization.py
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from collections import Counter
from pathlib import Path
from typing import Any, Iterable, Sequence

# Structural continuum: an amino acid, a dipeptide and a polypeptide differ by degree.
SUPER_EQUIVALENCE: tuple[frozenset[str], ...] = (frozenset({"amino acid", "peptide"}),)

# Reported as a bound only. Xenobiotics classifies by origin rather than structure, so
# it overlaps every structural class; collapsing it would inflate agreement meaninglessly.
MAXIMAL_COLLAPSE: tuple[frozenset[str], ...] = (
    frozenset({"amino acid", "peptide", "lipid", "xenobiotics"}),
)

UNSCOREABLE_SUB = re.compile(r"^(?:ask\b|unknown\b)", re.IGNORECASE)


def norm(value: str) -> str:
    return " ".join(str(value).strip().lower().split())


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


def canonical(label: str, groups: Sequence[frozenset[str]]) -> str:
    for group in groups:
        if label in group:
            return "|".join(sorted(group))
    return label


def super_agreement(rows: Iterable[dict[str, str]], groups: Sequence[frozenset[str]]) -> dict[str, Any]:
    scoreable = [row for row in rows if norm(row["true_super"])]
    correct = sum(
        1
        for row in scoreable
        if canonical(norm(row["true_super"]), groups) == canonical(norm(row["predicted_pathway"]), groups)
    )
    return wilson(correct, len(scoreable))


def sub_agreement(rows: Iterable[dict[str, str]], *, drop_unscoreable: bool) -> dict[str, Any]:
    scoreable = []
    for row in rows:
        label = row["true_sub"].strip()
        if not label:
            continue
        if drop_unscoreable and UNSCOREABLE_SUB.match(label):
            continue
        scoreable.append(row)
    correct = sum(1 for row in scoreable if row["correct_sub"].strip().lower() == "true")
    return wilson(correct, len(scoreable))


def share_of_ceiling(observed: dict[str, Any], ceiling: dict[str, Any]) -> float | None:
    """How much of the achievable agreement the predictor actually reaches."""

    if not observed["accuracy"] or not ceiling["accuracy"]:
        return None
    return round(observed["accuracy"] / ceiling["accuracy"], 4)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--pathway-accuracy",
        type=Path,
        default=Path("docs/paper/data/pathway/audit_pathway_accuracy.csv"),
    )
    parser.add_argument("--out", type=Path, default=Path("docs/paper/data/taxonomy_harmonization.json"))
    args = parser.parse_args()

    with args.pathway_accuracy.open(newline="") as handle:
        rows = list(csv.DictReader(handle))

    # The reference cohort is where the assignment copies the reference label, so its
    # agreement with the independent labels is the ceiling any predictor can reach.
    cohorts = {
        "all_assignments": rows,
        "inferred_only": [row for row in rows if row.get("source") == "neighbor_consensus"],
        "reference_ceiling": [row for row in rows if row.get("in_reference_db", "").strip().lower() == "true"],
    }

    def super_block(groups: Sequence[frozenset[str]]) -> dict[str, Any]:
        block = {name: super_agreement(cohort, groups) for name, cohort in cohorts.items()}
        block["inferred_share_of_ceiling"] = share_of_ceiling(
            block["inferred_only"], block["reference_ceiling"]
        )
        return block

    def sub_block(*, drop_unscoreable: bool) -> dict[str, Any]:
        block = {
            name: sub_agreement(cohort, drop_unscoreable=drop_unscoreable)
            for name, cohort in cohorts.items()
        }
        block["inferred_share_of_ceiling"] = share_of_ceiling(
            block["inferred_only"], block["reference_ceiling"]
        )
        return block

    excluded = Counter(
        row["true_sub"].strip()
        for row in rows
        if row["true_sub"].strip() and UNSCOREABLE_SUB.match(row["true_sub"].strip())
    )

    summary = {
        "method": (
            "Agreement recomputed under vocabulary corrections, applied identically to the "
            "predictor and to the annotation ceiling. The share of the ceiling reached is the "
            "figure to compare across scenarios; a correction that moves both equally has "
            "revealed a measurement artefact, not a better predictor."
        ),
        "super_pathway": {
            "equivalence_applied": [sorted(group) for group in SUPER_EQUIVALENCE],
            "raw": super_block(()),
            "harmonized": super_block(SUPER_EQUIVALENCE),
            "rejected_maximal_collapse": {
                "classes_merged": [sorted(group) for group in MAXIMAL_COLLAPSE],
                "reason": (
                    "Merges four of the eight super-pathways, including an origin-based class "
                    "that overlaps every structural one. Reported as an upper bound only."
                ),
                **super_block(MAXIMAL_COLLAPSE),
            },
        },
        "subpathway": {
            "raw": sub_block(drop_unscoreable=False),
            "scoreable_only": sub_block(drop_unscoreable=True),
            "excluded_labels": dict(excluded.most_common()),
            "excluded_rows": sum(excluded.values()),
            "excluded_reason": (
                "Internal review notes and explicit non-answers in the independent label "
                "column, against which no prediction can score."
            ),
        },
        "scope_note": (
            "Corrections are applied to the committed per-control comparison, not to the "
            "production predictor or to the reference database. A curated crosswalk between "
            "the two full subpathway vocabularies is not attempted here."
        ),
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2) + "\n")

    sup = summary["super_pathway"]
    sub = summary["subpathway"]
    print(
        f"[taxonomy-harmonization] super inferred {sup['raw']['inferred_only']['accuracy']:.1%} "
        f"-> {sup['harmonized']['inferred_only']['accuracy']:.1%}; share of ceiling "
        f"{sup['raw']['inferred_share_of_ceiling']:.1%} -> {sup['harmonized']['inferred_share_of_ceiling']:.1%}"
    )
    print(
        f"[taxonomy-harmonization] sub inferred {sub['raw']['inferred_only']['accuracy']:.1%} "
        f"-> {sub['scoreable_only']['inferred_only']['accuracy']:.1%}; share of ceiling "
        f"{sub['raw']['inferred_share_of_ceiling']:.1%} -> "
        f"{sub['scoreable_only']['inferred_share_of_ceiling']:.1%} "
        f"({sub['excluded_rows']} unscoreable rows excluded)"
    )
    print(f"[taxonomy-harmonization] wrote {args.out}")


if __name__ == "__main__":
    main()
