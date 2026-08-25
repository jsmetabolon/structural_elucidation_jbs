#!/usr/bin/env python3
"""Preliminary retention-screen measurement on the committed delivered predictions.

``tools/audit_topk_retention.py`` is the deciding audit for the per-candidate
retention screen, but it needs reranked top-10 outputs that a plain checkout
does not carry. This tool answers the same question on what *is* committed:
``data/integrated_best/`` records, per truth-verified control, the rank-1
candidate and the best candidate found in its top 10, each with its Tanimoto
similarity to the known structure. Joining those to the observed retention
indices in ``data/feature_rows/`` measures, on delivered predictions:

- how often the |z| >= 3 screen fires on a rank-1 candidate,
- whether the candidates it fires on are wrong (Tanimoto < 0.7) or correct
  (a false kill, which must not happen), and
- whether a flagged wrong rank-1 has a correct, unflagged alternative in its
  top 10 - the case where demotion visibly improves the internal report.

Two candidates per control instead of ten bounds the fire rate from below;
the full-list audit supersedes this one the moment outputs exist.

The retention math is imported from ``structural_elucidation.reports.confidence``
so this measures exactly the flag the report applies. A control's observed
index is its first RI-bearing feature row, which is also where the reranked
output's feature keys come from.

Usage:
    python3 tools/audit_integrated_best_retention.py
"""

from __future__ import annotations

import argparse
import csv
import glob
import json
import math
import sys
from pathlib import Path
from typing import Any, Callable, Sequence

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from structural_elucidation.reports.confidence import (  # noqa: E402
    RETENTION_IMPLAUSIBLE_Z,
    load_retention_calibration,
    retention_residual_z,
)
from structural_elucidation.utils.domain import ri_from_feature_key  # noqa: E402

# Same bars as tools/audit_topk_retention.py, read as a preliminary bound: two
# candidates per control cannot fire more often than ten would.
MIN_CONTROLS_FLAGGED_SHARE = 0.05
CORRECT_TANIMOTO = 0.7


def observed_ri_map(feature_rows_path: Path) -> dict[str, float]:
    """First finite observed retention index per control, from the feature keys."""
    out: dict[str, float] = {}
    with feature_rows_path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            control_id = str(row.get("control_id") or "").strip()
            ri = ri_from_feature_key(row.get("feature_key"))
            if control_id and math.isfinite(ri) and control_id not in out:
                out[control_id] = ri
    return out


def evaluate_controls(
    rows: Sequence[dict[str, Any]],
    observed: dict[str, float],
    lib_id: str,
    calibration: dict[str, dict],
    proxy_fn: Callable[[str], float],
    implausible_z: float = RETENTION_IMPLAUSIBLE_Z,
) -> dict[str, Any]:
    """Score each control's rank-1 and best-in-top-10 candidates and classify hits."""

    def z_for(observed_ri: float, smiles: str) -> float | None:
        proxy = proxy_fn(str(smiles or ""))
        return retention_residual_z(
            {
                "lib_id": lib_id,
                "observed_retention_index": observed_ri,
                "candidate_retention_proxy": proxy if math.isfinite(proxy) else None,
            },
            calibration,
        )

    scored = 0
    abs_z_rank1: list[float] = []
    rank1_flagged_wrong: list[dict[str, Any]] = []
    rank1_flagged_correct: list[dict[str, Any]] = []
    correct_best_flagged: list[dict[str, Any]] = []
    improvement_cases: list[dict[str, Any]] = []
    for row in rows:
        control_id = str(row.get("control_id") or "").strip()
        observed_ri = observed.get(control_id)
        if observed_ri is None:
            continue
        z_rank1 = z_for(observed_ri, row.get("rank1_smiles") or "")
        z_best = z_for(observed_ri, row.get("top10_smiles") or "")
        if z_rank1 is None and z_best is None:
            continue
        scored += 1
        rank1_tanimoto = float(row.get("rank1_best_tanimoto") or 0.0)
        best_tanimoto = float(row.get("top10_best_tanimoto") or 0.0)
        if z_rank1 is not None:
            abs_z_rank1.append(abs(z_rank1))
            if abs(z_rank1) >= implausible_z:
                hit = {
                    "control_id": control_id,
                    "z": round(z_rank1, 2),
                    "rank1_tanimoto": round(rank1_tanimoto, 4),
                }
                if rank1_tanimoto >= CORRECT_TANIMOTO:
                    rank1_flagged_correct.append(hit)
                else:
                    rank1_flagged_wrong.append(hit)
                    if (
                        best_tanimoto >= CORRECT_TANIMOTO
                        and z_best is not None
                        and abs(z_best) < implausible_z
                    ):
                        improvement_cases.append(
                            {
                                **hit,
                                "correct_alternative_z": round(z_best, 2),
                                "correct_alternative_tanimoto": round(best_tanimoto, 4),
                                "correct_alternative_rank": row.get("top10_rank") or None,
                            }
                        )
        if (
            z_best is not None
            and best_tanimoto >= CORRECT_TANIMOTO
            and abs(z_best) >= implausible_z
        ):
            correct_best_flagged.append(
                {
                    "control_id": control_id,
                    "z": round(z_best, 2),
                    "best_tanimoto": round(best_tanimoto, 4),
                }
            )
    abs_z_rank1.sort()
    p90 = (
        abs_z_rank1[max(0, math.ceil(0.9 * len(abs_z_rank1)) - 1)] if abs_z_rank1 else None
    )
    flagged_controls = len(rank1_flagged_wrong) + len(rank1_flagged_correct)
    return {
        "status": "evaluated",
        "n_controls": len(rows),
        "controls_scored": scored,
        "abs_z_rank1": {
            "p90": round(p90, 2) if p90 is not None else None,
            "max": round(abs_z_rank1[-1], 2) if abs_z_rank1 else None,
        },
        "rank1_flagged_controls": flagged_controls,
        "share_rank1_flagged": round(flagged_controls / scored, 4) if scored else None,
        "rank1_flagged_wrong": sorted(rank1_flagged_wrong, key=lambda item: item["control_id"]),
        "rank1_flagged_correct": sorted(
            rank1_flagged_correct, key=lambda item: item["control_id"]
        ),
        "correct_best_in_top10_flagged": sorted(
            correct_best_flagged, key=lambda item: item["control_id"]
        ),
        "improvement_cases": sorted(improvement_cases, key=lambda item: item["control_id"]),
    }


def preliminary_reading(libraries: dict[str, dict[str, Any]]) -> dict[str, Any]:
    """The same bars as the deciding audit, read as a lower bound."""
    evaluated = {
        lib: block for lib, block in libraries.items() if block.get("status") == "evaluated"
    }
    base = {
        "min_share_controls_with_flag": MIN_CONTROLS_FLAGGED_SHARE,
        "max_labeled_correct_flagged": 0,
    }
    if not evaluated:
        return {
            **base,
            "verdict": "NOT_EVALUATED",
            "reasoning": "No library produced scoreable controls.",
        }
    shares = [
        block["share_rank1_flagged"]
        for block in evaluated.values()
        if block["share_rank1_flagged"] is not None
    ]
    best_share = max(shares) if shares else None
    false_kills = sum(
        len(block["rank1_flagged_correct"]) + len(block["correct_best_in_top10_flagged"])
        for block in evaluated.values()
    )
    fires_enough = best_share is not None and best_share >= MIN_CONTROLS_FLAGGED_SHARE
    verdict = "GO" if fires_enough and false_kills == 0 else "NO_GO"
    if verdict == "GO":
        reasoning = (
            f"Scoring only two candidates per control the screen already touches "
            f"{best_share:.1%} of scoreable controls on its best library, every flagged "
            "candidate is wrong against the known structure, and no correct candidate is "
            "flagged. The full top-10 audit can only raise the fire rate; the false-kill "
            "bar it must still clear is the same."
        )
    else:
        reasons = []
        if not fires_enough:
            reasons.append(
                f"the screen touches {0.0 if best_share is None else best_share:.1%} of "
                f"scoreable controls at best, under the {MIN_CONTROLS_FLAGGED_SHARE:.0%} bar "
                "- though ten candidates may fire more often than these two"
            )
        if false_kills:
            reasons.append(
                f"{false_kills} correct candidate(s) would be flagged, which the deciding "
                "audit treats as disqualifying"
            )
        reasoning = "Preliminary reading is negative: " + "; ".join(reasons) + "."
    return {
        **base,
        "best_share_rank1_flagged": best_share,
        "correct_candidates_flagged": false_kills,
        "verdict": verdict,
        "reasoning": reasoning,
    }


def run_audit(
    feature_rows_dir: Path,
    integrated_dir: Path,
    calibration: dict[str, dict],
    proxy_fn: Callable[[str], float],
    implausible_z: float = RETENTION_IMPLAUSIBLE_Z,
) -> dict[str, Any]:
    libraries: dict[str, dict[str, Any]] = {}
    for lib in sorted(calibration):
        integrated_matches = sorted(glob.glob(str(integrated_dir / f"lib_{lib}_*.csv")))
        feature_matches = sorted(
            glob.glob(str(feature_rows_dir / f"lib{lib}_*_feature_rows.csv"))
        )
        if not integrated_matches or not feature_matches:
            libraries[lib] = {
                "status": "missing_inputs",
                "reason": "No committed integrated-best or feature-rows file for this library.",
            }
            continue
        with open(integrated_matches[0], newline="") as handle:
            rows = list(csv.DictReader(handle))
        libraries[lib] = evaluate_controls(
            rows,
            observed_ri_map(Path(feature_matches[0])),
            lib,
            calibration,
            proxy_fn,
            implausible_z,
        )
    return {
        "question": (
            "On the committed delivered predictions, does the |z| >= 3 retention screen "
            "fire on wrong candidates without ever firing on correct ones?"
        ),
        "cohort": (
            "Truth-verified controls in data/integrated_best (rank-1 candidate and best-in-"
            "top-10 candidate per control), observed retention indices joined from "
            "data/feature_rows feature keys."
        ),
        "screen": {
            "z": "structural_elucidation.reports.confidence.retention_residual_z",
            "implausible_abs_z": implausible_z,
            "correct_tanimoto": CORRECT_TANIMOTO,
        },
        "libraries": libraries,
        "preliminary_reading": preliminary_reading(libraries),
        "scope_note": (
            "A preliminary lower bound for tools/audit_topk_retention.py: it scores two "
            "candidates per control instead of the full reported list, on the cohort whose "
            "structures are known. The full-list audit on reranked outputs supersedes it."
        ),
    }


def render_markdown(summary: dict[str, Any]) -> str:
    lines = [
        "# Retention screen on the committed delivered predictions",
        "",
        "Generated by `tools/audit_integrated_best_retention.py`. A preliminary lower",
        "bound for the top-10 retention screen, measured on committed inputs only.",
        "",
        "| Library | Controls scored | Rank-1 flagged | ...of which wrong | False kills | Improvement cases |",
        "|---|---|---|---|---|---|",
    ]
    for lib_id, block in summary["libraries"].items():
        if block.get("status") != "evaluated":
            lines.append(f"| {lib_id} | - | - | - | - | {block['status']} |")
            continue
        false_kills = len(block["rank1_flagged_correct"]) + len(
            block["correct_best_in_top10_flagged"]
        )
        lines.append(
            f"| {lib_id} | {block['controls_scored']} | {block['rank1_flagged_controls']} | "
            f"{len(block['rank1_flagged_wrong'])} | {false_kills} | "
            f"{len(block['improvement_cases'])} |"
        )
    reading = summary["preliminary_reading"]
    lines += ["", "## Preliminary reading", "", f"**{reading['verdict']}**", "", reading["reasoning"], ""]
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--feature-rows-dir", type=Path, default=Path("data/feature_rows"))
    parser.add_argument("--integrated-dir", type=Path, default=Path("data/integrated_best"))
    parser.add_argument(
        "--calibration",
        type=Path,
        default=Path("docs/paper/data/retention/ri_calibration.json"),
    )
    parser.add_argument("--out-dir", type=Path, default=Path("docs/paper/data/retention"))
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    from structural_elucidation.nce.evidence import HAS_RDKIT, ri_proxy_from_smiles

    if not HAS_RDKIT:
        summary: dict[str, Any] = {
            "status": "not_evaluated",
            "reason": "RDKit is unavailable, so no candidate retention proxy can be computed.",
        }
    else:
        summary = run_audit(
            feature_rows_dir=args.feature_rows_dir,
            integrated_dir=args.integrated_dir,
            calibration=load_retention_calibration(args.calibration),
            proxy_fn=ri_proxy_from_smiles,
        )

    json_path = args.out_dir / "audit_integrated_best_retention.json"
    json_path.write_text(json.dumps(summary, indent=2) + "\n")
    markdown_path = args.out_dir / "audit_integrated_best_retention.md"
    if summary.get("status") == "not_evaluated":
        markdown_path.write_text(f"Not evaluated: {summary['reason']}\n")
        print(f"[integrated-retention] not evaluated: {summary['reason']}")
    else:
        markdown_path.write_text(render_markdown(summary) + "\n")
        for lib_id, block in summary["libraries"].items():
            if block.get("status") != "evaluated":
                print(f"[integrated-retention] lib {lib_id}: {block['status']}")
                continue
            print(
                f"[integrated-retention] lib {lib_id}: scored={block['controls_scored']} "
                f"rank1_flagged={block['rank1_flagged_controls']} "
                f"wrong={len(block['rank1_flagged_wrong'])} "
                f"false_kills={len(block['rank1_flagged_correct']) + len(block['correct_best_in_top10_flagged'])} "
                f"improvements={len(block['improvement_cases'])}"
            )
        print(f"[integrated-retention] preliminary -> {summary['preliminary_reading']['verdict']}")
    print(f"[integrated-retention] wrote {json_path} and {markdown_path}")


if __name__ == "__main__":
    main()
