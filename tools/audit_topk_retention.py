#!/usr/bin/env python3
"""Measure what a retention screen would do to the reported top-10 candidates.

The feasibility audit (``tools/audit_rt_feasibility.py``) established that the
recalibrated retention axis is soft and demote-only, and the confidence model
already applies it - to the rank-1 candidate alone. This audit asks the next
question: if the same standardized residual were computed for every reported
top-10 row, how often would the |z| >= 3 implausibility bar actually fire, whom
would it hit, and would it ever hit a candidate the truth set says is correct?

The answer decides the go / no-go gate for wiring a per-candidate flag into the
internal report:

- **GO** only if at least ``MIN_CONTROLS_FLAGGED_SHARE`` of scoreable controls
  have a flagged row on some library, and no ``label_ge_0p7`` row - a candidate
  known to be essentially correct - is ever flagged. A screen that fires almost
  never is not worth a rendering surface; one that fires on correct answers is
  worse than none.

Retention math is imported from ``structural_elucidation.reports.confidence``
(``retention_residual_z``) rather than reimplemented, so this audit measures
exactly the flag the report would apply. Libraries absent from
``ri_calibration.json`` did not clear the feasibility audit and are reported as
having no retention axis rather than a zero score.

The observed retention index is identical across a control's rows (it belongs
to the feature, not the candidate), so within-list discrimination comes
entirely from the per-candidate proxy.

Usage:
    python3 tools/audit_topk_retention.py \
        --rerank-root output/pipeline_equalization --libs 209,400
"""

from __future__ import annotations

import argparse
import csv
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

# Wire the per-row screen into the internal report only if it would actually
# touch a meaningful share of compounds. Below this it is a rendering surface
# nobody will ever see fire.
MIN_CONTROLS_FLAGGED_SHARE = 0.05
# ...and never fires on a candidate the truth set marks correct. The 3-sigma
# bar was chosen because the worst true structure in the feasibility audit sits
# at z ~ 2.9, so any correct candidate flagged here means the bar is wrong.
MAX_LABELED_CORRECT_FLAGGED = 0


def parse_label(value: Any) -> bool | None:
    """Interpret ``label_ge_0p7``: True/False when labeled, None when not."""
    text = str(value if value is not None else "").strip().lower()
    if text in {"", "nan", "none"}:
        return None
    if text == "true":
        return True
    if text == "false":
        return False
    try:
        return float(text) >= 1.0
    except ValueError:
        return None


def _int_or_none(value: Any) -> int | None:
    try:
        number = float(str(value).strip())
    except (TypeError, ValueError):
        return None
    return int(number) if math.isfinite(number) else None


def score_rows(
    rows: Sequence[dict[str, Any]],
    lib_id: Any,
    calibration: dict[str, dict],
    proxy_fn: Callable[[str], float],
    implausible_z: float = RETENTION_IMPLAUSIBLE_Z,
) -> list[dict[str, Any]]:
    """Attach the retention residual z and flag to each reported candidate row.

    A row scores only when everything the axis needs is there: an observed
    index embedded in its feature key, a finite proxy for its SMILES, and a
    library whose fit cleared the feasibility audit. Anything else leaves
    ``z`` as None, which downstream counts as unavailable rather than zero.
    """
    scored: list[dict[str, Any]] = []
    for row in rows:
        observed = ri_from_feature_key(row.get("feature_key"))
        smiles = str(row.get("predicted_smiles") or "").strip()
        proxy = proxy_fn(smiles) if smiles else math.nan
        row_lib = str(row.get("lib_id") or "").strip() or str(lib_id)
        z = retention_residual_z(
            {
                "lib_id": row_lib,
                "observed_retention_index": observed if math.isfinite(observed) else None,
                "candidate_retention_proxy": proxy if math.isfinite(proxy) else None,
            },
            calibration,
        )
        scored.append(
            {
                "control_id": str(row.get("control_id") or "").strip(),
                "rank": _int_or_none(row.get("rank")),
                "rerank_policy": str(row.get("rerank_policy") or "").strip(),
                "observed_ri": observed if math.isfinite(observed) else None,
                "proxy_ri": proxy if math.isfinite(proxy) else None,
                "z": z,
                "flagged": z is not None and abs(z) >= implausible_z,
                "label_correct": parse_label(row.get("label_ge_0p7")),
            }
        )
    return scored


def percentile(sorted_values: Sequence[float], share: float) -> float | None:
    """Nearest-rank percentile of an ascending sequence."""
    if not sorted_values:
        return None
    index = min(len(sorted_values) - 1, max(0, math.ceil(share * len(sorted_values)) - 1))
    return sorted_values[index]


def summarize_library(
    scored: Sequence[dict[str, Any]], max_examples: int = 15
) -> dict[str, Any]:
    """Aggregate one library's scored rows into the audit's per-library block."""
    controls = {row["control_id"] for row in scored if row["control_id"]}
    usable = [row for row in scored if row["z"] is not None]
    controls_scored = {row["control_id"] for row in usable}
    abs_z = sorted(abs(row["z"]) for row in usable)
    flagged = [row for row in usable if row["flagged"]]
    controls_with_flag = {row["control_id"] for row in flagged}
    rank1_flagged = {row["control_id"] for row in flagged if row["rank"] == 1}
    by_rank: dict[str, int] = {}
    for row in flagged:
        key = str(row["rank"]) if row["rank"] is not None else "unranked"
        by_rank[key] = by_rank.get(key, 0) + 1
    labeled = [row for row in usable if row["label_correct"] is not None]
    labeled_correct = [row for row in labeled if row["label_correct"]]
    correct_flagged = [row for row in labeled_correct if row["flagged"]]

    def _example(row: dict[str, Any]) -> dict[str, Any]:
        return {
            "control_id": row["control_id"],
            "rank": row["rank"],
            "rerank_policy": row["rerank_policy"],
            "z": round(row["z"], 2),
            "observed_ri": round(row["observed_ri"], 1),
            "proxy_ri": round(row["proxy_ri"], 1),
        }

    return {
        "status": "evaluated",
        "n_rows": len(scored),
        "n_controls": len(controls),
        "coverage": {
            "rows_with_observed_ri": sum(1 for row in scored if row["observed_ri"] is not None),
            "rows_with_proxy": sum(1 for row in scored if row["proxy_ri"] is not None),
            "rows_scored": len(usable),
            "share_rows_scored": round(len(usable) / len(scored), 4) if scored else None,
            "controls_scored": len(controls_scored),
        },
        "abs_z_distribution": {
            "p50": round(percentile(abs_z, 0.50), 3) if abs_z else None,
            "p90": round(percentile(abs_z, 0.90), 3) if abs_z else None,
            "p99": round(percentile(abs_z, 0.99), 3) if abs_z else None,
            "max": round(abs_z[-1], 3) if abs_z else None,
        },
        "flags": {
            "rows_flagged": len(flagged),
            "share_rows_flagged": round(len(flagged) / len(usable), 4) if usable else None,
            "controls_with_flag": len(controls_with_flag),
            "share_controls_with_flag": round(len(controls_with_flag) / len(controls_scored), 4)
            if controls_scored
            else None,
            "rank1_flagged_controls": len(rank1_flagged),
            "flagged_rows_by_rank": dict(sorted(by_rank.items(), key=lambda item: item[0])),
        },
        "false_kills": {
            "labeled_rows_scored": len(labeled),
            "labeled_correct_rows_scored": len(labeled_correct),
            "labeled_correct_rows_flagged": len(correct_flagged),
            "labeled_correct_flagged_examples": [_example(row) for row in correct_flagged],
            "max_abs_z_labeled_correct": round(max(abs(row["z"]) for row in labeled_correct), 3)
            if labeled_correct
            else None,
        },
        "flagged_examples": [
            _example(row)
            for row in sorted(flagged, key=lambda item: -abs(item["z"]))[: max(0, max_examples)]
        ],
    }


def gate_verdict(libraries: dict[str, dict[str, Any]]) -> dict[str, Any]:
    """GO / NO-GO for wiring the per-row flag into the internal report."""
    evaluated = {lib: block for lib, block in libraries.items() if block.get("status") == "evaluated"}
    base = {
        "min_share_controls_with_flag": MIN_CONTROLS_FLAGGED_SHARE,
        "max_labeled_correct_flagged": MAX_LABELED_CORRECT_FLAGGED,
    }
    if not evaluated:
        return {
            **base,
            "verdict": "NOT_EVALUATED",
            "reasoning": "No library produced scoreable rows, so there is nothing to decide on.",
        }
    shares = [
        block["flags"]["share_controls_with_flag"]
        for block in evaluated.values()
        if block["flags"]["share_controls_with_flag"] is not None
    ]
    best_share = max(shares) if shares else None
    correct_flagged = sum(
        block["false_kills"]["labeled_correct_rows_flagged"] for block in evaluated.values()
    )
    fires_enough = best_share is not None and best_share >= MIN_CONTROLS_FLAGGED_SHARE
    kills_none = correct_flagged <= MAX_LABELED_CORRECT_FLAGGED
    verdict = "GO" if fires_enough and kills_none else "NO_GO"
    if verdict == "GO":
        reasoning = (
            f"On its best library the screen touches {best_share:.1%} of scoreable controls "
            "and never fires on a candidate the truth set marks correct, so a demote-only "
            "flag in the internal report has something real to show without costing a "
            "known-correct answer."
        )
    else:
        reasons = []
        if not fires_enough:
            reasons.append(
                f"the screen touches {0.0 if best_share is None else best_share:.1%} of scoreable "
                f"controls at best, under the {MIN_CONTROLS_FLAGGED_SHARE:.0%} bar"
            )
        if not kills_none:
            reasons.append(
                f"{correct_flagged} labeled-correct row(s) would be flagged, and a screen that "
                "fires on known-correct answers is worse than none"
            )
        reasoning = "Do not wire the per-row flag: " + "; ".join(reasons) + "."
    return {
        **base,
        "best_share_controls_with_flag": best_share,
        "labeled_correct_rows_flagged": correct_flagged,
        "verdict": verdict,
        "reasoning": reasoning,
    }


def load_top10_rows(path: Path) -> list[dict[str, Any]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def run_audit(
    rerank_root: Path,
    libs: Sequence[str],
    top10_name: str,
    calibration: dict[str, dict],
    proxy_fn: Callable[[str], float],
    implausible_z: float = RETENTION_IMPLAUSIBLE_Z,
    max_examples: int = 15,
) -> dict[str, Any]:
    """Score every requested library's top-10 CSV and assemble the summary."""
    libraries: dict[str, dict[str, Any]] = {}
    for lib in libs:
        lib = str(lib).strip()
        if not lib:
            continue
        if lib not in calibration:
            libraries[lib] = {
                "status": "no_calibrated_axis",
                "reason": (
                    "Library did not clear tools/audit_rt_feasibility.py, so it has no "
                    "retention axis to screen with."
                ),
            }
            continue
        top10_path = rerank_root / f"rerank_lib{lib}" / top10_name
        if not top10_path.exists():
            libraries[lib] = {
                "status": "missing_top10",
                "reason": f"No reranked output at {top10_path}.",
            }
            continue
        scored = score_rows(
            load_top10_rows(top10_path), lib, calibration, proxy_fn, implausible_z
        )
        libraries[lib] = summarize_library(scored, max_examples=max_examples)
    return {
        "question": (
            "Would a demote-only retention flag on every reported top-10 candidate fire "
            "often enough to be worth rendering, without ever firing on a known-correct one?"
        ),
        "screen": {
            "z": "structural_elucidation.reports.confidence.retention_residual_z",
            "implausible_abs_z": implausible_z,
            "observed_ri": "parsed from each row's feature_key (identical across a control's rows)",
            "candidate_ri": "structural_elucidation.nce.evidence.ri_proxy_from_smiles per candidate",
        },
        "libraries": libraries,
        "gate": gate_verdict(libraries),
        "scope_note": (
            "Rows are counted as reported (unified_rerank_all_top10.csv includes the non-S0 "
            "policy lanes). label_ge_0p7 marks candidates essentially correct for the truth-"
            "verified controls only, so false-kill counts are bounded by that cohort, not the "
            "full X set."
        ),
    }


def render_markdown(summary: dict[str, Any]) -> str:
    lines = [
        "# Top-10 retention screen audit",
        "",
        "Generated by `tools/audit_topk_retention.py`. Whether a demote-only retention",
        "flag on every reported top-10 candidate is worth wiring into the internal",
        "report, measured on the reranked outputs before anything renders.",
        "",
    ]
    if summary.get("status") == "not_evaluated":
        return "\n".join(lines + [f"Not evaluated: {summary['reason']}", ""])
    lines += [
        "| Library | Rows | Scored | Controls with flag | Rows flagged | Rank-1 flagged | Correct rows flagged |",
        "|---|---|---|---|---|---|---|",
    ]
    for lib_id, block in summary["libraries"].items():
        if block.get("status") != "evaluated":
            lines.append(f"| {lib_id} | - | - | - | - | - | {block['status']} |")
            continue
        flags = block["flags"]
        share = flags["share_controls_with_flag"]
        lines.append(
            f"| {lib_id} | {block['n_rows']} | {block['coverage']['rows_scored']} "
            f"({block['coverage']['share_rows_scored']:.0%}) | "
            f"{flags['controls_with_flag']}"
            f"{'' if share is None else f' ({share:.1%})'} | "
            f"{flags['rows_flagged']} | {flags['rank1_flagged_controls']} | "
            f"{block['false_kills']['labeled_correct_rows_flagged']} |"
        )
    gate = summary["gate"]
    lines += [
        "",
        "## Gate",
        "",
        f"Fire on >= {gate['min_share_controls_with_flag']:.0%} of scoreable controls on some "
        f"library, and flag <= {gate['max_labeled_correct_flagged']} labeled-correct rows. "
        f"**{gate['verdict']}**",
        "",
        gate["reasoning"],
        "",
    ]
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rerank-root", type=Path, default=Path("output/pipeline_equalization"))
    parser.add_argument("--top10-name", default="unified_rerank_all_top10.csv")
    parser.add_argument("--libs", default="209,400")
    parser.add_argument(
        "--calibration",
        type=Path,
        default=Path("docs/paper/data/retention/ri_calibration.json"),
    )
    parser.add_argument("--out-dir", type=Path, default=Path("docs/paper/data/retention"))
    parser.add_argument("--implausible-z", type=float, default=RETENTION_IMPLAUSIBLE_Z)
    parser.add_argument("--max-examples", type=int, default=15)
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
            rerank_root=args.rerank_root,
            libs=[lib.strip() for lib in str(args.libs).split(",") if lib.strip()],
            top10_name=args.top10_name,
            calibration=load_retention_calibration(args.calibration),
            proxy_fn=ri_proxy_from_smiles,
            implausible_z=float(args.implausible_z),
            max_examples=int(args.max_examples),
        )

    json_path = args.out_dir / "audit_topk_retention.json"
    json_path.write_text(json.dumps(summary, indent=2) + "\n")
    markdown_path = args.out_dir / "audit_topk_retention.md"
    markdown_path.write_text(render_markdown(summary) + "\n")

    if summary.get("status") == "not_evaluated":
        print(f"[topk-retention] not evaluated: {summary['reason']}")
    else:
        for lib_id, block in summary["libraries"].items():
            if block.get("status") != "evaluated":
                print(f"[topk-retention] lib {lib_id}: {block['status']}")
                continue
            flags = block["flags"]
            print(
                f"[topk-retention] lib {lib_id}: rows={block['n_rows']} "
                f"scored={block['coverage']['rows_scored']} "
                f"flagged_rows={flags['rows_flagged']} "
                f"controls_with_flag={flags['controls_with_flag']} "
                f"correct_flagged={block['false_kills']['labeled_correct_rows_flagged']}"
            )
        print(f"[topk-retention] gate -> {summary['gate']['verdict']}")
    print(f"[topk-retention] wrote {json_path} and {markdown_path}")


if __name__ == "__main__":
    main()
