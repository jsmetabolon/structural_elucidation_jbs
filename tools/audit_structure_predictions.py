#!/usr/bin/env python3
"""Publication-grade structural recovery and confidence-calibration audit.

Outputs rank-k recovery at several Tanimoto thresholds, exact-source vs predicted-
source subgroups with Wilson intervals, leave-one-out cross-validated Platt
calibration (Brier score, calibration slope/intercept, reliability bins), and an
ablation that removes the pathway axis.

P0/P1 are evaluated only when their authoritative source files are supplied; the
92-control calibration set is never substituted for those cohorts.

Usage:
    python3 tools/audit_structure_predictions.py \
        --reports-root _audit_pull/reports \
        --out-dir docs/paper/data/structure
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import random
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "structural_elucidation" / "reports"))
import confidence as C  # noqa: E402

THRESHOLDS = (0.50, 0.70, 0.85, 0.999999)
TOP_K = (1, 3, 10)


def _num(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(number) else number


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


def bootstrap_median_ci(values: list[float], seed: int = 20260710, n_boot: int = 2000) -> list[float]:
    if not values:
        return [None, None]
    rng = random.Random(seed)
    medians = []
    for _ in range(n_boot):
        sample = [values[rng.randrange(len(values))] for _ in values]
        medians.append(float(np.median(sample)))
    return [round(float(np.quantile(medians, 0.025)), 4), round(float(np.quantile(medians, 0.975)), 4)]


def fit_logistic(x: np.ndarray, y: np.ndarray, iterations: int = 8000) -> tuple[float, float]:
    """Fit logistic P(y=1)=sigmoid(a*x+b) by stable gradient descent."""

    mu = float(x.mean())
    sd = float(x.std()) + 1e-9
    z = (x - mu) / sd
    a, b, lr = 0.0, math.log((float(y.mean()) + 1e-4) / (1.0 - float(y.mean()) + 1e-4)), 0.1
    for _ in range(iterations):
        linear = np.clip(a * z + b, -40, 40)
        p = 1.0 / (1.0 + np.exp(-linear))
        a -= lr * float(np.mean((p - y) * z))
        b -= lr * float(np.mean(p - y))
    return a / sd, b - a * mu / sd


def predict_logistic(x: float, a: float, b: float) -> float:
    linear = max(-40.0, min(40.0, a * x + b))
    return 1.0 / (1.0 + math.exp(-linear))


def loo_platt(raw: np.ndarray, y: np.ndarray) -> np.ndarray:
    out = np.zeros(len(raw), dtype=float)
    for index in range(len(raw)):
        keep = np.arange(len(raw)) != index
        a, b = fit_logistic(raw[keep], y[keep])
        out[index] = predict_logistic(float(raw[index]), a, b)
    return out


def reliability(probabilities: np.ndarray, y: np.ndarray, bins: int = 5) -> list[dict[str, Any]]:
    order = np.argsort(probabilities)
    groups = np.array_split(order, bins)
    rows = []
    for index, group in enumerate(groups, start=1):
        if len(group) == 0:
            continue
        correct = int(y[group].sum())
        metric = wilson(correct, len(group))
        rows.append(
            {
                "bin": index,
                "n": len(group),
                "mean_predicted": round(float(probabilities[group].mean()), 4),
                "observed": metric["accuracy"],
                "observed_ci95": metric["ci95"],
            }
        )
    return rows


def calibration_slope_intercept(probabilities: np.ndarray, y: np.ndarray) -> dict[str, float]:
    clipped = np.clip(probabilities, 1e-5, 1.0 - 1e-5)
    logits = np.log(clipped / (1.0 - clipped))
    slope, intercept = fit_logistic(logits, y)
    return {"slope": round(slope, 4), "intercept": round(intercept, 4)}


def _round_or_blank(value: float | None) -> str:
    """Render an axis score, leaving unavailable axes empty rather than zero.

    An unavailable axis is not a zero-scoring one, and writing 0 would make a
    missing MS2 spectrum look like a spectrum the candidate failed to explain.
    """

    return "" if value is None else f"{round(float(value), 6)}"


def raw_without_pathway(axes: list[C.Axis]) -> float:
    retained = [axis for axis in axes if axis.key != "pathway"]
    numerator = sum(axis.weight * float(axis.score) for axis in retained if axis.available and axis.score is not None)
    denominator = sum(axis.weight for axis in retained if axis.available and axis.score is not None)
    return numerator / denominator if denominator else 0.0


def ranked_candidates(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    return sorted(rows, key=lambda row: _num(row.get("rank")) or 1e9)


def load_records(root: Path) -> list[dict[str, Any]]:
    records = []
    for report_dir in sorted(root.glob("X-*")):
        if not report_dir.is_dir():
            continue
        slug = report_dir.name.lower()
        data = report_dir / "data"
        fields = C.read_fields_from_sidecars(data, slug)
        candidates = ranked_candidates(data / f"{slug}_top_candidates.csv")
        if not fields or not candidates:
            continue
        similarities = [_num(row.get("truth_tanimoto")) for row in candidates[:10]]
        rank1_truth = similarities[0] if similarities else None
        if rank1_truth is None or rank1_truth <= 0.0:
            continue
        similarities = [value if value is not None else 0.0 for value in similarities]
        axes = C.build_axes(fields)
        raw_full, completeness = C.raw_score(axes)
        by_key = {axis.key: axis for axis in axes}
        records.append(
            {
                "control_id": report_dir.name,
                "lib": str(candidates[0].get("lib_id") or "").removesuffix(".0"),
                "source": str(fields.get("source") or ""),
                "truth_tanimoto": similarities[0],
                "similarities": similarities,
                "raw_full": raw_full,
                "raw_without_pathway": raw_without_pathway(axes),
                "completeness": completeness,
                # Per-axis evidence, so tools/audit_thresholds.py can sweep the MSI
                # boundaries from the committed cohort instead of needing the bundle.
                "axis_structure": by_key["structure"].score,
                "axis_mass": by_key["mass"].score if by_key["mass"].available else None,
                "axis_ms2": by_key["ms2"].score if by_key["ms2"].available else None,
                "axis_pathway": by_key["pathway"].score if by_key["pathway"].available else None,
                "axis_retention": (
                    by_key["retention"].score if by_key["retention"].available else None
                ),
                "observed_mz_available": bool(fields.get("observed_mz_available")),
            }
        )
    return records


def rank_k_table(records: list[dict[str, Any]]) -> dict[str, Any]:
    table = {}
    for threshold in THRESHOLDS:
        threshold_key = "exact" if threshold > 0.99 else f"{threshold:.2f}"
        table[threshold_key] = {}
        for top_k in TOP_K:
            correct = sum(
                max(record["similarities"][:top_k], default=0.0) >= threshold for record in records
            )
            table[threshold_key][f"top{top_k}"] = wilson(correct, len(records))
    return table


def subgroup_table(records: list[dict[str, Any]]) -> dict[str, Any]:
    groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for record in records:
        groups[record["source"] or "other"].append(record)
    out = {}
    for source, group in sorted(groups.items()):
        similarities = [float(record["truth_tanimoto"]) for record in group]
        correct = sum(value >= 0.70 for value in similarities)
        out[source] = {
            "rank1_tanimoto_ge_0p70": wilson(correct, len(group)),
            "median_truth_tanimoto": round(float(np.median(similarities)), 4),
            "median_truth_tanimoto_ci95": bootstrap_median_ci(similarities),
        }
    return out


def calibration_result(records: list[dict[str, Any]], key: str) -> dict[str, Any]:
    raw = np.array([float(record[key]) for record in records])
    y = np.array([float(record["truth_tanimoto"] >= 0.70) for record in records])
    probabilities = loo_platt(raw, y)
    return {
        "method": "leave-one-out Platt scaling",
        "n": len(records),
        "brier_score": round(float(np.mean((probabilities - y) ** 2)), 4),
        "calibration": calibration_slope_intercept(probabilities, y),
        "reliability_bins": reliability(probabilities, y),
        "predictions": {
            record["control_id"]: round(float(probability), 6)
            for record, probability in zip(records, probabilities)
        },
    }


def p0_p1_status(args: argparse.Namespace) -> dict[str, Any]:
    p0 = args.p0_flags
    p1 = args.p1_truth_dir
    missing = [str(path) for path in (p0, p1) if path is None or not path.exists()]
    if missing:
        return {
            "status": "not_evaluated",
            "reason": "Authoritative P0/P1 cohort inputs are not available; the 92-control set is not substituted.",
            "required_inputs": [
                "models/v21_all_library_updated_summary/v20_all_libraries_top10_predictions_with_p0_flags.csv",
                "models/p1_reference_v1/p1_cohort.csv + p1_truth_audit.csv (and verified_plus_2c_resolved_v2)",
            ],
            "missing": missing,
        }
    return {
        "status": "inputs_present_but_not_yet_supported",
        "p0_flags": str(p0),
        "p1_truth_dir": str(p1),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reports-root", type=Path, default=Path("_audit_pull/reports"))
    parser.add_argument("--out-dir", type=Path, default=Path("docs/paper/data/structure"))
    parser.add_argument("--p0-flags", type=Path, default=None)
    parser.add_argument("--p1-truth-dir", type=Path, default=None)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    records = load_records(args.reports_root)
    if not records:
        raise SystemExit("No truth-verified structural records found")
    summary = {
        "truth_definition": "truth_tanimoto > 0 labels a control; correct if Tanimoto >= threshold",
        "n_truth_verified": len(records),
        "rank_k": rank_k_table(records),
        "source_subgroups": subgroup_table(records),
        "calibration_full": calibration_result(records, "raw_full"),
        "calibration_without_pathway": calibration_result(records, "raw_without_pathway"),
        "p0_p1": p0_p1_status(args),
    }
    summary_path = args.out_dir / "audit_structure_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")

    rows_path = args.out_dir / "audit_structure_controls.csv"
    full_predictions = summary["calibration_full"]["predictions"]
    ablated_predictions = summary["calibration_without_pathway"]["predictions"]
    with rows_path.open("w", newline="") as handle:
        fieldnames = [
            "control_id",
            "lib",
            "source",
            "truth_tanimoto",
            "top1_correct_0p70",
            "top3_correct_0p70",
            "top10_correct_0p70",
            "raw_full",
            "raw_without_pathway",
            "loo_probability_full",
            "loo_probability_without_pathway",
            "axis_structure",
            "axis_mass",
            "axis_ms2",
            "axis_pathway",
            "axis_retention",
            "observed_mz_available",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for record in records:
            similarities = record["similarities"]
            writer.writerow(
                {
                    "control_id": record["control_id"],
                    "lib": record["lib"],
                    "source": record["source"],
                    "truth_tanimoto": record["truth_tanimoto"],
                    "top1_correct_0p70": int(similarities[0] >= 0.70),
                    "top3_correct_0p70": int(max(similarities[:3]) >= 0.70),
                    "top10_correct_0p70": int(max(similarities[:10]) >= 0.70),
                    "raw_full": round(record["raw_full"], 6),
                    "raw_without_pathway": round(record["raw_without_pathway"], 6),
                    "loo_probability_full": full_predictions[record["control_id"]],
                    "loo_probability_without_pathway": ablated_predictions[record["control_id"]],
                    "axis_structure": _round_or_blank(record["axis_structure"]),
                    "axis_mass": _round_or_blank(record["axis_mass"]),
                    "axis_ms2": _round_or_blank(record["axis_ms2"]),
                    "axis_pathway": _round_or_blank(record["axis_pathway"]),
                    "axis_retention": _round_or_blank(record["axis_retention"]),
                    "observed_mz_available": int(record["observed_mz_available"]),
                }
            )
    print(f"[structure-audit] {len(records)} truth controls -> {summary_path} and {rows_path}")


if __name__ == "__main__":
    main()
