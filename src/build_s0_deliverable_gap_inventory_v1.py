#!/usr/bin/env python3
"""Inventory S0 controls that are absent from the current deliverable."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

import pandas as pd

from build_v14_all_library_prediction_review_v1 import clean_text
from build_v19b_high_level_accuracy_summary_v1 import LIB_CONFIG
from phase0_utils import read_csv_clean


DEFAULT_DELIVERABLE_CSV = Path("models/v20_all_library_updated_summary_v1/v20_updated_summary.csv")
DEFAULT_OUTPUT_ROOT = Path("models/v20_all_library_updated_summary_v1")
DEFAULT_KNOWN_X = Path("known_X.csv")
DEFAULT_CLASS_PANEL_ROOT = Path("models/v18_predicted_class_panel_v1")

EXPLICIT_BUCKETS = {
    (209, "X-26218"): "D feature_missing",
    (209, "X-26517"): "D feature_missing",
    (400, "X-14536"): "C DNI",
    (400, "X-14540"): "C DNI",
    (402, "X-14536"): "C DNI",
    (402, "X-14540"): "C DNI",
}

EXPLICIT_PREDICTION_STATUS = {
    (209, "X-26218"): "feature_missing",
    (209, "X-26517"): "feature_missing",
    (400, "X-14536"): "dni_proprietary",
    (400, "X-14540"): "dni_proprietary",
    (402, "X-14536"): "dni_proprietary",
    (402, "X-14540"): "dni_proprietary",
}

RECOMMENDED_ACTION = {
    "A specialist": "rebuild pool + rerank + backfill",
    "B audit-only": "merge audit specialist into pool + rerank + backfill",
    "C DNI": "attempt generation; document true-zero if still empty",
    "D feature_missing": "requires MS2 feature recovery before generation",
    "E no features": "document true-zero",
}


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--deliverable-csv", type=Path, default=DEFAULT_DELIVERABLE_CSV)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--known-x", type=Path, default=DEFAULT_KNOWN_X)
    parser.add_argument("--class-panel-root", type=Path, default=DEFAULT_CLASS_PANEL_ROOT)
    return parser.parse_args(argv)


def load_deliverable_keys(path: Path) -> set[tuple[int, str]]:
    df = pd.read_csv(path, low_memory=False)
    return {
        (int(lib_id), clean_text(control_id))
        for lib_id, control_id in zip(
            pd.to_numeric(df["library_id"], errors="coerce"),
            df["x_compound_id"].map(clean_text),
        )
        if pd.notna(lib_id) and clean_text(control_id)
    }


def load_known_lookup(path: Path) -> pd.DataFrame:
    known = read_csv_clean(path, dtype=str)
    out = known.rename(columns={"LIB_ID": "library_id", "Unknown": "control_id"}).copy()
    out["library_id"] = pd.to_numeric(out["library_id"], errors="coerce")
    out["control_id"] = out["control_id"].map(clean_text)
    out["known_level"] = out.get("Level", pd.Series([""] * len(out))).map(clean_text)
    out["known_sublevel"] = out.get("SubLevel", pd.Series([""] * len(out))).map(clean_text)
    return out[["library_id", "control_id", "known_level", "known_sublevel"]].drop_duplicates(
        subset=["library_id", "control_id"],
        keep="first",
    )


def load_prediction_status(scope: str, class_panel_root: Path) -> dict[str, str]:
    path = class_panel_root / scope / "per_control_class.csv"
    if not path.exists():
        return {}
    frame = pd.read_csv(path, low_memory=False)
    if "control_id" not in frame.columns:
        return {}
    status_col = "prediction_status" if "prediction_status" in frame.columns else ""
    if not status_col:
        return {}
    return {
        clean_text(control_id): clean_text(status)
        for control_id, status in zip(frame["control_id"], frame[status_col])
        if clean_text(control_id)
    }


def infer_bucket(*, library_id: int, control_id: str, n_candidates: int, source_label: str) -> str:
    explicit = EXPLICIT_BUCKETS.get((library_id, control_id))
    if explicit:
        return explicit
    if n_candidates > 0 and "nucleoside" in source_label.lower():
        return "A specialist"
    if n_candidates > 0 and "biotransformation" in source_label.lower():
        return "B audit-only"
    return "E no features"


def count_control_rows(path: Path, control_id: str) -> int:
    if not path.exists():
        return 0
    try:
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames is None or "control_id" not in reader.fieldnames:
                return 0
            return sum(1 for row in reader if clean_text(row.get("control_id")) == control_id)
    except OSError:
        return 0


def artifact_counts(row: pd.Series, control_id: str) -> dict[str, Any]:
    candidate_metrics_text = clean_text(row.get("candidate_metrics_csv"))
    denovo_text = clean_text(row.get("denovo_csv"))
    candidate_metrics_csv = Path(candidate_metrics_text) if candidate_metrics_text else None
    denovo_csv = Path(denovo_text) if denovo_text else None
    metrics_rows = count_control_rows(candidate_metrics_csv, control_id) if candidate_metrics_csv is not None else 0
    denovo_rows = count_control_rows(denovo_csv, control_id) if denovo_csv is not None else 0
    return {
        "audit_candidate_metrics_csv": candidate_metrics_text,
        "audit_denovo_csv": denovo_text,
        "audit_claimed_n_candidates": int(pd.to_numeric(row.get("n_candidates"), errors="coerce") or 0),
        "materialized_candidate_rows": metrics_rows,
        "materialized_denovo_rows": denovo_rows,
        "artifact_mismatch": bool(
            int(pd.to_numeric(row.get("n_candidates"), errors="coerce") or 0) > 0
            and metrics_rows == 0
            and denovo_rows == 0
        ),
    }


def build_inventory(args: argparse.Namespace) -> tuple[pd.DataFrame, dict[str, Any]]:
    deliverable_keys = load_deliverable_keys(args.deliverable_csv)
    known_lookup = load_known_lookup(args.known_x)

    records: list[dict[str, Any]] = []
    for lib_id, _ion_mode, scope, root in LIB_CONFIG:
        audit = pd.read_csv(root / "integrated_per_control_best.csv", low_memory=False)
        prediction_status_lookup = load_prediction_status(scope, args.class_panel_root)
        for _, row in audit.iterrows():
            control_id = clean_text(row.get("control_id"))
            if not control_id:
                continue
            key = (int(lib_id), control_id)
            in_deliverable = key in deliverable_keys
            if in_deliverable:
                continue
            n_candidates = int(pd.to_numeric(row.get("n_candidates"), errors="coerce") or 0)
            source_label = clean_text(row.get("source_label") or row.get("policy"))
            bucket = infer_bucket(
                library_id=int(lib_id),
                control_id=control_id,
                n_candidates=n_candidates,
                source_label=source_label,
            )
            prediction_status = (
                EXPLICIT_PREDICTION_STATUS.get(key)
                or prediction_status_lookup.get(control_id, "")
                or ("feature_missing" if bucket == "D feature_missing" else "")
                or ("no_feature_rows" if bucket in {"C DNI", "E no features"} else "")
            )
            actionable = bucket in {"A specialist", "B audit-only"}
            records.append(
                {
                    "library_id": int(lib_id),
                    "control_id": control_id,
                    "in_deliverable": False,
                    "bucket": bucket,
                    "n_candidates_integrated": n_candidates,
                    "rank1_tanimoto": pd.to_numeric(row.get("rank1_best_tanimoto"), errors="coerce"),
                    "best_source_label": source_label,
                    "prediction_status": prediction_status,
                    "actionable": actionable,
                    "recommended_action": RECOMMENDED_ACTION[bucket],
                    **artifact_counts(row, control_id),
                }
            )

    out = pd.DataFrame(records)
    out = out.merge(known_lookup, on=["library_id", "control_id"], how="left")
    out = out[
        [
            "library_id",
            "control_id",
            "in_deliverable",
            "bucket",
            "n_candidates_integrated",
            "rank1_tanimoto",
            "best_source_label",
            "known_level",
            "known_sublevel",
            "prediction_status",
            "actionable",
            "recommended_action",
            "audit_claimed_n_candidates",
            "materialized_candidate_rows",
            "materialized_denovo_rows",
            "artifact_mismatch",
            "audit_candidate_metrics_csv",
            "audit_denovo_csv",
        ]
    ].sort_values(["library_id", "control_id"], kind="mergesort").reset_index(drop=True)
    summary = {
        "n_missing_pairs": int(len(out)),
        "n_actionable": int(out["actionable"].fillna(False).sum()),
        "n_artifact_mismatch": int(out["artifact_mismatch"].fillna(False).sum()) if "artifact_mismatch" in out.columns else 0,
        "buckets": out["bucket"].value_counts().sort_index().to_dict(),
        "deliverable_csv": str(args.deliverable_csv),
    }
    return out, summary


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    args.output_root.mkdir(parents=True, exist_ok=True)
    inventory, summary = build_inventory(args)
    out_csv = args.output_root / "s0_deliverable_gap_inventory_v1.csv"
    inventory.to_csv(out_csv, index=False)
    (args.output_root / "manifest.json").write_text(
        json.dumps({"output_csv": str(out_csv), **summary}, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
