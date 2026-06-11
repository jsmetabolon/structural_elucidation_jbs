#!/usr/bin/env python3
"""Build a v20 deliverable by backfilling missing v20 recovery rows."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd

from build_s0_deliverable_gap_inventory_v1 import build_inventory
from build_v14_all_library_prediction_review_v1 import (
    DEFAULT_KNOWN_STRUCTURES,
    DEFAULT_KNOWN_X,
    DEFAULT_P0_THRESHOLD,
    DEFAULT_REFERENCE_DB,
    add_metadata_and_flags,
    add_source_columns,
    clean_text,
    enrich_pathway_annotations,
    load_known_metadata,
    load_pathway_assignment_map,
    load_reference_smiles,
)
from build_v19b_high_level_accuracy_summary_v1 import LIB_CONFIG
from build_v19b_updated_summary_v1 import (
    add_provisional_confidence,
    apply_amino_acid_pathway_overrides,
    ensure_full_columns,
    to_deliverable_format,
)
from phase0_utils import read_csv_clean


DEFAULT_BASE_FULL_CSV = Path(
    "models/v19b_all_library_updated_summary_v1/v19b_all_libraries_top10_predictions_with_p0_flags.csv"
)
DEFAULT_BASE_DELIVERABLE_CSV = Path("models/v19b_all_library_updated_summary_v1/v19b_updated_summary.csv")
DEFAULT_OUTPUT_ROOT = Path("models/v20_all_library_updated_summary_v1")
DEFAULT_BACKFILL_TOP10 = {
    "lib_209_negative": Path("models/rerank_v20_lib209_negative/unified_rerank_all_top10.csv"),
    "lib_400_positive": Path("models/rerank_v20_lib400_positive/unified_rerank_all_top10.csv"),
    "lib_402_positive": Path("models/rerank_v20_lib402_positive/unified_rerank_all_top10.csv"),
}


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-full-csv", type=Path, default=DEFAULT_BASE_FULL_CSV)
    parser.add_argument("--base-deliverable-csv", type=Path, default=DEFAULT_BASE_DELIVERABLE_CSV)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--known-x", type=Path, default=DEFAULT_KNOWN_X)
    parser.add_argument("--known-structures", type=Path, default=DEFAULT_KNOWN_STRUCTURES)
    parser.add_argument("--reference-db", type=Path, default=DEFAULT_REFERENCE_DB)
    parser.add_argument("--p0-threshold", type=float, default=DEFAULT_P0_THRESHOLD)
    parser.add_argument(
        "--backfill-top10",
        action="append",
        default=[],
        help="Override backfill top10 path as scope=/path/to/unified_rerank_top10.csv",
    )
    return parser.parse_args(argv)


def parse_backfill_map(values: list[str]) -> dict[str, Path]:
    mapping = {scope: path for scope, path in DEFAULT_BACKFILL_TOP10.items()}
    for value in values:
        scope, sep, raw_path = value.partition("=")
        if not sep or not clean_text(scope) or not clean_text(raw_path):
            raise ValueError(f"Invalid --backfill-top10 value: {value!r}")
        mapping[clean_text(scope)] = Path(clean_text(raw_path))
    return mapping


def scope_details(scope: str) -> tuple[int, str]:
    for lib_id, ion_mode, known_scope, _root in LIB_CONFIG:
        if known_scope == scope:
            return int(lib_id), ion_mode
    raise KeyError(scope)


def load_s0_keys() -> pd.DataFrame:
    records: list[dict[str, Any]] = []
    for lib_id, ion_mode, scope, root in LIB_CONFIG:
        audit = pd.read_csv(root / "integrated_per_control_best.csv", low_memory=False)
        cols = [c for c in ("control_id", "n_candidates", "rank1_best_tanimoto", "source_label", "policy") if c in audit.columns]
        scoped = audit[cols].copy()
        scoped["library_id"] = int(lib_id)
        scoped["ion_mode"] = ion_mode
        scoped["scope"] = scope
        records.append(scoped)
    out = pd.concat(records, ignore_index=True, sort=False)
    out["control_id"] = out["control_id"].map(clean_text)
    out["n_candidates"] = pd.to_numeric(out.get("n_candidates"), errors="coerce").fillna(0).astype(int)
    out["source_label"] = out.get("source_label", out.get("policy", pd.Series([""] * len(out)))).map(clean_text)
    return out[["scope", "library_id", "ion_mode", "control_id", "n_candidates", "rank1_best_tanimoto", "source_label"]]


def source_label_summary(frame: pd.DataFrame) -> str:
    for col in ("source_labels", "source_label", "candidate_source", "policy"):
        if col not in frame.columns:
            continue
        values = sorted({clean_text(value) for value in frame[col].dropna() if clean_text(value)})
        if values:
            return "|".join(values[:8])
    return ""


def load_backfill_keys(backfill_map: dict[str, Path]) -> pd.DataFrame:
    records: list[dict[str, Any]] = []
    seen: set[tuple[str, str]] = set()

    s0_keys = load_s0_keys()
    for row in s0_keys.itertuples(index=False):
        scope = clean_text(row.scope)
        control_id = clean_text(row.control_id)
        if not scope or not control_id:
            continue
        seen.add((scope, control_id))
        records.append(dict(row._asdict()))

    for scope, path in backfill_map.items():
        scope = clean_text(scope)
        if not scope or not path.exists():
            continue
        try:
            raw = read_csv_clean(path, low_memory=False)
        except pd.errors.EmptyDataError:
            continue
        if raw.empty or "control_id" not in raw.columns:
            continue
        lib_id, ion_mode = scope_details(scope)
        raw = raw.copy()
        raw["control_id"] = raw["control_id"].map(clean_text)
        for control_id, group in raw.loc[raw["control_id"].ne("")].groupby("control_id", sort=True):
            key = (scope, clean_text(control_id))
            if key in seen:
                continue
            seen.add(key)
            records.append(
                {
                    "scope": scope,
                    "library_id": int(lib_id),
                    "ion_mode": ion_mode,
                    "control_id": clean_text(control_id),
                    "n_candidates": int(len(group)),
                    "rank1_best_tanimoto": None,
                    "source_label": source_label_summary(group),
                }
            )

    if not records:
        return pd.DataFrame(
            columns=["scope", "library_id", "ion_mode", "control_id", "n_candidates", "rank1_best_tanimoto", "source_label"]
        )
    out = pd.DataFrame(records)
    return out.drop_duplicates(["scope", "control_id"], keep="first").reset_index(drop=True)


def load_backfill_rows(
    *,
    scope: str,
    control_id: str,
    path: Path,
    known_meta: pd.DataFrame,
    reference_smiles: dict[str, list[str]],
    p0_threshold: float,
    reference_db: Path,
) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    raw = read_csv_clean(path, low_memory=False)
    if "control_id" not in raw.columns:
        return pd.DataFrame()
    raw = raw.loc[raw["control_id"].map(clean_text).eq(control_id)].copy()
    if raw.empty:
        return raw
    raw = add_source_columns(
        raw,
        scope=scope,
        row_source_workflow="v20_gap_recovery_backfill",
        row_source_path=path,
        kpi_source_workflow="v20_gap_recovery_backfill",
        kpi_source_path=path,
    )
    if "predicted_pathway" not in raw.columns or "pathway_confidence_tier" not in raw.columns:
        raw = enrich_pathway_annotations(
            raw,
            load_pathway_assignment_map(scope=scope),
            reference_db=reference_db,
        )
    full = add_metadata_and_flags(
        raw,
        known_meta=known_meta,
        reference_smiles=reference_smiles,
        p0_threshold=p0_threshold,
    )
    full["row_source_workflow"] = "v20_gap_recovery_backfill"
    full["kpi_source_workflow"] = "v20_gap_recovery_backfill"
    full["metric_basis"] = "top10_rows"
    return full


def append_backfill_rows(args: argparse.Namespace) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    base_full = pd.read_csv(args.base_full_csv, low_memory=False)
    known_meta = load_known_metadata(args.known_x)
    reference_smiles = load_reference_smiles(args.known_structures)
    backfill_map = parse_backfill_map(args.backfill_top10)
    backfill_keys = load_backfill_keys(backfill_map)

    present = {
        (int(lib_id), clean_text(control_id))
        for lib_id, control_id in zip(
            pd.to_numeric(base_full["library_id"], errors="coerce"),
            base_full["control_id"].map(clean_text),
        )
        if pd.notna(lib_id) and clean_text(control_id)
    }

    backfilled: list[pd.DataFrame] = []
    changelog: list[dict[str, Any]] = []
    for row in backfill_keys.itertuples(index=False):
        key = (int(row.library_id), clean_text(row.control_id))
        if key in present:
            continue
        path = backfill_map.get(str(row.scope))
        if path is None:
            continue
        loaded = load_backfill_rows(
            scope=str(row.scope),
            control_id=clean_text(row.control_id),
            path=path,
            known_meta=known_meta,
            reference_smiles=reference_smiles,
            p0_threshold=float(args.p0_threshold),
            reference_db=args.reference_db,
        )
        if loaded.empty:
            continue
        backfilled.append(loaded)
        changelog.append(
            {
                "library_id": int(row.library_id),
                "scope": str(row.scope),
                "control_id": clean_text(row.control_id),
                "n_backfill_rows": int(len(loaded)),
                "backfill_top10_csv": str(path),
                "best_source_label": clean_text(row.source_label),
                "n_candidates_integrated": int(row.n_candidates),
            }
        )

    full = base_full.copy()
    if backfilled:
        full = pd.concat([base_full, *backfilled], ignore_index=True, sort=False)
        full = full.sort_values(["library_id", "control_id", "rank", "predicted_smiles"], kind="mergesort").reset_index(drop=True)
    full = ensure_full_columns(full)
    full = apply_amino_acid_pathway_overrides(full)
    full = add_provisional_confidence(full)
    deliverable = to_deliverable_format(full)
    return full, deliverable, pd.DataFrame(changelog)


def build_true_zero_report(
    *,
    inventory: pd.DataFrame,
    deliverable: pd.DataFrame,
) -> pd.DataFrame:
    present = {
        (int(lib_id), clean_text(control_id))
        for lib_id, control_id in zip(
            pd.to_numeric(deliverable["library_id"], errors="coerce"),
            deliverable["x_compound_id"].map(clean_text),
        )
        if pd.notna(lib_id) and clean_text(control_id)
    }
    zeros = inventory.loc[
        inventory.apply(
            lambda row: (int(row["library_id"]), clean_text(row["control_id"])) not in present,
            axis=1,
        )
    ].copy()
    reason_map = {
        "A specialist": "Specialist candidate source was identified, but no recovered rerank top-10 rows were available to backfill.",
        "B audit-only": "Integrated audit shows candidate coverage, but no materialized candidate-level rows were available to rebuild and rerank this control.",
        "C DNI": "No candidates after attempted generation; proprietary/DNI control remains undocumented in pool outputs.",
        "D feature_missing": "No recoverable MS2 feature rows; generation blocked upstream.",
        "E no features": "No candidate rows across integrated tracks after pool merge.",
    }
    zeros["reason"] = zeros["bucket"].map(reason_map).fillna("No candidate rows across integrated tracks.")
    return zeros[
        [
            "library_id",
            "control_id",
            "bucket",
            "known_level",
            "known_sublevel",
            "prediction_status",
            "reason",
            "n_candidates_integrated",
            "actionable",
            "recommended_action",
        ]
    ].rename(columns={"n_candidates_integrated": "n_candidates_all_tracks"}).sort_values(
        ["library_id", "control_id"],
        kind="mergesort",
    ).reset_index(drop=True)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    args.output_root.mkdir(parents=True, exist_ok=True)
    full, deliverable, changelog = append_backfill_rows(args)

    inventory, inventory_summary = build_inventory(
        argparse.Namespace(
            deliverable_csv=args.base_deliverable_csv,
            output_root=args.output_root,
            known_x=args.known_x,
            class_panel_root=Path("models/v18_predicted_class_panel_v1"),
        )
    )
    true_zero = build_true_zero_report(inventory=inventory, deliverable=deliverable)

    full_path = args.output_root / "v20_all_libraries_top10_predictions_with_p0_flags.csv"
    deliverable_path = args.output_root / "v20_updated_summary.csv"
    changelog_path = args.output_root / "v20_backfill_changelog.csv"
    zero_path = args.output_root / "s0_true_zero_predictions_v1.csv"
    inventory_path = args.output_root / "s0_deliverable_gap_inventory_v1.csv"

    full.to_csv(full_path, index=False)
    deliverable.to_csv(deliverable_path, index=False)
    changelog.to_csv(changelog_path, index=False)
    true_zero.to_csv(zero_path, index=False)
    inventory.to_csv(inventory_path, index=False)

    manifest = {
        "pipeline_version": "v20",
        "base_full_csv": str(args.base_full_csv),
        "base_deliverable_csv": str(args.base_deliverable_csv),
        "output_csv": str(deliverable_path),
        "full_csv": str(full_path),
        "backfill_changelog_csv": str(changelog_path),
        "true_zero_csv": str(zero_path),
        "gap_inventory_csv": str(inventory_path),
        "n_rows": int(len(deliverable)),
        "n_backfilled_controls": int(changelog["control_id"].nunique()) if not changelog.empty else 0,
        "n_backfilled_rows": int(changelog["n_backfill_rows"].sum()) if not changelog.empty else 0,
        "n_true_zero_pairs": int(len(true_zero)),
        "inventory_summary": inventory_summary,
    }
    (args.output_root / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
