#!/usr/bin/env python3
"""Build V19b-updated all-library summary CSV from V15 top-10 base with S0 rescue overlay."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd

from build_v14_all_library_prediction_review_v1 import (
    DEFAULT_KNOWN_STRUCTURES,
    DEFAULT_P0_THRESHOLD,
    DEFAULT_REFERENCE_DB,
    OUTPUT_COLUMNS,
    best_tanimoto_for_smiles,
    clean_text,
    enrich_pathway_annotations,
    load_pathway_assignment_map,
    load_reference_smiles,
)
from build_x12798_structure_resolution_report_v1 import (
    canonicalize_smiles,
    computed_exact_mass,
    computed_formula,
)
from assign_predicted_pathways import infer_pathway_override
from phase0_utils import canonicalize_smiles as x209_canonicalize_smiles

DEFAULT_V15_REVIEW = Path(
    "models/v15_all_library_prediction_review_v1/all_libraries_top10_predictions_with_p0_flags.csv"
)
DEFAULT_OUTPUT_ROOT = Path("models/v19b_all_library_updated_summary_v1")
DEFAULT_INTEGRATED_V19B_ROOT = Path("models/integrated_denovo_audit_v19b")

V19B_RESCUE_LIB_DIRS = (
    ("lib_305_negative", 305),
    ("lib_400_positive", 400),
    ("lib_402_positive", 402),
)

DELIVERABLE_COLUMNS = [
    "library_id",
    "ion_mode",
    "comp_id",
    "x_compound_id",
    "clean_external_info",
    "known_level",
    "known_sublevel",
    "has_control_reference",
    "rank",
    "control_similarity",
    "control_match",
    "control_match_status",
    "control_match_rank",
    "predicted_smiles",
    "predicted_formula",
    "predicted_exact_mass",
    "abs_mass_error_da",
    "primary_abs_mass_error_da",
    "mass_error_hypothesis",
    "true_pathway",
    "true_subpathway",
    "predicted_pathway",
    "predicted_subpathway",
    "pathway_match",
    "subpathway_match",
    "pathway_confidence_tier",
    "mass_match_score",
    "unified_score",
    "provisional_confidence",
    "candidate_name",
]

PATHWAY_TIER_WEIGHT = {"high": 1.0, "medium": 0.65, "low": 0.35}

DELIVERABLE_RENAME = {
    "has_p0_reference": "has_control_reference",
    "p0_similarity": "control_similarity",
    "p0_match_ge_0p7": "control_match",
    "p0_match_status": "control_match_status",
    "p0_match_rank": "control_match_rank",
    "unified_final_score": "unified_score",
}


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--v15-review-csv", type=Path, default=DEFAULT_V15_REVIEW)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--integrated-v19b-root", type=Path, default=DEFAULT_INTEGRATED_V19B_ROOT)
    parser.add_argument("--known-structures", type=Path, default=DEFAULT_KNOWN_STRUCTURES)
    parser.add_argument("--reference-db", type=Path, default=DEFAULT_REFERENCE_DB)
    parser.add_argument("--p0-threshold", type=float, default=DEFAULT_P0_THRESHOLD)
    return parser.parse_args(argv)


def _canonical(smiles: Any) -> str:
    text = clean_text(smiles)
    if not text:
        return ""
    return canonicalize_smiles(text) or x209_canonicalize_smiles(text) or ""


def load_v19b_rescue_table(integrated_root: Path) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for scope, lib_id in V19B_RESCUE_LIB_DIRS:
        path = integrated_root / scope / "integrated_per_control_best.csv"
        if not path.exists():
            raise FileNotFoundError(f"Missing V19b integrated per-control file: {path}")
        df = pd.read_csv(path, low_memory=False)
        df["scope"] = scope
        df["library_id"] = int(lib_id)
        frames.append(df)
    out = pd.concat(frames, ignore_index=True, sort=False)
    out["control_id"] = out["control_id"].map(clean_text)
    out["rank1_smiles"] = out["rank1_smiles"].map(clean_text)
    out["rank1_canonical_smiles"] = out["rank1_smiles"].map(_canonical)
    out["rank1_best_tanimoto"] = pd.to_numeric(out.get("rank1_best_tanimoto"), errors="coerce")
    out["v19b_source_label"] = out.get("source_label", pd.Series([""] * len(out))).map(clean_text)
    out["v19b_source_family"] = out.get("v19b_source_family", pd.Series([""] * len(out))).map(clean_text)
    return out.drop_duplicates(subset=["library_id", "control_id"], keep="first").reset_index(drop=True)


def recompute_structure_fields(smiles: str, target_neutral_mass: float | None) -> dict[str, Any]:
    formula = computed_formula(smiles)
    exact_mass = computed_exact_mass(smiles)
    abs_error = pd.NA
    if exact_mass is not None and target_neutral_mass is not None and not pd.isna(target_neutral_mass):
        abs_error = abs(float(exact_mass) - float(target_neutral_mass))
    return {
        "predicted_smiles": smiles,
        "predicted_formula": formula,
        "predicted_exact_mass": exact_mass,
        "abs_mass_error_da": abs_error,
        "primary_abs_mass_error_da": abs_error,
    }


def update_p0_flags(
    row: pd.Series,
    *,
    similarity: float,
    reference_smiles: dict[str, list[str]],
    p0_threshold: float,
) -> dict[str, Any]:
    control_id = clean_text(row.get("control_id"))
    has_ref = bool(row.get("has_p0_reference")) or control_id in reference_smiles
    sim = float(similarity)
    match_ge = has_ref and sim >= float(p0_threshold)
    status = "below_threshold"
    if match_ge:
        status = "matched_p0"
    elif not has_ref:
        status = "no_p0_reference"
    elif bool(row.get("is_true_unknown")):
        status = "true_unknown_no_p0_reference"
    match_rank = pd.NA
    if match_ge:
        match_rank = pd.to_numeric(row.get("rank"), errors="coerce")
    return {
        "has_p0_reference": has_ref,
        "p0_similarity": sim,
        "truth_tanimoto": sim,
        "p0_match_ge_0p7": match_ge,
        "p0_match_status": status,
        "p0_match_rank": match_rank,
        "label_ge_0p7": match_ge,
    }


def apply_v19b_overlay(
    base: pd.DataFrame,
    rescue: pd.DataFrame,
    *,
    reference_smiles: dict[str, list[str]],
    p0_threshold: float,
    reference_db: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    out = base.copy()
    out["rank"] = pd.to_numeric(out["rank"], errors="coerce")
    rescue_by_key = {
        (clean_text(row["control_id"]), int(row["library_id"])): row for _, row in rescue.iterrows()
    }

    changelog_records: list[dict[str, Any]] = []
    assignment_cache: dict[str, pd.DataFrame] = {}

    rank1_mask = out["rank"].eq(1)
    rank1 = out.loc[rank1_mask].copy()

    for idx, row in rank1.iterrows():
        control_id = clean_text(row.get("control_id"))
        lib_id = pd.to_numeric(row.get("library_id"), errors="coerce")
        if pd.isna(lib_id):
            continue
        rescue_row = rescue_by_key.get((control_id, int(lib_id)))
        if rescue_row is None:
            continue

        v19b_smiles = clean_text(rescue_row.get("rank1_smiles"))
        v19b_canon = clean_text(rescue_row.get("rank1_canonical_smiles")) or _canonical(v19b_smiles)
        v15_canon = _canonical(row.get("predicted_smiles"))
        if not v19b_canon or v19b_canon == v15_canon:
            continue

        target_mass = pd.to_numeric(row.get("target_neutral_mass"), errors="coerce")
        struct_fields = recompute_structure_fields(v19b_smiles, target_mass)
        integrated_tanimoto = pd.to_numeric(rescue_row.get("rank1_best_tanimoto"), errors="coerce")
        if pd.isna(integrated_tanimoto):
            similarity = best_tanimoto_for_smiles(v19b_smiles, reference_smiles.get(control_id, []))
        else:
            similarity = float(integrated_tanimoto)

        p0_fields = update_p0_flags(
            row,
            similarity=similarity,
            reference_smiles=reference_smiles,
            p0_threshold=p0_threshold,
        )

        scope = clean_text(rescue_row.get("scope")) or clean_text(row.get("scope"))
        if scope not in assignment_cache:
            assignment_cache[scope] = load_pathway_assignment_map(scope=scope)
        pathway_frame = pd.DataFrame([{**row.to_dict(), **struct_fields}])
        pathway_frame = enrich_pathway_annotations(
            pathway_frame,
            assignment_cache[scope],
            reference_db=reference_db,
        )
        pathway_row = pathway_frame.iloc[0]

        true_path = clean_text(row.get("true_pathway"))
        pred_path = clean_text(pathway_row.get("predicted_pathway"))
        true_sub = clean_text(row.get("true_subpathway"))
        pred_sub = clean_text(pathway_row.get("predicted_subpathway"))

        updates = {
            **struct_fields,
            **p0_fields,
            "predicted_pathway": clean_text(pathway_row.get("predicted_pathway")),
            "predicted_subpathway": clean_text(pathway_row.get("predicted_subpathway")),
            "pathway_confidence_tier": clean_text(pathway_row.get("pathway_confidence_tier")),
            "pathway_assignment_status": clean_text(pathway_row.get("pathway_assignment_status")),
            "pathway_vote_share": pathway_row.get("pathway_vote_share"),
            "pathway_best_reference_similarity": pathway_row.get("pathway_best_reference_similarity"),
            "pathway_match": bool(true_path and pred_path and true_path == pred_path),
            "subpathway_match": bool(true_sub and pred_sub and true_sub == pred_sub),
            "source_label": clean_text(rescue_row.get("v19b_source_label")) or "v19b_integrated_rescue",
            "source_labels": clean_text(rescue_row.get("v19b_source_label")) or "v19b_integrated_rescue",
            "row_source_workflow": "v19b_integrated_rescue_overlay",
            "kpi_source_workflow": "integrated_denovo_audit_v19b",
            "metric_basis": "v19b_overlay_rank1",
            "unified_final_score": similarity,
            "truth_tanimoto": similarity,
        }
        for column, value in updates.items():
            if column in out.columns:
                out.at[idx, column] = value

        changelog_records.append(
            {
                "control_id": control_id,
                "library_id": row.get("library_id"),
                "scope": scope,
                "rank": 1,
                "v15_rank1_smiles": clean_text(row.get("predicted_smiles")),
                "v19b_rank1_smiles": v19b_smiles,
                "v15_p0_similarity": pd.to_numeric(row.get("p0_similarity"), errors="coerce"),
                "v19b_p0_similarity": similarity,
                "v19b_source_label": clean_text(rescue_row.get("v19b_source_label")),
                "v19b_source_family": clean_text(rescue_row.get("v19b_source_family")),
            }
        )

    changelog = pd.DataFrame(changelog_records)
    return out, changelog


def _pathway_weight(tier: Any) -> float:
    return PATHWAY_TIER_WEIGHT.get(clean_text(tier).lower(), 0.5)


def _spectrum_component(row: pd.Series) -> float:
    values = []
    for column in ("ms2_match_score", "forward_spectrum_score", "nce_pipeline_score"):
        if column in row.index:
            value = pd.to_numeric(row.get(column), errors="coerce")
            if pd.notna(value):
                values.append(float(value))
    return max(values) if values else 0.0


def add_provisional_confidence(df: pd.DataFrame) -> pd.DataFrame:
    """Add row-level provisional_confidence in [0, 1].

    With P0 reference: structural Tanimoto (p0_similarity).
    Without reference: mass + within-compound rerank + pathway + spectrum, rank-discounted.
    """
    out = df.copy()
    out["rank"] = pd.to_numeric(out.get("rank"), errors="coerce")
    has_ref = out.get("has_p0_reference", pd.Series([False] * len(out))).fillna(False).astype(bool)
    p0_sim = pd.to_numeric(out.get("p0_similarity"), errors="coerce").fillna(0.0).clip(0.0, 1.0)
    confidences = pd.Series([0.0] * len(out), index=out.index, dtype=float)

    ref_mask = has_ref
    confidences.loc[ref_mask] = p0_sim.loc[ref_mask]

    non_ref = out.loc[~ref_mask].copy()
    if not non_ref.empty:
        group_cols = ["library_id", "control_id"]
        non_ref["_unified"] = pd.to_numeric(non_ref.get("unified_final_score"), errors="coerce").fillna(0.0)
        non_ref["_mass"] = pd.to_numeric(non_ref.get("mass_match_score"), errors="coerce").fillna(0.0).clip(0.0, 1.0)
        non_ref["_pathway"] = non_ref.get("pathway_confidence_tier", pd.Series([""] * len(non_ref))).map(_pathway_weight)
        non_ref["_spectrum"] = non_ref.apply(_spectrum_component, axis=1).clip(0.0, 1.0)
        non_ref["_rank"] = pd.to_numeric(non_ref["rank"], errors="coerce").fillna(10.0)
        non_ref["_rank_weight"] = 1.0 / non_ref["_rank"].clip(lower=1.0)

        max_unified = non_ref.groupby(group_cols, sort=False)["_unified"].transform("max").clip(lower=1e-12)
        non_ref["_rerank_rel"] = (non_ref["_unified"] / max_unified).clip(0.0, 1.0)

        rank1_map = (
            non_ref.loc[non_ref["_rank"].eq(1), group_cols + ["_unified"]]
            .drop_duplicates()
            .set_index(group_cols)["_unified"]
        )
        rank2_map = (
            non_ref.loc[non_ref["_rank"].eq(2), group_cols + ["_unified"]]
            .drop_duplicates()
            .set_index(group_cols)["_unified"]
        )

        def margin_factor(row: pd.Series) -> float:
            if float(row["_rank"]) != 1.0:
                return 1.0
            key = (row["library_id"], row["control_id"])
            u1 = float(rank1_map.get(key, row["_unified"]))
            u2 = float(rank2_map.get(key, 0.0))
            if u2 <= 0.0:
                return 1.0
            return float(min(1.0, max(0.0, (u1 - u2) / max(u1, 1e-12))))

        non_ref["_margin"] = non_ref.apply(margin_factor, axis=1)
        non_ref.loc[non_ref["_rank"].eq(1), "_rerank_rel"] *= non_ref.loc[non_ref["_rank"].eq(1), "_margin"]

        non_ref["_confidence"] = (
            0.35 * non_ref["_mass"]
            + 0.25 * non_ref["_rerank_rel"]
            + 0.15 * non_ref["_pathway"]
            + 0.15 * non_ref["_spectrum"]
            + 0.10 * non_ref["_rank_weight"]
        ).clip(0.0, 1.0)
        confidences.loc[non_ref.index] = non_ref["_confidence"]

    out["provisional_confidence"] = confidences.round(6)
    return out


def apply_amino_acid_pathway_overrides(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    predicted_path = out.get("predicted_pathway", pd.Series([""] * len(out), index=out.index)).map(clean_text)
    predicted_sub = out.get("predicted_subpathway", pd.Series([""] * len(out), index=out.index)).map(clean_text)
    candidate_name = out.get("candidate_name", pd.Series([""] * len(out), index=out.index)).map(clean_text)
    query_text = out.get("clean_external_info", pd.Series([""] * len(out), index=out.index)).map(clean_text)
    smiles = out.get("predicted_smiles", pd.Series([""] * len(out), index=out.index)).map(_canonical)

    override_path: list[str] = []
    override_sub: list[str] = []
    for idx in out.index:
        override = infer_pathway_override(
            smiles=smiles.loc[idx],
            query_name=query_text.loc[idx],
            candidate_name=candidate_name.loc[idx],
            current_pathway=predicted_path.loc[idx],
            current_subpathway=predicted_sub.loc[idx],
        )
        override_path.append(clean_text(override.get("predicted_pathway", predicted_path.loc[idx])))
        override_sub.append(clean_text(override.get("predicted_subpathway", predicted_sub.loc[idx])))

    out["predicted_pathway"] = override_path
    out["predicted_subpathway"] = override_sub
    true_path = out.get("true_pathway", pd.Series([""] * len(out), index=out.index)).map(clean_text)
    true_sub = out.get("true_subpathway", pd.Series([""] * len(out), index=out.index)).map(clean_text)
    out["pathway_match"] = true_path.ne("") & true_path.eq(out["predicted_pathway"].map(clean_text))
    out["subpathway_match"] = true_sub.ne("") & true_sub.eq(out["predicted_subpathway"].map(clean_text))
    return out


def to_deliverable_format(full: pd.DataFrame) -> pd.DataFrame:
    out = full.copy()
    for old, new in DELIVERABLE_RENAME.items():
        if old in out.columns:
            out[new] = out[old]
    present = [column for column in DELIVERABLE_COLUMNS if column in out.columns]
    missing = [column for column in DELIVERABLE_COLUMNS if column not in out.columns]
    for column in missing:
        out[column] = pd.NA
    deliverable = out[DELIVERABLE_COLUMNS].copy()
    for column in ("has_control_reference", "control_match", "pathway_match", "subpathway_match"):
        if column in deliverable.columns:
            deliverable[column] = deliverable[column].fillna(False).astype(bool)
    return deliverable.sort_values(["library_id", "x_compound_id", "rank"], kind="mergesort").reset_index(drop=True)


FULL_OUTPUT_COLUMNS = OUTPUT_COLUMNS + ["provisional_confidence"]


def ensure_full_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for column in FULL_OUTPUT_COLUMNS:
        if column not in out.columns:
            out[column] = pd.NA
    return out[FULL_OUTPUT_COLUMNS].copy()


def build_updated_summary(args: argparse.Namespace) -> dict[str, Any]:
    if not args.v15_review_csv.exists():
        raise FileNotFoundError(f"V15 review CSV not found: {args.v15_review_csv}")

    base = pd.read_csv(args.v15_review_csv, low_memory=False)
    rescue = load_v19b_rescue_table(args.integrated_v19b_root)
    reference_smiles = load_reference_smiles(args.known_structures)

    full, changelog = apply_v19b_overlay(
        base,
        rescue,
        reference_smiles=reference_smiles,
        p0_threshold=float(args.p0_threshold),
        reference_db=args.reference_db,
    )
    full = ensure_full_columns(full)
    full = apply_amino_acid_pathway_overrides(full)
    full = add_provisional_confidence(full)
    deliverable = to_deliverable_format(full)

    args.output_root.mkdir(parents=True, exist_ok=True)
    full_path = args.output_root / "v19b_all_libraries_top10_predictions_with_p0_flags.csv"
    deliverable_path = args.output_root / "v19b_updated_summary.csv"
    changelog_path = args.output_root / "v19b_overlay_changelog.csv"
    manifest_path = args.output_root / "manifest.json"

    full.to_csv(full_path, index=False)
    deliverable.to_csv(deliverable_path, index=False)
    changelog.to_csv(changelog_path, index=False)

    manifest = {
        "pipeline_version": "v19b",
        "base_source": str(args.v15_review_csv),
        "n_rows": int(len(deliverable)),
        "n_compounds": int(deliverable["x_compound_id"].nunique()),
        "n_overlay_changes": int(len(changelog)),
        "deliverable_csv": str(deliverable_path),
        "full_csv": str(full_path),
        "changelog_csv": str(changelog_path),
        "provisional_confidence_column": "provisional_confidence",
        "provisional_confidence_note": (
            "With P0 reference: equals control_similarity. "
            "Without reference: weighted mass/rerank/pathway/spectrum with rank discount."
        ),
        "v19b_rescue_sources": [
            str(args.integrated_v19b_root / scope / "integrated_per_control_best.csv")
            for scope, _ in V19B_RESCUE_LIB_DIRS
        ],
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")

    return {
        "manifest": manifest,
        "deliverable_path": deliverable_path,
        "full_path": full_path,
        "changelog_path": changelog_path,
        "n_overlay_changes": len(changelog),
    }


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    result = build_updated_summary(args)
    print(json.dumps(result["manifest"], indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
