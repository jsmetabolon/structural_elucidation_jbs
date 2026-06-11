#!/usr/bin/env python3
"""Build an all-library V14 top-10 prediction review workbook."""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path
from typing import Any

import pandas as pd

try:
    from rdkit import Chem, DataStructs
    from rdkit.Chem import AllChem

    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False

try:
    from openpyxl import Workbook
    from openpyxl.chart import BarChart, LineChart, Reference
    from openpyxl.formatting.rule import CellIsRule
    from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
    from openpyxl.utils import get_column_letter

    HAS_OPENPYXL = True
except ImportError:
    HAS_OPENPYXL = False

from phase0_utils import canonicalize_smiles, read_csv_clean


DEFAULT_SCOPES = (
    "lib_209_negative",
    "lib_305_negative",
    "lib_400_positive",
    "lib_402_positive",
)
DEFAULT_INPUT_ROOT = Path("models/v14_ranked_pathway_annotations_v1")
DEFAULT_OUTPUT_ROOT = Path("models/v14_all_library_prediction_review_v1")
DEFAULT_KNOWN_X = Path("known_X.csv")
DEFAULT_KNOWN_STRUCTURES = Path("known_X_all_available_structures_simple.csv")
DEFAULT_REFERENCE_DB = Path("Metabolite Synonym Database/metabolite_synonyms_merged_with_mtx_giga.sqlite")
DEFAULT_P0_THRESHOLD = 0.70
BEST_209_TOP10 = Path("models/rerank_lib209_step7_targeted_v2_auto/unified_rerank_top10.csv")
BEST_209_KPI_SUMMARY = Path(
    "models/integrated_denovo_audit_lib209_step7_targeted_v2_auto/lib_209_negative/integrated_audit_summary.csv"
)

OUTPUT_COLUMNS = [
    "scope",
    "library_id",
    "row_source_workflow",
    "row_source_path",
    "kpi_source_workflow",
    "kpi_source_path",
    "metric_basis",
    "ion_mode",
    "feature_key",
    "control_id",
    "comp_id",
    "x_compound_id",
    "known_compound_name",
    "clean_external_info",
    "known_level",
    "known_sublevel",
    "is_true_unknown",
    "has_p0_reference",
    "true_pathway",
    "true_subpathway",
    "known_mz",
    "parent_mz",
    "target_neutral_mass",
    "retention_index",
    "rank",
    "p0_similarity",
    "p0_match_ge_0p7",
    "p0_match_status",
    "p0_match_rank",
    "predicted_smiles",
    "predicted_formula",
    "predicted_exact_mass",
    "abs_mass_error_da",
    "primary_abs_mass_error_da",
    "mass_error_hypothesis",
    "strict_mass_pass",
    "predicted_pathway",
    "predicted_subpathway",
    "pathway_match",
    "subpathway_match",
    "pathway_confidence_tier",
    "pathway_assignment_status",
    "pathway_vote_share",
    "pathway_best_reference_similarity",
    "unified_final_score",
    "rerank_score",
    "mass_match_score",
    "formula_prior_score",
    "ms2_match_score",
    "forward_spectrum_score",
    "nce_pipeline_score",
    "unified_pool_rank",
    "rerank_policy",
    "source_label",
    "source_labels",
    "candidate_source",
    "candidate_name",
    "truth_tanimoto",
    "label_ge_0p7",
]

NUMERIC_COLUMNS = {
    "library_id",
    "known_mz",
    "parent_mz",
    "target_neutral_mass",
    "retention_index",
    "rank",
    "p0_similarity",
    "p0_match_rank",
    "predicted_exact_mass",
    "abs_mass_error_da",
    "primary_abs_mass_error_da",
    "pathway_vote_share",
    "pathway_best_reference_similarity",
    "unified_final_score",
    "rerank_score",
    "mass_match_score",
    "formula_prior_score",
    "ms2_match_score",
    "forward_spectrum_score",
    "nce_pipeline_score",
    "unified_pool_rank",
    "truth_tanimoto",
    "label_ge_0p7",
}
BOOLEAN_COLUMNS = {
    "is_true_unknown",
    "has_p0_reference",
    "p0_match_ge_0p7",
    "strict_mass_pass",
    "pathway_match",
    "subpathway_match",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build all-library V14 top-10 prediction review CSVs and workbook."
    )
    parser.add_argument("--input-root", type=Path, default=DEFAULT_INPUT_ROOT)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--known-x", type=Path, default=DEFAULT_KNOWN_X)
    parser.add_argument("--known-structures", type=Path, default=DEFAULT_KNOWN_STRUCTURES)
    parser.add_argument("--reference-db", type=Path, default=DEFAULT_REFERENCE_DB)
    parser.add_argument("--scopes", default=",".join(DEFAULT_SCOPES))
    parser.add_argument("--p0-threshold", type=float, default=DEFAULT_P0_THRESHOLD)
    parser.add_argument("--skip-workbook", action="store_true", help="Write CSV sidecars only.")
    return parser.parse_args()


def clean_text(value: object) -> str:
    if pd.isna(value):
        return ""
    text = str(value).strip()
    if text.lower() in {"nan", "none", "null", "<na>"}:
        return ""
    return text


def normalize_lib_id(value: object) -> str:
    text = clean_text(value)
    if not text:
        return ""
    numeric = pd.to_numeric(text, errors="coerce")
    if pd.notna(numeric):
        return str(int(float(numeric)))
    return text


def scope_to_library_id(scope: str) -> str:
    match = re.match(r"^lib_(\d+)_(negative|positive)$", clean_text(scope))
    if not match:
        raise ValueError(f"Unsupported scope format: {scope!r}")
    return match.group(1)


def parse_scope_list(scopes: str | list[str] | tuple[str, ...]) -> list[str]:
    if isinstance(scopes, str):
        values = scopes.split(",")
    else:
        values = list(scopes)
    parsed = [clean_text(scope) for scope in values if clean_text(scope)]
    if not parsed:
        raise ValueError("At least one scope is required.")
    return parsed


def default_source_config(scope: str, input_root: Path = DEFAULT_INPUT_ROOT) -> dict[str, str | Path]:
    if scope == "lib_209_negative":
        supplemental_path = input_root / scope / "top10_with_pathways.csv"
        return {
            "scope": scope,
            "row_source_workflow": "mixed_lib209_best_plus_v14_remaining_controls",
            "row_source_path": f"{BEST_209_TOP10};{supplemental_path}",
            "primary_row_source_workflow": "lib209_step7_targeted_v2_auto",
            "primary_row_source_path": BEST_209_TOP10,
            "supplemental_row_source_workflow": "v14_ranked_pathway_annotations_v1_remaining_209_controls",
            "supplemental_row_source_path": supplemental_path,
            "kpi_source_workflow": "integrated_denovo_audit_lib209_step7_targeted_v2_auto",
            "kpi_source_path": BEST_209_KPI_SUMMARY,
            "metric_basis": "integrated_audit",
        }
    return {
        "scope": scope,
        "row_source_workflow": "v14_ranked_pathway_annotations_v1",
        "row_source_path": input_root / scope / "top10_with_pathways.csv",
        "kpi_source_workflow": "v14_ranked_pathway_annotations_v1",
        "kpi_source_path": input_root / scope / "top10_with_pathways.csv",
        "metric_basis": "top10_rows",
    }


def source_config_for_scopes(scopes: list[str] | tuple[str, ...], input_root: Path = DEFAULT_INPUT_ROOT) -> dict[str, dict[str, str | Path]]:
    return {scope: default_source_config(scope, input_root=input_root) for scope in scopes}


def parse_feature_ri(feature_key: object) -> float | None:
    match = re.search(r"_ri([0-9.]+)", clean_text(feature_key))
    if not match:
        return None
    value = pd.to_numeric(match.group(1), errors="coerce")
    return None if pd.isna(value) else float(value)


def clean_known_name(value: object) -> str:
    text = clean_text(value)
    return re.sub(r"^KNOWN\s+", "", text, flags=re.IGNORECASE).strip()


def canonical_fingerprint(smiles: object):
    if not HAS_RDKIT:
        return None
    canonical = canonicalize_smiles(clean_text(smiles))
    if not canonical:
        return None
    mol = Chem.MolFromSmiles(canonical)
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)


def best_tanimoto_for_smiles(predicted_smiles: object, reference_smiles: list[str]) -> float:
    if not HAS_RDKIT or not reference_smiles:
        return 0.0
    predicted_fp = canonical_fingerprint(predicted_smiles)
    if predicted_fp is None:
        return 0.0
    best = 0.0
    for reference in reference_smiles:
        reference_fp = canonical_fingerprint(reference)
        if reference_fp is None:
            continue
        best = max(best, float(DataStructs.TanimotoSimilarity(predicted_fp, reference_fp)))
    return float(best)


def load_known_metadata(path: Path) -> pd.DataFrame:
    known = read_csv_clean(path, dtype=str)
    required = {"LIB_ID", "Unknown", "KNOWN", "Cleaned External Info", "Level", "SubLevel"}
    missing = sorted(required.difference(known.columns))
    if missing:
        raise ValueError(f"Missing required columns in {path}: {', '.join(missing)}")
    out = known.copy()
    out["library_id"] = out["LIB_ID"].map(normalize_lib_id)
    out["control_id"] = out["Unknown"].map(clean_text)
    out["known_compound_name"] = out["KNOWN"].map(clean_known_name)
    out["clean_external_info"] = out["Cleaned External Info"].map(clean_text)
    out["known_level"] = out["Level"].map(clean_text)
    out["known_sublevel"] = out["SubLevel"].map(clean_text)
    out["true_pathway"] = out.get("SuperPathway", pd.Series([""] * len(out))).map(clean_text)
    out["true_subpathway"] = out.get("Draft Subpathway", pd.Series([""] * len(out))).map(clean_text)
    out["known_mz"] = pd.to_numeric(out.get("m/z", pd.Series([pd.NA] * len(out))), errors="coerce")
    out = out.loc[out["library_id"].ne("") & out["control_id"].ne("")].copy()
    out = out.drop_duplicates(subset=["library_id", "control_id"], keep="first").reset_index(drop=True)
    return out[
        [
            "library_id",
            "control_id",
            "known_compound_name",
            "clean_external_info",
            "known_level",
            "known_sublevel",
            "true_pathway",
            "true_subpathway",
            "known_mz",
        ]
    ].copy()


def load_reference_smiles(path: Path) -> dict[str, list[str]]:
    references = read_csv_clean(path, dtype=str)
    required = {"x_compound_id", "smiles"}
    missing = sorted(required.difference(references.columns))
    if missing:
        raise ValueError(f"Missing required columns in {path}: {', '.join(missing)}")
    out: dict[str, list[str]] = {}
    for control_id, group in references.groupby(references["x_compound_id"].map(clean_text)):
        smiles = sorted(set(group["smiles"].map(clean_text).loc[lambda series: series.ne("")].tolist()))
        if control_id and smiles:
            out[control_id] = smiles
    return out


def first_present(df: pd.DataFrame, columns: list[str], default: object = pd.NA) -> pd.Series:
    out = pd.Series([default] * len(df), index=df.index, dtype=object)
    for column in columns:
        if column not in df.columns:
            continue
        values = df[column]
        if values.dtype == object:
            cleaned = values.map(clean_text)
            mask = out.map(clean_text).eq("") & cleaned.ne("")
        else:
            mask = pd.isna(out) & values.notna()
        out.loc[mask] = values.loc[mask]
    return out


def load_pathway_assignment_map(input_root: Path = DEFAULT_INPUT_ROOT, scope: str = "lib_209_negative") -> pd.DataFrame:
    path = input_root / scope / "unique_smiles_pathway_assignments.csv"
    if not path.exists():
        return pd.DataFrame(columns=["canonical_smiles"])
    assignments = read_csv_clean(path, low_memory=False)
    if "canonical_smiles" not in assignments.columns:
        return pd.DataFrame(columns=["canonical_smiles"])
    keep_columns = [
        "canonical_smiles",
        "predicted_pathway",
        "predicted_subpathway",
        "confidence_tier",
        "assignment_status",
        "vote_share",
        "best_reference_similarity",
    ]
    for column in keep_columns:
        if column not in assignments.columns:
            assignments[column] = pd.NA
    out = assignments[keep_columns].drop_duplicates(subset=["canonical_smiles"], keep="first").copy()
    return out.rename(
        columns={
            "confidence_tier": "pathway_confidence_tier",
            "assignment_status": "pathway_assignment_status",
            "vote_share": "pathway_vote_share",
            "best_reference_similarity": "pathway_best_reference_similarity",
        }
    )


def normalize_assignment_map(assignments: pd.DataFrame | None) -> pd.DataFrame:
    columns = [
        "canonical_smiles",
        "predicted_pathway",
        "predicted_subpathway",
        "pathway_confidence_tier",
        "pathway_assignment_status",
        "pathway_vote_share",
        "pathway_best_reference_similarity",
    ]
    if assignments is None or assignments.empty:
        return pd.DataFrame(columns=columns)
    out = assignments.copy()
    rename_map = {
        "confidence_tier": "pathway_confidence_tier",
        "assignment_status": "pathway_assignment_status",
        "vote_share": "pathway_vote_share",
        "best_reference_similarity": "pathway_best_reference_similarity",
    }
    out = out.rename(columns={old: new for old, new in rename_map.items() if old in out.columns})
    for column in columns:
        if column not in out.columns:
            out[column] = pd.NA
    out["canonical_smiles"] = out["canonical_smiles"].map(clean_text)
    return out.loc[out["canonical_smiles"].ne(""), columns].drop_duplicates(
        subset=["canonical_smiles"], keep="first"
    )


def generate_pathway_assignment_map(canonical_smiles: list[str], reference_db: Path) -> pd.DataFrame:
    unique_smiles = sorted(set(clean_text(smiles) for smiles in canonical_smiles if clean_text(smiles)))
    if not unique_smiles:
        return normalize_assignment_map(None)
    if not reference_db.exists():
        raise FileNotFoundError(reference_db)

    from assign_predicted_pathways import (
        DEFAULT_NEIGHBOR_SIMILARITY_FLOOR,
        DEFAULT_TOP_K_NEIGHBORS,
        DEFAULT_VOTE_SIMILARITY_POWER,
        build_pathway_mw_stats,
    )
    from assign_v14_ranked_pathways_v1 import (
        build_neighbor_table_parallel,
        load_v14_reference_library,
        summarize_unique_smiles_assignments,
    )

    reference_df = load_v14_reference_library(reference_db)
    neighbor_df = build_neighbor_table_parallel(
        unique_smiles,
        reference_df,
        top_k_neighbors=DEFAULT_TOP_K_NEIGHBORS,
        similarity_floor=DEFAULT_NEIGHBOR_SIMILARITY_FLOOR,
        jobs=1,
    )
    assignments, _ = summarize_unique_smiles_assignments(
        unique_smiles,
        neighbor_df,
        vote_similarity_power=DEFAULT_VOTE_SIMILARITY_POWER,
        pathway_mw_stats=build_pathway_mw_stats(reference_df),
    )
    return normalize_assignment_map(assignments)


def enrich_pathway_annotations(
    df: pd.DataFrame,
    assignment_map: pd.DataFrame | None,
    reference_db: Path | None = None,
) -> pd.DataFrame:
    out = df.copy()
    pathway_columns = [
        "predicted_pathway",
        "predicted_subpathway",
        "pathway_confidence_tier",
        "pathway_assignment_status",
        "pathway_vote_share",
        "pathway_best_reference_similarity",
    ]
    for column in pathway_columns:
        if column not in out.columns:
            out[column] = pd.NA

    out["_canonical_pathway_join_smiles"] = out.get("predicted_smiles", pd.Series([""] * len(out))).map(
        lambda value: clean_text(canonicalize_smiles(clean_text(value)))
    )
    assignment_map = normalize_assignment_map(assignment_map)

    rows_need_labels = (
        out["_canonical_pathway_join_smiles"].map(clean_text).ne("")
        & (
            out["predicted_pathway"].map(clean_text).eq("")
            | out["predicted_subpathway"].map(clean_text).eq("")
            | out["pathway_assignment_status"].map(clean_text).eq("")
        )
    )
    missing_canonicals = sorted(
        set(out.loc[rows_need_labels, "_canonical_pathway_join_smiles"].map(clean_text))
        - set(assignment_map["canonical_smiles"].map(clean_text))
    )
    if missing_canonicals and reference_db is not None:
        generated_map = generate_pathway_assignment_map(missing_canonicals, reference_db)
        assignment_map = (
            pd.concat([assignment_map, generated_map], ignore_index=True, sort=False)
            .drop_duplicates(subset=["canonical_smiles"], keep="first")
            .reset_index(drop=True)
        )

    if assignment_map.empty:
        for column in pathway_columns:
            out[column] = out[column].where(out[column].map(clean_text).ne(""), pd.NA)
        out.loc[out["pathway_assignment_status"].map(clean_text).eq(""), "pathway_assignment_status"] = (
            "missing_pathway_annotation"
        )
        return out.drop(columns=["_canonical_pathway_join_smiles"], errors="ignore")

    annotated = out.merge(
        assignment_map,
        left_on="_canonical_pathway_join_smiles",
        right_on="canonical_smiles",
        how="left",
        suffixes=("", "_assignment"),
        sort=False,
    )
    for column in pathway_columns:
        assignment_column = f"{column}_assignment"
        if assignment_column in annotated.columns:
            current = annotated[column].map(clean_text)
            incoming = annotated[assignment_column]
            annotated[column] = annotated[column].where(current.ne(""), incoming)
    annotated.loc[annotated["pathway_assignment_status"].map(clean_text).eq(""), "pathway_assignment_status"] = (
        "missing_pathway_annotation"
    )
    annotated = annotated.drop(
        columns=["_canonical_pathway_join_smiles", "canonical_smiles"]
        + [f"{column}_assignment" for column in pathway_columns if f"{column}_assignment" in annotated.columns],
        errors="ignore",
    )
    return annotated


def add_source_columns(
    df: pd.DataFrame,
    *,
    scope: str,
    row_source_workflow: str,
    row_source_path: Path,
    kpi_source_workflow: str,
    kpi_source_path: Path,
) -> pd.DataFrame:
    out = df.copy()
    out["scope"] = scope
    out["library_id"] = scope_to_library_id(scope)
    out["row_source_workflow"] = row_source_workflow
    out["row_source_path"] = str(row_source_path)
    out["kpi_source_workflow"] = kpi_source_workflow
    out["kpi_source_path"] = str(kpi_source_path)
    out["metric_basis"] = "top10_rows"
    return out


def load_scope_top10(
    scope: str,
    source_config: dict[str, str | Path],
    pathway_assignment_map: pd.DataFrame | None = None,
    reference_db: Path | None = None,
) -> pd.DataFrame:
    if scope == "lib_209_negative":
        primary_path = Path(source_config.get("primary_row_source_path", BEST_209_TOP10))
        supplemental_path = Path(
            source_config.get("supplemental_row_source_path", DEFAULT_INPUT_ROOT / scope / "top10_with_pathways.csv")
        )
        if not primary_path.exists():
            raise FileNotFoundError(primary_path)
        if not supplemental_path.exists():
            raise FileNotFoundError(supplemental_path)
        primary = add_source_columns(
            read_csv_clean(primary_path, low_memory=False),
            scope=scope,
            row_source_workflow=str(source_config.get("primary_row_source_workflow", "lib209_step7_targeted_v2_auto")),
            row_source_path=primary_path,
            kpi_source_workflow=str(source_config["kpi_source_workflow"]),
            kpi_source_path=Path(source_config["kpi_source_path"]),
        )
        primary = enrich_pathway_annotations(primary, pathway_assignment_map, reference_db=reference_db)
        supplemental = add_source_columns(
            read_csv_clean(supplemental_path, low_memory=False),
            scope=scope,
            row_source_workflow=str(
                source_config.get("supplemental_row_source_workflow", "v14_ranked_pathway_annotations_v1")
            ),
            row_source_path=supplemental_path,
            kpi_source_workflow=str(source_config["kpi_source_workflow"]),
            kpi_source_path=Path(source_config["kpi_source_path"]),
        )
        primary_controls = set(primary["control_id"].map(clean_text))
        supplemental = supplemental.loc[~supplemental["control_id"].map(clean_text).isin(primary_controls)].copy()
        return pd.concat([primary, supplemental], ignore_index=True, sort=False)

    path = Path(source_config["row_source_path"])
    if not path.exists():
        raise FileNotFoundError(path)
    df = add_source_columns(
        read_csv_clean(path, low_memory=False),
        scope=scope,
        row_source_workflow=str(source_config["row_source_workflow"]),
        row_source_path=path,
        kpi_source_workflow=str(source_config["kpi_source_workflow"]),
        kpi_source_path=Path(source_config["kpi_source_path"]),
    )
    return df


def add_metadata_and_flags(
    raw: pd.DataFrame,
    known_meta: pd.DataFrame,
    reference_smiles: dict[str, list[str]],
    p0_threshold: float = DEFAULT_P0_THRESHOLD,
) -> pd.DataFrame:
    out = raw.copy()
    out["library_id"] = out["library_id"].map(normalize_lib_id)
    out["control_id"] = out["control_id"].map(clean_text)
    meta = known_meta.copy()
    out = out.merge(meta, on=["library_id", "control_id"], how="left", sort=False).copy()
    out["x_compound_id"] = out["control_id"]
    out["is_true_unknown"] = out["known_level"].map(clean_text).eq("5") | out["known_sublevel"].map(clean_text).eq("5")
    out["has_p0_reference"] = out["control_id"].map(lambda value: clean_text(value) in reference_smiles)
    out["p0_similarity"] = pd.to_numeric(out.get("truth_tanimoto", pd.Series([pd.NA] * len(out))), errors="coerce")

    missing_similarity = out["p0_similarity"].isna() & out["has_p0_reference"]
    if missing_similarity.any():
        recomputed = []
        for row in out.loc[missing_similarity].itertuples(index=False):
            predicted = clean_text(getattr(row, "predicted_smiles", ""))
            control_id = clean_text(getattr(row, "control_id", ""))
            recomputed.append(best_tanimoto_for_smiles(predicted, reference_smiles.get(control_id, [])))
        out.loc[missing_similarity, "p0_similarity"] = recomputed
    out["p0_similarity"] = pd.to_numeric(out["p0_similarity"], errors="coerce").fillna(0.0)
    out["p0_match_ge_0p7"] = out["has_p0_reference"] & out["p0_similarity"].ge(float(p0_threshold))
    out["p0_match_status"] = "below_threshold"
    out.loc[out["p0_match_ge_0p7"], "p0_match_status"] = "matched_p0"
    out.loc[~out["has_p0_reference"], "p0_match_status"] = "no_p0_reference"
    out.loc[out["is_true_unknown"] & ~out["has_p0_reference"], "p0_match_status"] = "true_unknown_no_p0_reference"
    out["p0_match_rank"] = pd.NA
    out.loc[out["p0_match_ge_0p7"], "p0_match_rank"] = pd.to_numeric(out.loc[out["p0_match_ge_0p7"], "rank"], errors="coerce")

    true_path = out["true_pathway"].map(clean_text)
    pred_path = out.get("predicted_pathway", pd.Series([""] * len(out))).map(clean_text)
    true_sub = out["true_subpathway"].map(clean_text)
    pred_sub = out.get("predicted_subpathway", pd.Series([""] * len(out))).map(clean_text)
    out["pathway_match"] = true_path.ne("") & pred_path.ne("") & true_path.eq(pred_path)
    out["subpathway_match"] = true_sub.ne("") & pred_sub.ne("") & true_sub.eq(pred_sub)

    out["parent_mz"] = pd.to_numeric(first_present(out, ["parent_mz", "known_mz"], default=pd.NA), errors="coerce")
    out["target_neutral_mass"] = pd.to_numeric(out.get("target_neutral_mass", pd.Series([pd.NA] * len(out))), errors="coerce")
    out["retention_index"] = pd.to_numeric(
        first_present(out, ["retention_index", "mlims_ri"], default=pd.NA),
        errors="coerce",
    )
    out["retention_index"] = out["retention_index"].fillna(out["feature_key"].map(parse_feature_ri))

    for column in OUTPUT_COLUMNS:
        if column not in out.columns:
            out[column] = pd.NA
    result = out[OUTPUT_COLUMNS].copy()
    for column in NUMERIC_COLUMNS:
        if column in result.columns:
            result[column] = pd.to_numeric(result[column], errors="coerce")
    for column in BOOLEAN_COLUMNS:
        if column in result.columns:
            result[column] = result[column].fillna(False).astype(bool)
    return result.sort_values(["library_id", "control_id", "rank", "predicted_smiles"], kind="mergesort").reset_index(drop=True)


def join_unique_text(values: pd.Series) -> str:
    unique = sorted({clean_text(value) for value in values if clean_text(value)})
    return "; ".join(unique)


def load_integrated_kpi_counts(path: Path) -> dict[str, int]:
    if not path.exists():
        raise FileNotFoundError(path)
    summary = read_csv_clean(path)
    required = {"n_controls", "top10_ge_0p7", "rank1_ge_0p7"}
    missing = sorted(required.difference(summary.columns))
    if missing:
        raise ValueError(f"Missing required columns in {path}: {', '.join(missing)}")
    row = summary.iloc[0]
    return {
        "n_controls": int(pd.to_numeric(row["n_controls"], errors="raise")),
        "top10_ge_0p7": int(pd.to_numeric(row["top10_ge_0p7"], errors="raise")),
        "rank1_ge_0p7": int(pd.to_numeric(row["rank1_ge_0p7"], errors="raise")),
    }


def build_row_source_summary(combined: pd.DataFrame) -> pd.DataFrame:
    records: list[dict[str, Any]] = []
    for scope, group in combined.groupby("scope", sort=True):
        control_records: list[dict[str, Any]] = []
        for control_id, control_group in group.groupby("control_id", sort=False):
            rank_values = pd.to_numeric(control_group["rank"], errors="coerce")
            rank1_group = control_group.loc[rank_values.eq(1)]
            control_records.append(
                {
                    "control_id": control_id,
                    "has_p0_reference": bool(control_group["has_p0_reference"].any()),
                    "is_true_unknown": bool(control_group["is_true_unknown"].any()),
                    "top10_p0_hit": bool(control_group["p0_match_ge_0p7"].any()),
                    "rank1_p0_hit": bool(rank1_group["p0_match_ge_0p7"].any()) if not rank1_group.empty else False,
                }
            )
        control = pd.DataFrame(control_records)
        p0_controls = control.loc[control["has_p0_reference"]]
        rank1_rows = group.loc[pd.to_numeric(group["rank"], errors="coerce").eq(1)]
        records.append(
            {
                "scope": scope,
                "library_id": scope_to_library_id(scope),
                "metric_basis": "top10_rows",
                "row_source_workflow": join_unique_text(group["row_source_workflow"]),
                "row_source_path": join_unique_text(group["row_source_path"]),
                "kpi_source_workflow": join_unique_text(group["kpi_source_workflow"]),
                "kpi_source_path": join_unique_text(group["kpi_source_path"]),
                "n_top10_rows": int(group.shape[0]),
                "n_controls": int(control.shape[0]),
                "n_p0_reference_controls": int(p0_controls.shape[0]),
                "n_true_unknown_controls": int(control["is_true_unknown"].sum()),
                "n_true_unknown_rows": int(group["is_true_unknown"].sum()),
                "n_controls_with_top10_p0_match": int(p0_controls["top10_p0_hit"].sum()) if not p0_controls.empty else 0,
                "top10_p0_match_rate": (
                    float(p0_controls["top10_p0_hit"].mean()) if not p0_controls.empty else 0.0
                ),
                "n_controls_with_rank1_p0_match": int(p0_controls["rank1_p0_hit"].sum()) if not p0_controls.empty else 0,
                "rank1_p0_match_rate": (
                    float(p0_controls["rank1_p0_hit"].mean()) if not p0_controls.empty else 0.0
                ),
                "n_candidate_p0_hits": int(group["p0_match_ge_0p7"].sum()),
                "mean_rank1_unified_final_score": float(
                    pd.to_numeric(rank1_rows["unified_final_score"], errors="coerce").mean()
                )
                if not rank1_rows.empty
                else 0.0,
                "mean_top10_unified_final_score": float(
                    pd.to_numeric(group["unified_final_score"], errors="coerce").mean()
                )
                if not group.empty
                else 0.0,
            }
        )
    return pd.DataFrame(records).sort_values("library_id").reset_index(drop=True)


def build_summary_tables(
    combined: pd.DataFrame,
    source_configs: dict[str, dict[str, str | Path]] | None = None,
) -> dict[str, pd.DataFrame]:
    row_source_summary = build_row_source_summary(combined)
    library_summary = row_source_summary.copy()
    source_configs = source_configs or {}

    source_records: list[dict[str, Any]] = []
    reconciliation_records: list[dict[str, Any]] = []
    for idx, row in library_summary.iterrows():
        scope = clean_text(row["scope"])
        config = source_configs.get(scope, {})
        metric_basis = clean_text(config.get("metric_basis", "top10_rows")) or "top10_rows"
        row_source_workflow = clean_text(config.get("row_source_workflow", row["row_source_workflow"]))
        row_source_path = clean_text(config.get("row_source_path", row["row_source_path"]))
        kpi_source_workflow = clean_text(config.get("kpi_source_workflow", row["kpi_source_workflow"]))
        kpi_source_path = clean_text(config.get("kpi_source_path", row["kpi_source_path"]))

        library_summary.loc[idx, "metric_basis"] = metric_basis
        library_summary.loc[idx, "row_source_workflow"] = row_source_workflow
        library_summary.loc[idx, "row_source_path"] = row_source_path
        library_summary.loc[idx, "kpi_source_workflow"] = kpi_source_workflow
        library_summary.loc[idx, "kpi_source_path"] = kpi_source_path

        row_top10_hits = int(row["n_controls_with_top10_p0_match"])
        row_rank1_hits = int(row["n_controls_with_rank1_p0_match"])
        row_controls = int(row["n_p0_reference_controls"])
        headline_controls = row_controls
        headline_top10_hits = row_top10_hits
        headline_rank1_hits = row_rank1_hits

        if metric_basis == "integrated_audit":
            counts = load_integrated_kpi_counts(Path(kpi_source_path))
            headline_controls = counts["n_controls"]
            headline_top10_hits = counts["top10_ge_0p7"]
            headline_rank1_hits = counts["rank1_ge_0p7"]
            library_summary.loc[idx, "n_controls"] = headline_controls
            library_summary.loc[idx, "n_p0_reference_controls"] = headline_controls
            library_summary.loc[idx, "n_controls_with_top10_p0_match"] = headline_top10_hits
            library_summary.loc[idx, "top10_p0_match_rate"] = (
                float(headline_top10_hits / headline_controls) if headline_controls else 0.0
            )
            library_summary.loc[idx, "n_controls_with_rank1_p0_match"] = headline_rank1_hits
            library_summary.loc[idx, "rank1_p0_match_rate"] = (
                float(headline_rank1_hits / headline_controls) if headline_controls else 0.0
            )

        source_records.append(
            {
                "scope": scope,
                "library_id": row["library_id"],
                "row_source_workflow": row_source_workflow,
                "row_source_path": row_source_path,
                "kpi_source_workflow": kpi_source_workflow,
                "kpi_source_path": kpi_source_path,
                "metric_basis": metric_basis,
                "headline_n_controls": int(headline_controls),
                "headline_top10_ge_0p7": int(headline_top10_hits),
                "headline_rank1_ge_0p7": int(headline_rank1_hits),
                "headline_top10_rate": float(headline_top10_hits / headline_controls) if headline_controls else 0.0,
                "headline_rank1_rate": float(headline_rank1_hits / headline_controls) if headline_controls else 0.0,
            }
        )
        reconciliation_records.append(
            {
                "scope": scope,
                "library_id": row["library_id"],
                "row_source_metric_basis": "top10_rows",
                "headline_metric_basis": metric_basis,
                "row_source_controls": int(row_controls),
                "row_source_top10_hits": int(row_top10_hits),
                "row_source_rank1_hits": int(row_rank1_hits),
                "headline_controls": int(headline_controls),
                "headline_top10_hits": int(headline_top10_hits),
                "headline_rank1_hits": int(headline_rank1_hits),
                "top10_hit_delta": int(headline_top10_hits - row_top10_hits),
                "rank1_hit_delta": int(headline_rank1_hits - row_rank1_hits),
                "explanation": (
                    "Authoritative dashboard metrics use the integrated audit; candidate rows use selected top10 sources."
                    if metric_basis == "integrated_audit"
                    else "Headline metrics and candidate rows share the top10-row basis."
                ),
            }
        )

    rank_hits = combined.loc[combined["p0_match_ge_0p7"]].copy()
    rank_distribution = (
        rank_hits.pivot_table(index="scope", columns="rank", values="control_id", aggfunc="count", fill_value=0)
        .reindex(sorted(combined["scope"].dropna().unique()))
        .fillna(0)
    )
    for rank in range(1, 11):
        if rank not in rank_distribution.columns:
            rank_distribution[rank] = 0
    rank_distribution = rank_distribution[[rank for rank in range(1, 11)]].reset_index()
    rank_distribution.columns = ["scope"] + [f"rank_{rank}" for rank in range(1, 11)]

    known_unknown_counts = (
        combined.assign(control_type=combined["is_true_unknown"].map(lambda flag: "true_unknown" if flag else "known_or_review"))
        .pivot_table(index="scope", columns="control_type", values="control_id", aggfunc="count", fill_value=0)
        .reset_index()
    )
    for column in ("known_or_review", "true_unknown"):
        if column not in known_unknown_counts.columns:
            known_unknown_counts[column] = 0
    known_unknown_counts = known_unknown_counts[["scope", "known_or_review", "true_unknown"]]

    score_by_rank = (
        combined.groupby(["scope", "rank"], as_index=False)["unified_final_score"]
        .mean()
        .pivot(index="rank", columns="scope", values="unified_final_score")
        .reset_index()
        .sort_values("rank")
    )

    rank1 = combined.loc[pd.to_numeric(combined["rank"], errors="coerce").eq(1)].copy()
    tier_distribution = (
        rank1.pivot_table(index="scope", columns="pathway_confidence_tier", values="control_id", aggfunc="count", fill_value=0)
        .reset_index()
    )
    for tier in ("high", "medium", "low"):
        if tier not in tier_distribution.columns:
            tier_distribution[tier] = 0
    tier_distribution = tier_distribution[["scope", "high", "medium", "low"]]

    return {
        "library_summary": library_summary,
        "row_source_summary": row_source_summary,
        "source_provenance": pd.DataFrame(source_records).sort_values("library_id").reset_index(drop=True),
        "p0_metric_reconciliation": pd.DataFrame(reconciliation_records).sort_values("library_id").reset_index(drop=True),
        "rank_distribution": rank_distribution,
        "known_unknown_counts": known_unknown_counts,
        "score_by_rank": score_by_rank,
        "tier_distribution": tier_distribution,
    }


def write_csv_outputs(output_root: Path, combined: pd.DataFrame, summaries: dict[str, pd.DataFrame]) -> None:
    output_root.mkdir(parents=True, exist_ok=True)
    combined.to_csv(output_root / "all_libraries_top10_predictions_with_p0_flags.csv", index=False)
    combined.loc[combined["is_true_unknown"]].to_csv(output_root / "true_unknown_top10_predictions.csv", index=False)
    combined.loc[combined["p0_match_ge_0p7"]].to_csv(output_root / "p0_matches.csv", index=False)
    summaries["library_summary"].to_csv(output_root / "p0_match_summary_by_library.csv", index=False)
    summaries["row_source_summary"].to_csv(output_root / "p0_match_row_source_summary_by_library.csv", index=False)
    summaries["source_provenance"].to_csv(output_root / "source_provenance.csv", index=False)
    summaries["p0_metric_reconciliation"].to_csv(output_root / "p0_metric_reconciliation.csv", index=False)
    lib209_best = combined.loc[
        combined["scope"].eq("lib_209_negative")
        & combined["row_source_workflow"].eq("lib209_step7_targeted_v2_auto")
    ].copy()
    if not lib209_best.empty:
        lib209_best["canonical_smiles"] = lib209_best["predicted_smiles"].map(
            lambda value: clean_text(canonicalize_smiles(clean_text(value)))
        )
        assignment_columns = [
            "canonical_smiles",
            "predicted_smiles",
            "predicted_pathway",
            "predicted_subpathway",
            "pathway_confidence_tier",
            "pathway_assignment_status",
            "pathway_vote_share",
            "pathway_best_reference_similarity",
        ]
        lib209_best.loc[lib209_best["canonical_smiles"].ne(""), assignment_columns].drop_duplicates(
            subset=["canonical_smiles"], keep="first"
        ).sort_values("canonical_smiles").to_csv(
            output_root / "lib209_best_workflow_unique_smiles_pathway_assignments.csv", index=False
        )
    for name, df in summaries.items():
        if name not in {"library_summary", "row_source_summary", "source_provenance", "p0_metric_reconciliation"}:
            df.to_csv(output_root / f"{name}.csv", index=False)


def dataframe_to_rows(df: pd.DataFrame) -> list[list[Any]]:
    rows = [df.columns.tolist()]
    for row in df.itertuples(index=False, name=None):
        cleaned = []
        for value in row:
            if pd.isna(value):
                cleaned.append(None)
            elif isinstance(value, (pd.Timestamp,)):
                cleaned.append(value.to_pydatetime())
            elif isinstance(value, float) and not math.isfinite(value):
                cleaned.append(None)
            else:
                cleaned.append(value)
        rows.append(cleaned)
    return rows


def write_dataframe_sheet(wb: Workbook, sheet_name: str, df: pd.DataFrame, freeze: str = "A2") -> Any:
    ws = wb.create_sheet(sheet_name[:31])
    for row in dataframe_to_rows(df):
        ws.append(row)
    style_sheet(ws, freeze=freeze)
    return ws


def style_sheet(ws: Any, freeze: str = "A2") -> None:
    header_fill = PatternFill("solid", fgColor="1F4E78")
    header_font = Font(color="FFFFFF", bold=True)
    thin_gray = Side(style="thin", color="D9E2F3")
    ws.freeze_panes = freeze
    ws.auto_filter.ref = ws.dimensions
    for cell in ws[1]:
        cell.fill = header_fill
        cell.font = header_font
        cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
        cell.border = Border(bottom=thin_gray)
    max_widths: dict[int, int] = {}
    for row in ws.iter_rows(min_row=1, max_row=min(ws.max_row, 250), max_col=ws.max_column):
        for cell in row:
            text = "" if cell.value is None else str(cell.value)
            max_widths[cell.column] = min(max(max_widths.get(cell.column, 0), len(text) + 2), 48)
    for column_idx, width in max_widths.items():
        ws.column_dimensions[get_column_letter(column_idx)].width = max(10, width)
    for row in ws.iter_rows(min_row=2, max_row=ws.max_row):
        for cell in row:
            cell.alignment = Alignment(vertical="top", wrap_text=False)
    add_boolean_conditionals(ws)


def add_boolean_conditionals(ws: Any) -> None:
    headers = [cell.value for cell in ws[1]]
    green = PatternFill("solid", fgColor="C6EFCE")
    red = PatternFill("solid", fgColor="FFC7CE")
    gray = PatternFill("solid", fgColor="E7E6E6")
    for target in ("p0_match_ge_0p7", "pathway_match", "subpathway_match"):
        if target in headers and ws.max_row > 1:
            col = get_column_letter(headers.index(target) + 1)
            ref = f"{col}2:{col}{ws.max_row}"
            ws.conditional_formatting.add(ref, CellIsRule(operator="equal", formula=["TRUE"], fill=green))
            ws.conditional_formatting.add(ref, CellIsRule(operator="equal", formula=["FALSE"], fill=red))
    if "p0_match_status" in headers and ws.max_row > 1:
        col = get_column_letter(headers.index("p0_match_status") + 1)
        ref = f"{col}2:{col}{ws.max_row}"
        ws.conditional_formatting.add(ref, CellIsRule(operator="equal", formula=['"true_unknown_no_p0_reference"'], fill=gray))


def write_summary_sections(ws: Any, summaries: dict[str, pd.DataFrame]) -> dict[str, tuple[int, int]]:
    positions: dict[str, tuple[int, int]] = {}
    row = 1
    for name, df in summaries.items():
        positions[name] = (row + 1, row + len(df) + 1)
        ws.cell(row=row, column=1, value=name)
        ws.cell(row=row, column=1).font = Font(bold=True, size=13, color="1F4E78")
        for r_idx, values in enumerate(dataframe_to_rows(df), start=row + 1):
            for c_idx, value in enumerate(values, start=1):
                ws.cell(row=r_idx, column=c_idx, value=value)
        row += len(df) + 4
    style_sheet(ws, freeze="A1")
    return positions


def add_dashboard_chart(ws: Any, chart: Any, title: str, anchor: str) -> None:
    chart.title = title
    chart.width = 14
    chart.height = 7
    ws.add_chart(chart, anchor)


def build_dashboard(wb: Workbook, summaries: dict[str, pd.DataFrame], positions: dict[str, tuple[int, int]]) -> None:
    ws = wb.create_sheet("Dashboard", 0)
    ws["A1"] = "All-Library V14 Prediction Review"
    ws["A1"].font = Font(size=18, bold=True, color="1F4E78")
    ws["A2"] = "lib209 headline P0 metrics use best integrated audit; candidate rows use selected top10 source."
    ws["A2"].font = Font(italic=True, color="666666")

    summary = summaries["library_summary"]
    total_rows = int(summary["n_top10_rows"].sum())
    total_controls = int(summary["n_controls"].sum())
    total_hits = int(summary["n_controls_with_top10_p0_match"].sum())
    total_rank1_hits = int(summary["n_controls_with_rank1_p0_match"].sum())
    true_unknown_controls = int(summary["n_true_unknown_controls"].sum())
    kpis = [
        ("Top10 rows", total_rows),
        ("Controls", total_controls),
        ("P0 top10 hits", total_hits),
        ("P0 rank1 hits", total_rank1_hits),
        ("True unknown controls", true_unknown_controls),
    ]
    for idx, (label, value) in enumerate(kpis):
        col = 1 + idx * 3
        ws.cell(row=4, column=col, value=label)
        ws.cell(row=5, column=col, value=value)
        ws.cell(row=4, column=col).font = Font(bold=True, color="666666")
        ws.cell(row=5, column=col).font = Font(bold=True, size=16, color="1F4E78")

    summary_ws = wb["Summary Tables"]

    start, end = positions["library_summary"]
    summary_headers = [summary_ws.cell(row=start, column=col).value for col in range(1, summary_ws.max_column + 1)]
    chart = BarChart()
    chart.type = "bar"
    chart.style = 10
    for metric in ("top10_p0_match_rate", "rank1_p0_match_rate"):
        metric_col = summary_headers.index(metric) + 1
        data = Reference(summary_ws, min_col=metric_col, max_col=metric_col, min_row=start, max_row=end)
        chart.add_data(data, titles_from_data=True)
    cats = Reference(summary_ws, min_col=1, min_row=start + 1, max_row=end)
    chart.set_categories(cats)
    chart.y_axis.title = "Library"
    chart.x_axis.title = "Rate"
    add_dashboard_chart(ws, chart, "P0 Top10 and Rank1 Match Performance", "A8")

    start, end = positions["rank_distribution"]
    chart = BarChart()
    chart.type = "col"
    data = Reference(summary_ws, min_col=2, max_col=11, min_row=start, max_row=end)
    cats = Reference(summary_ws, min_col=1, min_row=start + 1, max_row=end)
    chart.add_data(data, titles_from_data=True)
    chart.set_categories(cats)
    chart.y_axis.title = "Candidate hits"
    add_dashboard_chart(ws, chart, "P0 Match Rank Distribution", "I8")

    start, end = positions["known_unknown_counts"]
    chart = BarChart()
    chart.type = "col"
    data = Reference(summary_ws, min_col=2, max_col=3, min_row=start, max_row=end)
    cats = Reference(summary_ws, min_col=1, min_row=start + 1, max_row=end)
    chart.add_data(data, titles_from_data=True)
    chart.set_categories(cats)
    chart.y_axis.title = "Rows"
    add_dashboard_chart(ws, chart, "Known/Review vs True Unknown Candidate Rows", "A24")

    start, end = positions["score_by_rank"]
    chart = LineChart()
    data = Reference(summary_ws, min_col=2, max_col=summary_ws.max_column, min_row=start, max_row=end)
    cats = Reference(summary_ws, min_col=1, min_row=start + 1, max_row=end)
    chart.add_data(data, titles_from_data=True)
    chart.set_categories(cats)
    chart.y_axis.title = "Mean unified_final_score"
    chart.x_axis.title = "Rank"
    add_dashboard_chart(ws, chart, "Mean V14 Score by Rank", "I24")

    start, end = positions["tier_distribution"]
    chart = BarChart()
    chart.type = "col"
    data = Reference(summary_ws, min_col=2, max_col=4, min_row=start, max_row=end)
    cats = Reference(summary_ws, min_col=1, min_row=start + 1, max_row=end)
    chart.add_data(data, titles_from_data=True)
    chart.set_categories(cats)
    chart.y_axis.title = "Rank1 controls"
    add_dashboard_chart(ws, chart, "Rank1 Pathway Confidence Tier", "A40")

    for col in range(1, 18):
        ws.column_dimensions[get_column_letter(col)].width = 14


def write_workbook(output_root: Path, combined: pd.DataFrame, summaries: dict[str, pd.DataFrame]) -> Path:
    if not HAS_OPENPYXL:
        raise RuntimeError("openpyxl is required to write the V14 review workbook")
    wb = Workbook()
    wb.remove(wb.active)
    summary_ws = wb.create_sheet("Summary Tables")
    positions = write_summary_sections(summary_ws, summaries)
    build_dashboard(wb, summaries, positions)

    write_dataframe_sheet(wb, "All Top10", combined)
    write_dataframe_sheet(wb, "True Unknown Top10", combined.loc[combined["is_true_unknown"]].copy())
    write_dataframe_sheet(wb, "Source Provenance", summaries["source_provenance"])
    for scope in sorted(combined["scope"].dropna().unique()):
        write_dataframe_sheet(wb, scope, combined.loc[combined["scope"].eq(scope)].copy())
    write_dataframe_sheet(wb, "P0 Matches", combined.loc[combined["p0_match_ge_0p7"]].copy())

    workbook_path = output_root / "v14_all_libraries_prediction_review.xlsx"
    wb.save(workbook_path)
    return workbook_path


def build_review_bundle(
    *,
    input_root: Path = DEFAULT_INPUT_ROOT,
    output_root: Path = DEFAULT_OUTPUT_ROOT,
    known_x: Path = DEFAULT_KNOWN_X,
    known_structures: Path = DEFAULT_KNOWN_STRUCTURES,
    reference_db: Path = DEFAULT_REFERENCE_DB,
    scopes: list[str] | tuple[str, ...] = DEFAULT_SCOPES,
    p0_threshold: float = DEFAULT_P0_THRESHOLD,
    skip_workbook: bool = False,
) -> dict[str, Path]:
    output_root.mkdir(parents=True, exist_ok=True)
    known_meta = load_known_metadata(known_x)
    reference_smiles = load_reference_smiles(known_structures)
    parsed_scopes = list(scopes)
    source_configs = source_config_for_scopes(parsed_scopes, input_root=input_root)
    pathway_assignment_map = load_pathway_assignment_map(input_root, "lib_209_negative")
    raw_frames = [
        load_scope_top10(
            scope,
            source_configs[scope],
            pathway_assignment_map=pathway_assignment_map if scope == "lib_209_negative" else None,
            reference_db=reference_db if scope == "lib_209_negative" else None,
        )
        for scope in parsed_scopes
    ]
    raw = pd.concat(raw_frames, ignore_index=True, sort=False)
    combined = add_metadata_and_flags(raw, known_meta, reference_smiles, p0_threshold=p0_threshold)
    summaries = build_summary_tables(combined, source_configs)
    write_csv_outputs(output_root, combined, summaries)
    paths = {
        "all_top10": output_root / "all_libraries_top10_predictions_with_p0_flags.csv",
        "true_unknown": output_root / "true_unknown_top10_predictions.csv",
        "summary": output_root / "p0_match_summary_by_library.csv",
        "row_source_summary": output_root / "p0_match_row_source_summary_by_library.csv",
        "source_provenance": output_root / "source_provenance.csv",
        "reconciliation": output_root / "p0_metric_reconciliation.csv",
        "lib209_best_pathway_assignments": output_root
        / "lib209_best_workflow_unique_smiles_pathway_assignments.csv",
        "p0_matches": output_root / "p0_matches.csv",
    }
    if not skip_workbook:
        paths["workbook"] = write_workbook(output_root, combined, summaries)
    return paths


def main() -> None:
    args = parse_args()
    paths = build_review_bundle(
        input_root=args.input_root,
        output_root=args.output_root,
        known_x=args.known_x,
        known_structures=args.known_structures,
        reference_db=args.reference_db,
        scopes=parse_scope_list(args.scopes),
        p0_threshold=float(args.p0_threshold),
        skip_workbook=bool(args.skip_workbook),
    )
    for label, path in paths.items():
        print(f"{label}: {path}")


if __name__ == "__main__":
    main()
