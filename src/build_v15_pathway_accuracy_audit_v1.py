#!/usr/bin/env python3
"""Audit V15 pathway/subpathway accuracy by library and P0/P1 cohort."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd

from build_v14_all_library_prediction_review_v1 import clean_text, scope_to_library_id
from phase0_utils import canonicalize_smiles
from x_domain_analog_utils import dataframe_to_markdown, write_csv


DEFAULT_REVIEW_CSV = Path("models/v15_all_library_prediction_review_v1/all_libraries_top10_predictions_with_p0_flags.csv")
DEFAULT_KNOWN_X = Path("known_X.csv")
DEFAULT_OUTPUT_ROOT = Path("models/v15_pathway_accuracy_audit_v1")
DEFAULT_V14_PATHWAY_ROOT = Path("models/v14_ranked_pathway_annotations_v1")
DEFAULT_INTEGRATED_ROOT = Path("models/integrated_denovo_audit_v15")
DEFAULT_PREDICTOR_ROOT = Path("models/v15_pathway_predictor_v2")
COHORTS = ("p0", "p1", "p0_and_p1", "all")
SURFACES = ("rerank", "integrated")
RANKS = ("1", "top10", "oracle")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--review-csv", type=Path, default=DEFAULT_REVIEW_CSV)
    parser.add_argument("--known-x", type=Path, default=DEFAULT_KNOWN_X)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--v14-pathway-root", type=Path, default=DEFAULT_V14_PATHWAY_ROOT)
    parser.add_argument("--integrated-root", type=Path, default=DEFAULT_INTEGRATED_ROOT)
    parser.add_argument("--predicted-source", choices=("v1", "v2"), default="v1")
    parser.add_argument("--predictor-root", type=Path, default=DEFAULT_PREDICTOR_ROOT)
    return parser.parse_args()


def norm_label(value: object) -> str:
    return clean_text(value).lower()


def bool_series(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False).astype(bool)
    text = series.fillna("").astype(str).str.strip().str.lower()
    return text.isin({"true", "1", "yes", "y", "t"})


def canonical(value: object) -> str:
    return clean_text(canonicalize_smiles(clean_text(value)))


def normalize_level(value: object) -> str:
    text = clean_text(value)
    numeric = pd.to_numeric(text, errors="coerce")
    if pd.notna(numeric):
        return str(int(float(numeric)))
    return text


def load_review(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, low_memory=False)
    required = {
        "control_id",
        "library_id",
        "scope",
        "predicted_smiles",
        "rank",
        "pathway_match",
        "subpathway_match",
        "true_pathway",
        "true_subpathway",
        "has_p0_reference",
        "p0_match_ge_0p7",
        "known_level",
        "known_sublevel",
    }
    missing = sorted(required.difference(df.columns))
    if missing:
        raise ValueError(f"Missing required V15 review columns: {', '.join(missing)}")
    df = df.copy()
    df["control_id"] = df["control_id"].map(clean_text)
    df["scope"] = df["scope"].map(clean_text)
    df["library_id"] = pd.to_numeric(df["library_id"], errors="coerce").astype("Int64")
    df["rank"] = pd.to_numeric(df["rank"], errors="coerce")
    df["has_p0_reference"] = bool_series(df["has_p0_reference"])
    df["p0_match_ge_0p7"] = bool_series(df["p0_match_ge_0p7"])
    df["pathway_match"] = bool_series(df["pathway_match"])
    df["subpathway_match"] = bool_series(df["subpathway_match"])
    df["known_level_norm"] = df["known_level"].map(normalize_level)
    df["known_sublevel_norm"] = df["known_sublevel"].map(normalize_level)
    df["true_pathway_norm"] = df["true_pathway"].map(norm_label)
    df["true_subpathway_norm"] = df["true_subpathway"].map(norm_label)
    return df


def control_table(review: pd.DataFrame) -> pd.DataFrame:
    cols = [
        "scope",
        "library_id",
        "control_id",
        "has_p0_reference",
        "known_level_norm",
        "known_sublevel_norm",
        "true_pathway",
        "true_subpathway",
        "true_pathway_norm",
        "true_subpathway_norm",
    ]
    out = review[cols].drop_duplicates(subset=["scope", "control_id"], keep="first").copy()
    out["is_p0"] = out["has_p0_reference"].astype(bool)
    out["is_p1"] = out["known_level_norm"].isin({"1", "2"})
    return out


def cohort_controls(controls: pd.DataFrame, cohort: str) -> pd.DataFrame:
    if cohort == "p0":
        return controls.loc[controls["is_p0"]].copy()
    if cohort == "p1":
        return controls.loc[controls["is_p1"]].copy()
    if cohort == "p0_and_p1":
        return controls.loc[controls["is_p0"] & controls["is_p1"]].copy()
    if cohort == "all":
        return controls.copy()
    raise ValueError(f"Unsupported cohort: {cohort}")


def matches_from_predictions(pred: pd.DataFrame, controls: pd.DataFrame) -> pd.DataFrame:
    if pred.empty:
        return pd.DataFrame(columns=["scope", "control_id", "pathway_match", "subpathway_match", "both_match"])
    keep = controls[
        ["scope", "control_id", "true_pathway_norm", "true_subpathway_norm"]
    ].drop_duplicates(subset=["scope", "control_id"])
    out = pred.merge(keep, on=["scope", "control_id"], how="left", sort=False)
    out["predicted_pathway_norm"] = out["predicted_pathway"].map(norm_label)
    out["predicted_subpathway_norm"] = out["predicted_subpathway"].map(norm_label)
    out["pathway_match"] = out["true_pathway_norm"].ne("") & out["predicted_pathway_norm"].eq(out["true_pathway_norm"])
    out["subpathway_match"] = out["true_subpathway_norm"].ne("") & out["predicted_subpathway_norm"].eq(
        out["true_subpathway_norm"]
    )
    out["both_match"] = out["pathway_match"] & out["subpathway_match"]
    return out[["scope", "control_id", "pathway_match", "subpathway_match", "both_match"]].copy()


def load_assignment_map(scope: str, root: Path) -> pd.DataFrame:
    path = root / scope / "unique_smiles_pathway_assignments.csv"
    if not path.exists():
        return pd.DataFrame(columns=["canonical_smiles", "predicted_pathway", "predicted_subpathway"])
    df = pd.read_csv(path, low_memory=False)
    for col in ("canonical_smiles", "predicted_pathway", "predicted_subpathway"):
        if col not in df.columns:
            df[col] = ""
    return df[["canonical_smiles", "predicted_pathway", "predicted_subpathway"]].drop_duplicates(
        subset=["canonical_smiles"], keep="first"
    )


def predictions_from_smiles(scope: str, smiles_df: pd.DataFrame, assignment_map: pd.DataFrame) -> pd.DataFrame:
    if smiles_df.empty:
        return pd.DataFrame(columns=["scope", "control_id", "predicted_pathway", "predicted_subpathway"])
    out = smiles_df.copy()
    out["canonical_smiles"] = out["predicted_smiles"].map(canonical)
    out = out.merge(assignment_map, on="canonical_smiles", how="left", sort=False)
    out["scope"] = scope
    out["predicted_pathway"] = out["predicted_pathway"].fillna("")
    out["predicted_subpathway"] = out["predicted_subpathway"].fillna("")
    return out[["scope", "control_id", "predicted_pathway", "predicted_subpathway"]].copy()


def rerank_v1_matches(review: pd.DataFrame, scope: str, rank: str, v14_root: Path, controls: pd.DataFrame) -> pd.DataFrame:
    scoped = review.loc[review["scope"].eq(scope)].copy()
    if rank == "1":
        rows = scoped.loc[scoped["rank"].eq(1)].copy()
        out = rows[["scope", "control_id", "pathway_match", "subpathway_match"]].copy()
        out["both_match"] = out["pathway_match"] & out["subpathway_match"]
        return out
    if rank == "top10":
        rows = scoped.loc[scoped["rank"].between(1, 10, inclusive="both")].copy()
        return (
            rows.groupby(["scope", "control_id"], as_index=False)
            .agg(pathway_match=("pathway_match", "max"), subpathway_match=("subpathway_match", "max"))
            .assign(both_match=lambda df: df["pathway_match"] & df["subpathway_match"])
        )
    if rank == "oracle":
        path = v14_root / scope / "all_scored_candidates_with_pathways.csv"
        if not path.exists():
            return pd.DataFrame(columns=["scope", "control_id", "pathway_match", "subpathway_match", "both_match"])
        usecols = [
            "control_id",
            "predicted_pathway",
            "predicted_subpathway",
        ]
        rows = pd.read_csv(path, usecols=lambda col: col in usecols, low_memory=False)
        rows["scope"] = scope
        matches = matches_from_predictions(rows, controls)
        return (
            matches.groupby(["scope", "control_id"], as_index=False)
            .agg(pathway_match=("pathway_match", "max"), subpathway_match=("subpathway_match", "max"))
            .assign(both_match=lambda df: df["pathway_match"] & df["subpathway_match"])
        )
    raise ValueError(rank)


def integrated_v1_matches(scope: str, rank: str, integrated_root: Path, v14_root: Path, controls: pd.DataFrame) -> pd.DataFrame:
    path = integrated_root / scope / "integrated_per_control_best.csv"
    if not path.exists():
        return pd.DataFrame(columns=["scope", "control_id", "pathway_match", "subpathway_match", "both_match"])
    df = pd.read_csv(path, low_memory=False)
    smiles_col = {"1": "rank1_smiles", "top10": "top10_smiles", "oracle": "full_pool_oracle_smiles"}[rank]
    if smiles_col not in df.columns:
        return pd.DataFrame(columns=["scope", "control_id", "pathway_match", "subpathway_match", "both_match"])
    smiles_df = df[["control_id", smiles_col]].rename(columns={smiles_col: "predicted_smiles"}).copy()
    smiles_df = smiles_df.loc[smiles_df["predicted_smiles"].map(clean_text).ne("")]
    pred = predictions_from_smiles(scope, smiles_df, load_assignment_map(scope, v14_root))
    return matches_from_predictions(pred, controls)


def v2_matches(scope: str, rank: str, predictor_root: Path, controls: pd.DataFrame) -> pd.DataFrame:
    path = predictor_root / "per_control_predictions.csv"
    if not path.exists():
        return pd.DataFrame(columns=["scope", "control_id", "pathway_match", "subpathway_match", "both_match"])
    df = pd.read_csv(path, low_memory=False)
    df["scope"] = df["scope"].map(clean_text) if "scope" in df.columns else ""
    scoped = df.loc[df["scope"].eq(scope)].copy()
    if scoped.empty:
        return pd.DataFrame(columns=["scope", "control_id", "pathway_match", "subpathway_match", "both_match"])
    pred = scoped.rename(
        columns={
            "predicted_pathway_v2": "predicted_pathway",
            "predicted_subpathway_v2": "predicted_subpathway",
        }
    )
    pred["scope"] = scope
    return matches_from_predictions(pred, controls)


def summarize_matches(
    *,
    controls: pd.DataFrame,
    matches: pd.DataFrame,
    scope: str,
    surface: str,
    rank: str,
    cohort: str,
) -> dict[str, Any]:
    cohort_df = cohort_controls(controls.loc[controls["scope"].eq(scope)].copy(), cohort)
    pred = cohort_df.merge(matches, on=["scope", "control_id"], how="left", sort=False)
    for col in ("pathway_match", "subpathway_match", "both_match"):
        pred[col] = pred[col].fillna(False).astype(bool)
    path_den = pred.loc[pred["true_pathway_norm"].ne("")]
    sub_den = pred.loc[pred["true_subpathway_norm"].ne("")]
    both_den = pred.loc[pred["true_pathway_norm"].ne("") & pred["true_subpathway_norm"].ne("")]
    pathway_correct = int(path_den["pathway_match"].sum())
    subpathway_correct = int(sub_den["subpathway_match"].sum())
    both_correct = int(both_den["both_match"].sum())
    return {
        "lib_id": int(scope_to_library_id(scope)),
        "ion_mode": scope.rsplit("_", 1)[1],
        "scope": scope,
        "cohort": cohort,
        "surface": surface,
        "rank": rank,
        "n_controls": int(cohort_df["control_id"].nunique()),
        "n_with_truth_pathway": int(path_den["control_id"].nunique()),
        "n_with_truth_subpathway": int(sub_den["control_id"].nunique()),
        "pathway_correct": pathway_correct,
        "subpathway_correct": subpathway_correct,
        "both_correct": both_correct,
        "pathway_accuracy": float(pathway_correct / len(path_den)) if len(path_den) else 0.0,
        "subpathway_accuracy": float(subpathway_correct / len(sub_den)) if len(sub_den) else 0.0,
        "both_accuracy": float(both_correct / len(both_den)) if len(both_den) else 0.0,
    }


def build_audit(
    review_csv: Path,
    output_root: Path,
    v14_pathway_root: Path,
    integrated_root: Path,
    predicted_source: str,
    predictor_root: Path,
) -> dict[str, Path]:
    output_root.mkdir(parents=True, exist_ok=True)
    review = load_review(review_csv)
    controls = control_table(review)
    records = []
    for scope in sorted(review["scope"].dropna().unique()):
        for surface in SURFACES:
            for rank in RANKS:
                if predicted_source == "v2":
                    matches = v2_matches(scope, rank, predictor_root, controls)
                elif surface == "rerank":
                    matches = rerank_v1_matches(review, scope, rank, v14_pathway_root, controls)
                else:
                    matches = integrated_v1_matches(scope, rank, integrated_root, v14_pathway_root, controls)
                for cohort in COHORTS:
                    records.append(
                        summarize_matches(
                            controls=controls,
                            matches=matches,
                            scope=scope,
                            surface=surface,
                            rank=rank,
                            cohort=cohort,
                        )
                    )
    breakdown = pd.DataFrame(records).sort_values(["lib_id", "cohort", "surface", "rank"]).reset_index(drop=True)
    write_csv(output_root / "per_lib_p0_p1_breakdown.csv", breakdown)
    summary = (
        breakdown.groupby(["cohort", "surface", "rank"], as_index=False)
        .agg(
            n_controls=("n_controls", "sum"),
            n_with_truth_pathway=("n_with_truth_pathway", "sum"),
            n_with_truth_subpathway=("n_with_truth_subpathway", "sum"),
            pathway_correct=("pathway_correct", "sum"),
            subpathway_correct=("subpathway_correct", "sum"),
            both_correct=("both_correct", "sum"),
        )
        .copy()
    )
    summary["pathway_accuracy"] = summary["pathway_correct"] / summary["n_with_truth_pathway"].replace({0: pd.NA})
    summary["subpathway_accuracy"] = summary["subpathway_correct"] / summary["n_with_truth_subpathway"].replace({0: pd.NA})
    both_den = summary[["n_with_truth_pathway", "n_with_truth_subpathway"]].min(axis=1).replace({0: pd.NA})
    summary["both_accuracy"] = summary["both_correct"] / both_den
    write_csv(output_root / "cross_lib_summary.csv", summary)
    report = render_report(breakdown, summary, predicted_source)
    (output_root / "summary.md").write_text(report)
    return {
        "breakdown": output_root / "per_lib_p0_p1_breakdown.csv",
        "summary": output_root / "cross_lib_summary.csv",
        "report": output_root / "summary.md",
    }


def render_report(breakdown: pd.DataFrame, summary: pd.DataFrame, predicted_source: str) -> str:
    focus = breakdown.loc[
        breakdown["cohort"].isin(["p0", "p1"])
        & breakdown["surface"].eq("rerank")
        & breakdown["rank"].isin(["1", "top10"])
    ].copy()
    return "\n".join(
        [
            f"# V15 Pathway Accuracy Audit ({predicted_source.upper()})",
            "",
            "## Per-Library P0/P1 Focus",
            "",
            dataframe_to_markdown(
                focus[
                    [
                        "scope",
                        "cohort",
                        "rank",
                        "n_with_truth_pathway",
                        "pathway_correct",
                        "pathway_accuracy",
                        "n_with_truth_subpathway",
                        "subpathway_correct",
                        "subpathway_accuracy",
                        "both_accuracy",
                    ]
                ]
            ),
            "",
            "## Cross-Lib Summary",
            "",
            dataframe_to_markdown(summary),
            "",
        ]
    )


def main() -> None:
    args = parse_args()
    paths = build_audit(
        review_csv=args.review_csv,
        output_root=args.output_root,
        v14_pathway_root=args.v14_pathway_root,
        integrated_root=args.integrated_root,
        predicted_source=args.predicted_source,
        predictor_root=args.predictor_root,
    )
    print(json.dumps({key: str(value) for key, value in paths.items()}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
