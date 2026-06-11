#!/usr/bin/env python3
"""Build V19b high-level accuracy summary CSV (P0/P1 structure + pathway)."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd

from build_v14_all_library_prediction_review_v1 import clean_text, scope_to_library_id
from build_v15_pathway_accuracy_audit_v1 import bool_series, norm_label, normalize_level

DEFAULT_OUTPUT_ROOT = Path("models/v19b_high_level_accuracy_summary_v1")
DEFAULT_REVIEW_CSV = Path(
    "models/v19b_all_library_updated_summary_v1/v19b_all_libraries_top10_predictions_with_p0_flags.csv"
)
DEFAULT_DELIVERABLE_CSV = Path("models/v19b_all_library_updated_summary_v1/v19b_updated_summary.csv")
DEFAULT_P1_ROOT = Path("models/p1_reference_v1")
DEFAULT_KNOWN_X = Path("known_X.csv")
DEFAULT_TAU = 0.7
PLACEHOLDER_PATHWAY_LABELS = frozenset(
    {
        "unknown",
        "do not include",
        "can add this",
        "ask dts",
        "ask annie/matt/brian i",
    }
)

LIB_CONFIG = (
    (209, "negative", "lib_209_negative", Path("models/integrated_denovo_audit_v15/lib_209_negative")),
    (305, "negative", "lib_305_negative", Path("models/integrated_denovo_audit_v19b/lib_305_negative")),
    (400, "positive", "lib_400_positive", Path("models/integrated_denovo_audit_v19b/lib_400_positive")),
    (402, "positive", "lib_402_positive", Path("models/integrated_denovo_audit_v19b/lib_402_positive")),
)

P1_LIB_DIRS = (
    ("lib_209_negative_plus_2c_v19", 209, "negative"),
    ("lib_305_negative_plus_2c_v19", 305, "negative"),
    ("lib_400_positive_plus_2c_v19", 400, "positive"),
    ("lib_402_positive_plus_2c_v19", 402, "positive"),
)
P1_PREFERRED_DIR_SUFFIX = "_plus_2c_v20_pool"
P1_LEGACY_DIR_SUFFIX = "_plus_2c_v19"

SUMMARY_COLUMNS = [
    "metric_group",
    "cohort",
    "accuracy_target",
    "rank_view",
    "lib_id",
    "ion_mode",
    "breakdown",
    "eligible_scope",
    "n_eligible",
    "n_correct",
    "pct_correct",
    "notes",
]

ROLLUP_SPECS = (
    ("negative_libs", "negative", {209, 305}),
    ("positive_libs", "positive", {400, 402}),
    ("all_libs", "", {209, 305, 400, 402}),
)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--review-csv", type=Path, default=DEFAULT_REVIEW_CSV)
    parser.add_argument("--deliverable-csv", type=Path, default=DEFAULT_DELIVERABLE_CSV)
    parser.add_argument("--p1-root", type=Path, default=DEFAULT_P1_ROOT)
    parser.add_argument("--known-x", type=Path, default=DEFAULT_KNOWN_X)
    parser.add_argument("--tau", type=float, default=DEFAULT_TAU)
    parser.add_argument(
        "--s0-cohort-mode",
        choices=("full_benchmark", "deliverable_intersect"),
        default="full_benchmark",
        help="Whether the P0 pathway cohort is all 111 benchmark controls or only those present in the review CSV.",
    )
    return parser.parse_args(argv)


def is_scorable_pathway_label(series: pd.Series) -> pd.Series:
    normed = series.map(norm_label)
    return normed.ne("") & ~normed.isin(PLACEHOLDER_PATHWAY_LABELS)


def has_known_pathway_label(series: pd.Series) -> pd.Series:
    return is_scorable_pathway_label(series)


def has_known_subpathway_label(series: pd.Series) -> pd.Series:
    return is_scorable_pathway_label(series)


def pct(n_correct: int, n_eligible: int) -> float:
    if n_eligible <= 0:
        return 0.0
    return round(100.0 * float(n_correct) / float(n_eligible), 4)


def summary_row(
    *,
    metric_group: str,
    cohort: str,
    accuracy_target: str,
    rank_view: str,
    lib_id: int | str,
    ion_mode: str,
    breakdown: str,
    eligible_scope: str,
    n_eligible: int,
    n_correct: int,
    notes: str,
) -> dict[str, Any]:
    lib_text = "" if breakdown != "per_library" else str(int(lib_id))
    return {
        "metric_group": metric_group,
        "cohort": cohort,
        "accuracy_target": accuracy_target,
        "rank_view": rank_view,
        "lib_id": lib_text,
        "ion_mode": ion_mode if breakdown == "per_library" else breakdown.replace("_libs", ""),
        "breakdown": breakdown,
        "eligible_scope": eligible_scope,
        "n_eligible": int(n_eligible),
        "n_correct": int(n_correct),
        "pct_correct": pct(n_correct, n_eligible),
        "notes": notes,
    }


def rollup_rows(
    per_lib_rows: list[dict[str, Any]],
    *,
    metric_group: str,
    cohort: str,
    accuracy_target: str,
    rank_view: str,
    eligible_scope: str,
    notes: str,
) -> list[dict[str, Any]]:
    lib_rows = [
        r
        for r in per_lib_rows
        if r["breakdown"] == "per_library"
        and r["rank_view"] == rank_view
        and r["eligible_scope"] == eligible_scope
        and (accuracy_target == "" or r["accuracy_target"] == accuracy_target)
    ]
    out: list[dict[str, Any]] = []
    for breakdown, ion_mode, lib_ids in ROLLUP_SPECS:
        subset = [r for r in lib_rows if int(r["lib_id"]) in lib_ids]
        if not subset:
            continue
        n_eligible = int(sum(r["n_eligible"] for r in subset))
        n_correct = int(sum(r["n_correct"] for r in subset))
        out.append(
            summary_row(
                metric_group=metric_group,
                cohort=cohort,
                accuracy_target=accuracy_target,
                rank_view=rank_view,
                lib_id="",
                ion_mode=ion_mode,
                breakdown=breakdown,
                eligible_scope=eligible_scope,
                n_eligible=n_eligible,
                n_correct=n_correct,
                notes=notes,
            )
        )
    return out


def load_p0_lib_counts(root: Path, tau: float) -> tuple[int, int, int, str]:
    summary_path = root / "integrated_audit_summary.csv"
    if summary_path.exists():
        row = pd.read_csv(summary_path, low_memory=False).iloc[0]
        notes = f"integrated_audit_summary.csv (authoritative gate counts)"
        return (
            int(row["n_controls"]),
            int(row["rank1_ge_0p7"]),
            int(row["top10_ge_0p7"]),
            notes,
        )
    df = pd.read_csv(root / "integrated_per_control_best.csv", low_memory=False)
    notes = f"integrated_per_control_best.csv (tau>={tau})"
    return (
        len(df),
        int(pd.to_numeric(df["rank1_best_tanimoto"], errors="coerce").ge(tau).sum()),
        int(pd.to_numeric(df["top10_best_tanimoto"], errors="coerce").ge(tau).sum()),
        notes,
    )


def load_p0_structure_rows(
    tau: float,
    *,
    deliverable_csv: Path | None = None,
    s0_cohort_mode: str = "full_benchmark",
) -> list[dict[str, Any]]:
    if s0_cohort_mode == "deliverable_intersect":
        if deliverable_csv is None:
            raise ValueError("deliverable_csv is required for deliverable_intersect S0 structure rows")
        return load_p0_structure_rows_from_deliverable(deliverable_csv, tau)

    per_lib: list[dict[str, Any]] = []
    for lib_id, ion_mode, scope, root in LIB_CONFIG:
        n_eligible, rank1_hits, top10_hits, notes = load_p0_lib_counts(root, tau)
        version = "v19b" if "v19b" in str(root) else "v15"
        notes = f"{notes} ({version})"
        per_lib.append(
            summary_row(
                metric_group="P0_structure",
                cohort="p0",
                accuracy_target="structure",
                rank_view="rank1",
                lib_id=lib_id,
                ion_mode=ion_mode,
                breakdown="per_library",
                eligible_scope="s0_benchmark_controls",
                n_eligible=n_eligible,
                n_correct=rank1_hits,
                notes=notes,
            )
        )
        per_lib.append(
            summary_row(
                metric_group="P0_structure",
                cohort="p0",
                accuracy_target="structure",
                rank_view="top10",
                lib_id=lib_id,
                ion_mode=ion_mode,
                breakdown="per_library",
                eligible_scope="s0_benchmark_controls",
                n_eligible=n_eligible,
                n_correct=top10_hits,
                notes=notes,
            )
        )
    rows = list(per_lib)
    rows.extend(
        rollup_rows(
            per_lib,
            metric_group="P0_structure",
            cohort="p0",
            accuracy_target="structure",
            rank_view="rank1",
            eligible_scope="s0_benchmark_controls",
            notes="S0 integrated audit rollup",
        )
    )
    rows.extend(
        rollup_rows(
            per_lib,
            metric_group="P0_structure",
            cohort="p0",
            accuracy_target="structure",
            rank_view="top10",
            eligible_scope="s0_benchmark_controls",
            notes="S0 integrated audit rollup",
        )
    )
    return rows


def load_p0_structure_rows_from_deliverable(deliverable_csv: Path, tau: float) -> list[dict[str, Any]]:
    df = pd.read_csv(deliverable_csv, low_memory=False)
    df["rank"] = pd.to_numeric(df["rank"], errors="coerce")
    df["library_id"] = pd.to_numeric(df["library_id"], errors="coerce")
    df["has_control_reference"] = bool_series(df["has_control_reference"])
    df["control_similarity"] = pd.to_numeric(df["control_similarity"], errors="coerce").fillna(0.0)

    per_lib: list[dict[str, Any]] = []
    for lib_id, ion_mode, _, _ in LIB_CONFIG:
        lib_df = df.loc[df["library_id"].eq(lib_id) & df["has_control_reference"]].copy()
        controls = lib_df.drop_duplicates(subset=["x_compound_id"])
        n_eligible = int(len(controls))
        rank1_hits = int(lib_df.loc[lib_df["rank"].eq(1), "control_similarity"].ge(tau).sum())
        top10 = (
            lib_df.loc[lib_df["rank"].between(1, 10, inclusive="both")]
            .groupby("x_compound_id", as_index=False)["control_similarity"]
            .max()
        )
        top10_hits = int(top10["control_similarity"].ge(tau).sum())
        notes = f"{deliverable_csv.name} (deliverable-intersect S0 controls; tau>={tau})"
        per_lib.append(
            summary_row(
                metric_group="P0_structure",
                cohort="p0",
                accuracy_target="structure",
                rank_view="rank1",
                lib_id=lib_id,
                ion_mode=ion_mode,
                breakdown="per_library",
                eligible_scope="s0_deliverable_intersect_controls",
                n_eligible=n_eligible,
                n_correct=rank1_hits,
                notes=notes,
            )
        )
        per_lib.append(
            summary_row(
                metric_group="P0_structure",
                cohort="p0",
                accuracy_target="structure",
                rank_view="top10",
                lib_id=lib_id,
                ion_mode=ion_mode,
                breakdown="per_library",
                eligible_scope="s0_deliverable_intersect_controls",
                n_eligible=n_eligible,
                n_correct=top10_hits,
                notes=notes,
            )
        )

    rows = list(per_lib)
    for rank_view in ("rank1", "top10"):
        rows.extend(
            rollup_rows(
                per_lib,
                metric_group="P0_structure",
                cohort="p0",
                accuracy_target="structure",
                rank_view=rank_view,
                eligible_scope="s0_deliverable_intersect_controls",
                notes="Deliverable-intersect S0 structure rollup",
            )
        )
    return rows


P1_STRUCTURE_SCOPE = "p1_with_reference"
P1_TRUTH_SAFE_SUFFIX = "_audit_v2_truth_safe"


def resolve_p1_summary_dir(p1_root: Path, legacy_dirname: str) -> str:
    preferred = legacy_dirname.replace(P1_LEGACY_DIR_SUFFIX, P1_PREFERRED_DIR_SUFFIX)
    for dirname in (preferred, legacy_dirname):
        if (p1_root / dirname / "p1_summary.json").exists():
            return dirname
    return preferred


def p1_summary_rows_for_scope(
    *,
    path: Path,
    lib_id: int,
    ion_mode: str,
    eligible_scope: str,
    scope_note: str,
) -> list[dict[str, Any]]:
    summary = json.loads(path.read_text())
    correctness = summary.get("correctness", {})
    if eligible_scope == "p1_with_reference":
        n_eligible = int(summary.get("p1_with_reference", 0))
    else:
        n_eligible = int(correctness.get("p1_with_pool", 0))
    correct_rank1 = int(correctness.get("correct_rank1", 0))
    correct_top10 = int(correctness.get("correct_top10", 0))
    notes = f"{path} ({scope_note})"
    return [
        summary_row(
            metric_group="P1_structure",
            cohort="p1",
            accuracy_target="structure",
            rank_view="rank1",
            lib_id=lib_id,
            ion_mode=ion_mode,
            breakdown="per_library",
            eligible_scope=eligible_scope,
            n_eligible=n_eligible,
            n_correct=correct_rank1,
            notes=notes,
        ),
        summary_row(
            metric_group="P1_structure",
            cohort="p1",
            accuracy_target="structure",
            rank_view="top10",
            lib_id=lib_id,
            ion_mode=ion_mode,
            breakdown="per_library",
            eligible_scope=eligible_scope,
            n_eligible=n_eligible,
            n_correct=correct_top10,
            notes=notes,
        ),
    ]


def load_p1_structure_rows(p1_root: Path) -> list[dict[str, Any]]:
    per_lib: list[dict[str, Any]] = []
    for dirname, lib_id, ion_mode in P1_LIB_DIRS:
        dirname = resolve_p1_summary_dir(p1_root, dirname)
        path = p1_root / dirname / "p1_summary.json"
        per_lib.extend(
            p1_summary_rows_for_scope(
                path=path,
                lib_id=lib_id,
                ion_mode=ion_mode,
                eligible_scope=P1_STRUCTURE_SCOPE,
                scope_note="P1 compounds with verified or resolved-2C reference SMILES; without pool count as not correct; verified_plus_2c_resolved_v2",
            )
        )
        per_lib.extend(
            p1_summary_rows_for_scope(
                path=path,
                lib_id=lib_id,
                ion_mode=ion_mode,
                eligible_scope="p1_with_pool",
                scope_note="P1 compounds with reference SMILES and non-empty candidate pool; verified_plus_2c_resolved_v2",
            )
        )
        safe_path = p1_root / f"{dirname}{P1_TRUTH_SAFE_SUFFIX}" / "p1_summary.json"
        if safe_path.exists():
            per_lib.extend(
                p1_summary_rows_for_scope(
                    path=safe_path,
                    lib_id=lib_id,
                    ion_mode=ion_mode,
                    eligible_scope="p1_with_pool_and_truth_safe",
                    scope_note="P1 compounds with non-empty pool using verified plus safe/review 2C truth only",
                )
            )
    rows = list(per_lib)
    for eligible_scope in sorted({row["eligible_scope"] for row in per_lib}):
        rollup_note = {
            "p1_with_reference": "P1 rollup (sum correct / sum with reference SMILES)",
            "p1_with_pool": "P1 rollup (sum correct / sum with non-empty candidate pool)",
            "p1_with_pool_and_truth_safe": "P1 truth-safe rollup (safe/review 2C only; non-empty candidate pool)",
        }.get(eligible_scope, "P1 rollup")
        rows.extend(
            rollup_rows(
                per_lib,
                metric_group="P1_structure",
                cohort="p1",
                accuracy_target="structure",
                rank_view="rank1",
                eligible_scope=eligible_scope,
                notes=rollup_note,
            )
        )
        rows.extend(
            rollup_rows(
                per_lib,
                metric_group="P1_structure",
                cohort="p1",
                accuracy_target="structure",
                rank_view="top10",
                eligible_scope=eligible_scope,
                notes=rollup_note,
            )
        )
    return rows


def load_all_with_reference_rows(deliverable_csv: Path, tau: float) -> list[dict[str, Any]]:
    df = pd.read_csv(deliverable_csv, low_memory=False)
    df["rank"] = pd.to_numeric(df["rank"], errors="coerce")
    df["library_id"] = pd.to_numeric(df["library_id"], errors="coerce")
    df["has_control_reference"] = bool_series(df["has_control_reference"])
    df["control_similarity"] = pd.to_numeric(df["control_similarity"], errors="coerce").fillna(0.0)

    per_lib: list[dict[str, Any]] = []
    for lib_id, ion_mode, _, _ in LIB_CONFIG:
        lib_df = df.loc[df["library_id"].eq(lib_id) & df["has_control_reference"]]
        controls = lib_df.drop_duplicates(subset=["x_compound_id"])
        n_eligible = len(controls)
        rank1 = lib_df.loc[lib_df["rank"].eq(1)]
        rank1_hits = int(rank1["control_similarity"].ge(tau).sum())
        top10 = (
            lib_df.loc[lib_df["rank"].between(1, 10, inclusive="both")]
            .groupby("x_compound_id", as_index=False)["control_similarity"]
            .max()
        )
        top10_hits = int(top10["control_similarity"].ge(tau).sum())
        notes = f"{deliverable_csv.name} (all compounds with P0 reference)"
        per_lib.append(
            summary_row(
                metric_group="P0_structure",
                cohort="all_with_reference",
                accuracy_target="structure",
                rank_view="rank1",
                lib_id=lib_id,
                ion_mode=ion_mode,
                breakdown="per_library",
                eligible_scope="compounds_with_p0_reference",
                n_eligible=n_eligible,
                n_correct=rank1_hits,
                notes=notes,
            )
        )
        per_lib.append(
            summary_row(
                metric_group="P0_structure",
                cohort="all_with_reference",
                accuracy_target="structure",
                rank_view="top10",
                lib_id=lib_id,
                ion_mode=ion_mode,
                breakdown="per_library",
                eligible_scope="compounds_with_p0_reference",
                n_eligible=n_eligible,
                n_correct=top10_hits,
                notes=notes,
            )
        )
    rows = list(per_lib)
    rows.extend(
        rollup_rows(
            per_lib,
            metric_group="P0_structure",
            cohort="all_with_reference",
            accuracy_target="structure",
            rank_view="rank1",
            eligible_scope="compounds_with_p0_reference",
            notes="Deliverable rollup",
        )
    )
    rows.extend(
        rollup_rows(
            per_lib,
            metric_group="P0_structure",
            cohort="all_with_reference",
            accuracy_target="structure",
            rank_view="top10",
            eligible_scope="compounds_with_p0_reference",
            notes="Deliverable rollup",
        )
    )
    return rows


def load_x_compound_coverage_rows(deliverable_csv: Path) -> list[dict[str, Any]]:
    df = pd.read_csv(deliverable_csv, low_memory=False)
    df["rank"] = pd.to_numeric(df["rank"], errors="coerce")
    df["library_id"] = pd.to_numeric(df["library_id"], errors="coerce")

    per_lib: list[dict[str, Any]] = []
    for lib_id, ion_mode, _, _ in LIB_CONFIG:
        lib_df = df.loc[df["library_id"].eq(lib_id)].copy()
        controls = lib_df.drop_duplicates(subset=["x_compound_id"])
        n_eligible = int(len(controls))
        n_ranked = int(
            lib_df.loc[lib_df["rank"].between(1, 10, inclusive="both"), "x_compound_id"].nunique()
        )
        notes = f"{deliverable_csv.name} (all available X compounds in the full v20 list)"
        per_lib.append(
            summary_row(
                metric_group="X_compound_coverage",
                cohort="all_available",
                accuracy_target="ranked_candidate_rows",
                rank_view="top10",
                lib_id=lib_id,
                ion_mode=ion_mode,
                breakdown="per_library",
                eligible_scope="full_list_x_compounds",
                n_eligible=n_eligible,
                n_correct=n_ranked,
                notes=notes,
            )
        )
    rows = list(per_lib)
    rows.extend(
        rollup_rows(
            per_lib,
            metric_group="X_compound_coverage",
            cohort="all_available",
            accuracy_target="ranked_candidate_rows",
            rank_view="top10",
            eligible_scope="full_list_x_compounds",
            notes="Full-list X-compound coverage rollup",
        )
    )
    return rows


def load_s0_benchmark_pathway_truth(known_x: Path) -> pd.DataFrame:
    """S0 integrated-audit controls with pathway truth from known_X."""
    rows: list[dict[str, Any]] = []
    kx = pd.read_csv(known_x, low_memory=False)
    required = {"LIB_ID", "Unknown", "SuperPathway", "Draft Subpathway", "Level", "SubLevel"}
    missing = sorted(required.difference(kx.columns))
    if missing:
        raise ValueError(f"{known_x} missing columns: {missing}")
    kx = kx.rename(columns={"LIB_ID": "library_id", "Unknown": "control_id"})
    kx["library_id"] = pd.to_numeric(kx["library_id"], errors="coerce")
    kx["control_id"] = kx["control_id"].map(clean_text)
    kx["true_pathway"] = kx["SuperPathway"].map(clean_text)
    kx["true_subpathway"] = kx["Draft Subpathway"].map(clean_text)
    kx["known_level"] = kx["Level"].map(normalize_level)
    kx["known_sublevel"] = kx["SubLevel"].map(clean_text)
    kx_lookup = kx.drop_duplicates(subset=["library_id", "control_id"], keep="first")

    for lib_id, ion_mode, scope, root in LIB_CONFIG:
        audit = pd.read_csv(root / "integrated_per_control_best.csv", low_memory=False)
        for control_id in audit["control_id"].astype(str):
            rows.append(
                {
                    "scope": scope,
                    "library_id": int(lib_id),
                    "ion_mode": ion_mode,
                    "control_id": clean_text(control_id),
                }
            )
    if len(rows) != 111:
        raise ValueError(f"Expected 111 S0 benchmark controls, found {len(rows)}")
    out = pd.DataFrame(rows)
    out = out.merge(
        kx_lookup[
            ["library_id", "control_id", "true_pathway", "true_subpathway", "known_level", "known_sublevel"]
        ],
        on=["library_id", "control_id"],
        how="left",
    )
    out["is_s0_benchmark"] = True
    return out


def control_table_from_review(
    review: pd.DataFrame,
    known_x: Path,
    *,
    s0_cohort_mode: str = "full_benchmark",
) -> pd.DataFrame:
    cols = [
        "scope",
        "library_id",
        "control_id",
        "has_p0_reference",
        "known_level",
        "known_sublevel",
        "true_pathway",
        "true_subpathway",
    ]
    review_ctrl = review[cols].drop_duplicates(subset=["scope", "control_id"], keep="first").copy()
    review_ctrl["has_p0_reference"] = bool_series(review_ctrl["has_p0_reference"])
    s0_truth = load_s0_benchmark_pathway_truth(known_x)

    merge_how = "outer" if s0_cohort_mode == "full_benchmark" else "left"
    out = review_ctrl.merge(
        s0_truth,
        on=["scope", "library_id", "control_id"],
        how=merge_how,
        suffixes=("_review", "_s0"),
    )
    out["library_id"] = pd.to_numeric(out["library_id"], errors="coerce")
    out["control_id"] = out["control_id"].map(clean_text)
    out["has_p0_reference"] = out["has_p0_reference"].fillna(False).astype(bool)
    out["is_s0_benchmark"] = out["is_s0_benchmark"].fillna(False).astype(bool)
    out["is_p0"] = out["is_s0_benchmark"]
    out["known_level"] = out["known_level_s0"].combine_first(out["known_level_review"])
    out["known_sublevel"] = out["known_sublevel_s0"].combine_first(out["known_sublevel_review"])
    out["true_pathway"] = out.apply(
        lambda r: clean_text(r["true_pathway_s0"])
        if bool(r["is_s0_benchmark"])
        else clean_text(r.get("true_pathway_review", r.get("true_pathway", ""))),
        axis=1,
    )
    out["true_subpathway"] = out.apply(
        lambda r: clean_text(r["true_subpathway_s0"])
        if bool(r["is_s0_benchmark"])
        else clean_text(r.get("true_subpathway_review", r.get("true_subpathway", ""))),
        axis=1,
    )
    out["known_level_norm"] = out["known_level"].map(normalize_level)
    out["is_p1"] = out["known_level_norm"].isin({"1", "2"})
    out["has_known_pathway"] = has_known_pathway_label(out["true_pathway"])
    out["has_known_subpathway"] = has_known_subpathway_label(out["true_subpathway"])
    return out[
        [
            "scope",
            "library_id",
            "control_id",
            "has_p0_reference",
            "is_s0_benchmark",
            "known_level",
            "known_sublevel",
            "known_level_norm",
            "is_p0",
            "is_p1",
            "true_pathway",
            "true_subpathway",
            "has_known_pathway",
            "has_known_subpathway",
        ]
    ]


def cohort_controls(controls: pd.DataFrame, cohort: str) -> pd.DataFrame:
    if cohort == "p0":
        return controls.loc[controls["is_p0"]].copy()
    if cohort == "p1":
        return controls.loc[controls["is_p1"]].copy()
    if cohort == "all":
        return controls.copy()
    raise ValueError(f"Unsupported cohort: {cohort}")


def rank_matches(review: pd.DataFrame, scope: str, rank_view: str) -> pd.DataFrame:
    scoped = review.loc[review["scope"].astype(str).eq(scope)].copy()
    scoped["pathway_match"] = bool_series(scoped["pathway_match"])
    scoped["subpathway_match"] = bool_series(scoped["subpathway_match"])
    scoped["both_match"] = scoped["pathway_match"] & scoped["subpathway_match"]
    if rank_view == "rank1":
        rows = scoped.loc[pd.to_numeric(scoped["rank"], errors="coerce").eq(1)]
        return rows[["scope", "control_id", "pathway_match", "subpathway_match", "both_match"]].copy()
    if rank_view == "top10":
        rows = scoped.loc[pd.to_numeric(scoped["rank"], errors="coerce").between(1, 10, inclusive="both")]
        return (
            rows.groupby(["scope", "control_id"], as_index=False)
            .agg(
                pathway_match=("pathway_match", "max"),
                subpathway_match=("subpathway_match", "max"),
                both_match=("both_match", "max"),
            )
            .copy()
        )
    raise ValueError(rank_view)


def pathway_metric_rows(
    controls: pd.DataFrame,
    matches: pd.DataFrame,
    *,
    scope: str,
    lib_id: int,
    ion_mode: str,
    cohort: str,
    rank_view: str,
    p0_eligible_scope: str,
) -> list[dict[str, Any]]:
    cohort_df = cohort_controls(controls.loc[controls["scope"].eq(scope)], cohort)
    merged = cohort_df.merge(matches, on=["scope", "control_id"], how="left")
    for col in ("pathway_match", "subpathway_match", "both_match"):
        merged[col] = merged[col].fillna(False).astype(bool)

    rows: list[dict[str, Any]] = []
    if cohort == "p0":
        specs = (
            ("pathway", merged, "pathway_match"),
            ("subpathway", merged, "subpathway_match"),
            ("pathway_and_subpathway", merged, "both_match"),
        )
    else:
        specs = (
            ("pathway", merged.loc[merged["has_known_pathway"]], "pathway_match"),
            ("subpathway", merged.loc[merged["has_known_subpathway"]], "subpathway_match"),
            (
                "pathway_and_subpathway",
                merged.loc[merged["has_known_pathway"] & merged["has_known_subpathway"]],
                "both_match",
            ),
        )
    scope_by_cohort = {
        "p0": p0_eligible_scope,
        "p1": "p1_cohort_scorable_pathway",
        "all": "all_controls_scorable_pathway",
    }
    subpathway_scope_by_cohort = {
        "p0": p0_eligible_scope,
        "p1": "p1_cohort_scorable_subpathway",
        "all": "all_controls_scorable_subpathway",
    }
    both_scope_by_cohort = {
        "p0": p0_eligible_scope,
        "p1": "p1_cohort_scorable_pathway_and_subpathway",
        "all": "all_controls_scorable_pathway_and_subpathway",
    }
    scope_for_target = {
        "pathway": scope_by_cohort,
        "subpathway": subpathway_scope_by_cohort,
        "pathway_and_subpathway": both_scope_by_cohort,
    }
    for accuracy_target, denom_df, match_col in specs:
        eligible_scope = scope_for_target[accuracy_target][cohort]
        notes = (
            f"v19b_all_libraries_top10_predictions_with_p0_flags.csv "
            f"({cohort}, {rank_view}; {eligible_scope}; "
            f"p0 cohort={p0_eligible_scope}; placeholder pathway labels count as incorrect)"
        )
        n_eligible = int(denom_df["control_id"].nunique())
        n_correct = int(denom_df[match_col].sum()) if n_eligible else 0
        rows.append(
            summary_row(
                metric_group="pathway",
                cohort=cohort,
                accuracy_target=accuracy_target,
                rank_view=rank_view,
                lib_id=lib_id,
                ion_mode=ion_mode,
                breakdown="per_library",
                eligible_scope=eligible_scope,
                n_eligible=n_eligible,
                n_correct=n_correct,
                notes=notes,
            )
        )
    return rows


def load_pathway_rows(
    review_csv: Path,
    known_x: Path,
    *,
    s0_cohort_mode: str = "full_benchmark",
) -> list[dict[str, Any]]:
    review = pd.read_csv(review_csv, low_memory=False)
    required = {
        "scope",
        "library_id",
        "control_id",
        "rank",
        "pathway_match",
        "subpathway_match",
        "true_pathway",
        "true_subpathway",
        "has_p0_reference",
        "known_level",
    }
    missing = sorted(required.difference(review.columns))
    if missing:
        raise ValueError(f"Missing review columns: {', '.join(missing)}")

    controls = control_table_from_review(review, known_x, s0_cohort_mode=s0_cohort_mode)
    p0_eligible_scope = (
        "s0_deliverable_intersect_controls"
        if s0_cohort_mode == "deliverable_intersect"
        else "s0_benchmark_controls"
    )
    per_lib: list[dict[str, Any]] = []
    for lib_id, ion_mode, scope, _ in LIB_CONFIG:
        for cohort in ("p0", "p1", "all"):
            for rank_view in ("rank1", "top10"):
                matches = rank_matches(review, scope, rank_view)
                per_lib.extend(
                    pathway_metric_rows(
                        controls,
                        matches,
                        scope=scope,
                        lib_id=lib_id,
                        ion_mode=ion_mode,
                        cohort=cohort,
                        rank_view=rank_view,
                        p0_eligible_scope=p0_eligible_scope,
                    )
                )

    rows = list(per_lib)
    for cohort in ("p0", "p1", "all"):
        for rank_view in ("rank1", "top10"):
            for accuracy_target in ("pathway", "subpathway", "pathway_and_subpathway"):
                subset = [r for r in per_lib if r["cohort"] == cohort and r["accuracy_target"] == accuracy_target]
                scope_for_target = {
                    "pathway": {
                        "p0": p0_eligible_scope,
                        "p1": "p1_cohort_scorable_pathway",
                        "all": "all_controls_scorable_pathway",
                    },
                    "subpathway": {
                        "p0": p0_eligible_scope,
                        "p1": "p1_cohort_scorable_subpathway",
                        "all": "all_controls_scorable_subpathway",
                    },
                    "pathway_and_subpathway": {
                        "p0": p0_eligible_scope,
                        "p1": "p1_cohort_scorable_pathway_and_subpathway",
                        "all": "all_controls_scorable_pathway_and_subpathway",
                    },
                }
                rows.extend(
                    rollup_rows(
                        subset,
                        metric_group="pathway",
                        cohort=cohort,
                        accuracy_target=accuracy_target,
                        rank_view=rank_view,
                        eligible_scope=scope_for_target[accuracy_target][cohort],
                        notes="Pathway rollup",
                    )
                )
    return rows


def build_accuracy_summary(args: argparse.Namespace) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    rows.extend(
        load_p0_structure_rows(
            float(args.tau),
            deliverable_csv=args.deliverable_csv,
            s0_cohort_mode=args.s0_cohort_mode,
        )
    )
    rows.extend(load_p1_structure_rows(args.p1_root))
    rows.extend(load_all_with_reference_rows(args.deliverable_csv, float(args.tau)))
    rows.extend(load_x_compound_coverage_rows(args.deliverable_csv))
    rows.extend(load_pathway_rows(args.review_csv, args.known_x, s0_cohort_mode=args.s0_cohort_mode))
    return pd.DataFrame(rows, columns=SUMMARY_COLUMNS)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    args.output_root.mkdir(parents=True, exist_ok=True)
    summary = build_accuracy_summary(args)
    out_path = args.output_root / "v19b_high_level_accuracy_summary.csv"
    summary.to_csv(out_path, index=False)
    manifest = {
        "pipeline_version": "v19b",
        "tau": float(args.tau),
        "n_rows": int(len(summary)),
        "output_csv": str(out_path),
        "sources": {
            "review_csv": str(args.review_csv),
            "deliverable_csv": str(args.deliverable_csv),
            "p1_root": str(args.p1_root),
            "p0_integrated": [str(root / "integrated_per_control_best.csv") for _, _, _, root in LIB_CONFIG],
        },
    }
    (args.output_root / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
