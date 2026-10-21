#!/usr/bin/env python3
"""Add pathway/subpathway predictions to equalization rerank top-10 CSVs.

Joins Metabolon super/sub pathway labels from the metabolite synonym reference
DB (exact or neighbor-based assignment) onto per-library rerank deliverables.

Pipeline context:
    First enrichment stage after equalization rerank. Reads
    ``unified_rerank_all_top10.csv`` per lib under ``pipeline_equalization/``.

Key I/O:
    Inputs: per-lib top-10 CSV, metabolite synonym SQLite DB, optional assignment cache.
    Outputs: ``*_with_pathways.csv``, ``pathway_assignment_map.csv``, summary JSON.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sqlite3
from collections import Counter
from pathlib import Path
from typing import Any

import pandas as pd

from structural_elucidation.paths import REPO_ROOT, resolve
from structural_elucidation.utils.chemistry import canonicalize_smiles


# ── Defaults and output column schema ────────────────────────────────
DEFAULT_EQUALIZATION_ROOT: Path = Path("models/pipeline_equalization")
DEFAULT_REFERENCE_DB: Path = Path("Metabolite Synonym Database/metabolite_synonyms_merged_with_mtx_giga.sqlite")
PATHWAY_COLUMNS: list[str] = [
    "predicted_pathway",
    "predicted_subpathway",
    "pathway_confidence_tier",
    "pathway_assignment_status",
    "pathway_vote_share",
    "pathway_best_reference_similarity",
]
COMPOUND_PATHWAY_COLUMNS: list[str] = [
    "compound_predicted_pathway",
    "compound_predicted_subpathway",
    "compound_predicted_pathway_raw",
    "compound_pathway_source",
    "compound_pathway_vote_share",
    "compound_pathway_n_voters",
    "compound_pathway_agrees_rank1",
    "compound_pathway_status",
    "compound_pathway_name_boosted",
]
VOTE_POLICIES = ("exact_first", "all_eligible", "exact_if_neighbors_agree")
DEFAULT_VOTE_POLICY = "exact_if_neighbors_agree"
AA_PEPTIDE_LABELS = {"Amino Acid", "Peptide"}
AA_PEPTIDE_MERGED = "Amino Acid / Peptide"
NAME_COLUMNS = ("candidate_name", "query_name", "library_name")
SUBPATHWAY_NAME_BOOST = 0.15
NAME_PRIOR_MIN_COUNT = 4
NAME_PRIOR_MIN_PRECISION = 0.90
NAME_PRIOR_HOLDOUT_FRAC = 0.25
NAME_PRIOR_HOLDOUT_SEED = 0
TOKEN_STOPWORDS = {
    "acid",
    "and",
    "beta",
    "compound",
    "containing",
    "derivative",
    "drug",
    "ester",
    "form",
    "isomer",
    "metabolism",
    "mono",
    "none",
    "other",
    "pathway",
    "with",
}
MASS_SCORE_FLOOR = 0.90
MASS_ERROR_DA_MAX = 0.02
NUMERIC_MASS_COLUMNS = ("formula_prior_score", "mass_match_score", "abs_mass_error_da")
MASS_FLAG_COLUMNS = ("strict_mass_pass", "gate_mass_within_1da")
CANDIDATE_SCORE_COLUMNS = ("rerank_score_calibrated", "rerank_score", "unified_final_score")
EXACT_STATUS = "exact_reference_match"
NO_PATHWAY_STATUSES = {"", "no_reference_match", "missing_pathway_annotation"}
EMPTY_COMPOUND_ASSIGNMENT: dict[str, Any] = {
    "compound_predicted_pathway": "",
    "compound_predicted_subpathway": "",
    "compound_predicted_pathway_raw": "",
    "compound_pathway_source": "none",
    "compound_pathway_vote_share": pd.NA,
    "compound_pathway_n_voters": 0,
    "compound_pathway_agrees_rank1": pd.NA,
    "compound_pathway_status": "no_eligible_candidates",
    "compound_pathway_name_boosted": 0,
}


# ── CLI and parsing helpers ──────────────────────────────────────────
def clean_text(value: object) -> str:
    """Normalize a cell value to a stripped string, mapping null sentinels to empty."""
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except Exception:
        pass
    text = str(value).strip()
    return "" if text.lower() in {"nan", "none", "null", "<na>"} else text


def parse_args() -> argparse.Namespace:
    """Build the pathway enrichment CLI argument parser."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--equalization-root", type=Path, default=DEFAULT_EQUALIZATION_ROOT)
    parser.add_argument("--reference-db", type=Path, default=DEFAULT_REFERENCE_DB)
    parser.add_argument("--libs", nargs="+", default=["209,305,400,402"])
    parser.add_argument("--input-name", default="unified_rerank_all_top10.csv")
    parser.add_argument("--output-name", default="unified_rerank_all_top10_with_pathways.csv")
    parser.add_argument("--summary-json", type=Path, default=None)
    parser.add_argument("--assignment-cache", type=Path, default=None)
    parser.add_argument("--allow-missing", action="store_true")
    parser.add_argument(
        "--compound-vote-policy",
        choices=list(VOTE_POLICIES),
        default=DEFAULT_VOTE_POLICY,
        help="How mass-consistent exact catalog hits enter the compound vote.",
    )
    parser.add_argument(
        "--collapse-amino-acid-peptide",
        action="store_true",
        default=True,
        help="Tally Amino Acid and Peptide as one super-pathway (raw winner kept in sidecar).",
    )
    parser.add_argument(
        "--no-collapse-amino-acid-peptide",
        action="store_false",
        dest="collapse_amino_acid_peptide",
    )
    parser.add_argument(
        "--subpathway-name-prior",
        action="store_true",
        default=True,
        help="Boost subpathway votes with reference-DB name tokens (held-out controls skipped).",
    )
    parser.add_argument(
        "--no-subpathway-name-prior",
        action="store_false",
        dest="subpathway_name_prior",
    )
    return parser.parse_args()


def parse_libs(values: list[str]) -> list[int]:
    """Parse comma-separated library IDs from CLI nargs list."""
    libs: list[int] = []
    for value in values:
        for token in str(value).split(","):
            token = token.strip()
            if token:
                libs.append(int(token))
    return libs


def canonicalize_for_pathway(value: object) -> str:
    """Canonicalize SMILES for pathway join keys."""
    text = clean_text(value)
    if not text:
        return ""
    return clean_text(canonicalize_smiles(text))


def rerank_path(root: Path, lib_id: int, name: str) -> Path:
    """Path to a named CSV under ``rerank_lib{lib_id}/``."""
    return root / f"rerank_lib{int(lib_id)}" / name


def load_inputs(root: Path, libs: list[int], input_name: str, allow_missing: bool) -> dict[int, pd.DataFrame]:
    """Load per-library rerank top-10 CSVs from the equalization tree."""
    frames: dict[int, pd.DataFrame] = {}
    missing: list[Path] = []
    for lib_id in libs:
        path = rerank_path(root, lib_id, input_name)
        if not path.exists():
            missing.append(path)
            continue
        frame = pd.read_csv(path, low_memory=False)
        if "predicted_smiles" not in frame.columns:
            raise ValueError(f"{path} missing required column: predicted_smiles")
        frames[int(lib_id)] = frame
    if missing and not allow_missing:
        joined = "\n".join(f"  {path}" for path in missing)
        raise FileNotFoundError(f"Missing rerank top10 input CSVs:\n{joined}")
    return frames


# ── Pathway assignment map ─────────────────────────────────────────
def load_assignment_cache(path: Path | None) -> pd.DataFrame:
    """Load cached SMILES→pathway assignments if present."""
    if path is None or not path.exists():
        return normalize_assignment_map(None)
    return normalize_assignment_map(pd.read_csv(path, low_memory=False))


def normalize_assignment_map(assignments: pd.DataFrame | None) -> pd.DataFrame:
    """Normalize assignment-map columns to the pathway enrichment schema."""
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


def load_reference_pathways(reference_db: Path) -> pd.DataFrame:
    """Load canonical SMILES with Metabolon super/sub labels from SQLite."""
    query = """
        SELECT
            smiles,
            metabolon_super,
            metabolon_sub
        FROM metabolites
        WHERE TRIM(COALESCE(smiles, '')) <> ''
          AND TRIM(COALESCE(metabolon_super, '')) <> ''
          AND TRIM(COALESCE(metabolon_sub, '')) <> ''
    """
    connection = sqlite3.connect(str(reference_db))
    try:
        reference = pd.read_sql_query(query, connection)
    finally:
        connection.close()
    if reference.empty:
        return pd.DataFrame(columns=["canonical_smiles", "metabolon_super", "metabolon_sub"])
    reference["canonical_smiles"] = reference["smiles"].map(canonicalize_for_pathway)
    reference["metabolon_super"] = reference["metabolon_super"].map(clean_text)
    reference["metabolon_sub"] = reference["metabolon_sub"].map(clean_text)
    return reference.loc[
        reference["canonical_smiles"].ne("") & reference["metabolon_super"].ne("") & reference["metabolon_sub"].ne(""),
        ["canonical_smiles", "metabolon_super", "metabolon_sub"],
    ].drop_duplicates()


def generate_exact_reference_assignment_map(canonical_smiles: list[str], reference_db: Path) -> pd.DataFrame:
    """Assign pathways by exact canonical-SMILES match in the reference DB."""
    unique_smiles = sorted(set(clean_text(smiles) for smiles in canonical_smiles if clean_text(smiles)))
    if not unique_smiles:
        return normalize_assignment_map(None)
    reference = load_reference_pathways(reference_db)
    rows: list[dict[str, Any]] = []
    for smiles in unique_smiles:
        matches = reference.loc[reference["canonical_smiles"].eq(smiles)]
        if matches.empty:
            rows.append(
                {
                    "canonical_smiles": smiles,
                    "predicted_pathway": "",
                    "predicted_subpathway": "",
                    "pathway_confidence_tier": "low",
                    "pathway_assignment_status": "no_reference_match",
                    "pathway_vote_share": 0.0,
                    "pathway_best_reference_similarity": 0.0,
                }
            )
            continue
        pathway = clean_text(matches["metabolon_super"].mode().iloc[0])
        subpathway = clean_text(matches.loc[matches["metabolon_super"].eq(pathway), "metabolon_sub"].mode().iloc[0])
        rows.append(
            {
                "canonical_smiles": smiles,
                "predicted_pathway": pathway,
                "predicted_subpathway": subpathway,
                "pathway_confidence_tier": "high",
                "pathway_assignment_status": "exact_reference_match",
                "pathway_vote_share": 1.0,
                "pathway_best_reference_similarity": 1.0,
            }
        )
    return normalize_assignment_map(pd.DataFrame(rows))


def generate_pathway_assignment_map(canonical_smiles: list[str], reference_db: Path) -> pd.DataFrame:
    """Generate assignments via neighbor search when legacy builder is importable."""
    try:
        from build_v14_all_library_prediction_review_v1 import (  # type: ignore[import-not-found]
            generate_pathway_assignment_map as generate_with_neighbors,
        )
    except ImportError:
        return generate_exact_reference_assignment_map(canonical_smiles, reference_db)
    return normalize_assignment_map(generate_with_neighbors(canonical_smiles, reference_db))


def enrich_pathway_annotations(df: pd.DataFrame, assignment_map: pd.DataFrame | None) -> pd.DataFrame:
    """Left-join pathway columns onto rerank rows by canonical SMILES."""
    out = df.copy()
    for column in PATHWAY_COLUMNS:
        if column not in out.columns:
            out[column] = pd.NA

    out["_canonical_pathway_join_smiles"] = out.get("predicted_smiles", pd.Series([""] * len(out))).map(
        canonicalize_for_pathway
    )
    assignment_map = normalize_assignment_map(assignment_map)

    if assignment_map.empty:
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
    for column in PATHWAY_COLUMNS:
        assignment_column = f"{column}_assignment"
        if assignment_column in annotated.columns:
            current = annotated[column].map(clean_text)
            annotated[column] = annotated[column].where(current.ne(""), annotated[assignment_column])
    annotated.loc[annotated["pathway_assignment_status"].map(clean_text).eq(""), "pathway_assignment_status"] = (
        "missing_pathway_annotation"
    )
    return annotated.drop(
        columns=["_canonical_pathway_join_smiles", "canonical_smiles"]
        + [f"{column}_assignment" for column in PATHWAY_COLUMNS if f"{column}_assignment" in annotated.columns],
        errors="ignore",
    )


def apply_amino_acid_pathway_overrides(df: pd.DataFrame) -> pd.DataFrame:
    """Apply legacy amino-acid pathway overrides when the builder module exists."""
    try:
        from build_v19b_updated_summary_v1 import (  # type: ignore[import-not-found]
            apply_amino_acid_pathway_overrides as apply_legacy_overrides,
        )
    except ImportError:
        apply_legacy_overrides = None
    if apply_legacy_overrides is not None:
        return apply_legacy_overrides(df)

    out = df.copy()
    true_path = out.get("true_pathway", pd.Series([""] * len(out), index=out.index)).map(clean_text)
    true_sub = out.get("true_subpathway", pd.Series([""] * len(out), index=out.index)).map(clean_text)
    predicted_path = out.get("predicted_pathway", pd.Series([""] * len(out), index=out.index)).map(clean_text)
    predicted_sub = out.get("predicted_subpathway", pd.Series([""] * len(out), index=out.index)).map(clean_text)
    if "pathway_match" in out.columns or true_path.ne("").any():
        out["pathway_match"] = true_path.ne("") & true_path.eq(predicted_path)
    if "subpathway_match" in out.columns or true_sub.ne("").any():
        out["subpathway_match"] = true_sub.ne("") & true_sub.eq(predicted_sub)
    return out


def _finite_number(value: object) -> float | None:
    """Parse a numeric cell, returning None for missing/non-finite values."""
    if value is None:
        return None
    try:
        if pd.isna(value):
            return None
    except Exception:
        pass
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if number != number or number in {float("inf"), float("-inf")}:
        return None
    return number


def _truthy_flag(value: object) -> bool:
    """Interpret a gate/pass column as a boolean."""
    text = clean_text(value).lower()
    if text in {"1", "true", "yes", "t"}:
        return True
    if text in {"", "0", "false", "no", "f"}:
        return False
    number = _finite_number(value)
    return bool(number)


def control_id_column(df: pd.DataFrame) -> str | None:
    """Preferred X-compound identifier column on a rerank frame."""
    for column in ("control_id", "x_compound_id"):
        if column in df.columns:
            return column
    return None


def mass_columns_present(df: pd.DataFrame) -> bool:
    """True when any numeric mass-consistency column exists on the frame."""
    return any(column in df.columns for column in NUMERIC_MASS_COLUMNS)


def is_mass_consistent(row: pd.Series, *, numeric_mass_columns_present: bool) -> bool:
    """Whether a candidate is allowed to vote on the compound pathway."""
    if numeric_mass_columns_present:
        formula = _finite_number(row.get("formula_prior_score")) if "formula_prior_score" in row.index else None
        mass = _finite_number(row.get("mass_match_score")) if "mass_match_score" in row.index else None
        error = _finite_number(row.get("abs_mass_error_da")) if "abs_mass_error_da" in row.index else None
        if formula is not None and formula >= MASS_SCORE_FLOOR:
            return True
        if mass is not None and mass >= MASS_SCORE_FLOOR:
            return True
        if error is not None and error <= MASS_ERROR_DA_MAX:
            return True
        return False
    for column in MASS_FLAG_COLUMNS:
        if column in row.index and _truthy_flag(row.get(column)):
            return True
    return False


def is_pathway_eligible(row: pd.Series, *, numeric_mass_columns_present: bool) -> bool:
    """Eligible voters have a pathway label, a catalog status, and a mass fit."""
    if not clean_text(row.get("predicted_pathway")):
        return False
    if clean_text(row.get("pathway_assignment_status")) in NO_PATHWAY_STATUSES:
        return False
    return is_mass_consistent(row, numeric_mass_columns_present=numeric_mass_columns_present)


def eligible_mask(df: pd.DataFrame, *, numeric_mass_columns_present: bool) -> pd.Series:
    """Vectorized eligibility for the compound-level vote."""
    path_ok = df["predicted_pathway"].map(clean_text).ne("") if "predicted_pathway" in df.columns else pd.Series(False, index=df.index)
    status = df["pathway_assignment_status"].map(clean_text) if "pathway_assignment_status" in df.columns else pd.Series("", index=df.index)
    status_ok = ~status.isin(NO_PATHWAY_STATUSES)
    if numeric_mass_columns_present:
        formula = pd.to_numeric(df["formula_prior_score"], errors="coerce") if "formula_prior_score" in df.columns else pd.Series(float("nan"), index=df.index)
        mass = pd.to_numeric(df["mass_match_score"], errors="coerce") if "mass_match_score" in df.columns else pd.Series(float("nan"), index=df.index)
        error = pd.to_numeric(df["abs_mass_error_da"], errors="coerce") if "abs_mass_error_da" in df.columns else pd.Series(float("nan"), index=df.index)
        mass_ok = ((formula >= MASS_SCORE_FLOOR) | (mass >= MASS_SCORE_FLOOR) | (error <= MASS_ERROR_DA_MAX)).fillna(False)
    else:
        mass_ok = pd.Series(False, index=df.index)
        for column in MASS_FLAG_COLUMNS:
            if column in df.columns:
                mass_ok = mass_ok | df[column].map(_truthy_flag)
    return path_ok & status_ok & mass_ok


def candidate_score_series(df: pd.DataFrame) -> pd.Series:
    """Vectorized candidate strength: calibrated rerank, then raw, then 1/rank."""
    score = pd.Series(0.0, index=df.index, dtype="float64")
    filled = pd.Series(False, index=df.index)
    for column in CANDIDATE_SCORE_COLUMNS:
        if column not in df.columns:
            continue
        values = pd.to_numeric(df[column], errors="coerce")
        use = values.notna() & ~filled
        if bool(use.any()):
            score.loc[use] = values.loc[use]
            filled = filled | use
    if "rank" in df.columns:
        ranks = pd.to_numeric(df["rank"], errors="coerce")
        use = ~filled & ranks.notna() & (ranks > 0)
        if bool(use.any()):
            score.loc[use] = 1.0 / ranks.loc[use]
    return score


def catalog_weight_series(df: pd.DataFrame) -> pd.Series:
    """Vectorized catalog evidence weight for each candidate row."""
    status = (
        df["pathway_assignment_status"].map(clean_text)
        if "pathway_assignment_status" in df.columns
        else pd.Series("", index=df.index)
    )
    weight = pd.Series(0.0, index=df.index, dtype="float64")
    exact = status.eq(EXACT_STATUS)
    weight.loc[exact] = 1.0
    share = (
        pd.to_numeric(df["pathway_vote_share"], errors="coerce")
        if "pathway_vote_share" in df.columns
        else pd.Series(float("nan"), index=df.index)
    )
    use_share = ~exact & share.notna() & (share > 0) & (share <= 1)
    weight.loc[use_share] = share.loc[use_share]
    similarity = (
        pd.to_numeric(df["pathway_best_reference_similarity"], errors="coerce")
        if "pathway_best_reference_similarity" in df.columns
        else pd.Series(0.0, index=df.index)
    ).fillna(0.0).clip(lower=0)
    rest = ~exact & ~use_share
    weight.loc[rest] = similarity.loc[rest] ** 2
    return weight


def candidate_score(row: pd.Series) -> float:
    """Relative candidate strength: calibrated rerank, then raw, then 1/rank."""
    return float(candidate_score_series(row.to_frame().T).iloc[0])


def catalog_weight(row: pd.Series) -> float:
    """How strongly this structure's catalog evidence should count."""
    return float(catalog_weight_series(row.to_frame().T).iloc[0])


def _ranked_tally(weights: dict[str, float], first_seen: dict[str, int]) -> list[tuple[str, float]]:
    """Sort a tally by weight; ties keep the class whose first voter ranked higher."""
    return sorted(weights.items(), key=lambda item: (-item[1], first_seen[item[0]]))


def collapse_super_label(name: str, collapse_amino_acid_peptide: bool) -> str:
    """Map Amino Acid and Peptide onto one tally bucket when requested."""
    if collapse_amino_acid_peptide and name in AA_PEPTIDE_LABELS:
        return AA_PEPTIDE_MERGED
    return name


def supers_agree(left: str, right: str, collapse_amino_acid_peptide: bool) -> bool:
    """True when two super-pathway labels count as the same class."""
    a = collapse_super_label(clean_text(left), collapse_amino_acid_peptide)
    b = collapse_super_label(clean_text(right), collapse_amino_acid_peptide)
    return bool(a) and a == b


def tokenize_name(text: str) -> set[str]:
    """Tokens used by the subpathway name prior (copied from the Spec2Mol assigner)."""
    normalized = clean_text(text).lower()
    if not normalized:
        return set()
    tokens: set[str] = set()
    for token in re.findall(r"[a-z0-9][a-z0-9'\-/]{2,}", normalized):
        token = token.strip("-_/'")
        if not token or token in TOKEN_STOPWORDS or token.isdigit():
            continue
        tokens.add(token)
    return tokens


_NAME_PRIOR_CACHE: dict[str, dict[tuple[str, str], dict[str, float | str]]] = {}


def load_subpathway_name_priors(reference_db: Path) -> dict[tuple[str, str], dict[str, float | str]]:
    """High-precision name-token → subpathway map from the reference catalog."""
    cache_key = str(Path(reference_db).resolve()) if Path(reference_db).exists() else ""
    if cache_key and cache_key in _NAME_PRIOR_CACHE:
        return _NAME_PRIOR_CACHE[cache_key]
    if not reference_db.exists():
        return {}
    query = """
        SELECT canonical_name, metabolon_super, metabolon_sub
        FROM metabolites
        WHERE TRIM(COALESCE(canonical_name, '')) <> ''
          AND TRIM(COALESCE(metabolon_super, '')) <> ''
          AND TRIM(COALESCE(metabolon_sub, '')) <> ''
    """
    connection = sqlite3.connect(str(reference_db))
    try:
        reference = pd.read_sql_query(query, connection)
    finally:
        connection.close()
    if reference.empty:
        return {}
    token_counts: dict[tuple[str, str], Counter[str]] = {}
    for row in reference.itertuples(index=False):
        pathway = clean_text(row.metabolon_super)
        subpathway = clean_text(row.metabolon_sub)
        if not pathway or not subpathway:
            continue
        for token in tokenize_name(clean_text(row.canonical_name)):
            key = (pathway, token)
            token_counts.setdefault(key, Counter())[subpathway] += 1
    prior_map: dict[tuple[str, str], dict[str, float | str]] = {}
    for key, counter in token_counts.items():
        total = int(sum(counter.values()))
        if total < NAME_PRIOR_MIN_COUNT:
            continue
        subpathway, count = sorted(counter.items(), key=lambda item: (-item[1], item[0]))[0]
        precision = float(count) / float(total)
        if precision < NAME_PRIOR_MIN_PRECISION:
            continue
        prior_map[key] = {"subpathway": subpathway, "precision": precision, "count": float(total)}
    if cache_key:
        _NAME_PRIOR_CACHE[cache_key] = prior_map
    return prior_map


def in_name_prior_holdout(
    control_id: object,
    frac: float = NAME_PRIOR_HOLDOUT_FRAC,
    seed: int = NAME_PRIOR_HOLDOUT_SEED,
) -> bool:
    """Deterministic hold-out so name-boost evaluation is not fit on every known_X id."""
    text = clean_text(control_id)
    if not text or frac <= 0:
        return False
    digest = hashlib.md5(f"{int(seed)}:{text}".encode("utf-8")).hexdigest()
    return (int(digest[:8], 16) / 0xFFFFFFFF) < float(frac)


def _row_query_name(row: dict[str, Any]) -> str:
    for column in NAME_COLUMNS:
        name = clean_text(row.get(column))
        if name:
            return name
    return ""


def _subpathway_name_boosts(
    rows: list[dict[str, Any]],
    winning_super: str,
    priors: dict[tuple[str, str], dict[str, float | str]],
) -> dict[str, float]:
    boosts: dict[str, float] = {}
    if not priors or not winning_super:
        return boosts
    names = {_row_query_name(row) for row in rows}
    names.discard("")
    if not names:
        return boosts
    super_keys = {winning_super}
    if winning_super == AA_PEPTIDE_MERGED:
        super_keys = set(AA_PEPTIDE_LABELS)
    for name in names:
        for token in tokenize_name(name):
            for super_label in super_keys:
                payload = priors.get((super_label, token))
                if not payload:
                    continue
                subpathway = clean_text(payload.get("subpathway", ""))
                if not subpathway:
                    continue
                boosts[subpathway] = boosts.get(subpathway, 0.0) + float(payload.get("precision", 0.0))
    return boosts


def vote_compound_pathway(
    voters: pd.DataFrame,
    *,
    collapse_amino_acid_peptide: bool = False,
    name_priors: dict[tuple[str, str], dict[str, float | str]] | None = None,
    apply_name_boost: bool = False,
) -> dict[str, Any]:
    """Weighted super-then-sub vote over eligible candidate rows."""
    if voters.empty:
        return dict(EMPTY_COMPOUND_ASSIGNMENT)

    scored = voters.copy()
    if "_candidate_score" not in scored.columns:
        scored["_candidate_score"] = candidate_score_series(scored)
    if "_catalog_weight" not in scored.columns:
        scored["_catalog_weight"] = catalog_weight_series(scored)
    max_score = float(scored["_candidate_score"].max())
    if max_score > 0:
        scored["_candidate_rel"] = scored["_candidate_score"] / max_score
    else:
        scored["_candidate_rel"] = 1.0
    scored["_vote_weight"] = scored["_candidate_rel"] * scored["_catalog_weight"]
    scored = scored.sort_values(
        by=["_vote_weight", "rank"] if "rank" in scored.columns else ["_vote_weight"],
        ascending=[False, True] if "rank" in scored.columns else [False],
        kind="mergesort",
    )
    raw_weights: dict[str, float] = {}
    raw_first: dict[str, int] = {}
    super_weights: dict[str, float] = {}
    super_first: dict[str, int] = {}
    rows = scored.to_dict("records")
    for index, row in enumerate(rows):
        super_pathway = clean_text(row.get("predicted_pathway"))
        weight = float(row.get("_vote_weight") or 0.0)
        if not super_pathway or weight <= 0:
            continue
        raw_weights[super_pathway] = raw_weights.get(super_pathway, 0.0) + weight
        raw_first.setdefault(super_pathway, index)
        bucket = collapse_super_label(super_pathway, collapse_amino_acid_peptide)
        super_weights[bucket] = super_weights.get(bucket, 0.0) + weight
        super_first.setdefault(bucket, index)
    if not super_weights:
        return dict(EMPTY_COMPOUND_ASSIGNMENT)

    ranked_super = _ranked_tally(super_weights, super_first)
    winning_bucket, winning_weight = ranked_super[0]
    total_weight = sum(super_weights.values())
    ranked_raw = _ranked_tally(raw_weights, raw_first)
    raw_super = ranked_raw[0][0]
    if winning_bucket == AA_PEPTIDE_MERGED and collapse_amino_acid_peptide:
        displayed_super = max(
            ((name, raw_weights.get(name, 0.0), raw_first.get(name, 10**9)) for name in AA_PEPTIDE_LABELS if name in raw_weights),
            key=lambda item: (item[1], -item[2]),
            default=(raw_super, 0.0, 0),
        )[0]
    else:
        displayed_super = winning_bucket

    def _backs_winner(name: str) -> bool:
        return collapse_super_label(name, collapse_amino_acid_peptide) == winning_bucket

    sub_weights: dict[str, float] = {}
    sub_first: dict[str, int] = {}
    winner_rows = [row for row in rows if _backs_winner(clean_text(row.get("predicted_pathway")))]
    for index, row in enumerate(winner_rows):
        subpathway = clean_text(row.get("predicted_subpathway"))
        weight = float(row.get("_vote_weight") or 0.0)
        if not subpathway or weight <= 0:
            continue
        sub_weights[subpathway] = sub_weights.get(subpathway, 0.0) + weight
        sub_first.setdefault(subpathway, index)

    name_boosted = 0
    if apply_name_boost and name_priors and sub_weights:
        boosts = _subpathway_name_boosts(winner_rows, winning_bucket, name_priors)
        if boosts:
            winner_total = sum(sub_weights.values())
            for subpathway, boost in boosts.items():
                if subpathway not in sub_weights:
                    continue
                capped = min(0.30 * winner_total, float(boost) * SUBPATHWAY_NAME_BOOST * winner_total)
                sub_weights[subpathway] = sub_weights[subpathway] + capped
                name_boosted = 1

    ranked_sub = _ranked_tally(sub_weights, sub_first)
    winning_sub = ranked_sub[0][0] if ranked_sub else ""

    exact = all(clean_text(row.get("pathway_assignment_status")) == EXACT_STATUS for row in rows)
    source = "exact_mass_consistent" if exact else "weighted_candidate_consensus"
    return {
        "compound_predicted_pathway": displayed_super,
        "compound_predicted_subpathway": winning_sub,
        "compound_predicted_pathway_raw": raw_super,
        "compound_pathway_source": source,
        "compound_pathway_vote_share": round(winning_weight / total_weight, 6) if total_weight else pd.NA,
        "compound_pathway_n_voters": int(len(rows)),
        "compound_pathway_agrees_rank1": pd.NA,
        "compound_pathway_status": source,
        "compound_pathway_name_boosted": name_boosted,
    }


def select_compound_voters(
    eligible: pd.DataFrame,
    *,
    policy: str,
    collapse_amino_acid_peptide: bool,
) -> pd.DataFrame:
    """Choose which eligible rows vote, given the exact-hit policy."""
    if eligible.empty or policy == "all_eligible":
        return eligible
    if "pathway_assignment_status" not in eligible.columns:
        return eligible
    exact = eligible.loc[eligible["pathway_assignment_status"].map(clean_text).eq(EXACT_STATUS)]
    if exact.empty or policy == "exact_first":
        return exact if not exact.empty else eligible
    neighbors = eligible.loc[eligible["pathway_assignment_status"].map(clean_text).ne(EXACT_STATUS)]
    if neighbors.empty:
        return exact
    exact_vote = vote_compound_pathway(exact, collapse_amino_acid_peptide=collapse_amino_acid_peptide)
    neighbor_vote = vote_compound_pathway(neighbors, collapse_amino_acid_peptide=collapse_amino_acid_peptide)
    if supers_agree(
        exact_vote["compound_predicted_pathway"],
        neighbor_vote["compound_predicted_pathway"],
        collapse_amino_acid_peptide,
    ):
        return exact
    return eligible


def assign_compound_pathways(
    df: pd.DataFrame,
    *,
    policy: str = DEFAULT_VOTE_POLICY,
    collapse_amino_acid_peptide: bool = True,
    name_priors: dict[tuple[str, str], dict[str, float | str]] | None = None,
    apply_name_boost: bool = True,
    name_prior_holdout_frac: float = NAME_PRIOR_HOLDOUT_FRAC,
    name_prior_holdout_seed: int = NAME_PRIOR_HOLDOUT_SEED,
) -> pd.DataFrame:
    """Broadcast one mass-gated, score-weighted pathway per X compound."""
    if policy not in VOTE_POLICIES:
        raise ValueError(f"Unknown compound vote policy: {policy}")
    out = df.copy()
    empty = dict(EMPTY_COMPOUND_ASSIGNMENT)
    id_col = control_id_column(out)
    if id_col is None:
        for column, value in empty.items():
            out[column] = value
        return out

    numeric_mass = mass_columns_present(out)
    eligible_rows = eligible_mask(out, numeric_mass_columns_present=numeric_mass)
    out["_candidate_score"] = candidate_score_series(out)
    out["_catalog_weight"] = catalog_weight_series(out)
    by_control: dict[Any, dict[str, Any]] = {}
    for control_id, group in out.groupby(id_col, sort=False, dropna=False):
        if not clean_text(control_id):
            by_control[control_id] = dict(empty)
            continue
        eligible = group.loc[eligible_rows.loc[group.index]]
        voters = select_compound_voters(
            eligible,
            policy=policy,
            collapse_amino_acid_peptide=collapse_amino_acid_peptide,
        )
        holdout = in_name_prior_holdout(control_id, name_prior_holdout_frac, name_prior_holdout_seed)
        payload = vote_compound_pathway(
            voters,
            collapse_amino_acid_peptide=collapse_amino_acid_peptide,
            name_priors=name_priors,
            apply_name_boost=bool(apply_name_boost and name_priors and not holdout),
        )
        ranked = group.sort_values("rank", kind="mergesort") if "rank" in group.columns else group
        rank1_pathway = clean_text(ranked.iloc[0].get("predicted_pathway")) if len(ranked) else ""
        assigned = clean_text(payload["compound_predicted_pathway"])
        if assigned:
            payload["compound_pathway_agrees_rank1"] = int(
                supers_agree(assigned, rank1_pathway, collapse_amino_acid_peptide)
            )
        by_control[control_id] = payload

    keys = out[id_col]
    for column in COMPOUND_PATHWAY_COLUMNS:
        out[column] = keys.map(lambda cid, col=column: by_control.get(cid, empty).get(col))
    return out.drop(columns=["_candidate_score", "_catalog_weight"], errors="ignore")


def build_assignment_map(
    *,
    frames: dict[int, pd.DataFrame],
    reference_db: Path,
    assignment_cache: Path | None,
) -> tuple[pd.DataFrame, dict[str, int]]:
    """Build or extend the SMILES→pathway cache for all unique input SMILES."""
    unique_smiles = sorted(
        {
            canonical
            for frame in frames.values()
            for canonical in frame["predicted_smiles"].map(canonicalize_for_pathway)
            if canonical
        }
    )
    cached = load_assignment_cache(assignment_cache)
    cached_smiles = set(cached["canonical_smiles"].map(clean_text)) if not cached.empty else set()
    missing = [smiles for smiles in unique_smiles if smiles not in cached_smiles]
    generated = normalize_assignment_map(None)
    if missing:
        if not reference_db.exists():
            raise FileNotFoundError(reference_db)
        generated = generate_pathway_assignment_map(missing, reference_db)
    assignment_map = (
        pd.concat([cached, generated], ignore_index=True, sort=False)
        .pipe(normalize_assignment_map)
        .drop_duplicates(subset=["canonical_smiles"], keep="first")
        .reset_index(drop=True)
    )
    if assignment_cache is not None:
        assignment_cache.parent.mkdir(parents=True, exist_ok=True)
        assignment_map.to_csv(assignment_cache, index=False)
    stats = {
        "n_unique_input_smiles": int(len(unique_smiles)),
        "n_cached_assignment_rows": int(len(cached)),
        "n_missing_assignment_rows_requested": int(len(missing)),
        "n_generated_assignment_rows": int(len(generated)),
        "n_assignment_rows": int(len(assignment_map)),
    }
    return assignment_map, stats


def ordered_columns(original_columns: list[str], enriched: pd.DataFrame) -> list[str]:
    """Preserve original column order while appending new pathway columns."""
    columns = list(original_columns)
    for column in PATHWAY_COLUMNS + COMPOUND_PATHWAY_COLUMNS:
        if column not in columns and column in enriched.columns:
            columns.append(column)
    return columns + [column for column in enriched.columns if column not in columns and column not in {"pathway_match", "subpathway_match"}]


def enrich_one(
    *,
    frame: pd.DataFrame,
    assignment_map: pd.DataFrame,
    policy: str = DEFAULT_VOTE_POLICY,
    collapse_amino_acid_peptide: bool = True,
    name_priors: dict[tuple[str, str], dict[str, float | str]] | None = None,
    apply_name_boost: bool = True,
) -> pd.DataFrame:
    """Enrich one library frame and restore original column ordering."""
    original_columns = list(frame.columns)
    out = enrich_pathway_annotations(frame, assignment_map=assignment_map)
    out = apply_amino_acid_pathway_overrides(out)
    out = assign_compound_pathways(
        out,
        policy=policy,
        collapse_amino_acid_peptide=collapse_amino_acid_peptide,
        name_priors=name_priors,
        apply_name_boost=apply_name_boost,
    )
    if "pathway_match" not in original_columns:
        out = out.drop(columns=["pathway_match"], errors="ignore")
    if "subpathway_match" not in original_columns:
        out = out.drop(columns=["subpathway_match"], errors="ignore")
    return out.loc[:, ordered_columns(original_columns, out)]


def summarize_lib(lib_id: int, input_path: Path, output_path: Path, frame: pd.DataFrame) -> dict[str, Any]:
    """Summarize pathway assignment coverage for one library output."""
    status = frame["pathway_assignment_status"].map(clean_text) if "pathway_assignment_status" in frame.columns else pd.Series([], dtype=object)
    tier = frame["pathway_confidence_tier"].map(clean_text) if "pathway_confidence_tier" in frame.columns else pd.Series([], dtype=object)
    path = frame["predicted_pathway"].map(clean_text) if "predicted_pathway" in frame.columns else pd.Series([], dtype=object)
    sub = frame["predicted_subpathway"].map(clean_text) if "predicted_subpathway" in frame.columns else pd.Series([], dtype=object)
    smiles = frame["predicted_smiles"].map(clean_text) if "predicted_smiles" in frame.columns else pd.Series([], dtype=object)
    id_col = control_id_column(frame)
    n_controls = 0
    n_controls_with_compound_pathway = 0
    n_rank1_agree = 0
    compound_source_counts: dict[str, int] = {}
    if id_col and not frame.empty:
        per_control = frame.drop_duplicates(subset=[id_col], keep="first")
        n_controls = int(len(per_control))
        if "compound_predicted_pathway" in per_control.columns:
            assigned = per_control["compound_predicted_pathway"].map(clean_text).ne("")
            n_controls_with_compound_pathway = int(assigned.sum())
            if "compound_pathway_agrees_rank1" in per_control.columns:
                agrees = pd.to_numeric(per_control["compound_pathway_agrees_rank1"], errors="coerce")
                n_rank1_agree = int(agrees.loc[assigned].fillna(0).sum())
        if "compound_pathway_source" in per_control.columns:
            # Counted raw: clean_text would fold the literal "none" source into "".
            sources = per_control["compound_pathway_source"].fillna("").astype(str).str.strip()
            compound_source_counts = {
                str(k): int(v) for k, v in sources.value_counts(dropna=False).sort_index().items()
            }
    return {
        "lib_id": int(lib_id),
        "input_csv": str(input_path),
        "output_csv": str(output_path),
        "n_rows": int(len(frame)),
        "n_rows_with_smiles": int(smiles.ne("").sum()),
        "n_unique_smiles": int(len({canonicalize_for_pathway(value) for value in smiles if clean_text(value)} - {""})),
        "n_rows_with_predicted_pathway": int(path.ne("").sum()),
        "n_rows_with_predicted_subpathway": int(sub.ne("").sum()),
        "n_controls": n_controls,
        "n_controls_with_compound_pathway": n_controls_with_compound_pathway,
        "n_controls_compound_agrees_rank1": n_rank1_agree,
        "compound_pathway_source_counts": compound_source_counts,
        "pathway_assignment_status_counts": {str(k): int(v) for k, v in status.value_counts(dropna=False).sort_index().items()},
        "pathway_confidence_tier_counts": {str(k): int(v) for k, v in tier.value_counts(dropna=False).sort_index().items()},
    }


# ── CLI entry point ──────────────────────────────────────────────────
def run(args: argparse.Namespace) -> dict[str, Any]:
    """Run pathway enrichment for all requested libraries."""
    repo_root = REPO_ROOT
    root = resolve(args.equalization_root, repo_root)
    reference_db = resolve(args.reference_db, repo_root)
    libs = parse_libs(args.libs)
    if not libs:
        raise ValueError("No libraries requested")
    frames = load_inputs(root, libs, args.input_name, bool(args.allow_missing))
    if not frames:
        raise FileNotFoundError(f"No rerank inputs found under {root}")
    assignment_cache = args.assignment_cache
    if assignment_cache is None:
        assignment_cache = root / "pathway_assignment_map.csv"
    else:
        assignment_cache = resolve(assignment_cache, repo_root)
    assignment_map, assignment_stats = build_assignment_map(
        frames=frames,
        reference_db=reference_db,
        assignment_cache=assignment_cache,
    )
    name_priors = load_subpathway_name_priors(reference_db) if args.subpathway_name_prior else None

    lib_summaries: list[dict[str, Any]] = []
    for lib_id in libs:
        if lib_id not in frames:
            continue
        input_path = rerank_path(root, lib_id, args.input_name)
        output_path = rerank_path(root, lib_id, args.output_name)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        enriched = enrich_one(
            frame=frames[lib_id],
            assignment_map=assignment_map,
            policy=args.compound_vote_policy,
            collapse_amino_acid_peptide=bool(args.collapse_amino_acid_peptide),
            name_priors=name_priors,
            apply_name_boost=bool(args.subpathway_name_prior),
        )
        enriched.to_csv(output_path, index=False)
        lib_summaries.append(summarize_lib(lib_id, input_path, output_path, enriched))

    summary_path = args.summary_json
    if summary_path is None:
        summary_path = root / "pathway_enrichment_summary.json"
    else:
        summary_path = resolve(summary_path, repo_root)
    summary = {
        "equalization_root": str(root),
        "reference_db": str(reference_db),
        "assignment_cache_csv": str(assignment_cache),
        "input_name": args.input_name,
        "output_name": args.output_name,
        "compound_vote_policy": args.compound_vote_policy,
        "collapse_amino_acid_peptide": bool(args.collapse_amino_acid_peptide),
        "subpathway_name_prior": bool(args.subpathway_name_prior),
        "n_subpathway_name_priors": int(len(name_priors or {})),
        **assignment_stats,
        "libs": lib_summaries,
    }
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return summary


def main() -> int:
    """CLI entry point."""
    run(parse_args())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
