#!/usr/bin/env python3
"""Generate retrieval-anchored analog candidate pools for X and hard holdouts."""

from __future__ import annotations

from structural_elucidation.chemistry.adduct_ladder import (
    hypotheses as adduct_ladder_hypotheses,
    include_polar_for_lib,
)
from structural_elucidation.paths import REPO_ROOT

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from structural_elucidation.utils.denovo import exact_mass_from_smiles, formula_from_smiles, read_table
from structural_elucidation.evaluation.controls import HAS_RDKIT
from structural_elucidation.utils.domain import (
    domain_distance,
    feature_summary,
    lib_mode,
    load_known_x_comp_ids,
    load_training_structures,
    parse_libs,
    safe_float,
    safe_int,
    write_csv,
    write_json,
)

if HAS_RDKIT:
    from rdkit import Chem, DataStructs
    from rdkit.Chem import AllChem


PROFILE_TOLERANCE = {"strict_1da": 1.0, "balanced": 5.0}


def fp_from_smiles(smiles: Any):
    """Compute a Morgan fingerprint for one SMILES string."""
    if not HAS_RDKIT:
        return None
    text = str(smiles or "").strip()
    if not text or text.lower() == "nan":
        return None
    mol = Chem.MolFromSmiles(text)
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)


def max_tanimoto(fp, anchor_fps: list[Any]) -> float:
    """Return the maximum Tanimoto vs a list of anchor fingerprints."""
    if fp is None or not anchor_fps:
        return 0.0
    values = [float(DataStructs.TanimotoSimilarity(fp, anchor_fp)) for anchor_fp in anchor_fps if anchor_fp is not None]
    return max(values) if values else 0.0


def mass_match(
    pred_mass: float,
    parent_mz: float,
    ion_mode: str,
    *,
    include_polar: bool = True,
) -> dict[str, Any]:
    """Evaluate precursor mass agreement for a predicted neutral mass."""
    rows = []
    for adduct, neutral in adduct_ladder_hypotheses(
        float(parent_mz), ion_mode, include_polar=include_polar
    ):
        rows.append(
            {
                "adduct": adduct,
                "neutral_mass": float(neutral),
                "abs_mass_error_da": abs(float(pred_mass) - float(neutral)),
                "hypothesis": adduct,
                "isotope_shift_count": 0,
            }
        )
    best = min(rows, key=lambda row: row["abs_mass_error_da"]) if rows else {}
    return best


def mass_score(abs_err: float, tolerance_da: float) -> float:
    """Convert absolute mass error to a clipped [0, 1] score."""
    value = safe_float(abs_err, float("inf"))
    if value <= 0.001:
        return 1.0
    if value >= float(tolerance_da):
        return 0.0
    return float(max(0.0, 1.0 - value / float(tolerance_da)))


def mass_error_bin(abs_err: float, strict_pass: bool) -> int:
    """Bucket absolute mass error for reporting bins."""
    value = safe_float(abs_err, float("inf"))
    if bool(strict_pass):
        return 0
    if value <= 0.5:
        return 1
    if value <= 1.0:
        return 2
    if value <= 2.0:
        return 3
    return 4


def load_retrieval_report(baseline_root: Path, lib_id: int) -> pd.DataFrame:
    """Load the baseline retrieval report for one library."""
    path = baseline_root / f"lib_{int(lib_id)}_{lib_mode(int(lib_id))}" / "retrieval" / "retrieval_report.csv"
    if not path.exists():
        return pd.DataFrame()
    df = read_table(path)
    if df.empty:
        return df
    df["rank"] = pd.to_numeric(df.get("rank"), errors="coerce").fillna(999999)
    return df.sort_values(["feature_key", "rank"], kind="mergesort")


def load_x_targets(unknown_feature_table: Path, libs: list[int], max_features: int) -> pd.DataFrame:
    """Load X-domain feature targets for analog generation."""
    targets = feature_summary(read_table(unknown_feature_table), default_scope="x_control")
    targets = targets[targets["lib_id"].astype("Int64").isin(libs) & targets["xid"].astype(str).ne("")].copy()
    if max_features > 0:
        targets = targets.groupby("lib_id", group_keys=False).head(max_features).copy()
    return targets.reset_index(drop=True)


def load_hard_targets(hard_holdout_dir: Path, libs: list[int], max_features: int) -> pd.DataFrame:
    """Load hard holdout targets for analog generation."""
    rows = []
    for lib_id in libs:
        path = hard_holdout_dir / f"x_like_hard_holdout_{int(lib_id)}.csv"
        if not path.exists():
            path = hard_holdout_dir / f"non_x_holdout_{int(lib_id)}.csv"
        if not path.exists():
            continue
        df = read_table(path)
        if df.empty:
            continue
        work = df.copy()
        work["holdout_scope"] = "non_x_known"
        work["xid"] = ""
        work["lib_id"] = int(lib_id)
        work["ion_mode"] = lib_mode(int(lib_id))
        work["peak_count"] = pd.to_numeric(work.get("peak_count"), errors="coerce").fillna(0)
        rows.append(work)
    out = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    if max_features > 0 and not out.empty:
        out = out.groupby("lib_id", group_keys=False).head(max_features).copy()
    return out.reset_index(drop=True)


def prepare_candidate_source(
    *,
    repo_root: Path,
    libs: list[int],
    known_x: Path,
    hard_holdout_dir: Path,
) -> pd.DataFrame:
    """Build the merged candidate source table for analog runs."""
    candidates = load_training_structures(repo_root, libs)
    if candidates.empty:
        return candidates
    exclude_comp_ids = set(load_known_x_comp_ids(known_x))
    hard_all = hard_holdout_dir / "x_like_hard_holdout_all.csv"
    if hard_all.exists():
        hard = read_table(hard_all)
        if "comp_id" in hard.columns:
            exclude_comp_ids.update(
                int(value) for value in pd.to_numeric(hard["comp_id"], errors="coerce").dropna().astype(int)
            )
    candidates = candidates[~candidates["comp_id"].astype(int).isin(exclude_comp_ids)].copy()
    candidates = candidates[
        candidates["canonical_smiles_key"].astype(str).str.strip().ne("")
        & ~candidates["canonical_smiles_key"].astype(str).str.lower().eq("nan")
    ].copy()
    candidates["exact_mass"] = candidates["canonical_smiles_key"].map(exact_mass_from_smiles)
    candidates["candidate_formula"] = candidates["canonical_smiles_key"].map(formula_from_smiles)
    candidates["fp"] = candidates["canonical_smiles_key"].map(fp_from_smiles)
    candidates["exact_mass"] = pd.to_numeric(candidates["exact_mass"], errors="coerce")
    return candidates.dropna(subset=["exact_mass", "fp"]).copy()


def anchors_for_x(target: pd.Series, retrieval_df: pd.DataFrame, anchor_top_n: int) -> list[dict[str, Any]]:
    """Select retrieval anchors for an X-domain target feature."""
    if retrieval_df.empty:
        return []
    group = retrieval_df[retrieval_df["feature_key"].astype(str).eq(str(target.get("feature_key", "")))].head(
        int(anchor_top_n)
    )
    rows = []
    for _, row in group.iterrows():
        smiles = str(row.get("smiles", "") or "")
        fp = fp_from_smiles(smiles)
        if fp is None:
            continue
        rows.append(
            {
                "anchor_source": "retrieval",
                "anchor_smiles": smiles,
                "anchor_rank": safe_int(row.get("rank"), 999999),
                "anchor_score": safe_float(row.get("retrieval_score"), 0.0),
                "anchor_formula": str(row.get("candidate_formula", "") or ""),
                "anchor_fragment_score": safe_float(row.get("fragment_match_count_score"), 0.0),
                "anchor_fp": fp,
            }
        )
    return rows


def anchors_for_hard(target: pd.Series, candidate_source: pd.DataFrame, anchor_top_n: int) -> list[dict[str, Any]]:
    """Select candidate-source anchors for a hard holdout target."""
    scoped = candidate_source[candidate_source["lib_id"].astype(int).eq(int(target.get("lib_id")))].copy()
    if scoped.empty:
        return []
    rows = []
    for _, row in scoped.iterrows():
        parts = domain_distance(
            target_parent_mz=safe_float(target.get("parent_mz"), 0.0),
            candidate_parent_mz=safe_float(row.get("parent_mz"), 0.0),
            target_ri=safe_float(target.get("ri"), 0.0),
            candidate_ri=safe_float(row.get("ri"), 0.0),
            target_peak_count=safe_float(target.get("peak_count"), 0.0),
            candidate_peak_count=0.0,
        )
        rows.append({**row.to_dict(), **parts})
    scored = pd.DataFrame(rows).sort_values(["domain_distance", "parent_mz"], kind="mergesort").head(int(anchor_top_n))
    anchors = []
    for rank, (_, row) in enumerate(scored.iterrows(), start=1):
        fp = row.get("fp")
        if fp is None:
            continue
        anchors.append(
            {
                "anchor_source": "nearest_known_feature_proxy",
                "anchor_smiles": str(row.get("canonical_smiles_key", "")),
                "anchor_rank": int(rank),
                "anchor_score": float(1.0 / max(rank, 1)),
                "anchor_formula": str(row.get("candidate_formula", "") or row.get("formula", "") or ""),
                "anchor_fragment_score": 0.0,
                "anchor_fp": fp,
            }
        )
    return anchors


def score_feature_candidates(
    *,
    target: pd.Series,
    candidate_source: pd.DataFrame,
    anchors: list[dict[str, Any]],
    tolerance_da: float,
    max_candidates_per_feature: int,
    include_polar: bool = True,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Score and rank analog candidates for one target feature."""
    lib_id = int(target.get("lib_id"))
    ion_mode = str(target.get("ion_mode") or lib_mode(lib_id))
    parent_mz = safe_float(target.get("parent_mz"), 0.0)
    anchor_fps = [anchor["anchor_fp"] for anchor in anchors if anchor.get("anchor_fp") is not None]
    anchor_formulas = {
        str(anchor.get("anchor_formula", "")) for anchor in anchors if str(anchor.get("anchor_formula", ""))
    }
    best_anchor = anchors[0] if anchors else {}
    rows: list[dict[str, Any]] = []
    scoped = candidate_source[candidate_source["lib_id"].astype(int).eq(lib_id)].copy()
    if scoped.empty:
        return [], [
            {
                "feature_key": str(target.get("feature_key", "")),
                "comp_id": safe_int(target.get("comp_id"), -1),
                "holdout_scope": str(target.get("holdout_scope", "")),
                "anchor_source": str(anchor.get("anchor_source", "")),
                "anchor_rank": safe_int(anchor.get("anchor_rank"), 0),
                "anchor_smiles": str(anchor.get("anchor_smiles", "")),
                "anchor_formula": str(anchor.get("anchor_formula", "")),
                "anchor_score": safe_float(anchor.get("anchor_score"), 0.0),
            }
            for anchor in anchors
        ]
    if anchor_fps:
        scoped_fps = list(scoped["fp"])
        tanimoto_matrix = [DataStructs.BulkTanimotoSimilarity(anchor_fp, scoped_fps) for anchor_fp in anchor_fps]
        anchor_tanimotos = np.max(np.asarray(tanimoto_matrix, dtype=float), axis=0)
    else:
        anchor_tanimotos = np.zeros(scoped.shape[0], dtype=float)

    prefilter_size = max(int(max_candidates_per_feature) * 8, 200)
    exact_masses = pd.to_numeric(scoped["exact_mass"], errors="coerce").to_numpy(dtype=float)
    neutral_masses = np.asarray(
        [
            float(neutral)
            for _, neutral in adduct_ladder_hypotheses(
                parent_mz, ion_mode, include_polar=include_polar
            )
        ],
        dtype=float,
    )
    rough_mass_errors = np.min(np.abs(exact_masses.reshape(-1, 1) - neutral_masses.reshape(1, -1)), axis=1)
    anchor_order = np.argsort(-anchor_tanimotos)[:prefilter_size]
    mass_order = np.argsort(rough_mass_errors)[:prefilter_size]
    candidate_positions = sorted(set(int(pos) for pos in np.concatenate([anchor_order, mass_order])))

    for pos in candidate_positions:
        cand = scoped.iloc[int(pos)]
        pred_mass = safe_float(cand.get("exact_mass"), math.nan)
        if not np.isfinite(pred_mass):
            continue
        match = mass_match(pred_mass, parent_mz, ion_mode, include_polar=include_polar)
        abs_err = safe_float(match.get("abs_mass_error_da"), math.inf)
        strict_pass = abs_err <= 1.0
        anchor_tan = float(anchor_tanimotos[pos])
        candidate_formula = str(cand.get("candidate_formula", "") or cand.get("formula", "") or "")
        formula_anchor_match = bool(candidate_formula and candidate_formula in anchor_formulas)
        score = (
            0.42 * mass_score(abs_err, max(float(tolerance_da), 1.0))
            + 0.36 * anchor_tan
            + 0.12 * safe_float(best_anchor.get("anchor_score"), 0.0)
            + 0.06 * safe_float(best_anchor.get("anchor_fragment_score"), 0.0)
            + (0.04 if formula_anchor_match else 0.0)
        )
        rows.append(
            {
                "feature_key": str(target.get("feature_key", "")),
                "comp_id": safe_int(target.get("comp_id"), -1),
                "parent_mz": parent_mz,
                "predicted_smiles": str(cand.get("canonical_smiles_key", "")),
                "predicted_exact_mass": pred_mass,
                "predicted_formula": candidate_formula,
                "formula_condition_source": str(best_anchor.get("anchor_source", "none")),
                "formula_condition_formula": str(best_anchor.get("anchor_formula", "")),
                "formula_condition_nonzero": bool(best_anchor.get("anchor_formula", "")),
                "abs_mass_error_da": abs_err,
                "primary_abs_mass_error_da": abs_err,
                "mass_error_bin": mass_error_bin(abs_err, strict_pass),
                "mass_error_hypothesis": str(match.get("hypothesis", "")),
                "mass_hypothesis_adduct": str(match.get("adduct", "")),
                "mass_hypothesis_neutral_mass": safe_float(match.get("neutral_mass")),
                "mass_hypothesis_isotope_shift_count": safe_int(match.get("isotope_shift_count"), 0),
                "strict_mass_pass": bool(strict_pass),
                "mass_match_score": mass_score(abs_err, max(float(tolerance_da), 1.0)),
                "mean_log_prob": float(score),
                "retrieval_tanimoto": float(anchor_tan),
                "substructure_filter_pass": True,
                "generation_support_count": int(max(1, round(anchor_tan * 10))),
                "generation_support_fraction": float(anchor_tan),
                "analog_score": float(score),
                "anchor_source": str(best_anchor.get("anchor_source", "")),
                "anchor_rank": safe_int(best_anchor.get("anchor_rank"), 0),
                "anchor_smiles": str(best_anchor.get("anchor_smiles", "")),
                "anchor_formula": str(best_anchor.get("anchor_formula", "")),
                "anchor_score": safe_float(best_anchor.get("anchor_score"), 0.0),
                "anchor_fragment_score": safe_float(best_anchor.get("anchor_fragment_score"), 0.0),
                "anchor_tanimoto": float(anchor_tan),
                "formula_anchor_match": bool(formula_anchor_match),
            }
        )
    rows = sorted(
        rows,
        key=lambda row: (
            -float(row["analog_score"]),
            float(row["abs_mass_error_da"]),
            -float(row["anchor_tanimoto"]),
            str(row["predicted_smiles"]),
        ),
    )[: int(max_candidates_per_feature)]
    for rank, row in enumerate(rows, start=1):
        row["candidate_file_rank"] = int(rank)
        row["rank"] = int(rank)
    anchor_rows = [
        {
            "feature_key": str(target.get("feature_key", "")),
            "comp_id": safe_int(target.get("comp_id"), -1),
            "holdout_scope": str(target.get("holdout_scope", "")),
            "anchor_source": str(anchor.get("anchor_source", "")),
            "anchor_rank": safe_int(anchor.get("anchor_rank"), 0),
            "anchor_smiles": str(anchor.get("anchor_smiles", "")),
            "anchor_formula": str(anchor.get("anchor_formula", "")),
            "anchor_score": safe_float(anchor.get("anchor_score"), 0.0),
        }
        for anchor in anchors
    ]
    return rows, anchor_rows


def write_profile(
    *,
    profile_dir: Path,
    lib_id: int,
    profile: str,
    candidate_rows: list[dict[str, Any]],
    anchor_rows: list[dict[str, Any]],
    top_k: int,
) -> dict[str, str]:
    """Write analog profile CSV and denovo report artifacts."""
    profile_dir.mkdir(parents=True, exist_ok=True)
    tolerance = PROFILE_TOLERANCE.get(profile.replace("non_x_holdout_", ""), 5.0)
    candidate_df = pd.DataFrame(candidate_rows)
    if not candidate_df.empty:
        candidate_df = candidate_df.sort_values(["feature_key", "candidate_file_rank"], kind="mergesort")
    denovo_rows = []
    if not candidate_df.empty:
        kept = candidate_df[
            pd.to_numeric(candidate_df["abs_mass_error_da"], errors="coerce").le(float(tolerance))
        ].copy()
        kept = kept.sort_values(["feature_key", "candidate_file_rank"], kind="mergesort")
        for feature_key, group in kept.groupby("feature_key", sort=False):
            for rank, (_, row) in enumerate(group.head(int(top_k)).iterrows(), start=1):
                payload = row.to_dict()
                payload["rank"] = int(rank)
                payload["feature_label"] = str(feature_key)
                payload["is_valid"] = True
                payload["candidate_latent_cosine"] = float("nan")
                payload["samples_drawn"] = int(len(group))
                denovo_rows.append(payload)
    denovo_df = pd.DataFrame(denovo_rows)
    anchor_df = pd.DataFrame(anchor_rows)
    outputs = {
        "candidate_metrics_csv": str(profile_dir / f"merged_{int(lib_id)}_candidate_metrics_all.csv"),
        "denovo_csv": str(profile_dir / f"merged_{int(lib_id)}_denovo_report.csv"),
        "anchor_audit_csv": str(profile_dir / "analog_anchor_audit.csv"),
        "summary_json": str(profile_dir / "analog_generation_summary.json"),
    }
    write_csv(Path(outputs["candidate_metrics_csv"]), candidate_df)
    write_csv(profile_dir / "candidate_metrics_all.csv", candidate_df)
    write_csv(Path(outputs["denovo_csv"]), denovo_df)
    write_csv(profile_dir / "denovo_report.csv", denovo_df)
    write_csv(Path(outputs["anchor_audit_csv"]), anchor_df)
    summary = {
        "profile": str(profile),
        "lib_id": int(lib_id),
        "mass_tolerance_da": float(tolerance),
        "n_candidate_rows": int(candidate_df.shape[0]),
        "n_denovo_rows": int(denovo_df.shape[0]),
        "n_features_with_candidates": int(candidate_df["feature_key"].nunique()) if not candidate_df.empty else 0,
        "n_features_with_denovo_rows": int(denovo_df["feature_key"].nunique()) if not denovo_df.empty else 0,
        "outputs": outputs,
    }
    write_json(Path(outputs["summary_json"]), summary)
    (profile_dir / f"merged_{int(lib_id)}_selfies_run_summary.json").write_text(
        json.dumps({"mass_tolerance_da": float(tolerance), "analog_generation": summary}, indent=2) + "\n"
    )
    (profile_dir / "selfies_run_summary.json").write_text(
        json.dumps({"mass_tolerance_da": float(tolerance), "analog_generation": summary}, indent=2) + "\n"
    )
    return outputs


def run_analog_generation(
    *,
    repo_root: Path,
    output_root: Path,
    baseline_root: Path,
    hard_holdout_dir: Path,
    known_x: Path,
    unknown_feature_table: Path,
    libs: list[int],
    profiles: list[str],
    scopes: list[str],
    max_features_per_scope: int,
    anchor_top_n: int,
    max_candidates_per_feature: int,
    top_k: int,
    drop_polar_adducts: bool = False,
) -> dict[str, Any]:
    """Run retrieval-anchored analog generation across targets."""
    if not HAS_RDKIT:
        raise RuntimeError("RDKit is required for retrieval-anchored analog generation.")
    run_root = output_root / "runs" / "retrieval_anchored_analog_v1"
    candidate_source = prepare_candidate_source(
        repo_root=repo_root,
        libs=libs,
        known_x=known_x,
        hard_holdout_dir=hard_holdout_dir,
    )
    outputs: list[dict[str, Any]] = []
    x_targets = load_x_targets(unknown_feature_table, libs, max_features_per_scope)
    hard_targets = load_hard_targets(hard_holdout_dir, libs, max_features_per_scope)
    for lib_id in libs:
        mode = lib_mode(int(lib_id))
        retrieval_df = load_retrieval_report(baseline_root, int(lib_id))
        scoped_targets = []
        if "x_control" in scopes and not x_targets.empty:
            scoped_targets.append(("x_control", x_targets[x_targets["lib_id"].astype("Int64").eq(int(lib_id))].copy()))
        if "hard_holdout" in scopes and not hard_targets.empty:
            scoped_targets.append(
                ("hard_holdout", hard_targets[hard_targets["lib_id"].astype("Int64").eq(int(lib_id))].copy())
            )
        for scope, targets in scoped_targets:
            if targets.empty:
                continue
            include_polar = include_polar_for_lib(int(lib_id)) if drop_polar_adducts else True
            for profile in profiles:
                output_profile = profile if scope == "x_control" else f"non_x_holdout_{profile}"
                profile_dir = run_root / f"lib_{int(lib_id)}_{mode}" / "selfies_profiles" / output_profile
                candidate_rows: list[dict[str, Any]] = []
                anchor_rows: list[dict[str, Any]] = []
                for _, target in targets.iterrows():
                    target = target.copy()
                    target["holdout_scope"] = "x_control" if scope == "x_control" else "non_x_known"
                    anchors = (
                        anchors_for_x(target, retrieval_df, anchor_top_n)
                        if scope == "x_control"
                        else anchors_for_hard(target, candidate_source, anchor_top_n)
                    )
                    rows, anchors_out = score_feature_candidates(
                        target=target,
                        candidate_source=candidate_source,
                        anchors=anchors,
                        tolerance_da=PROFILE_TOLERANCE.get(profile, 5.0),
                        max_candidates_per_feature=max_candidates_per_feature,
                        include_polar=include_polar,
                    )
                    for row in rows:
                        row["holdout_scope"] = str(target["holdout_scope"])
                        row["analog_profile_scope"] = scope
                    candidate_rows.extend(rows)
                    anchor_rows.extend(anchors_out)
                profile_outputs = write_profile(
                    profile_dir=profile_dir,
                    lib_id=int(lib_id),
                    profile=output_profile,
                    candidate_rows=candidate_rows,
                    anchor_rows=anchor_rows,
                    top_k=int(top_k),
                )
                outputs.append(
                    {
                        "lib_id": int(lib_id),
                        "ion_mode": mode,
                        "scope": scope,
                        "profile": output_profile,
                        "profile_dir": str(profile_dir),
                        **profile_outputs,
                    }
                )
    payload = {
        "output_root": str(output_root),
        "run_root": str(run_root),
        "libs": [int(lib_id) for lib_id in libs],
        "profiles": profiles,
        "scopes": scopes,
        "candidate_source_rows": int(candidate_source.shape[0]),
        "outputs": outputs,
    }
    write_json(output_root / "analog_generation_summary.json", payload)
    return payload


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments for analog generation."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", type=Path, default=Path("models/denovo_x_domain_analog_v1"))
    parser.add_argument(
        "--baseline-root",
        type=Path,
        default=Path("models/merged_structuremasst_bundle_excl_x_unresolved_2026-04-17_metabolon"),
    )
    parser.add_argument("--hard-holdout-dir", type=Path, default=Path("models/denovo_x_domain_analog_v1/hard_holdouts"))
    parser.add_argument("--known-x", type=Path, default=Path("known_X.csv"))
    parser.add_argument(
        "--unknown-feature-table", type=Path, default=Path("merged_encoder/data/ions_unknown_merged.csv")
    )
    parser.add_argument("--libs", type=str, default="209,305")
    parser.add_argument("--profiles", type=str, default="strict_1da,balanced")
    parser.add_argument("--scopes", type=str, default="x_control,hard_holdout")
    parser.add_argument("--max-features-per-scope", type=int, default=0)
    parser.add_argument("--anchor-top-n", type=int, default=5)
    parser.add_argument("--max-candidates-per-feature", type=int, default=50)
    parser.add_argument("--top-k", type=int, default=10)
    parser.add_argument(
        "--drop-polar-adducts",
        action="store_true",
        help="On NEG-only libs 209/305, omit platform_hint=polar ladder rows (Cl/formate).",
    )
    return parser.parse_args()


def resolve(path: Path, repo_root: Path) -> Path:
    """Resolve a repo-relative or absolute path against the repo root."""
    path = path.expanduser()
    return path if path.is_absolute() else repo_root / path


def main() -> None:
    """CLI entry point for X-domain analog generation."""
    args = parse_args()
    repo_root = REPO_ROOT
    payload = run_analog_generation(
        repo_root=repo_root,
        output_root=resolve(args.output_root, repo_root),
        baseline_root=resolve(args.baseline_root, repo_root),
        hard_holdout_dir=resolve(args.hard_holdout_dir, repo_root),
        known_x=resolve(args.known_x, repo_root),
        unknown_feature_table=resolve(args.unknown_feature_table, repo_root),
        libs=parse_libs(args.libs),
        profiles=[part.strip() for part in str(args.profiles).split(",") if part.strip()],
        scopes=[part.strip() for part in str(args.scopes).split(",") if part.strip()],
        max_features_per_scope=int(args.max_features_per_scope),
        anchor_top_n=int(args.anchor_top_n),
        max_candidates_per_feature=int(args.max_candidates_per_feature),
        top_k=int(args.top_k),
        drop_polar_adducts=bool(args.drop_polar_adducts),
    )
    print(json.dumps({"run_root": payload["run_root"], "outputs": len(payload["outputs"])}, indent=2))


if __name__ == "__main__":
    main()
