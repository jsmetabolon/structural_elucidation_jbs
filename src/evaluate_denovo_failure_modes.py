#!/usr/bin/env python3
"""Failure-mode audit for de novo SELFIES candidate pools."""

from __future__ import annotations

import argparse
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd

try:
    from rdkit import Chem
    from rdkit.Chem.MolStandardize import rdMolStandardize
except ImportError:
    Chem = None
    rdMolStandardize = None

from build_post_run_decision_package import (
    LibContext,
    dataframe_to_markdown,
    mass_coverage_row,
)
from denovo_improvement_utils import (
    canonicalize_smiles,
    formula_from_smiles,
    read_table,
)
from evaluate_denovo_controls import HAS_RDKIT, load_controls, tanimoto
from evaluate_selfies_candidate_pools import (
    discover_profile_dirs,
    find_artifact,
    load_mass_tolerance,
)


GOOD_TANIMOTO = 0.7
_UNCHARGER = rdMolStandardize.Uncharger() if rdMolStandardize is not None else None


@dataclass(frozen=True)
class AuditControl:
    control_id: str
    holdout_scope: str
    lib_id: int
    ion_mode: str
    comp_id: int | None
    match_key: str
    smiles: list[str]
    formulas: set[str]


def extract_xid(feature_key: Any) -> str | None:
    match = re.search(r"(X-\d+)", str(feature_key or ""))
    return match.group(1) if match else None


def load_x_controls(
    *,
    known_x: Path,
    known_x_all: Path,
    lib_id: int,
    ion_mode: str,
) -> list[AuditControl]:
    controls, _summary = load_controls(
        known_x,
        known_x_all,
        control_lib_id=int(lib_id),
    )
    out: list[AuditControl] = []
    for xid, smiles_values in sorted(controls.items()):
        formulas = {formula_from_smiles(smiles) for smiles in smiles_values}
        out.append(
            AuditControl(
                control_id=str(xid),
                holdout_scope="x_control",
                lib_id=int(lib_id),
                ion_mode=str(ion_mode),
                comp_id=None,
                match_key=str(xid),
                smiles=list(smiles_values),
                formulas={formula for formula in formulas if formula},
            )
        )
    return out


def load_non_x_holdouts(holdout_dir: Path, lib_id: int, ion_mode: str) -> list[AuditControl]:
    path = holdout_dir / f"non_x_holdout_{int(lib_id)}.csv"
    if not path.exists():
        return []
    df = read_table(path)
    required = {"comp_id", "canonical_smiles"}
    if not required.issubset(set(df.columns)):
        return []
    out: list[AuditControl] = []
    for _, row in df.iterrows():
        try:
            comp_id = int(row["comp_id"])
        except Exception:
            continue
        canonical = canonicalize_smiles(row.get("canonical_smiles"))
        if canonical is None:
            continue
        formula = str(row.get("formula", "") or "").strip() or formula_from_smiles(canonical)
        out.append(
            AuditControl(
                control_id=f"non_x_{int(lib_id)}_{comp_id}",
                holdout_scope="non_x_known",
                lib_id=int(lib_id),
                ion_mode=str(ion_mode),
                comp_id=comp_id,
                match_key=str(comp_id),
                smiles=[canonical],
                formulas={formula} if formula else set(),
            )
        )
    return out


def normalize_prediction_pool(df: pd.DataFrame, source_path: Path | None = None) -> pd.DataFrame:
    out = df.copy()
    if out.empty:
        return out
    if "predicted_smiles" not in out.columns and "smiles" in out.columns:
        out["predicted_smiles"] = out["smiles"]
    if "rank" not in out.columns:
        out["rank"] = np.arange(1, len(out) + 1)
    if "feature_key" not in out.columns:
        out["feature_key"] = ""
    if "comp_id" not in out.columns:
        out["comp_id"] = np.nan
    out["_feature_key_str"] = out["feature_key"].astype(str)
    out["_xid"] = out["_feature_key_str"].map(extract_xid)
    out["_comp_id"] = pd.to_numeric(out["comp_id"], errors="coerce")
    out["_rank"] = pd.to_numeric(out["rank"], errors="coerce").fillna(999999).astype(int)
    out["_predicted_smiles"] = out["predicted_smiles"].fillna("").astype(str)
    out["_predicted_canonical"] = out["_predicted_smiles"].map(lambda value: canonicalize_smiles(value) or "")
    if "predicted_formula" in out.columns:
        out["_predicted_formula"] = out["predicted_formula"].fillna("").astype(str)
    elif "candidate_formula" in out.columns:
        out["_predicted_formula"] = out["candidate_formula"].fillna("").astype(str)
    else:
        out["_predicted_formula"] = out["_predicted_canonical"].map(formula_from_smiles)
    out["_abs_mass_error_da"] = pd.to_numeric(out.get("abs_mass_error_da"), errors="coerce")
    out["_source_path"] = str(source_path) if source_path is not None else ""
    return out


def pool_for_control(pool: pd.DataFrame, control: AuditControl) -> pd.DataFrame:
    if pool.empty:
        return pool
    if control.holdout_scope == "x_control":
        return pool[pool["_xid"].astype(str).eq(str(control.match_key))].copy()
    comp_mask = pool["_comp_id"].eq(float(control.comp_id)) if control.comp_id is not None else False
    key_pattern = rf"(?:^|[^A-Za-z0-9])(?:comp(?:ound)?[_-]?id|comp)[_-]?{re.escape(str(control.match_key))}(?:$|[^A-Za-z0-9])"
    key_mask = pool["_feature_key_str"].str.contains(key_pattern, regex=True, na=False)
    return pool[comp_mask | key_mask].copy()


def neutralize_smiles_for_tanimoto(smiles: str) -> str:
    canonical = canonicalize_smiles(smiles) or ""
    if not HAS_RDKIT or Chem is None or _UNCHARGER is None or not canonical:
        return canonical
    try:
        mol = Chem.MolFromSmiles(canonical)
        if mol is None:
            return canonical
        neutral = _UNCHARGER.uncharge(mol)
        Chem.SanitizeMol(neutral)
        return Chem.MolToSmiles(neutral, isomericSmiles=True)
    except Exception:
        return canonical


def best_tanimoto_for_smiles(predicted: str, controls: Iterable[str]) -> tuple[bool, float]:
    if not predicted:
        return False, 0.0
    exact = False
    best = 0.0
    control_values = list(controls)
    for control_smiles in control_values:
        if predicted == control_smiles:
            exact = True
        best = max(best, tanimoto(predicted, control_smiles))
    if best < GOOD_TANIMOTO:
        predicted_neutral = neutralize_smiles_for_tanimoto(predicted)
        for control_smiles in control_values:
            control_neutral = neutralize_smiles_for_tanimoto(control_smiles)
            if predicted_neutral and predicted_neutral == control_neutral:
                exact = True
            if predicted_neutral and control_neutral:
                best = max(best, tanimoto(predicted_neutral, control_neutral))
    return exact, float(best)


def annotate_pool(pool: pd.DataFrame, control: AuditControl) -> pd.DataFrame:
    out = pool.copy()
    if out.empty:
        return out
    exact_values: list[bool] = []
    tan_values: list[float] = []
    formula_matches: list[bool] = []
    for _, row in out.iterrows():
        exact, best = best_tanimoto_for_smiles(
            str(row.get("_predicted_canonical", "")),
            control.smiles,
        )
        exact_values.append(bool(exact))
        tan_values.append(float(best))
        formula_matches.append(str(row.get("_predicted_formula", "")) in control.formulas)
    out["_exact_match"] = exact_values
    out["_best_tanimoto"] = tan_values
    out["_formula_match"] = formula_matches
    return out


def classify_failure(
    *,
    annotated_pool: pd.DataFrame,
    top_pool: pd.DataFrame,
    mass_tolerance_da: float,
) -> str:
    if annotated_pool.empty:
        return "no_candidate"
    mass_pool = annotated_pool[
        pd.to_numeric(annotated_pool["_abs_mass_error_da"], errors="coerce").le(float(mass_tolerance_da))
    ].copy()
    if mass_pool.empty:
        return "no_mass_match"
    if bool(top_pool.get("_best_tanimoto", pd.Series(dtype=float)).ge(GOOD_TANIMOTO).any()):
        return "good_candidate_top10"
    if bool(annotated_pool["_best_tanimoto"].ge(GOOD_TANIMOTO).any()):
        return "good_candidate_misranked"
    if not bool(mass_pool["_formula_match"].any()):
        return "mass_match_wrong_formula"
    return "formula_match_wrong_structure"


def control_audit_row(
    *,
    control: AuditControl,
    full_pool: pd.DataFrame,
    top_pool: pd.DataFrame,
    mass_tolerance_da: float,
    experiment: str,
    profile: str,
    denovo_path: Path | None,
    candidate_path: Path | None,
) -> dict[str, Any]:
    full_annotated = annotate_pool(pool_for_control(full_pool, control), control)
    top_annotated = annotate_pool(pool_for_control(top_pool, control), control)
    failure_mode = classify_failure(
        annotated_pool=full_annotated,
        top_pool=top_annotated,
        mass_tolerance_da=float(mass_tolerance_da),
    )
    best_full = (
        full_annotated.sort_values(
            ["_best_tanimoto", "_exact_match", "_abs_mass_error_da"],
            ascending=[False, False, True],
            na_position="last",
        ).head(1)
        if not full_annotated.empty
        else pd.DataFrame()
    )
    best_top = (
        top_annotated.sort_values(
            ["_best_tanimoto", "_exact_match", "_rank"],
            ascending=[False, False, True],
            na_position="last",
        ).head(1)
        if not top_annotated.empty
        else pd.DataFrame()
    )

    def first_value(df: pd.DataFrame, col: str, default: Any = "") -> Any:
        if df.empty or col not in df.columns:
            return default
        value = df.iloc[0].get(col, default)
        if pd.isna(value):
            return default
        return value

    mass_series = pd.to_numeric(full_annotated.get("_abs_mass_error_da"), errors="coerce")
    row = {
        "experiment": str(experiment),
        "lib_id": int(control.lib_id),
        "ion_mode": str(control.ion_mode),
        "profile": str(profile),
        "holdout_scope": str(control.holdout_scope),
        "control_id": str(control.control_id),
        "comp_id": "" if control.comp_id is None else int(control.comp_id),
        "n_candidates": int(full_annotated.shape[0]),
        "n_topk": int(top_annotated.shape[0]),
        "n_mass_match": int(mass_series.le(float(mass_tolerance_da)).sum()) if not full_annotated.empty else 0,
        "n_formula_match": int(full_annotated.get("_formula_match", pd.Series(dtype=bool)).sum()) if not full_annotated.empty else 0,
        "failure_mode": failure_mode,
        "top10_best_tanimoto": round(float(first_value(best_top, "_best_tanimoto", 0.0)), 4),
        "top10_best_rank": first_value(best_top, "_rank", ""),
        "full_pool_best_tanimoto": round(float(first_value(best_full, "_best_tanimoto", 0.0)), 4),
        "full_pool_best_rank": first_value(best_full, "_rank", ""),
        "full_pool_best_abs_mass_error_da": first_value(best_full, "_abs_mass_error_da", ""),
        "full_pool_best_formula_match": bool(first_value(best_full, "_formula_match", False)),
        "full_pool_best_smiles": first_value(best_full, "_predicted_smiles", ""),
        "denovo_csv": str(denovo_path) if denovo_path is not None else "",
        "candidate_metrics_csv": str(candidate_path) if candidate_path is not None else "",
    }
    for window in [0.5, 1.0, 2.0, 5.0]:
        key = str(window).replace(".", "p")
        row[f"has_candidate_within_{key}_da"] = bool(mass_series.le(float(window)).any()) if not full_annotated.empty else False
    return row


def metric_summary(rows: list[dict[str, Any]]) -> dict[str, Any]:
    n_controls = len(rows)
    covered = [row for row in rows if int(row.get("n_candidates", 0) or 0) > 0]
    top_values = [float(row.get("top10_best_tanimoto", 0.0) or 0.0) for row in rows]
    oracle_values = [float(row.get("full_pool_best_tanimoto", 0.0) or 0.0) for row in rows]
    return {
        "n_controls": int(n_controls),
        "n_controls_covered": int(len(covered)),
        "rank1_or_top10_mean_best_tanimoto": round(float(np.mean(top_values)) if top_values else 0.0, 4),
        "full_pool_oracle_mean_best_tanimoto": round(float(np.mean(oracle_values)) if oracle_values else 0.0, 4),
        "top10_tanimoto_ge_0p7": int(sum(value >= GOOD_TANIMOTO for value in top_values)),
        "full_pool_oracle_tanimoto_ge_0p7": int(sum(value >= GOOD_TANIMOTO for value in oracle_values)),
    }


def evaluate_failure_modes(
    *,
    run_root: Path,
    known_x: Path,
    known_x_all: Path,
    holdout_dir: Path | None,
    output_dir: Path,
    top_k: int,
) -> dict[str, Any]:
    if not HAS_RDKIT:
        raise RuntimeError("RDKit is required for de novo failure-mode evaluation.")
    output_dir.mkdir(parents=True, exist_ok=True)

    per_control_rows: list[dict[str, Any]] = []
    summary_rows: list[dict[str, Any]] = []
    mass_rows: list[dict[str, Any]] = []
    for experiment, lib_dir, lib_id, ion_mode, profile, profile_dir in discover_profile_dirs(run_root):
        controls = load_x_controls(
            known_x=known_x,
            known_x_all=known_x_all,
            lib_id=int(lib_id),
            ion_mode=str(ion_mode),
        )
        if holdout_dir is not None:
            controls.extend(load_non_x_holdouts(holdout_dir, int(lib_id), str(ion_mode)))
        denovo_path = find_artifact(profile_dir, int(lib_id), "denovo_report.csv")
        candidate_path = find_artifact(profile_dir, int(lib_id), "candidate_metrics_all.csv")
        if candidate_path is None or not candidate_path.exists():
            continue
        candidate_df = pd.read_csv(candidate_path)
        candidate_pool = normalize_prediction_pool(candidate_df, candidate_path)
        denovo_pool = (
            normalize_prediction_pool(pd.read_csv(denovo_path), denovo_path)
            if denovo_path is not None and denovo_path.exists()
            else pd.DataFrame()
        )
        if not denovo_pool.empty:
            denovo_pool = denovo_pool[denovo_pool["_rank"].le(int(top_k))].copy()
        tolerance = load_mass_tolerance(profile, profile_dir, int(lib_id))
        ctx = LibContext(
            lib_dir=lib_dir,
            lib_id=int(lib_id),
            ion_mode=str(ion_mode),
            manifest={"lib_id": int(lib_id), "ion_mode": str(ion_mode)},
        )
        mass_row = mass_coverage_row(ctx, str(profile), profile_dir, candidate_df)
        mass_row["experiment"] = str(experiment)
        mass_rows.append(mass_row)
        lane_rows = [
            control_audit_row(
                control=control,
                full_pool=candidate_pool,
                top_pool=denovo_pool,
                mass_tolerance_da=float(tolerance),
                experiment=str(experiment),
                profile=str(profile),
                denovo_path=denovo_path,
                candidate_path=candidate_path,
            )
            for control in controls
        ]
        per_control_rows.extend(lane_rows)
        for scope in ["x_control", "non_x_known"]:
            scoped_rows = [row for row in lane_rows if row["holdout_scope"] == scope]
            if not scoped_rows:
                continue
            mode_counts = pd.Series([row["failure_mode"] for row in scoped_rows]).value_counts().to_dict()
            summary_rows.append(
                {
                    "experiment": str(experiment),
                    "lib_id": int(lib_id),
                    "ion_mode": str(ion_mode),
                    "profile": str(profile),
                    "holdout_scope": scope,
                    **metric_summary(scoped_rows),
                    **{f"failure_{key}": int(value) for key, value in mode_counts.items()},
                }
            )

    per_control_df = pd.DataFrame(per_control_rows)
    summary_df = pd.DataFrame(summary_rows)
    mass_df = pd.DataFrame(mass_rows)
    per_control_csv = output_dir / "failure_mode_per_control.csv"
    summary_csv = output_dir / "failure_mode_summary.csv"
    mass_csv = output_dir / "failure_mode_mass_coverage.csv"
    per_control_df.to_csv(per_control_csv, index=False)
    summary_df.to_csv(summary_csv, index=False)
    mass_df.to_csv(mass_csv, index=False)
    report_md = output_dir / "candidate_pool_improvement_targets.md"
    report_md.write_text(render_report(summary_df, output_dir))
    payload = {
        "run_root": str(run_root),
        "known_x": str(known_x),
        "known_x_all": str(known_x_all),
        "holdout_dir": str(holdout_dir) if holdout_dir is not None else None,
        "top_k": int(top_k),
        "outputs": {
            "per_control_csv": str(per_control_csv),
            "summary_csv": str(summary_csv),
            "mass_csv": str(mass_csv),
            "report_md": str(report_md),
            "summary_json": str(output_dir / "failure_mode_audit_summary.json"),
        },
    }
    (output_dir / "failure_mode_audit_summary.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n"
    )
    return payload


def render_report(summary_df: pd.DataFrame, output_dir: Path) -> str:
    lines = [
        "# De Novo Candidate-Pool Failure Audit",
        "",
        "This audit separates candidate generation failures from ranking failures for X controls and non-X holdouts when matching feature rows exist.",
        "",
    ]
    if summary_df.empty:
        lines.append("No candidate-pool rows were available.")
    else:
        display_cols = [
            "lib_id",
            "profile",
            "experiment",
            "holdout_scope",
            "n_controls",
            "n_controls_covered",
            "top10_tanimoto_ge_0p7",
            "full_pool_oracle_tanimoto_ge_0p7",
            "full_pool_oracle_mean_best_tanimoto",
        ]
        available = [col for col in display_cols if col in summary_df.columns]
        lines.append(dataframe_to_markdown(summary_df[available]))
    lines.extend(
        [
            "",
            "## Outputs",
            f"- Per-control rows: `{output_dir / 'failure_mode_per_control.csv'}`",
            f"- Summary rows: `{output_dir / 'failure_mode_summary.csv'}`",
            f"- Mass coverage rows: `{output_dir / 'failure_mode_mass_coverage.csv'}`",
            "",
        ]
    )
    return "\n".join(lines)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--known-x", type=Path, default=Path("known_X.csv"))
    parser.add_argument(
        "--known-x-all",
        type=Path,
        default=Path("known_X_all_available_structures_simple.csv"),
    )
    parser.add_argument("--holdout-dir", type=Path, default=Path(""))
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("models/denovo_generation_improvement_v1/failure_audit"),
    )
    parser.add_argument("--top-k", type=int, default=10)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run_root = Path(args.run_root).expanduser().resolve()
    known_x = Path(args.known_x).expanduser()
    known_x_all = Path(args.known_x_all).expanduser()
    holdout_dir = Path(args.holdout_dir).expanduser() if str(args.holdout_dir).strip() else None
    output_dir = Path(args.output_dir).expanduser()
    repo_root = Path(__file__).resolve().parent
    if not known_x.is_absolute():
        known_x = repo_root / known_x
    if not known_x_all.is_absolute():
        known_x_all = repo_root / known_x_all
    if holdout_dir is not None and not holdout_dir.is_absolute():
        holdout_dir = repo_root / holdout_dir
    if not output_dir.is_absolute():
        output_dir = repo_root / output_dir
    payload = evaluate_failure_modes(
        run_root=run_root,
        known_x=known_x,
        known_x_all=known_x_all,
        holdout_dir=holdout_dir,
        output_dir=output_dir,
        top_k=int(args.top_k),
    )
    print(json.dumps(payload["outputs"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
