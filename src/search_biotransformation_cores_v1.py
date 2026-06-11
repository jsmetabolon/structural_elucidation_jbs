#!/usr/bin/env python3
"""Search candidate core structures for biotransformation-guided generation."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from biotransformation_guided_v1_utils import (
    TRACK_ANNOTATION_ASSISTED,
    TRACK_CONTROL_BLIND,
    TRANSFORMS,
    formula_distance_score,
    peptide_candidates_from_composition,
)
from denovo_improvement_utils import exact_mass_from_smiles, formula_from_smiles
from nce_evidence_utils import load_excluded_structure_keys
from run_retrieval_anchored_denovo_assembly_v4 import resolve
from x_domain_analog_utils import (
    dataframe_to_markdown,
    load_known_x_comp_ids,
    load_training_structures,
    safe_float,
    write_csv,
    write_json,
)


DEFAULT_OUTPUT_ROOT = Path("models/biotransformation_guided_denovo_v1")


def parse_csv(text: Any) -> list[str]:
    return [part.strip() for part in str(text or "").split(",") if part.strip()]


def inchikey_connectivity(value: Any) -> str:
    return str(value or "").split("-", 1)[0]


def load_exclusion_keys(repo_root: Path) -> tuple[set[str], set[str], set[str], set[int]]:
    paths = [
        repo_root / "known_X_all_available_structures_simple.csv",
        repo_root / "models/x_known_curation_v1/curated_x_structures_broad.csv",
        repo_root / "models/x_known_curation_v1/known_x_evidence_catalog.csv",
        repo_root / "models/denovo_generation_improvement_v1/holdouts/non_x_holdout_209.csv",
        repo_root / "models/denovo_generation_improvement_v1/holdouts/non_x_holdout_305.csv",
        repo_root / "models/denovo_x_domain_analog_v1/hard_holdouts/x_like_hard_holdout_all.csv",
    ]
    excluded = load_excluded_structure_keys(paths)
    canonical = {str(value) for value in excluded.get("canonical_smiles", pd.Series(dtype=str)).astype(str) if str(value)}
    inchikey = {str(value) for value in excluded.get("inchikey", pd.Series(dtype=str)).astype(str) if str(value)}
    connectivity = {inchikey_connectivity(value) for value in inchikey if inchikey_connectivity(value)}
    comp_ids = load_known_x_comp_ids(repo_root / "known_X.csv")
    return canonical, inchikey, connectivity, comp_ids


def build_source_pool(repo_root: Path, libs: list[int]) -> pd.DataFrame:
    source = load_training_structures(repo_root, libs)
    if source.empty:
        return pd.DataFrame()
    canonical_excl, inchikey_excl, connectivity_excl, comp_id_excl = load_exclusion_keys(repo_root)
    out = source.copy()
    out["comp_id"] = pd.to_numeric(out.get("comp_id"), errors="coerce")
    out = out[~out["comp_id"].isin(comp_id_excl)].copy()
    out = out[~out["canonical_smiles_key"].astype(str).isin(canonical_excl)].copy()
    out = out[~out["inchikey_key"].astype(str).isin(inchikey_excl)].copy()
    out = out[~out["inchikey_key"].map(inchikey_connectivity).isin(connectivity_excl)].copy()
    out["source_exact_mass"] = out["canonical_smiles_key"].map(exact_mass_from_smiles)
    out["source_formula"] = out["formula"].where(out["formula"].astype(str).ne(""), out["canonical_smiles_key"].map(formula_from_smiles))
    out = out[pd.to_numeric(out["source_exact_mass"], errors="coerce").notna()].copy()
    out["source_exact_mass"] = pd.to_numeric(out["source_exact_mass"], errors="coerce")
    return out.drop_duplicates(["lib_id", "canonical_smiles_key"]).reset_index(drop=True)


def source_score(row: pd.Series, source: pd.Series, mass_tolerance_da: float) -> float:
    mass_error = abs(safe_float(source.get("source_exact_mass"), math.inf) - safe_float(row.get("core_mass"), math.nan))
    mass_score = max(0.0, 1.0 - mass_error / max(float(mass_tolerance_da), 1e-9))
    formula_score = formula_distance_score(source.get("source_formula", ""), row.get("core_formula", ""))
    return float(np.clip(0.58 * mass_score + 0.24 * formula_score + 0.18 * safe_float(row.get("evidence_confidence"), 0.0), 0.0, 1.0))


def search_cores(args: argparse.Namespace) -> dict[str, Any]:
    repo_root = Path(__file__).resolve().parent
    output_root = resolve(args.output_root, repo_root)
    output_root.mkdir(parents=True, exist_ok=True)
    evidence_path = resolve(args.evidence_csv, repo_root)
    if not evidence_path.exists():
        raise FileNotFoundError(f"Missing evidence CSV: {evidence_path}")
    evidence = pd.read_csv(evidence_path)
    tracks = set(parse_csv(args.tracks))
    if tracks:
        evidence = evidence[evidence["track"].astype(str).isin(tracks)].copy()
    controls = set(parse_csv(args.controls))
    if controls:
        evidence = evidence[evidence["control_id"].astype(str).isin(controls)].copy()
    source_pool = build_source_pool(repo_root, [209])
    rows: list[dict[str, Any]] = []
    for _, hyp in evidence.iterrows():
        transform = str(hyp.get("transform_name", ""))
        track = str(hyp.get("track", TRACK_CONTROL_BLIND))
        core_mass = safe_float(hyp.get("core_mass"), math.nan)
        if transform == "peptide":
            composition = str(hyp.get("peptide_composition", "") or "")
            if not composition:
                continue
            smiles_values = peptide_candidates_from_composition(composition, max_permutations=int(args.max_peptide_permutations))
            for rank, smiles in enumerate(smiles_values, start=1):
                rows.append(
                    {
                        **hyp.to_dict(),
                        "core_rank": int(rank),
                        "core_source_type": "peptide_composition",
                        "source_comp_id": "",
                        "core_smiles": smiles,
                        "source_exact_mass": exact_mass_from_smiles(smiles),
                        "source_formula": formula_from_smiles(smiles),
                        "core_mass_error_da": abs(exact_mass_from_smiles(smiles) - core_mass) if np.isfinite(core_mass) else math.nan,
                        "core_formula_score": 1.0,
                        "core_source_score": 0.95 if track == TRACK_ANNOTATION_ASSISTED else 0.75,
                    }
                )
            continue
        if transform not in TRANSFORMS:
            continue
        if not np.isfinite(core_mass) or source_pool.empty:
            continue
        scoped = source_pool[source_pool["lib_id"].astype(int).eq(209)].copy()
        scoped["core_mass_error_da"] = (pd.to_numeric(scoped["source_exact_mass"], errors="coerce") - core_mass).abs()
        scoped = scoped[scoped["core_mass_error_da"].le(float(args.core_tolerance_da))].copy()
        if scoped.empty:
            continue
        scoped["core_formula_score"] = scoped["source_formula"].map(lambda formula: formula_distance_score(formula, hyp.get("core_formula", "")))
        scoped["core_source_score"] = scoped.apply(lambda source: source_score(hyp, source, float(args.core_tolerance_da)), axis=1)
        scoped = scoped.sort_values(
            ["core_source_score", "core_mass_error_da", "core_formula_score", "canonical_smiles_key"],
            ascending=[False, True, False, True],
            kind="mergesort",
        ).head(int(args.max_cores_per_hypothesis))
        for rank, source in enumerate(scoped.itertuples(index=False), start=1):
            row = hyp.to_dict()
            source_dict = source._asdict()
            rows.append(
                {
                    **row,
                    "core_rank": int(rank),
                    "core_source_type": "known_non_x_core_mass_search",
                    "source_comp_id": source_dict.get("comp_id", ""),
                    "core_smiles": source_dict.get("canonical_smiles_key", ""),
                    "source_exact_mass": source_dict.get("source_exact_mass", math.nan),
                    "source_formula": source_dict.get("source_formula", ""),
                    "core_mass_error_da": source_dict.get("core_mass_error_da", math.nan),
                    "core_formula_score": source_dict.get("core_formula_score", 0.0),
                    "core_source_score": source_dict.get("core_source_score", 0.0),
                    "source_asset_path": source_dict.get("asset_path", ""),
                }
            )
    cores = pd.DataFrame(rows)
    if not cores.empty:
        cores = cores.sort_values(
            ["track", "control_id", "evidence_confidence", "core_source_score", "core_mass_error_da"],
            ascending=[True, True, False, False, True],
            kind="mergesort",
        ).drop_duplicates(["track", "control_id", "mass_hypothesis_adduct", "transform_name", "core_smiles"], keep="first")
    summary = (
        cores.groupby(["track", "transform_name"], dropna=False)
        .agg(n_core_rows=("core_smiles", "size"), n_controls=("control_id", "nunique"))
        .reset_index()
        if not cores.empty
        else pd.DataFrame(columns=["track", "transform_name", "n_core_rows", "n_controls"])
    )
    write_csv(output_root / "core_search_candidates.csv", cores)
    write_csv(output_root / "core_search_summary.csv", summary)
    report = "\n".join(
        [
            "# Biotransformation Core Search V1",
            "",
            f"Source pool rows after leakage exclusions: {int(source_pool.shape[0])}",
            "",
            "## Core Search Summary",
            dataframe_to_markdown(summary),
            "",
        ]
    )
    (output_root / "core_search_report.md").write_text(report)
    payload = {
        "output_root": str(output_root),
        "evidence_csv": str(evidence_path),
        "n_source_pool_rows": int(source_pool.shape[0]),
        "n_core_rows": int(cores.shape[0]),
        "n_core_controls": int(cores["control_id"].nunique()) if not cores.empty else 0,
        "core_candidates_csv": str(output_root / "core_search_candidates.csv"),
    }
    write_json(output_root / "core_search_summary.json", payload)
    print(json.dumps(payload, indent=2, sort_keys=True))
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--evidence-csv", type=Path, default=DEFAULT_OUTPUT_ROOT / "biotransformation_evidence.csv")
    parser.add_argument("--tracks", default=f"{TRACK_CONTROL_BLIND},{TRACK_ANNOTATION_ASSISTED}")
    parser.add_argument("--controls", default="")
    parser.add_argument("--core-tolerance-da", type=float, default=1.5)
    parser.add_argument("--max-cores-per-hypothesis", type=int, default=40)
    parser.add_argument("--max-peptide-permutations", type=int, default=24)
    search_cores(parser.parse_args())


if __name__ == "__main__":
    main()
