#!/usr/bin/env python3
# flake8: noqa
"""Generate multi-conjugate candidates in the biotransformation candidate shape."""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from denovo_improvement_utils import canonicalize_smiles, exact_mass_from_smiles, formula_from_smiles
from evaluate_denovo_failure_modes import load_x_controls
from nce_evidence_utils import HAS_RDKIT, TARGETED_REACTION_SMARTS, AllChem, Chem
from rerank_denovo_assembly_v3 import eval_reranked, write_denovo_report
from run_retrieval_anchored_denovo_assembly_v4 import resolve
from x_domain_analog_utils import dataframe_to_markdown, safe_float, safe_int, write_csv, write_json


DEFAULT_OUTPUT_ROOT = Path("models/multi_conjugate_candidates_v1")
RUN_NAME = "multi_conjugate_v1"


def lib_scope(args: argparse.Namespace) -> str:
    return f"lib_{int(args.lib_id)}_{str(args.ion_mode)}"


@dataclass(frozen=True)
class BaseMolecule:
    name: str
    smiles: str
    family: str


BASE_LIBRARY: tuple[BaseMolecule, ...] = (
    BaseMolecule("glycine", "NCC(=O)O", "free_amino_acid"),
    BaseMolecule("alanine", "CC(N)C(=O)O", "free_amino_acid"),
    BaseMolecule("serine", "NC(CO)C(=O)O", "free_amino_acid"),
    BaseMolecule("threonine", "CC(O)C(N)C(=O)O", "free_amino_acid"),
    BaseMolecule("glutamic_acid", "NC(CCC(=O)O)C(=O)O", "free_amino_acid"),
    BaseMolecule("phenylalanine", "NC(Cc1ccccc1)C(=O)O", "free_amino_acid"),
    BaseMolecule("tyrosine", "NC(Cc1ccc(O)cc1)C(=O)O", "free_amino_acid"),
    BaseMolecule("tryptophan", "NC(Cc1c[nH]c2ccccc12)C(=O)O", "free_amino_acid"),
    BaseMolecule("histidine", "NC(Cc1c[nH]cn1)C(=O)O", "free_amino_acid"),
    BaseMolecule("indole_acrylic_acid", "O=C(O)C=Cc1c[nH]c2ccccc12", "aromatic_acid_core"),
    BaseMolecule("benzoic_acid", "O=C(O)c1ccccc1", "aromatic_acid_core"),
    BaseMolecule("salicylic_acid", "O=C(O)c1ccccc1O", "phenolic_acid_core"),
    BaseMolecule("cholic_acid_core", "CC(C)C1CCC2C3CCC4CC(O)CCC4(C)C3CCC12C", "bile_acid_core"),
    BaseMolecule("chenodeoxycholic_acid_core", "CC(C)C1CCC2C3CCC4CC(O)CCC4(C)C3CCC12", "bile_acid_core"),
)


CONJUGATE_REACTIONS: dict[str, tuple[str, ...]] = {
    "lactoyl": ("lactoyl_amine_conjugation", "lactoyl_hydroxyl_conjugation"),
    "acetyl": ("amine_acetylation", "hydroxyl_acetylation"),
    "glycine": ("amine_glycine_conjugation",),
    "glucuronide": ("glucuronide_hydroxyl_conjugation", "glucuronide_carboxyl_conjugation", "glucuronide_amine_conjugation"),
    "sulfate": ("hydroxyl_sulfation",),
    "phosphate": ("hydroxyl_phosphorylation",),
    "taurine": ("carboxyl_taurine_amide",),
    "carnitine": ("carboxyl_carnitine_ester",),
    "methyl": ("carboxyl_methylation",),
}


REQUIRED_COLUMNS = [
    "feature_key",
    "control_id",
    "lib_id",
    "ion_mode",
    "profile",
    "track",
    "predicted_smiles",
    "predicted_exact_mass",
    "predicted_formula",
    "target_neutral_mass",
    "abs_mass_error_da",
    "mass_match_score",
    "evidence_confidence",
    "forward_spectrum_score",
    "operator_family",
    "generation_score",
    "novelty_score",
    "rank",
]


def reaction_for_name(reaction_name: str) -> Any:
    smarts = TARGETED_REACTION_SMARTS.get(reaction_name)
    if not smarts or AllChem is None:
        return None
    try:
        return AllChem.ReactionFromSmarts(smarts)
    except Exception:
        return None


def canonical_product_smiles(product: Any) -> str:
    if Chem is None:
        return ""
    try:
        Chem.SanitizeMol(product)
        smiles = Chem.MolToSmiles(product, canonical=True)
    except Exception:
        return ""
    return canonicalize_smiles(smiles) or ""


def apply_reaction(base_smiles: str, reaction_name: str, max_products: int) -> list[str]:
    if not HAS_RDKIT or Chem is None:
        return []
    reaction = reaction_for_name(reaction_name)
    if reaction is None:
        return []
    mol = Chem.MolFromSmiles(str(base_smiles))
    if mol is None:
        return []
    products: list[str] = []
    seen: set[str] = set()
    try:
        product_sets = reaction.RunReactants((mol,))
    except Exception:
        return []
    for product_set in product_sets:
        for product in product_set:
            smiles = canonical_product_smiles(product)
            if smiles and smiles not in seen:
                seen.add(smiles)
                products.append(smiles)
            if len(products) >= int(max_products):
                return products
    return products


def mass_match_score(abs_error: float, tolerance_da: float) -> float:
    if not np.isfinite(abs_error):
        return 0.0
    return float(np.clip(1.0 - abs(float(abs_error)) / max(float(tolerance_da), 1e-9), 0.0, 1.0))


def seed_reactions(conjugate_name: Any) -> tuple[str, ...]:
    text = str(conjugate_name or "").strip().lower()
    if text in CONJUGATE_REACTIONS:
        return CONJUGATE_REACTIONS[text]
    return tuple(name for key, names in CONJUGATE_REACTIONS.items() if key in text for name in names)


def build_rows(seed_queue: pd.DataFrame, args: argparse.Namespace) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows: list[dict[str, Any]] = []
    skipped: list[dict[str, Any]] = []
    if not HAS_RDKIT:
        return pd.DataFrame(), pd.DataFrame([{"reason": "rdkit_unavailable", "n_seed_rows": int(seed_queue.shape[0])}])
    for _, seed in seed_queue.iterrows():
        control_id = str(seed.get("control_id", "") or "")
        feature_key = str(seed.get("feature_key", "") or "")
        conjugate_name = str(seed.get("conjugate_name", "") or "")
        target_neutral_mass = safe_float(seed.get("neutral_mass"), math.nan)
        if not control_id or not feature_key or not np.isfinite(target_neutral_mass):
            skipped.append({"control_id": control_id, "conjugate_name": conjugate_name, "reason": "missing_control_feature_or_mass"})
            continue
        reactions = seed_reactions(conjugate_name)
        if not reactions:
            skipped.append({"control_id": control_id, "conjugate_name": conjugate_name, "reason": "no_targeted_reaction_smarts"})
            continue
        for reaction_name in reactions:
            for base in BASE_LIBRARY:
                products = apply_reaction(base.smiles, reaction_name, int(args.max_products_per_base_reaction))
                if not products:
                    continue
                for product_smiles in products:
                    exact_mass = exact_mass_from_smiles(product_smiles)
                    if not np.isfinite(exact_mass):
                        continue
                    abs_error = abs(float(exact_mass) - float(target_neutral_mass))
                    score = mass_match_score(abs_error, float(args.mass_score_tolerance_da))
                    formula_prior = safe_float(seed.get("formula_prior_score"), 0.0)
                    ms1_supported = str(seed.get("ms1_buddy_supported", False)).lower() in {"true", "1", "yes", "y"}
                    evidence = float(np.clip(0.35 + 0.35 * formula_prior + (0.15 if ms1_supported else 0.0), 0.0, 1.0))
                    forward = float(np.clip(0.35 + 0.25 * score + (0.15 if ms1_supported else 0.0), 0.0, 1.0))
                    generation = float(np.clip(0.45 * score + 0.25 * evidence + 0.20 * forward + 0.10 * formula_prior, 0.0, 1.0))
                    rows.append(
                        {
                            "feature_key": feature_key,
                            "control_id": control_id,
                            "comp_id": safe_int(seed.get("comp_id"), -1),
                            "lib_id": safe_int(seed.get("lib_id"), int(args.lib_id)),
                            "ion_mode": str(seed.get("ion_mode", str(args.ion_mode)) or str(args.ion_mode)),
                            "profile": str(args.profile),
                            "track": str(args.track),
                            "predicted_smiles": product_smiles,
                            "predicted_exact_mass": round(float(exact_mass), 7),
                            "predicted_formula": formula_from_smiles(product_smiles),
                            "target_neutral_mass": round(float(target_neutral_mass), 7),
                            "abs_mass_error_da": round(float(abs_error), 7),
                            "mass_match_score": round(score, 6),
                            "evidence_confidence": round(evidence, 6),
                            "forward_spectrum_score": round(forward, 6),
                            "operator_family": f"multi_conjugate:{reaction_name}",
                            "generation_score": round(generation, 6),
                            "novelty_score": 1.0,
                            "base_name": base.name,
                            "base_family": base.family,
                            "conjugate_name": conjugate_name,
                            "source_functional_group": seed.get("source_functional_group", ""),
                            "expected_evidence": seed.get("expected_evidence", ""),
                            "adduct_label": seed.get("adduct_label", ""),
                            "formula_prior_score": formula_prior,
                            "ms1_buddy_supported": ms1_supported,
                            "candidate_source": "multi_conjugate_candidates_v1",
                        }
                    )
    candidates = pd.DataFrame(rows)
    if not candidates.empty:
        candidates = candidates.sort_values(
            ["feature_key", "mass_match_score", "generation_score", "evidence_confidence", "abs_mass_error_da", "predicted_smiles"],
            ascending=[True, False, False, False, True, True],
            kind="mergesort",
        ).drop_duplicates(["feature_key", "predicted_smiles"], keep="first")
        candidates["rank"] = candidates.groupby("feature_key", dropna=False).cumcount() + 1
        candidates = candidates[pd.to_numeric(candidates["rank"], errors="coerce").le(int(args.max_candidates_per_feature))].reset_index(drop=True)
        for column in REQUIRED_COLUMNS:
            if column not in candidates.columns:
                candidates[column] = ""
        candidates = candidates[REQUIRED_COLUMNS + [col for col in candidates.columns if col not in REQUIRED_COLUMNS]]
    return candidates, pd.DataFrame(skipped)


def run(args: argparse.Namespace) -> dict[str, Any]:
    repo_root = Path(__file__).resolve().parent
    output_root = resolve(args.output_root, repo_root)
    seed_queue_path = resolve(args.seed_queue, repo_root)
    seed_queue = pd.read_csv(seed_queue_path, low_memory=False) if seed_queue_path.exists() else pd.DataFrame()
    candidates, skipped = build_rows(seed_queue, args)
    out_dir = (
        output_root
        / "runs"
        / RUN_NAME
        / f"track_{args.track}"
        / lib_scope(args)
        / "selfies_profiles"
        / str(args.profile)
    )
    candidate_path = out_dir / "candidate_metrics_all.csv"
    denovo_path = out_dir / "denovo_report.csv"
    write_csv(candidate_path, candidates)
    ranked = candidates[candidates.get("rank", pd.Series(dtype=int)).le(10)].copy() if not candidates.empty else candidates
    write_denovo_report(ranked, denovo_path)
    controls = load_x_controls(
        known_x=repo_root / "known_X.csv",
        known_x_all=repo_root / "known_X_all_available_structures_simple.csv",
        lib_id=int(args.lib_id),
        ion_mode=str(args.ion_mode),
    )
    metrics, per_control = eval_reranked(
        candidate_df=candidates,
        ranked_df=ranked,
        controls=controls,
        lib_id=int(args.lib_id),
        ion_mode=str(args.ion_mode),
        profile=str(args.profile),
        policy=RUN_NAME,
        top_k=10,
        candidate_path=candidate_path,
        denovo_path=denovo_path,
    )
    write_csv(output_root / "multi_conjugate_metrics.csv", pd.DataFrame(metrics))
    write_csv(output_root / "multi_conjugate_per_control.csv", pd.DataFrame(per_control))
    write_csv(output_root / "multi_conjugate_candidate_skips.csv", skipped)
    summary = (
        candidates.groupby(["conjugate_name", "operator_family"], dropna=False)
        .agg(
            n_controls=("control_id", "nunique"),
            n_candidates=("predicted_smiles", "count"),
            best_mass_match_score=("mass_match_score", "max"),
        )
        .reset_index()
        if not candidates.empty
        else pd.DataFrame()
    )
    write_csv(output_root / "multi_conjugate_candidate_summary.csv", summary)
    report = [
        "# Multi-Conjugate Candidates V1",
        "",
        f"Seed queue: `{seed_queue_path}`",
        "",
        "## Summary",
        dataframe_to_markdown(summary),
        "",
        "## Skipped Seeds",
        dataframe_to_markdown(skipped.head(80) if not skipped.empty else skipped),
        "",
    ]
    (output_root / "multi_conjugate_candidate_report.md").write_text("\n".join(report))
    payload = {
        "output_root": str(output_root),
        "candidate_metrics_all": str(candidate_path),
        "n_seed_rows": int(seed_queue.shape[0]),
        "n_candidate_rows": int(candidates.shape[0]),
        "n_controls": int(candidates["control_id"].nunique()) if not candidates.empty else 0,
        "rdkit_available": bool(HAS_RDKIT),
    }
    write_json(output_root / "multi_conjugate_candidate_summary.json", payload)
    print(json.dumps(payload, indent=2, sort_keys=True))
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--seed-queue", type=Path, default=Path("models/multi_conjugate_generation_v1/multi_conjugate_seed_queue.csv"))
    parser.add_argument("--lib-id", type=int, default=209)
    parser.add_argument("--ion-mode", default="negative")
    parser.add_argument("--profile", default="strict_1da")
    parser.add_argument("--track", default="control_blind")
    parser.add_argument("--mass-score-tolerance-da", type=float, default=0.15)
    parser.add_argument("--max-products-per-base-reaction", type=int, default=4)
    parser.add_argument("--max-candidates-per-feature", type=int, default=200)
    run(parser.parse_args())


if __name__ == "__main__":
    main()
