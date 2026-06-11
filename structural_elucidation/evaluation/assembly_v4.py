#!/usr/bin/env python3
"""Targeted retrieval-anchored de novo assembly V4 hybrid.

V4 keeps V2's failure-aware, mass-filtered generation path as the base and
adds V1-style SELFIES/BRICS diversity only after profile mass filtering. It
writes multiple predeclared ranking variants over one physical candidate pool.
"""

from __future__ import annotations

import argparse
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from structural_elucidation.paths import REPO_ROOT, resolve

import structural_elucidation.evaluation.assembly_v1 as v1
from structural_elucidation.utils.denovo import (
    exact_mass_from_smiles,
    formula_from_smiles,
    inchikey_from_smiles,
    read_table,
)
from structural_elucidation.evaluation.controls import HAS_RDKIT
from structural_elucidation.evaluation.failure_modes import load_x_controls
from structural_elucidation.evaluation.analog import (
    PROFILE_TOLERANCE,
    anchors_for_hard,
    anchors_for_x,
    fp_from_smiles,
    load_hard_targets,
    load_retrieval_report,
    load_x_targets,
    mass_error_bin,
    mass_match,
    mass_score,
)
from structural_elucidation.utils.domain import (
    dataframe_to_markdown,
    lib_mode,
    parse_libs,
    safe_float,
    safe_int,
    write_csv,
    write_json,
)

if HAS_RDKIT:
    from rdkit import Chem, DataStructs
    from rdkit.Chem import rdFMCS


DEFAULT_OUTPUT_ROOT = Path("models/denovo_pipeline_v4/assembly_v4")
DEFAULT_BASELINE_ROOT = Path("models/merged_structuremasst_bundle_excl_x_unresolved_2026-04-17_metabolon")
DEFAULT_ANALOG_ROOT = Path("models/denovo_x_domain_analog_v1")
DEFAULT_ABLATION_ROOT = Path("models/decoder_ablation_experiments")
DEFAULT_HARD_HOLDOUT_DIR = Path("models/denovo_x_domain_analog_v1/hard_holdouts")
DEFAULT_FAILURE_AUDIT = Path("models/denovo_pipeline_v2/audit/per_control_failure_audit.csv")
GOOD_TANIMOTO = 0.7
RUN_NAME = "retrieval_anchored_denovo_assembly_v4"
RANKING_VARIANTS = ("v4_balanced", "v4_anchor_heavy", "v4_mass_first", "v4_no_qed_novelty")
DEFAULT_VARIANT_BY_PROFILE = {"balanced": "v4_balanced", "strict_1da": "v4_mass_first"}


def parse_csv(text: str) -> list[str]:
    return [part.strip() for part in str(text).split(",") if part.strip()]


def parse_targets(text: str) -> tuple[list[int], list[str]]:
    libs: list[int] = []
    profiles: list[str] = []
    for item in parse_csv(text):
        if ":" not in item:
            continue
        lib_text, profile = item.split(":", 1)
        lib_id = safe_int(lib_text, -1)
        if lib_id > 0 and lib_id not in libs:
            libs.append(lib_id)
        if profile and profile not in profiles:
            profiles.append(profile)
    return libs, profiles


def load_failure_map(path: Path) -> dict[tuple[str, int, str, str], str]:
    if not path.exists():
        return {}
    df = read_table(path)
    if df.empty:
        return {}
    out: dict[tuple[str, int, str, str], str] = {}
    for _, row in df.iterrows():
        scope = str(row.get("holdout_scope", ""))
        lib_id = safe_int(row.get("lib_id"), -1)
        profile = str(row.get("base_profile", "") or row.get("profile", "")).replace("non_x_holdout_", "")
        control_id = str(row.get("control_id", "") or "")
        comp_id = str(row.get("comp_id", "") or "")
        bucket = str(row.get("failure_bucket", "") or "")
        if control_id:
            out[(scope, lib_id, profile, control_id)] = bucket
        if comp_id and comp_id.lower() != "nan":
            out[(scope, lib_id, profile, comp_id)] = bucket
    return out


def target_control_key(target: pd.Series, scope: str) -> str:
    if scope == "x_control":
        return str(target.get("xid", "") or "")
    comp_id = target.get("comp_id", "")
    return str(safe_int(comp_id, -1)) if not pd.isna(comp_id) else ""


def failure_bucket_for_target(
    failure_map: dict[tuple[str, int, str, str], str],
    target: pd.Series,
    scope: str,
    profile: str,
) -> str:
    key = (
        "x_control" if scope == "x_control" else "non_x_known",
        safe_int(target.get("lib_id"), -1),
        str(profile).replace("non_x_holdout_", ""),
        target_control_key(target, scope),
    )
    return failure_map.get(key, "")


def heavy_atom_count(smiles: str) -> int:
    mol = v1.mol_from_smiles(smiles)
    return int(mol.GetNumHeavyAtoms()) if mol is not None else 0


def mass_delta_directed_candidates(
    *,
    anchor_smiles: str,
    parent_mz: float,
    ion_mode: str,
    tolerance_da: float,
    limit: int,
) -> dict[str, set[str]]:
    anchor = v1.canonical_or_empty(anchor_smiles)
    if not anchor:
        return {}
    anchor_mass = safe_float(exact_mass_from_smiles(anchor), math.nan)
    anchor_error = safe_float(mass_match(anchor_mass, parent_mz, ion_mode).get("abs_mass_error_da"), math.inf)
    raw: dict[str, set[str]] = defaultdict(set)
    for smiles in v1.atom_substitution_candidates(anchor, limit=max(20, int(limit))):
        raw[smiles].add("mass_delta_directed_atom_isostere")
    for smiles in v1.bond_order_candidates(anchor, limit=max(12, int(limit) // 2)):
        raw[smiles].add("formula_preserving_bond_order_edit")
    for smiles in v1.tautomer_candidates(anchor, limit=max(12, int(limit) // 2)):
        raw[smiles].add("formula_preserving_tautomer_protomer")
    out: dict[str, set[str]] = defaultdict(set)
    for smiles, modes in raw.items():
        mass = safe_float(exact_mass_from_smiles(smiles), math.nan)
        if not np.isfinite(mass):
            continue
        abs_err = safe_float(mass_match(mass, parent_mz, ion_mode).get("abs_mass_error_da"), math.inf)
        if abs_err <= float(tolerance_da) or abs_err < anchor_error:
            out[smiles].update(modes)
        if len(out) >= int(limit):
            break
    return out


def mcs_rgroup_edit_candidates(anchor_smiles: str, seed_smiles: list[str], limit: int = 120) -> dict[str, set[str]]:
    if not HAS_RDKIT:
        return {}
    anchor = v1.canonical_or_empty(anchor_smiles)
    anchor_mol = v1.mol_from_smiles(anchor)
    if anchor_mol is None:
        return {}
    out: dict[str, set[str]] = defaultdict(set)
    min_atoms = max(3, int(0.35 * max(heavy_atom_count(anchor), 1)))
    effective_limit = min(int(limit), 48)
    for seed in seed_smiles[:4]:
        seed = v1.canonical_or_empty(seed)
        seed_mol = v1.mol_from_smiles(seed)
        if seed_mol is None:
            continue
        try:
            mcs = rdFMCS.FindMCS([anchor_mol, seed_mol], timeout=2)
        except Exception:
            continue
        if int(mcs.numAtoms) < min_atoms or not mcs.smartsString:
            continue
        patt = Chem.MolFromSmarts(mcs.smartsString)
        if patt is None:
            continue
        raw = set()
        raw.update(v1.tautomer_candidates(seed, limit=8))
        raw.update(v1.bond_order_candidates(seed, limit=8))
        raw.update(v1.atom_substitution_candidates(seed, limit=16))
        for candidate in raw:
            mol = v1.mol_from_smiles(candidate)
            if mol is None:
                continue
            try:
                if mol.HasSubstructMatch(patt):
                    out[candidate].add("mcs_rgroup_edit")
            except Exception:
                continue
            if len(out) >= effective_limit:
                return out
    return out


def mass_constrained_fragment_replacements(
    *,
    anchor_smiles: str,
    source_smiles: list[str],
    parent_mz: float,
    ion_mode: str,
    tolerance_da: float,
    limit: int,
) -> dict[str, set[str]]:
    out: dict[str, set[str]] = defaultdict(set)
    effective_limit = min(int(limit), 40)
    for smiles in v1.brics_recombination_candidates(anchor_smiles, source_smiles[:4], limit=effective_limit):
        mass = safe_float(exact_mass_from_smiles(smiles), math.nan)
        if not np.isfinite(mass):
            continue
        abs_err = safe_float(mass_match(mass, parent_mz, ion_mode).get("abs_mass_error_da"), math.inf)
        if abs_err <= float(tolerance_da):
            out[smiles].add("mass_constrained_fragment_replacement")
        if len(out) >= effective_limit:
            break
    return out


def should_expand_fragments(failure_bucket: str) -> bool:
    return failure_bucket in {"low_structural_similarity", "x_domain_specific_gap", "no_candidate", "exact_known_only"}


def should_focus_mass(failure_bucket: str) -> bool:
    return failure_bucket in {"no_mass_compatible_novel", "no_candidate", "exact_known_only"}


def merge_generated(target: dict[str, set[str]], source: dict[str, set[str]]) -> None:
    for smiles, modes in source.items():
        target[smiles].update(modes)


def broad_mass_filtered_candidates(
    *,
    anchor_smiles: str,
    fragment_source_smiles: list[str],
    token_pool: list[str],
    parent_mz: float,
    ion_mode: str,
    tolerance_da: float,
    limit: int,
) -> dict[str, set[str]]:
    """Run broad V1-style operators but only retain profile mass-compatible rows."""
    anchor = v1.canonical_or_empty(anchor_smiles)
    if not anchor:
        return {}
    raw: dict[str, set[str]] = defaultdict(set)
    for smiles in v1.selfies_mutation_candidates(anchor, token_pool, limit=max(30, int(limit))):
        raw[smiles].add("broad_mass_filtered_selfies_mutation")
    for smiles in v1.brics_recombination_candidates(anchor, fragment_source_smiles, limit=max(40, int(limit))):
        raw[smiles].add("broad_mass_filtered_brics_recombination")
    out: dict[str, set[str]] = defaultdict(set)
    for smiles, modes in raw.items():
        mass = safe_float(exact_mass_from_smiles(smiles), math.nan)
        if not np.isfinite(mass):
            continue
        abs_err = safe_float(mass_match(mass, parent_mz, ion_mode).get("abs_mass_error_da"), math.inf)
        if abs_err <= float(tolerance_da):
            out[smiles].update(modes)
        if len(out) >= int(limit):
            break
    return out


def generate_for_anchor_v4(
    *,
    anchor: dict[str, Any],
    fragment_source_smiles: list[str],
    token_pool: list[str],
    parent_mz: float,
    ion_mode: str,
    tolerance_da: float,
    failure_bucket: str,
    per_anchor_limit: int,
) -> dict[str, set[str]]:
    anchor_smiles = v1.canonical_or_empty(anchor.get("anchor_smiles", ""))
    if not anchor_smiles:
        return {}
    generated: dict[str, set[str]] = defaultdict(set)
    generated[anchor_smiles].add("anchor_identity_benchmark_rejected_if_known")
    taut_limit = max(20, int(per_anchor_limit) // 3)
    for smiles in v1.tautomer_candidates(anchor_smiles, limit=taut_limit):
        generated[smiles].add("formula_preserving_tautomer_protomer")
    for smiles in v1.bond_order_candidates(anchor_smiles, limit=max(12, int(per_anchor_limit) // 4)):
        generated[smiles].add("formula_preserving_bond_order_edit")
    merge_generated(
        generated,
        mass_delta_directed_candidates(
            anchor_smiles=anchor_smiles,
            parent_mz=parent_mz,
            ion_mode=ion_mode,
            tolerance_da=tolerance_da,
            limit=max(40, int(per_anchor_limit) // 2),
        ),
    )
    if should_expand_fragments(failure_bucket):
        merge_generated(
            generated,
            mcs_rgroup_edit_candidates(
                anchor_smiles,
                fragment_source_smiles,
                limit=max(40, int(per_anchor_limit) // 2),
            ),
        )
        merge_generated(
            generated,
            mass_constrained_fragment_replacements(
                anchor_smiles=anchor_smiles,
                source_smiles=fragment_source_smiles,
                parent_mz=parent_mz,
                ion_mode=ion_mode,
                tolerance_da=tolerance_da,
                limit=max(40, int(per_anchor_limit) // 2),
            ),
        )
    anchor_mass = safe_float(exact_mass_from_smiles(anchor_smiles), math.nan)
    anchor_abs_err = safe_float(mass_match(anchor_mass, parent_mz, ion_mode).get("abs_mass_error_da"), math.inf)
    if anchor_abs_err <= max(float(tolerance_da) * 2.0, 2.0) or should_focus_mass(failure_bucket):
        for smiles in v1.selfies_mutation_candidates(
            anchor_smiles, token_pool, limit=max(30, int(per_anchor_limit) // 3)
        ):
            generated[smiles].add("selfies_seed_mutation_after_mass_seed")
    merge_generated(
        generated,
        broad_mass_filtered_candidates(
            anchor_smiles=anchor_smiles,
            fragment_source_smiles=fragment_source_smiles,
            token_pool=token_pool,
            parent_mz=parent_mz,
            ion_mode=ion_mode,
            tolerance_da=tolerance_da,
            limit=max(50, int(per_anchor_limit) // 2),
        ),
    )
    return generated


def primary_operator(modes: set[str]) -> str:
    priority = [
        "mass_delta_directed_atom_isostere",
        "mcs_rgroup_edit",
        "mass_constrained_fragment_replacement",
        "broad_mass_filtered_brics_recombination",
        "broad_mass_filtered_selfies_mutation",
        "formula_preserving_tautomer_protomer",
        "formula_preserving_bond_order_edit",
        "selfies_seed_mutation_after_mass_seed",
        "anchor_identity_benchmark_rejected_if_known",
    ]
    for item in priority:
        if item in modes:
            return item
    return sorted(modes)[0] if modes else ""


def is_broad_operator(operator_set: str) -> bool:
    return "broad_mass_filtered_" in str(operator_set)


def broad_cap_for_lib(lib_id: int) -> int:
    if int(lib_id) == 305:
        return 25
    if int(lib_id) == 209:
        return 75
    return 50


def score_components(
    *,
    abs_err: float,
    tolerance: float,
    anchor_tan: float,
    anchor_score: float,
    support_count: int,
    plausibility: float,
    novelty_score: float,
    formula_anchor_match: bool,
    operator_family: str,
) -> dict[str, float]:
    mass_component = mass_score(abs_err, max(float(tolerance), 1.0))
    support_component = min(1.0, safe_float(support_count, 0.0) / 3.0)
    formula_component = 1.0 if formula_anchor_match else 0.0
    operator_bonus = 0.05 if operator_family in {"mass_delta_directed_atom_isostere", "mcs_rgroup_edit"} else 0.0
    if str(operator_family).startswith("broad_mass_filtered_"):
        operator_bonus = 0.03
    scores = {
        "score_component_mass": float(mass_component),
        "score_component_anchor_tanimoto": float(anchor_tan),
        "score_component_anchor_score": float(anchor_score),
        "score_component_operator_support": float(support_component),
        "score_component_qed": float(plausibility),
        "score_component_novelty": float(novelty_score),
        "score_component_formula_match": float(formula_component),
        "score_component_operator_bonus": float(operator_bonus),
    }
    scores["v4_balanced_score"] = float(
        0.42 * mass_component
        + 0.22 * anchor_tan
        + 0.08 * anchor_score
        + 0.08 * support_component
        + 0.04 * plausibility
        + 0.06 * formula_component
        + 0.02 * novelty_score
        + operator_bonus
    )
    scores["v4_anchor_heavy_score"] = float(
        0.32 * mass_component
        + 0.35 * anchor_tan
        + 0.10 * anchor_score
        + 0.06 * support_component
        + 0.03 * plausibility
        + 0.08 * formula_component
        + 0.01 * novelty_score
        + operator_bonus
    )
    scores["v4_mass_first_score"] = float(
        0.58 * mass_component
        + 0.15 * anchor_tan
        + 0.07 * anchor_score
        + 0.06 * support_component
        + 0.02 * plausibility
        + 0.05 * formula_component
        + 0.01 * novelty_score
        + operator_bonus
    )
    scores["v4_no_qed_novelty_score"] = float(
        0.46 * mass_component
        + 0.28 * anchor_tan
        + 0.10 * anchor_score
        + 0.08 * support_component
        + 0.06 * formula_component
        + operator_bonus
    )
    return scores


def score_generated_candidates_v4(
    *,
    target: pd.Series,
    anchors: list[dict[str, Any]],
    candidate_source: pd.DataFrame,
    exclusion_canonicals: set[str],
    exclusion_inchikeys: set[str],
    profile: str,
    max_candidates_per_feature: int,
    failure_bucket: str,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    lib_id = int(target.get("lib_id"))
    ion_mode = str(target.get("ion_mode") or lib_mode(lib_id))
    parent_mz = safe_float(target.get("parent_mz"), 0.0)
    tolerance = PROFILE_TOLERANCE.get(str(profile).replace("non_x_holdout_", ""), 5.0)
    seed_limit = 10 if should_focus_mass(failure_bucket) or should_expand_fragments(failure_bucket) else 5
    seed_anchors = v1.mass_near_known_seed_anchors(
        candidate_source=candidate_source,
        lib_id=lib_id,
        parent_mz=parent_mz,
        ion_mode=ion_mode,
        tolerance_da=max(float(tolerance), 5.0),
        limit=seed_limit,
    )
    anchor_smiles_seen = {str(anchor.get("anchor_smiles", "")) for anchor in anchors}
    anchors = list(anchors) + [
        seed
        for seed in seed_anchors
        if str(seed.get("anchor_smiles", "")) and str(seed.get("anchor_smiles", "")) not in anchor_smiles_seen
    ]
    db_fps = list(candidate_source[candidate_source["lib_id"].astype(int).eq(lib_id)]["fp"])
    anchor_fps = [anchor.get("anchor_fp") for anchor in anchors if anchor.get("anchor_fp") is not None]
    anchor_formulas = {
        str(anchor.get("anchor_formula", "")) for anchor in anchors if str(anchor.get("anchor_formula", ""))
    }
    fragment_limit = 40 if should_expand_fragments(failure_bucket) else 24
    fragment_source_smiles = v1.select_fragment_source_smiles(
        candidate_source=candidate_source,
        lib_id=lib_id,
        parent_mz=parent_mz,
        ion_mode=ion_mode,
        anchor_fps=anchor_fps,
        limit=fragment_limit,
    )
    token_pool = v1.collect_selfies_tokens(
        [str(anchor.get("anchor_smiles", "")) for anchor in anchors] + fragment_source_smiles
    )
    generated_modes: dict[str, set[str]] = defaultdict(set)
    anchor_by_candidate: dict[str, dict[str, Any]] = {}
    anchor_cap = 6 if should_expand_fragments(failure_bucket) or should_focus_mass(failure_bucket) else 4
    anchors = anchors[:anchor_cap]
    per_anchor_limit = min(max(int(max_candidates_per_feature), 120), 160)
    for anchor in anchors:
        for smiles, modes in generate_for_anchor_v4(
            anchor=anchor,
            fragment_source_smiles=fragment_source_smiles,
            token_pool=token_pool,
            parent_mz=parent_mz,
            ion_mode=ion_mode,
            tolerance_da=tolerance,
            failure_bucket=failure_bucket,
            per_anchor_limit=per_anchor_limit,
        ).items():
            generated_modes[smiles].update(modes)
            anchor_by_candidate.setdefault(smiles, anchor)
    rows: list[dict[str, Any]] = []
    structure_rows: list[dict[str, Any]] = []
    novelty_rows: list[dict[str, Any]] = []
    for smiles, modes in generated_modes.items():
        canonical = v1.canonical_or_empty(smiles)
        if not canonical:
            continue
        key = inchikey_from_smiles(canonical)
        exact_known = canonical in exclusion_canonicals or (key and key in exclusion_inchikeys)
        fp = fp_from_smiles(canonical)
        if fp is None:
            continue
        nearest_db_tan = v1.nearest_db_similarity(fp, db_fps)
        novelty_rows.append(
            {
                "feature_key": str(target.get("feature_key", "")),
                "comp_id": safe_int(target.get("comp_id"), -1),
                "candidate_smiles": canonical,
                "inchikey": key,
                "exact_known_structure_rejected": bool(exact_known),
                "nearest_db_tanimoto": nearest_db_tan,
                "generation_modes": "|".join(sorted(modes)),
                "operator_family": primary_operator(modes),
                "assembly_version": "v4",
                "failure_bucket": failure_bucket,
            }
        )
        if exact_known:
            continue
        exact_mass = exact_mass_from_smiles(canonical)
        if not np.isfinite(safe_float(exact_mass)):
            continue
        match = mass_match(float(exact_mass), parent_mz, ion_mode)
        abs_err = safe_float(match.get("abs_mass_error_da"), math.inf)
        anchor = anchor_by_candidate.get(smiles, anchors[0] if anchors else {})
        anchor_tan = v1.max_anchor_similarity(fp, anchor_fps)
        formula = formula_from_smiles(canonical)
        formula_anchor_match = bool(formula and formula in anchor_formulas)
        plausibility = v1.qed_score(canonical)
        novelty_score = max(0.0, 1.0 - min(nearest_db_tan, 1.0))
        support_count = len(modes)
        operator_family = primary_operator(modes)
        components = score_components(
            abs_err=abs_err,
            tolerance=tolerance,
            anchor_tan=anchor_tan,
            anchor_score=safe_float(anchor.get("anchor_score"), 0.0),
            support_count=support_count,
            plausibility=plausibility,
            novelty_score=novelty_score,
            formula_anchor_match=formula_anchor_match,
            operator_family=operator_family,
        )
        score = components["v4_balanced_score"]
        payload = {
            "feature_key": str(target.get("feature_key", "")),
            "comp_id": safe_int(target.get("comp_id"), -1),
            "parent_mz": parent_mz,
            "predicted_smiles": canonical,
            "predicted_exact_mass": float(exact_mass),
            "predicted_formula": formula,
            "formula_condition_source": str(anchor.get("anchor_source", "none")),
            "formula_condition_formula": str(anchor.get("anchor_formula", "")),
            "formula_condition_nonzero": bool(anchor.get("anchor_formula", "")),
            "abs_mass_error_da": abs_err,
            "primary_abs_mass_error_da": abs_err,
            "mass_error_bin": mass_error_bin(abs_err, abs_err <= 1.0),
            "mass_error_hypothesis": str(match.get("hypothesis", "")),
            "mass_hypothesis_adduct": str(match.get("adduct", "")),
            "mass_hypothesis_neutral_mass": safe_float(match.get("neutral_mass")),
            "mass_hypothesis_isotope_shift_count": safe_int(match.get("isotope_shift_count"), 0),
            "strict_mass_pass": bool(abs_err <= 1.0),
            "strict_1da_pass": bool(abs_err <= 1.0),
            "profile_mass_pass": bool(abs_err <= float(tolerance)),
            "mass_match_score": mass_score(abs_err, max(float(tolerance), 1.0)),
            "mean_log_prob": float(score),
            "retrieval_tanimoto": float(anchor_tan),
            "anchor_tanimoto": float(anchor_tan),
            "substructure_filter_pass": True,
            "generation_support_count": int(support_count),
            "generation_support_fraction": min(1.0, support_count / 4.0),
            "assembly_score": float(score),
            "analog_score": float(score),
            "assembly_generation_modes": "|".join(sorted(modes)),
            "operator_set": "|".join(sorted(modes)),
            "operator_family": operator_family,
            "operator_count": int(support_count),
            "operator_selection_failure_bucket": str(failure_bucket),
            "generation_track": "assembly_de_novo",
            "assembly_track": "novel_denovo_assembly_v4",
            "assembly_version": "v4",
            "anchor_source": str(anchor.get("anchor_source", "")),
            "anchor_rank": safe_int(anchor.get("anchor_rank"), 0),
            "anchor_smiles": str(anchor.get("anchor_smiles", "")),
            "anchor_formula": str(anchor.get("anchor_formula", "")),
            "anchor_score": safe_float(anchor.get("anchor_score"), 0.0),
            "anchor_fragment_score": safe_float(anchor.get("anchor_fragment_score"), 0.0),
            "formula_anchor_match": bool(formula_anchor_match),
            "nearest_db_tanimoto": float(nearest_db_tan),
            "novelty_score": float(novelty_score),
            "qed_score": float(plausibility),
            "exact_known_structure": False,
            **components,
        }
        structure_rows.append(
            {
                **payload,
                "structure_only_reason": "profile_mass_pass" if abs_err <= float(tolerance) else "mass_incompatible",
            }
        )
        novelty_rows[-1]["abs_mass_error_da"] = abs_err
        novelty_rows[-1]["profile_mass_pass"] = bool(abs_err <= float(tolerance))
        novelty_rows[-1]["strict_1da_pass"] = bool(abs_err <= 1.0)
        if abs_err > float(tolerance):
            novelty_rows[-1]["rejected_reason"] = "mass_incompatible_for_profile"
            continue
        rows.append(payload)
    broad_cap = broad_cap_for_lib(lib_id)
    non_broad_rows = [row for row in rows if not is_broad_operator(str(row.get("operator_set", "")))]
    broad_rows = [row for row in rows if is_broad_operator(str(row.get("operator_set", "")))]
    broad_rows = sorted(
        broad_rows,
        key=lambda row: (
            -safe_float(row.get("v4_balanced_score"), 0.0),
            safe_float(row.get("abs_mass_error_da"), 9999.0),
            -safe_float(row.get("anchor_tanimoto"), 0.0),
            str(row.get("predicted_smiles", "")),
        ),
    )[:broad_cap]
    rows = non_broad_rows + broad_rows
    rows = sorted(
        rows,
        key=lambda row: (
            -safe_float(row.get("v4_balanced_score"), 0.0),
            safe_float(row.get("abs_mass_error_da"), 9999.0),
            -safe_float(row.get("anchor_tanimoto"), 0.0),
            str(row.get("predicted_smiles", "")),
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
            "assembly_version": "v4",
            "failure_bucket": str(failure_bucket),
        }
        for anchor in anchors
    ]
    return rows, anchor_rows, novelty_rows, structure_rows


def base_profile_name(profile: str) -> str:
    return str(profile).replace("non_x_holdout_", "")


def default_variant_for_profile(profile: str) -> str:
    return DEFAULT_VARIANT_BY_PROFILE.get(base_profile_name(profile), "v4_balanced")


def sort_for_variant(df: pd.DataFrame, variant: str) -> pd.DataFrame:
    if df.empty:
        return df.copy()
    score_col = f"{variant}_score"
    work = df.copy()
    if score_col not in work.columns:
        work[score_col] = pd.to_numeric(work.get("assembly_score", 0.0), errors="coerce").fillna(0.0)
    return work.sort_values(
        [score_col, "abs_mass_error_da", "anchor_tanimoto", "predicted_smiles"],
        ascending=[False, True, False, True],
        kind="mergesort",
    )


def build_denovo_report(candidate_df: pd.DataFrame, *, variant: str, top_k: int) -> pd.DataFrame:
    if candidate_df.empty:
        return pd.DataFrame()
    rows: list[dict[str, Any]] = []
    for feature_key, group in candidate_df.groupby("feature_key", sort=False):
        ranked = sort_for_variant(group, variant)
        for rank, (_, row) in enumerate(ranked.head(int(top_k)).iterrows(), start=1):
            payload = row.to_dict()
            payload["rank"] = int(rank)
            payload["candidate_file_rank"] = int(rank)
            payload["score_variant"] = str(variant)
            payload["assembly_score"] = safe_float(
                payload.get(f"{variant}_score"), safe_float(payload.get("assembly_score"), 0.0)
            )
            payload["mean_log_prob"] = payload["assembly_score"]
            payload["feature_label"] = str(feature_key)
            payload["is_valid"] = True
            payload["candidate_latent_cosine"] = float("nan")
            payload["samples_drawn"] = int(group.shape[0])
            rows.append(payload)
    return pd.DataFrame(rows)


def write_profile_outputs_v4(
    *,
    profile_dir: Path,
    lib_id: int,
    profile: str,
    candidate_rows: list[dict[str, Any]],
    anchor_rows: list[dict[str, Any]],
    novelty_rows: list[dict[str, Any]],
    structure_rows: list[dict[str, Any]],
    top_k: int,
) -> dict[str, str]:
    profile_dir.mkdir(parents=True, exist_ok=True)
    default_variant = default_variant_for_profile(profile)
    candidate_df = pd.DataFrame(candidate_rows)
    if not candidate_df.empty:
        candidate_df = sort_for_variant(candidate_df, default_variant)
        ranked_rows = []
        for feature_key, group in candidate_df.groupby("feature_key", sort=False):
            for rank, (_, row) in enumerate(group.iterrows(), start=1):
                payload = row.to_dict()
                payload["rank"] = int(rank)
                payload["candidate_file_rank"] = int(rank)
                payload["score_variant"] = default_variant
                payload["assembly_score"] = safe_float(
                    payload.get(f"{default_variant}_score"), safe_float(payload.get("assembly_score"), 0.0)
                )
                payload["mean_log_prob"] = payload["assembly_score"]
                ranked_rows.append(payload)
        candidate_df = pd.DataFrame(ranked_rows)
    denovo_df = build_denovo_report(candidate_df, variant=default_variant, top_k=top_k)
    anchor_df = pd.DataFrame(anchor_rows)
    novelty_df = pd.DataFrame(novelty_rows)
    structure_df = pd.DataFrame(structure_rows)
    outputs = {
        "candidate_metrics_csv": str(profile_dir / "candidate_metrics_all.csv"),
        "denovo_csv": str(profile_dir / "denovo_report.csv"),
        "structure_only_pool_csv": str(profile_dir / "structure_only_pool.csv"),
        "anchor_fragment_audit_csv": str(profile_dir / "anchor_fragment_audit.csv"),
        "novelty_audit_csv": str(profile_dir / "novelty_audit.csv"),
        "summary_json": str(profile_dir / "assembly_generation_summary.json"),
    }
    write_csv(Path(outputs["candidate_metrics_csv"]), candidate_df)
    write_csv(profile_dir / f"merged_{int(lib_id)}_candidate_metrics_all.csv", candidate_df)
    write_csv(Path(outputs["denovo_csv"]), denovo_df)
    write_csv(profile_dir / f"merged_{int(lib_id)}_denovo_report.csv", denovo_df)
    write_csv(Path(outputs["structure_only_pool_csv"]), structure_df)
    write_csv(Path(outputs["anchor_fragment_audit_csv"]), anchor_df)
    write_csv(Path(outputs["novelty_audit_csv"]), novelty_df)
    variant_outputs: dict[str, str] = {}
    for variant in RANKING_VARIANTS:
        variant_dir = profile_dir / "ranking_variants" / variant
        variant_dir.mkdir(parents=True, exist_ok=True)
        variant_report = build_denovo_report(candidate_df, variant=variant, top_k=top_k)
        write_csv(variant_dir / "denovo_report.csv", variant_report)
        write_csv(variant_dir / f"merged_{int(lib_id)}_denovo_report.csv", variant_report)
        variant_outputs[variant] = str(variant_dir / "denovo_report.csv")
    tolerance = PROFILE_TOLERANCE.get(base_profile_name(profile), 5.0)
    summary = {
        "profile": str(profile),
        "lib_id": int(lib_id),
        "assembly_version": "v4",
        "default_score_variant": default_variant,
        "ranking_variants": list(RANKING_VARIANTS),
        "mass_tolerance_da": float(tolerance),
        "n_candidate_rows": int(candidate_df.shape[0]),
        "n_structure_only_rows": int(structure_df.shape[0]),
        "n_denovo_rows": int(denovo_df.shape[0]),
        "n_features_with_candidates": int(candidate_df["feature_key"].nunique()) if not candidate_df.empty else 0,
        "n_features_with_denovo_rows": int(denovo_df["feature_key"].nunique()) if not denovo_df.empty else 0,
        "n_novelty_rows_considered": int(novelty_df.shape[0]),
        "n_exact_known_rejected": (
            int(novelty_df["exact_known_structure_rejected"].sum())
            if not novelty_df.empty and "exact_known_structure_rejected" in novelty_df.columns
            else 0
        ),
        "n_mass_incompatible_rejected": (
            int(
                novelty_df.get("rejected_reason", pd.Series(dtype=str))
                .astype(str)
                .eq("mass_incompatible_for_profile")
                .sum()
            )
            if not novelty_df.empty
            else 0
        ),
        "outputs": {**outputs, "ranking_variant_denovo_reports": variant_outputs},
    }
    write_json(Path(outputs["summary_json"]), summary)
    (profile_dir / "selfies_run_summary.json").write_text(
        json.dumps({"mass_tolerance_da": summary["mass_tolerance_da"], "denovo_assembly": summary}, indent=2) + "\n"
    )
    (profile_dir / f"merged_{int(lib_id)}_selfies_run_summary.json").write_text(
        json.dumps({"mass_tolerance_da": summary["mass_tolerance_da"], "denovo_assembly": summary}, indent=2) + "\n"
    )
    return outputs


def run_denovo_assembly_v4(
    *,
    repo_root: Path,
    output_root: Path,
    baseline_root: Path,
    hard_holdout_dir: Path,
    known_x: Path,
    known_x_all: Path,
    unknown_feature_table: Path,
    failure_audit_csv: Path,
    libs: list[int],
    profiles: list[str],
    scopes: list[str],
    max_features_per_scope: int,
    anchor_top_n: int,
    max_candidates_per_feature: int,
    top_k: int,
    all_x_features: bool = False,
) -> dict[str, Any]:
    if not HAS_RDKIT:
        raise RuntimeError("RDKit is required for retrieval-anchored de novo assembly V4.")
    run_root = output_root / "runs" / RUN_NAME
    candidate_source = v1.prepare_fragment_source(repo_root=repo_root, libs=libs)
    exclusion_canonicals, exclusion_inchikeys = v1.load_structure_exclusion_keys(
        candidate_source=candidate_source,
        known_x_all=known_x_all,
        hard_holdout_dir=hard_holdout_dir,
    )
    failure_map = load_failure_map(failure_audit_csv)
    x_targets = load_x_targets(unknown_feature_table, libs, 0)
    if not all_x_features and not x_targets.empty:
        control_xids: set[str] = set()
        for lib_id in libs:
            for control in load_x_controls(
                known_x=known_x,
                known_x_all=known_x_all,
                lib_id=int(lib_id),
                ion_mode=lib_mode(int(lib_id)),
            ):
                control_xids.add(str(control.control_id))
        x_targets = x_targets[x_targets["xid"].astype(str).isin(control_xids)].copy()
    if max_features_per_scope > 0 and not x_targets.empty:
        x_targets = x_targets.groupby("lib_id", group_keys=False).head(int(max_features_per_scope)).copy()
    hard_targets = load_hard_targets(hard_holdout_dir, libs, max_features_per_scope)
    outputs: list[dict[str, Any]] = []
    for lib_id in libs:
        mode = lib_mode(int(lib_id))
        retrieval_df = load_retrieval_report(baseline_root, int(lib_id))
        scoped_targets: list[tuple[str, pd.DataFrame]] = []
        if "x_control" in scopes and not x_targets.empty:
            scoped_targets.append(("x_control", x_targets[x_targets["lib_id"].astype("Int64").eq(int(lib_id))].copy()))
        if "hard_holdout" in scopes and not hard_targets.empty:
            scoped_targets.append(
                ("hard_holdout", hard_targets[hard_targets["lib_id"].astype("Int64").eq(int(lib_id))].copy())
            )
        for scope, targets in scoped_targets:
            if targets.empty:
                continue
            for profile in profiles:
                output_profile = profile if scope == "x_control" else f"non_x_holdout_{profile}"
                profile_dir = run_root / f"lib_{int(lib_id)}_{mode}" / "selfies_profiles" / output_profile
                candidate_rows: list[dict[str, Any]] = []
                anchor_rows: list[dict[str, Any]] = []
                novelty_rows: list[dict[str, Any]] = []
                structure_rows: list[dict[str, Any]] = []
                for _, row in targets.iterrows():
                    target = row.copy()
                    target["holdout_scope"] = "x_control" if scope == "x_control" else "non_x_known"
                    anchors = (
                        anchors_for_x(target, retrieval_df, anchor_top_n)
                        if scope == "x_control"
                        else anchors_for_hard(target, candidate_source, anchor_top_n)
                    )
                    failure_bucket = failure_bucket_for_target(failure_map, target, scope, profile)
                    rows, anchors_out, novelty_out, structure_out = score_generated_candidates_v4(
                        target=target,
                        anchors=anchors,
                        candidate_source=candidate_source,
                        exclusion_canonicals=exclusion_canonicals,
                        exclusion_inchikeys=exclusion_inchikeys,
                        profile=profile,
                        max_candidates_per_feature=max_candidates_per_feature,
                        failure_bucket=failure_bucket,
                    )
                    for candidate in rows:
                        candidate["holdout_scope"] = str(target["holdout_scope"])
                        candidate["assembly_profile_scope"] = scope
                    for candidate in structure_out:
                        candidate["holdout_scope"] = str(target["holdout_scope"])
                        candidate["assembly_profile_scope"] = scope
                    candidate_rows.extend(rows)
                    anchor_rows.extend(anchors_out)
                    novelty_rows.extend(novelty_out)
                    structure_rows.extend(structure_out)
                profile_outputs = write_profile_outputs_v4(
                    profile_dir=profile_dir,
                    lib_id=int(lib_id),
                    profile=output_profile,
                    candidate_rows=candidate_rows,
                    anchor_rows=anchor_rows,
                    novelty_rows=novelty_rows,
                    structure_rows=structure_rows,
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
        "all_x_features": bool(all_x_features),
        "x_target_rows": int(x_targets.shape[0]),
        "hard_target_rows": int(hard_targets.shape[0]),
        "known_structure_canonical_exclusions": int(len(exclusion_canonicals)),
        "known_structure_inchikey_exclusions": int(len(exclusion_inchikeys)),
        "failure_audit_csv": str(failure_audit_csv),
        "failure_map_rows": int(len(failure_map)),
        "has_selfies": bool(v1.HAS_SELFIES),
        "outputs": outputs,
    }
    write_json(output_root / "assembly_v4_generation_summary.json", payload)
    return payload


def v4_acceptance(metrics: pd.DataFrame) -> dict[str, Any]:
    rows = []
    if metrics.empty:
        return {"accepted_rows": [], "next_branch": "No V4 metrics were produced."}
    assembly = metrics[metrics.get("track", pd.Series(dtype=str)).astype(str).eq("novel_denovo_assembly_v4")].copy()
    if assembly.empty:
        assembly = metrics[metrics.get("track", pd.Series(dtype=str)).astype(str).eq("novel_denovo_assembly")].copy()
    for _, row in assembly.iterrows():
        if str(row.get("holdout_scope")) != "x_control" or str(row.get("metric_scope")) not in {
            "physical_full_pool_oracle",
            "top10",
        }:
            continue
        lib_id = safe_int(row.get("lib_id"), -1)
        profile = str(row.get("base_profile", "") or row.get("profile", ""))
        mean = safe_float(row.get("mean_best_tanimoto"), 0.0)
        ge07 = safe_int(row.get("n_tanimoto_ge_0p7"), 0)
        accepted = False
        reason = ""
        if lib_id == 305 and profile == "strict_1da":
            accepted = mean > 0.4370 and ge07 >= 2
            reason = "beat_305_strict_analog_v1_oracle_reference"
        elif lib_id == 305 and profile == "balanced":
            accepted = mean > 0.4486 and ge07 >= 1
            reason = "beat_305_balanced_analog_v1_oracle_reference_without_high_hit_loss"
        elif lib_id == 209 and profile == "strict_1da":
            accepted = mean >= 0.3631
            reason = "meets_209_strict_baseline_floor"
        else:
            accepted = ge07 >= 1 and mean >= 0.3
            reason = "diagnostic_or_secondary_gate"
        rows.append({**row.to_dict(), "accepted_against_reference": bool(accepted), "reference_reason": reason})
    accepted_rows = [row for row in rows if row.get("accepted_against_reference")]
    next_branch = (
        "V4 improved against references; inspect score variants and promote only if top-10 clears retrieval/analog gates."
        if accepted_rows
        else "V4 did not beat historical oracle references; inspect operator-level failures before more generation."
    )
    return {"accepted_rows": accepted_rows, "reference_rows": rows, "next_branch": next_branch}


def write_v4_aggregate_aliases(output_root: Path, aggregate_payload: dict[str, Any]) -> dict[str, Any]:
    aggregate_dir = output_root / "aggregate"
    metrics_path = aggregate_dir / "denovo_assembly_v4_metrics.csv"
    deltas_path = aggregate_dir / "denovo_assembly_v4_deltas.csv"
    per_path = aggregate_dir / "denovo_assembly_v4_per_control.csv"
    metrics = read_table(metrics_path) if metrics_path.exists() else pd.DataFrame()
    if not metrics.empty and "track" in metrics.columns:
        metrics.loc[metrics["track"].astype(str).eq("novel_denovo_assembly"), "track"] = "novel_denovo_assembly_v4"
    deltas = read_table(deltas_path) if deltas_path.exists() else pd.DataFrame()
    per_control = read_table(per_path) if per_path.exists() else pd.DataFrame()
    recommendations = v4_acceptance(metrics)
    write_csv(aggregate_dir / "denovo_assembly_v4_metrics.csv", metrics)
    write_csv(aggregate_dir / "denovo_assembly_v4_deltas.csv", deltas)
    write_csv(aggregate_dir / "denovo_assembly_v4_per_control.csv", per_control)
    write_json(aggregate_dir / "next_branch_recommendations.json", recommendations)
    display_cols = [
        col
        for col in [
            "track",
            "lib_id",
            "profile",
            "holdout_scope",
            "metric_scope",
            "n_controls_covered",
            "mean_best_tanimoto",
            "n_tanimoto_ge_0p7",
        ]
        if col in metrics.columns
    ]
    report = "\n".join(
        [
            "# Retrieval-Anchored De Novo Assembly V4",
            "",
            f"Recommendation: {recommendations.get('next_branch', '')}",
            "",
            "## Track Metrics",
            dataframe_to_markdown(metrics[display_cols] if display_cols else metrics),
            "",
            "## Acceptance Rows",
            dataframe_to_markdown(pd.DataFrame(recommendations.get("reference_rows", []))),
            "",
            "## V4 Aggregator Payload",
            "```json",
            json.dumps(aggregate_payload, indent=2, sort_keys=True),
            "```",
            "",
        ]
    )
    (aggregate_dir / "denovo_assembly_v4_report.md").write_text(report)
    summary = {
        "metrics_csv": str(aggregate_dir / "denovo_assembly_v4_metrics.csv"),
        "deltas_csv": str(aggregate_dir / "denovo_assembly_v4_deltas.csv"),
        "per_control_csv": str(aggregate_dir / "denovo_assembly_v4_per_control.csv"),
        "report_md": str(aggregate_dir / "denovo_assembly_v4_report.md"),
        "recommendations_json": str(aggregate_dir / "next_branch_recommendations.json"),
        "recommendations": recommendations,
    }
    write_json(aggregate_dir / "denovo_assembly_v4_aggregate_summary.json", summary)
    return summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--baseline-root", type=Path, default=DEFAULT_BASELINE_ROOT)
    parser.add_argument("--analog-root", type=Path, default=DEFAULT_ANALOG_ROOT)
    parser.add_argument("--ablation-root", type=Path, default=DEFAULT_ABLATION_ROOT)
    parser.add_argument("--hard-holdout-dir", type=Path, default=DEFAULT_HARD_HOLDOUT_DIR)
    parser.add_argument("--failure-audit-csv", type=Path, default=DEFAULT_FAILURE_AUDIT)
    parser.add_argument("--known-x", type=Path, default=Path("known_X.csv"))
    parser.add_argument("--known-x-all", type=Path, default=Path("known_X_all_available_structures_simple.csv"))
    parser.add_argument(
        "--unknown-feature-table", type=Path, default=Path("merged_encoder/data/ions_unknown_merged.csv")
    )
    parser.add_argument("--targets", type=str, default="")
    parser.add_argument("--libs", type=str, default="305,209")
    parser.add_argument("--profiles", type=str, default="strict_1da,balanced")
    parser.add_argument("--scopes", type=str, default="x_control,hard_holdout")
    parser.add_argument("--max-features-per-scope", type=int, default=0)
    parser.add_argument("--anchor-top-n", type=int, default=5)
    parser.add_argument("--max-candidates-per-feature", type=int, default=300)
    parser.add_argument("--top-k", type=int, default=10)
    parser.add_argument("--all-x-features", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    repo_root = REPO_ROOT
    output_root = resolve(args.output_root, repo_root)
    baseline_root = resolve(args.baseline_root, repo_root)
    analog_root = resolve(args.analog_root, repo_root)
    ablation_root = resolve(args.ablation_root, repo_root)
    hard_holdout_dir = resolve(args.hard_holdout_dir, repo_root)
    known_x = resolve(args.known_x, repo_root)
    known_x_all = resolve(args.known_x_all, repo_root)
    unknown_feature_table = resolve(args.unknown_feature_table, repo_root)
    failure_audit_csv = resolve(args.failure_audit_csv, repo_root)
    target_libs, target_profiles = parse_targets(args.targets)
    libs = target_libs or parse_libs(args.libs)
    profiles = target_profiles or parse_csv(args.profiles)
    generation = run_denovo_assembly_v4(
        repo_root=repo_root,
        output_root=output_root,
        baseline_root=baseline_root,
        hard_holdout_dir=hard_holdout_dir,
        known_x=known_x,
        known_x_all=known_x_all,
        unknown_feature_table=unknown_feature_table,
        failure_audit_csv=failure_audit_csv,
        libs=libs,
        profiles=profiles,
        scopes=parse_csv(args.scopes),
        max_features_per_scope=int(args.max_features_per_scope),
        anchor_top_n=int(args.anchor_top_n),
        max_candidates_per_feature=int(args.max_candidates_per_feature),
        top_k=int(args.top_k),
        all_x_features=bool(args.all_x_features),
    )
    print(json.dumps({"run_root": generation["run_root"], "outputs": len(generation["outputs"])}, indent=2))


if __name__ == "__main__":
    main()
