#!/usr/bin/env python3
"""Retrieval-anchored de novo assembly for X-domain candidates.

This branch uses retrieval hits as constraints and fragment/scaffold priors,
but final candidates must be novel with respect to the known-structure assets
and holdout/control structures.
"""

from __future__ import annotations

from structural_elucidation.paths import REPO_ROOT

import argparse
import json
import math
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd

from structural_elucidation.evaluation.aggregate import (
    build_deltas,
    build_hard_x_gap,
    evaluate_profile,
    load_baseline_metrics,
)
from structural_elucidation.utils.denovo import (
    canonicalize_smiles,
    exact_mass_from_smiles,
    formula_from_smiles,
    inchikey_from_smiles,
    read_table,
)
from structural_elucidation.evaluation.controls import HAS_RDKIT
from structural_elucidation.evaluation.failure_modes import load_non_x_holdouts, load_x_controls
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
    load_training_structures,
    parse_libs,
    safe_float,
    safe_int,
    write_csv,
    write_json,
)

if HAS_RDKIT:
    from rdkit import Chem, DataStructs, RDLogger
    from rdkit.Chem import AllChem, BRICS, QED
    from rdkit.Chem.MolStandardize import rdMolStandardize

    RDLogger.DisableLog("rdApp.*")

try:  # SELFIES is optional in this workspace.
    import selfies as sf

    HAS_SELFIES = True
except Exception:  # pragma: no cover - optional dependency
    sf = None
    HAS_SELFIES = False


DEFAULT_OUTPUT_ROOT = Path("models/denovo_assembly_v1")
DEFAULT_BASELINE_ROOT = Path("models/merged_structuremasst_bundle_excl_x_unresolved_2026-04-17_metabolon")
DEFAULT_ANALOG_ROOT = Path("models/denovo_x_domain_analog_v1")
DEFAULT_HARD_HOLDOUT_DIR = Path("models/denovo_x_domain_analog_v1/hard_holdouts")
DEFAULT_ABLATION_ROOT = Path("models/decoder_ablation_experiments")
GOOD_TANIMOTO = 0.7


@dataclass(frozen=True)
class GeneratedCandidate:
    smiles: str
    canonical_smiles: str
    inchikey: str
    generation_modes: tuple[str, ...]
    anchor_smiles: str
    anchor_rank: int
    anchor_score: float
    anchor_formula: str
    anchor_source: str


def canonical_or_empty(smiles: Any) -> str:
    return canonicalize_smiles(smiles) or ""


def mol_from_smiles(smiles: Any):
    if not HAS_RDKIT:
        return None
    canonical = canonical_or_empty(smiles)
    if not canonical:
        return None
    return Chem.MolFromSmiles(canonical)


def canonical_from_mol(mol: Any) -> str:
    if not HAS_RDKIT or mol is None:
        return ""
    try:
        Chem.SanitizeMol(mol)
        return Chem.MolToSmiles(mol, isomericSmiles=True)
    except Exception:
        return ""


def atom_substitution_candidates(anchor_smiles: str, limit: int = 80) -> set[str]:
    """Generate small valid atom substitutions around an anchor."""
    mol = mol_from_smiles(anchor_smiles)
    if mol is None:
        return set()
    replacements = {
        6: [7, 8, 16],
        7: [6, 8],
        8: [6, 7, 16],
        16: [6, 8],
        9: [17, 35],
        17: [9, 35],
        35: [9, 17],
    }
    out: set[str] = set()
    for atom in mol.GetAtoms():
        choices = replacements.get(atom.GetAtomicNum(), [])
        for atomic_num in choices:
            rw = Chem.RWMol(mol)
            rw_atom = rw.GetAtomWithIdx(atom.GetIdx())
            rw_atom.SetAtomicNum(int(atomic_num))
            rw_atom.SetFormalCharge(0)
            canonical = canonical_from_mol(rw.GetMol())
            if canonical and canonical != canonical_or_empty(anchor_smiles):
                out.add(canonical)
            if len(out) >= int(limit):
                return out
    return out


def bond_order_candidates(anchor_smiles: str, limit: int = 60) -> set[str]:
    """Generate conservative single/double bond perturbations."""
    mol = mol_from_smiles(anchor_smiles)
    if mol is None:
        return set()
    out: set[str] = set()
    for bond in mol.GetBonds():
        if bond.IsInRing():
            continue
        current = bond.GetBondType()
        if current not in {Chem.BondType.SINGLE, Chem.BondType.DOUBLE}:
            continue
        new_type = Chem.BondType.DOUBLE if current == Chem.BondType.SINGLE else Chem.BondType.SINGLE
        rw = Chem.RWMol(mol)
        rw.RemoveBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
        rw.AddBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), new_type)
        canonical = canonical_from_mol(rw.GetMol())
        if canonical and canonical != canonical_or_empty(anchor_smiles):
            out.add(canonical)
        if len(out) >= int(limit):
            return out
    return out


def tautomer_candidates(anchor_smiles: str, limit: int = 80) -> set[str]:
    """Enumerate formula-preserving tautomers where RDKit can do so."""
    mol = mol_from_smiles(anchor_smiles)
    if mol is None:
        return set()
    out: set[str] = set()
    try:
        enumerator = rdMolStandardize.TautomerEnumerator()
        for tautomer in enumerator.Enumerate(mol):
            canonical = canonical_from_mol(tautomer)
            if canonical and canonical != canonical_or_empty(anchor_smiles):
                out.add(canonical)
            if len(out) >= int(limit):
                break
    except Exception:
        return out
    return out


def selfies_mutation_candidates(anchor_smiles: str, token_pool: list[str], limit: int = 120) -> set[str]:
    """Use SELFIES token edits when the optional package is available."""
    if not HAS_SELFIES:
        return set()
    try:
        encoded = sf.encoder(anchor_smiles)
        tokens = list(sf.split_selfies(encoded))
    except Exception:
        return set()
    if not tokens:
        return set()
    pool = [token for token in token_pool if token]
    if not pool:
        pool = sorted(set(tokens))
    out: set[str] = set()
    for idx in range(len(tokens)):
        for token in pool[:24]:
            if token == tokens[idx]:
                continue
            edited = list(tokens)
            edited[idx] = token
            try:
                canonical = canonical_or_empty(sf.decoder("".join(edited)))
            except Exception:
                canonical = ""
            if canonical and canonical != canonical_or_empty(anchor_smiles):
                out.add(canonical)
            if len(out) >= int(limit):
                return out
    for idx in range(len(tokens) + 1):
        for token in pool[:12]:
            edited = list(tokens)
            edited.insert(idx, token)
            try:
                canonical = canonical_or_empty(sf.decoder("".join(edited)))
            except Exception:
                canonical = ""
            if canonical and canonical != canonical_or_empty(anchor_smiles):
                out.add(canonical)
            if len(out) >= int(limit):
                return out
    return out


def brics_recombination_candidates(
    anchor_smiles: str,
    source_smiles: Iterable[str],
    limit: int = 160,
) -> set[str]:
    """Build novel molecules from anchor and same-lib BRICS fragments."""
    if not HAS_RDKIT:
        return set()
    fragments: set[str] = set()
    effective_limit = min(int(limit), 80)
    for smiles in [anchor_smiles, *list(source_smiles)[:8]]:
        mol = mol_from_smiles(smiles)
        if mol is None:
            continue
        try:
            fragments.update(BRICS.BRICSDecompose(mol))
        except Exception:
            continue
    if len(fragments) < 2:
        return set()
    out: set[str] = set()
    try:
        for mol in BRICS.BRICSBuild(list(sorted(fragments))[:40], maxDepth=2):
            canonical = canonical_from_mol(mol)
            if canonical and canonical != canonical_or_empty(anchor_smiles):
                out.add(canonical)
            if len(out) >= effective_limit:
                break
    except Exception:
        return out
    return out


def collect_selfies_tokens(smiles_values: Iterable[str], limit: int = 64) -> list[str]:
    if not HAS_SELFIES:
        return []
    tokens: list[str] = []
    seen: set[str] = set()
    for smiles in smiles_values:
        try:
            encoded = sf.encoder(smiles)
            parts = list(sf.split_selfies(encoded))
        except Exception:
            continue
        for token in parts:
            if token not in seen:
                seen.add(token)
                tokens.append(token)
            if len(tokens) >= int(limit):
                return tokens
    return tokens


def load_structure_exclusion_keys(
    *,
    candidate_source: pd.DataFrame,
    known_x_all: Path,
    hard_holdout_dir: Path,
) -> tuple[set[str], set[str]]:
    canonical_keys = {
        str(value)
        for value in candidate_source.get("canonical_smiles_key", pd.Series(dtype=str)).astype(str)
        if str(value)
    }
    inchikeys = {
        str(value) for value in candidate_source.get("inchikey_key", pd.Series(dtype=str)).astype(str) if str(value)
    }
    if known_x_all.exists():
        truth = read_table(known_x_all)
        for col in ["smiles", "SMILES", "canonical_smiles"]:
            if col in truth.columns:
                for value in truth[col].dropna().tolist():
                    canonical = canonical_or_empty(value)
                    if canonical:
                        canonical_keys.add(canonical)
                        key = inchikey_from_smiles(canonical)
                        if key:
                            inchikeys.add(key)
                break
    for path in sorted(hard_holdout_dir.glob("*.csv")):
        try:
            holdouts = read_table(path)
        except Exception:
            continue
        for col in ["canonical_smiles", "smiles", "SMILES"]:
            if col in holdouts.columns:
                for value in holdouts[col].dropna().tolist():
                    canonical = canonical_or_empty(value)
                    if canonical:
                        canonical_keys.add(canonical)
                        key = inchikey_from_smiles(canonical)
                        if key:
                            inchikeys.add(key)
                break
        if "inchikey" in holdouts.columns:
            inchikeys.update(
                str(value)
                for value in holdouts["inchikey"].dropna().astype(str)
                if str(value) and str(value).lower() != "nan"
            )
    return canonical_keys, inchikeys


def prepare_fragment_source(
    *,
    repo_root: Path,
    libs: list[int],
) -> pd.DataFrame:
    source = load_training_structures(repo_root, libs)
    if source.empty:
        return source
    work = source.copy()
    work["canonical_smiles_key"] = work["canonical_smiles_key"].fillna("").astype(str)
    work = work[work["canonical_smiles_key"].str.strip().ne("")].copy()
    work["exact_mass"] = work["canonical_smiles_key"].map(exact_mass_from_smiles)
    work["candidate_formula"] = work["canonical_smiles_key"].map(formula_from_smiles)
    work["fp"] = work["canonical_smiles_key"].map(fp_from_smiles)
    work["exact_mass"] = pd.to_numeric(work["exact_mass"], errors="coerce")
    return work.dropna(subset=["exact_mass", "fp"]).copy()


def nearest_db_similarity(fp: Any, db_fps: list[Any]) -> float:
    if fp is None or not db_fps:
        return 0.0
    try:
        return float(max(DataStructs.BulkTanimotoSimilarity(fp, db_fps)))
    except Exception:
        return 0.0


def max_anchor_similarity(fp: Any, anchor_fps: list[Any]) -> float:
    if fp is None or not anchor_fps:
        return 0.0
    try:
        return float(max(DataStructs.BulkTanimotoSimilarity(fp, anchor_fps)))
    except Exception:
        return 0.0


def qed_score(smiles: str) -> float:
    mol = mol_from_smiles(smiles)
    if mol is None:
        return 0.0
    try:
        return float(QED.qed(mol))
    except Exception:
        return 0.0


def select_fragment_source_smiles(
    *,
    candidate_source: pd.DataFrame,
    lib_id: int,
    parent_mz: float,
    ion_mode: str,
    anchor_fps: list[Any],
    limit: int = 24,
) -> list[str]:
    scoped = candidate_source[candidate_source["lib_id"].astype(int).eq(int(lib_id))].copy()
    if scoped.empty:
        return []
    exact_masses = pd.to_numeric(scoped["exact_mass"], errors="coerce").to_numpy(dtype=float)
    mass_errors = []
    for mass in exact_masses:
        match = mass_match(float(mass), float(parent_mz), ion_mode)
        mass_errors.append(safe_float(match.get("abs_mass_error_da"), 9999.0))
    if anchor_fps:
        source_fps = list(scoped["fp"])
        matrix = [DataStructs.BulkTanimotoSimilarity(anchor_fp, source_fps) for anchor_fp in anchor_fps]
        anchor_scores = np.max(np.asarray(matrix, dtype=float), axis=0)
    else:
        anchor_scores = np.zeros(scoped.shape[0], dtype=float)
    work = scoped.copy()
    work["_mass_error"] = mass_errors
    work["_anchor_similarity"] = anchor_scores
    work["_fragment_source_score"] = -0.05 * work["_mass_error"] + work["_anchor_similarity"]
    selected = work.sort_values(
        ["_fragment_source_score", "_anchor_similarity", "_mass_error"],
        ascending=[False, False, True],
        kind="mergesort",
    ).head(int(limit))
    return selected["canonical_smiles_key"].astype(str).tolist()


def mass_near_known_seed_anchors(
    *,
    candidate_source: pd.DataFrame,
    lib_id: int,
    parent_mz: float,
    ion_mode: str,
    tolerance_da: float,
    limit: int = 5,
) -> list[dict[str, Any]]:
    """Use known structures as generation seeds, never as final candidates."""
    scoped = candidate_source[candidate_source["lib_id"].astype(int).eq(int(lib_id))].copy()
    if scoped.empty:
        return []
    rows: list[dict[str, Any]] = []
    for _, row in scoped.iterrows():
        mass = safe_float(row.get("exact_mass"), math.nan)
        if not np.isfinite(mass):
            continue
        match = mass_match(mass, parent_mz, ion_mode)
        abs_err = safe_float(match.get("abs_mass_error_da"), math.inf)
        if abs_err > float(tolerance_da):
            continue
        rows.append(
            {
                "anchor_source": "mass_near_known_seed",
                "anchor_smiles": str(row.get("canonical_smiles_key", "")),
                "anchor_rank": 0,
                "anchor_score": mass_score(abs_err, max(float(tolerance_da), 1.0)),
                "anchor_formula": str(row.get("candidate_formula", "") or row.get("formula", "")),
                "anchor_fragment_score": 0.0,
                "anchor_fp": row.get("fp"),
                "_abs_mass_error_da": abs_err,
            }
        )
    rows = sorted(
        rows,
        key=lambda item: (
            safe_float(item.get("_abs_mass_error_da"), math.inf),
            str(item.get("anchor_smiles", "")),
        ),
    )[: int(limit)]
    for rank, row in enumerate(rows, start=1):
        row["anchor_rank"] = int(rank)
    return rows


def generate_for_anchor(
    *,
    anchor: dict[str, Any],
    fragment_source_smiles: list[str],
    token_pool: list[str],
    per_anchor_limit: int,
) -> dict[str, set[str]]:
    anchor_smiles = canonical_or_empty(anchor.get("anchor_smiles", ""))
    if not anchor_smiles:
        return {}
    generated: dict[str, set[str]] = defaultdict(set)
    generated[anchor_smiles].add("anchor_identity")
    for smiles in tautomer_candidates(anchor_smiles, limit=max(20, per_anchor_limit // 2)):
        generated[smiles].add("tautomer_enumeration")
    for smiles in atom_substitution_candidates(anchor_smiles, limit=max(20, per_anchor_limit // 2)):
        generated[smiles].add("atom_substitution")
    for smiles in bond_order_candidates(anchor_smiles, limit=max(12, per_anchor_limit // 4)):
        generated[smiles].add("bond_order_mutation")
    for smiles in selfies_mutation_candidates(anchor_smiles, token_pool, limit=max(30, per_anchor_limit)):
        generated[smiles].add("selfies_mutation")
    for smiles in brics_recombination_candidates(
        anchor_smiles, fragment_source_smiles, limit=max(40, per_anchor_limit)
    ):
        generated[smiles].add("brics_recombination")
    return generated


def score_generated_candidates(
    *,
    target: pd.Series,
    anchors: list[dict[str, Any]],
    candidate_source: pd.DataFrame,
    exclusion_canonicals: set[str],
    exclusion_inchikeys: set[str],
    profile: str,
    max_candidates_per_feature: int,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    lib_id = int(target.get("lib_id"))
    ion_mode = str(target.get("ion_mode") or lib_mode(lib_id))
    parent_mz = safe_float(target.get("parent_mz"), 0.0)
    tolerance = PROFILE_TOLERANCE.get(str(profile).replace("non_x_holdout_", ""), 5.0)
    seed_anchors = mass_near_known_seed_anchors(
        candidate_source=candidate_source,
        lib_id=lib_id,
        parent_mz=parent_mz,
        ion_mode=ion_mode,
        tolerance_da=tolerance,
        limit=5,
    )
    anchors = list(anchors) + [
        seed
        for seed in seed_anchors
        if str(seed.get("anchor_smiles", ""))
        and str(seed.get("anchor_smiles", "")) not in {str(anchor.get("anchor_smiles", "")) for anchor in anchors}
    ]
    db_fps = list(candidate_source[candidate_source["lib_id"].astype(int).eq(lib_id)]["fp"])
    anchor_fps = [anchor.get("anchor_fp") for anchor in anchors if anchor.get("anchor_fp") is not None]
    anchor_formulas = {
        str(anchor.get("anchor_formula", "")) for anchor in anchors if str(anchor.get("anchor_formula", ""))
    }
    fragment_source_smiles = select_fragment_source_smiles(
        candidate_source=candidate_source,
        lib_id=lib_id,
        parent_mz=parent_mz,
        ion_mode=ion_mode,
        anchor_fps=anchor_fps,
    )
    token_pool = collect_selfies_tokens(
        [str(anchor.get("anchor_smiles", "")) for anchor in anchors] + fragment_source_smiles
    )
    generated_modes: dict[str, set[str]] = defaultdict(set)
    anchor_by_candidate: dict[str, dict[str, Any]] = {}
    per_anchor_limit = max(int(max_candidates_per_feature), 80)
    for anchor in anchors:
        for smiles, modes in generate_for_anchor(
            anchor=anchor,
            fragment_source_smiles=fragment_source_smiles,
            token_pool=token_pool,
            per_anchor_limit=per_anchor_limit,
        ).items():
            generated_modes[smiles].update(modes)
            anchor_by_candidate.setdefault(smiles, anchor)
    rows: list[dict[str, Any]] = []
    novelty_rows: list[dict[str, Any]] = []
    for smiles, modes in generated_modes.items():
        canonical = canonical_or_empty(smiles)
        if not canonical:
            continue
        key = inchikey_from_smiles(canonical)
        exact_known = canonical in exclusion_canonicals or (key and key in exclusion_inchikeys)
        fp = fp_from_smiles(canonical)
        if fp is None:
            continue
        nearest_db_tan = nearest_db_similarity(fp, db_fps)
        novelty_rows.append(
            {
                "feature_key": str(target.get("feature_key", "")),
                "comp_id": safe_int(target.get("comp_id"), -1),
                "candidate_smiles": canonical,
                "inchikey": key,
                "exact_known_structure_rejected": bool(exact_known),
                "nearest_db_tanimoto": nearest_db_tan,
                "generation_modes": "|".join(sorted(modes)),
            }
        )
        if exact_known:
            continue
        exact_mass = exact_mass_from_smiles(canonical)
        if not np.isfinite(safe_float(exact_mass)):
            continue
        match = mass_match(float(exact_mass), parent_mz, ion_mode)
        abs_err = safe_float(match.get("abs_mass_error_da"), math.inf)
        if abs_err > float(tolerance):
            continue
        anchor = anchor_by_candidate.get(smiles, anchors[0] if anchors else {})
        anchor_tan = max_anchor_similarity(fp, anchor_fps)
        formula = formula_from_smiles(canonical)
        formula_anchor_match = bool(formula and formula in anchor_formulas)
        plausibility = qed_score(canonical)
        novelty_score = max(0.0, 1.0 - min(nearest_db_tan, 1.0))
        support_count = len(modes)
        score = (
            0.40 * mass_score(abs_err, max(float(tolerance), 1.0))
            + 0.24 * anchor_tan
            + 0.10 * safe_float(anchor.get("anchor_score"), 0.0)
            + 0.08 * min(1.0, support_count / 3.0)
            + 0.08 * plausibility
            + (0.06 if formula_anchor_match else 0.0)
            + 0.04 * novelty_score
        )
        rows.append(
            {
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
                "assembly_track": "novel_denovo_assembly",
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
            }
        )
    rows = sorted(
        rows,
        key=lambda row: (
            -safe_float(row.get("assembly_score"), 0.0),
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
        }
        for anchor in anchors
    ]
    return rows, anchor_rows, novelty_rows


def write_profile_outputs(
    *,
    profile_dir: Path,
    lib_id: int,
    profile: str,
    candidate_rows: list[dict[str, Any]],
    anchor_rows: list[dict[str, Any]],
    novelty_rows: list[dict[str, Any]],
    top_k: int,
) -> dict[str, str]:
    profile_dir.mkdir(parents=True, exist_ok=True)
    candidate_df = pd.DataFrame(candidate_rows)
    if not candidate_df.empty:
        candidate_df = candidate_df.sort_values(["feature_key", "candidate_file_rank"], kind="mergesort")
    denovo_rows: list[dict[str, Any]] = []
    if not candidate_df.empty:
        for feature_key, group in candidate_df.groupby("feature_key", sort=False):
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
    novelty_df = pd.DataFrame(novelty_rows)
    outputs = {
        "candidate_metrics_csv": str(profile_dir / "candidate_metrics_all.csv"),
        "denovo_csv": str(profile_dir / "denovo_report.csv"),
        "anchor_fragment_audit_csv": str(profile_dir / "anchor_fragment_audit.csv"),
        "novelty_audit_csv": str(profile_dir / "novelty_audit.csv"),
        "summary_json": str(profile_dir / "assembly_generation_summary.json"),
    }
    write_csv(Path(outputs["candidate_metrics_csv"]), candidate_df)
    write_csv(profile_dir / f"merged_{int(lib_id)}_candidate_metrics_all.csv", candidate_df)
    write_csv(Path(outputs["denovo_csv"]), denovo_df)
    write_csv(profile_dir / f"merged_{int(lib_id)}_denovo_report.csv", denovo_df)
    write_csv(Path(outputs["anchor_fragment_audit_csv"]), anchor_df)
    write_csv(Path(outputs["novelty_audit_csv"]), novelty_df)
    summary = {
        "profile": str(profile),
        "lib_id": int(lib_id),
        "mass_tolerance_da": float(PROFILE_TOLERANCE.get(str(profile).replace("non_x_holdout_", ""), 5.0)),
        "n_candidate_rows": int(candidate_df.shape[0]),
        "n_denovo_rows": int(denovo_df.shape[0]),
        "n_features_with_candidates": int(candidate_df["feature_key"].nunique()) if not candidate_df.empty else 0,
        "n_features_with_denovo_rows": int(denovo_df["feature_key"].nunique()) if not denovo_df.empty else 0,
        "n_novelty_rows_considered": int(novelty_df.shape[0]),
        "n_exact_known_rejected": (
            int(novelty_df["exact_known_structure_rejected"].sum()) if not novelty_df.empty else 0
        ),
        "outputs": outputs,
    }
    write_json(Path(outputs["summary_json"]), summary)
    (profile_dir / "selfies_run_summary.json").write_text(
        json.dumps({"mass_tolerance_da": summary["mass_tolerance_da"], "denovo_assembly": summary}, indent=2) + "\n"
    )
    (profile_dir / f"merged_{int(lib_id)}_selfies_run_summary.json").write_text(
        json.dumps({"mass_tolerance_da": summary["mass_tolerance_da"], "denovo_assembly": summary}, indent=2) + "\n"
    )
    return outputs


def run_denovo_assembly(
    *,
    repo_root: Path,
    output_root: Path,
    baseline_root: Path,
    hard_holdout_dir: Path,
    known_x: Path,
    known_x_all: Path,
    unknown_feature_table: Path,
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
        raise RuntimeError("RDKit is required for retrieval-anchored de novo assembly.")
    run_root = output_root / "runs" / "retrieval_anchored_denovo_assembly_v1"
    candidate_source = prepare_fragment_source(repo_root=repo_root, libs=libs)
    exclusion_canonicals, exclusion_inchikeys = load_structure_exclusion_keys(
        candidate_source=candidate_source,
        known_x_all=known_x_all,
        hard_holdout_dir=hard_holdout_dir,
    )
    x_targets = load_x_targets(unknown_feature_table, libs, max_features_per_scope)
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
                for _, row in targets.iterrows():
                    target = row.copy()
                    target["holdout_scope"] = "x_control" if scope == "x_control" else "non_x_known"
                    anchors = (
                        anchors_for_x(target, retrieval_df, anchor_top_n)
                        if scope == "x_control"
                        else anchors_for_hard(target, candidate_source, anchor_top_n)
                    )
                    rows, anchors_out, novelty_out = score_generated_candidates(
                        target=target,
                        anchors=anchors,
                        candidate_source=candidate_source,
                        exclusion_canonicals=exclusion_canonicals,
                        exclusion_inchikeys=exclusion_inchikeys,
                        profile=profile,
                        max_candidates_per_feature=max_candidates_per_feature,
                    )
                    for candidate in rows:
                        candidate["holdout_scope"] = str(target["holdout_scope"])
                        candidate["assembly_profile_scope"] = scope
                    candidate_rows.extend(rows)
                    anchor_rows.extend(anchors_out)
                    novelty_rows.extend(novelty_out)
                profile_outputs = write_profile_outputs(
                    profile_dir=profile_dir,
                    lib_id=int(lib_id),
                    profile=output_profile,
                    candidate_rows=candidate_rows,
                    anchor_rows=anchor_rows,
                    novelty_rows=novelty_rows,
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
        "has_selfies": bool(HAS_SELFIES),
        "outputs": outputs,
    }
    write_json(output_root / "assembly_generation_summary.json", payload)
    return payload


def evaluate_retrieval_track(
    *,
    baseline_root: Path,
    known_x: Path,
    known_x_all: Path,
    libs: list[int],
    top_k: int,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for lib_id in libs:
        mode = lib_mode(int(lib_id))
        retrieval_df = load_retrieval_report(baseline_root, int(lib_id))
        if retrieval_df.empty:
            continue
        pool = retrieval_df.copy()
        if "predicted_smiles" not in pool.columns and "smiles" in pool.columns:
            pool["predicted_smiles"] = pool["smiles"]
        temp = DEFAULT_OUTPUT_ROOT / "__unused__"
        controls = load_x_controls(known_x=known_x, known_x_all=known_x_all, lib_id=int(lib_id), ion_mode=mode)
        from structural_elucidation.evaluation.failure_modes import (
            annotate_pool,
            normalize_prediction_pool,
            pool_for_control,
        )

        normalized = normalize_prediction_pool(pool, baseline_root)
        for metric_scope, cutoff in [("rank1", 1), ("top10", int(top_k)), ("full_pool_oracle", 999999)]:
            metric_rows = []
            scoped_pool = normalized[normalized["_rank"].le(cutoff)] if cutoff < 999999 else normalized
            for control in controls:
                annotated = annotate_pool(pool_for_control(scoped_pool, control), control)
                if annotated.empty:
                    metric_rows.append({"best_tanimoto": 0.0, "exact_match": False, "n_candidates": 0})
                    continue
                best = annotated.sort_values(["_best_tanimoto", "_rank"], ascending=[False, True]).head(1)
                metric_rows.append(
                    {
                        "best_tanimoto": safe_float(best.iloc[0].get("_best_tanimoto"), 0.0),
                        "exact_match": bool(best.iloc[0].get("_exact_match", False)),
                        "n_candidates": int(annotated.shape[0]),
                    }
                )
            values = [safe_float(row["best_tanimoto"], 0.0) for row in metric_rows]
            covered = [row for row in metric_rows if safe_int(row.get("n_candidates"), 0) > 0]
            rows.append(
                {
                    "track": "retrieval_only",
                    "lib_id": int(lib_id),
                    "ion_mode": mode,
                    "profile": "retrieval",
                    "base_profile": "retrieval",
                    "holdout_scope": "x_control",
                    "metric_scope": metric_scope,
                    "n_controls": len(metric_rows),
                    "n_controls_covered": len(covered),
                    "n_exact_match": int(sum(bool(row["exact_match"]) for row in metric_rows)),
                    "n_tanimoto_ge_0p7": int(sum(value >= GOOD_TANIMOTO for value in values)),
                    "mean_best_tanimoto": round(float(np.mean(values)) if values else 0.0, 4),
                    "median_best_tanimoto": round(float(np.median(values)) if values else 0.0, 4),
                }
            )
    return pd.DataFrame(rows)


def load_database_analog_metrics(analog_root: Path, libs: list[int]) -> pd.DataFrame:
    path = analog_root / "aggregate" / "x_domain_analog_v1_metrics.csv"
    if not path.exists():
        return pd.DataFrame()
    df = read_table(path)
    if df.empty or "lib_id" not in df.columns:
        return pd.DataFrame()
    df = df[df["lib_id"].astype("Int64").isin(libs)].copy()
    if df.empty:
        return df
    df["track"] = "database_analog_policy"
    return df


def aggregate_denovo_assembly(
    *,
    output_root: Path,
    run_root: Path,
    baseline_root: Path,
    analog_root: Path,
    ablation_root: Path,
    hard_holdout_dir: Path,
    known_x: Path,
    known_x_all: Path,
    libs: list[int],
    top_k: int,
) -> dict[str, Any]:
    metric_rows: list[dict[str, Any]] = []
    per_control_rows: list[dict[str, Any]] = []
    for lib_dir in sorted(run_root.glob("lib_*_*")):
        if not lib_dir.is_dir():
            continue
        parts = lib_dir.name.split("_")
        if len(parts) < 3:
            continue
        lib_id = int(parts[1])
        if lib_id not in libs:
            continue
        ion_mode = "_".join(parts[2:])
        for profile_dir in sorted((lib_dir / "selfies_profiles").glob("*")):
            if not profile_dir.is_dir():
                continue
            profile = profile_dir.name
            holdout_scope = "non_x_known" if profile.startswith("non_x_holdout_") else "x_control"
            controls = (
                load_non_x_holdouts(hard_holdout_dir, lib_id, ion_mode)
                if holdout_scope == "non_x_known"
                else load_x_controls(known_x=known_x, known_x_all=known_x_all, lib_id=lib_id, ion_mode=ion_mode)
            )
            metrics, per_rows = evaluate_profile(
                profile_dir=profile_dir,
                lib_id=lib_id,
                ion_mode=ion_mode,
                profile=profile,
                holdout_scope=holdout_scope,
                controls=controls,
                top_k=top_k,
            )
            for row in metrics:
                row["track"] = "novel_denovo_assembly"
            for row in per_rows:
                row["track"] = "novel_denovo_assembly"
            metric_rows.extend(metrics)
            per_control_rows.extend(per_rows)

    assembly_metrics = pd.DataFrame(metric_rows)
    if not assembly_metrics.empty:
        ordered_cols = ["track"] + [col for col in assembly_metrics.columns if col != "track"]
        assembly_metrics = assembly_metrics[ordered_cols]
    retrieval_metrics = evaluate_retrieval_track(
        baseline_root=baseline_root,
        known_x=known_x,
        known_x_all=known_x_all,
        libs=libs,
        top_k=top_k,
    )
    analog_metrics = load_database_analog_metrics(analog_root, libs)
    all_metrics = pd.concat(
        [df for df in [retrieval_metrics, analog_metrics, assembly_metrics] if not df.empty],
        ignore_index=True,
        sort=False,
    )
    baseline = load_baseline_metrics(ablation_root)
    assembly_deltas = build_deltas(assembly_metrics, baseline) if not assembly_metrics.empty else pd.DataFrame()
    if not assembly_deltas.empty:
        assembly_deltas.insert(0, "track", "novel_denovo_assembly")
    hard_gap = build_hard_x_gap(assembly_metrics) if not assembly_metrics.empty else pd.DataFrame()
    novelty_rows: list[pd.DataFrame] = []
    for path in sorted(run_root.glob("lib_*_*/selfies_profiles/*/novelty_audit.csv")):
        try:
            df = read_table(path)
        except Exception:
            continue
        if df.empty:
            continue
        lib_id = int(path.parts[-4].split("_")[1])
        profile = path.parts[-2]
        df["lib_id"] = lib_id
        df["profile"] = profile
        novelty_rows.append(df)
    novelty = pd.concat(novelty_rows, ignore_index=True) if novelty_rows else pd.DataFrame()
    if not novelty.empty:
        novelty_summary = (
            novelty.groupby(["lib_id", "profile"], dropna=False)
            .agg(
                n_candidates_considered=("candidate_smiles", "count"),
                n_exact_known_rejected=("exact_known_structure_rejected", "sum"),
                mean_nearest_db_tanimoto=("nearest_db_tanimoto", "mean"),
                median_nearest_db_tanimoto=("nearest_db_tanimoto", "median"),
            )
            .reset_index()
        )
    else:
        novelty_summary = pd.DataFrame()
    accepted = []
    if not assembly_deltas.empty:
        accepted = assembly_deltas[
            assembly_deltas["metric_scope"].eq("full_pool_oracle") & assembly_deltas["acceptance_pass"].astype(bool)
        ].to_dict(orient="records")
    high_hit_safe = [row for row in accepted if safe_int(row.get("delta_ge07"), 0) >= 0]
    if high_hit_safe:
        next_branch = "De novo assembly improved full-pool oracle without losing high-similarity hits; prioritize reranking of novel assembly candidates."
    elif accepted:
        next_branch = "De novo assembly produced a mean-only oracle gain but lost high-similarity hits; inspect per-control failures and improve assembly operators before reranking."
    else:
        next_branch = "Novel assembly did not clear X-control oracle gates; improve assembly operators and feature constraints before reranking."
    recommendations = {
        "accepted_candidate_pool_rows": accepted,
        "high_hit_safe_candidate_pool_rows": high_hit_safe,
        "next_branch": next_branch,
    }
    aggregate_dir = output_root / "aggregate"
    aggregate_dir.mkdir(parents=True, exist_ok=True)
    write_csv(aggregate_dir / "denovo_assembly_v1_metrics.csv", all_metrics)
    write_csv(aggregate_dir / "denovo_assembly_v1_per_control.csv", pd.DataFrame(per_control_rows))
    write_csv(aggregate_dir / "denovo_assembly_v1_deltas.csv", assembly_deltas)
    write_csv(aggregate_dir / "hard_holdout_vs_x_gap.csv", hard_gap)
    write_csv(aggregate_dir / "novelty_vs_quality.csv", novelty_summary)
    write_json(aggregate_dir / "next_branch_recommendations.json", recommendations)
    report = render_report(all_metrics, assembly_deltas, novelty_summary, hard_gap, recommendations)
    (aggregate_dir / "denovo_assembly_v1_report.md").write_text(report)
    payload = {
        "metrics_csv": str(aggregate_dir / "denovo_assembly_v1_metrics.csv"),
        "deltas_csv": str(aggregate_dir / "denovo_assembly_v1_deltas.csv"),
        "novelty_vs_quality_csv": str(aggregate_dir / "novelty_vs_quality.csv"),
        "recommendations_json": str(aggregate_dir / "next_branch_recommendations.json"),
        "report_md": str(aggregate_dir / "denovo_assembly_v1_report.md"),
        "n_metric_rows": int(all_metrics.shape[0]),
        "n_per_control_rows": int(len(per_control_rows)),
        "recommendations": recommendations,
    }
    write_json(aggregate_dir / "denovo_assembly_v1_aggregate_summary.json", payload)
    return payload


def render_report(
    metrics: pd.DataFrame,
    deltas: pd.DataFrame,
    novelty: pd.DataFrame,
    gap: pd.DataFrame,
    recommendations: dict[str, Any],
) -> str:
    metric_cols = [
        col
        for col in [
            "track",
            "lib_id",
            "profile",
            "holdout_scope",
            "metric_scope",
            "n_controls",
            "n_controls_covered",
            "mean_best_tanimoto",
            "n_tanimoto_ge_0p7",
        ]
        if col in metrics.columns
    ]
    lines = [
        "# Retrieval-Anchored De Novo Assembly V1",
        "",
        f"Recommendation: {recommendations.get('next_branch', '')}",
        "",
        "## Track Metrics",
        dataframe_to_markdown(metrics[metric_cols] if metric_cols else metrics),
        "",
        "## Novel Assembly Deltas Versus Baseline SELFIES",
        dataframe_to_markdown(deltas),
        "",
        "## Novelty Versus Quality",
        dataframe_to_markdown(novelty),
        "",
        "## Hard Holdout Versus X Gap",
        dataframe_to_markdown(gap),
        "",
    ]
    return "\n".join(lines)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--baseline-root", type=Path, default=DEFAULT_BASELINE_ROOT)
    parser.add_argument("--analog-root", type=Path, default=DEFAULT_ANALOG_ROOT)
    parser.add_argument("--ablation-root", type=Path, default=DEFAULT_ABLATION_ROOT)
    parser.add_argument("--hard-holdout-dir", type=Path, default=DEFAULT_HARD_HOLDOUT_DIR)
    parser.add_argument("--known-x", type=Path, default=Path("known_X.csv"))
    parser.add_argument("--known-x-all", type=Path, default=Path("known_X_all_available_structures_simple.csv"))
    parser.add_argument(
        "--unknown-feature-table", type=Path, default=Path("merged_encoder/data/ions_unknown_merged.csv")
    )
    parser.add_argument("--libs", type=str, default="305,209")
    parser.add_argument("--profiles", type=str, default="strict_1da,balanced")
    parser.add_argument("--scopes", type=str, default="x_control,hard_holdout")
    parser.add_argument("--max-features-per-scope", type=int, default=0)
    parser.add_argument("--anchor-top-n", type=int, default=5)
    parser.add_argument("--max-candidates-per-feature", type=int, default=250)
    parser.add_argument("--top-k", type=int, default=10)
    parser.add_argument(
        "--all-x-features",
        action="store_true",
        help="Generate for every X feature. Default is structured X controls only for benchmark runs.",
    )
    return parser.parse_args()


def resolve(path: Path, repo_root: Path) -> Path:
    path = path.expanduser()
    return path if path.is_absolute() else repo_root / path


def parse_csv(text: str) -> list[str]:
    return [part.strip() for part in str(text).split(",") if part.strip()]


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
    libs = parse_libs(args.libs)
    profiles = parse_csv(args.profiles)
    scopes = parse_csv(args.scopes)
    generation = run_denovo_assembly(
        repo_root=repo_root,
        output_root=output_root,
        baseline_root=baseline_root,
        hard_holdout_dir=hard_holdout_dir,
        known_x=known_x,
        known_x_all=known_x_all,
        unknown_feature_table=unknown_feature_table,
        libs=libs,
        profiles=profiles,
        scopes=scopes,
        max_features_per_scope=int(args.max_features_per_scope),
        anchor_top_n=int(args.anchor_top_n),
        max_candidates_per_feature=int(args.max_candidates_per_feature),
        top_k=int(args.top_k),
        all_x_features=bool(args.all_x_features),
    )
    aggregate = aggregate_denovo_assembly(
        output_root=output_root,
        run_root=Path(generation["run_root"]),
        baseline_root=baseline_root,
        analog_root=analog_root,
        ablation_root=ablation_root,
        hard_holdout_dir=hard_holdout_dir,
        known_x=known_x,
        known_x_all=known_x_all,
        libs=libs,
        top_k=int(args.top_k),
    )
    print(
        json.dumps(
            {"run_root": generation["run_root"], "outputs": len(generation["outputs"]), "aggregate": aggregate},
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
