#!/usr/bin/env python3
"""Utilities for peptide-first de novo generation.

This ports the useful, dependency-light peptide pieces from the old NCE
``aminoAcids.py`` flow into the current spec2mol workflow.
"""

from __future__ import annotations

import itertools
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from denovo_improvement_utils import canonicalize_smiles, exact_mass_from_smiles, formula_counts, formula_from_smiles, inchikey_from_smiles
from nce_evidence_utils import ion_and_loss_bins, neutral_mass_from_precursor
from peptide_composition_v1_utils import normalize_peptide_composition, parse_residue_tokens
from peptide_permutation_v1_utils import normalize_ordered_peptide_sequence
from peptide_spectrum_supervision_v1_utils import calibrated_spectrum_score, load_calibration
from x_domain_analog_utils import safe_float, safe_int, xid_from_feature_key


DEFAULT_OUTPUT_ROOT = Path("models/peptide_first_denovo_v1")
POLICY = "protected_champion_peptide_first_guarded"
GOOD_TANIMOTO = 0.7
PROTECTED_209_RECOVERED = {
    "X-15853",
    "X-17855",
    "X-17861",
    "X-17919",
    "X-18345",
    "X-18887",
    "X-21410",
    "X-21821",
    "X-24137",
    "X-24972",
    "X-25459",
    "X-25523",
    "X-25737",
    "X-25744",
}

PROTON_MASS = 1.007276466812
WATER_FORMULA = {"H": 2, "O": 1}
WATER_MASS = 18.010564684


@dataclass(frozen=True)
class ResidueSpec:
    token: str
    formula: str
    mass: float
    sidechain_smiles: str
    first_only: bool = False


RESIDUES: dict[str, ResidueSpec] = {
    "Ala": ResidueSpec("Ala", "C3H5NO", 71.037114, "(C)"),
    "Arg": ResidueSpec("Arg", "C6H12N4O", 156.101111, "(CCCNC(=N)N)"),
    "Asn": ResidueSpec("Asn", "C4H6N2O2", 114.042927, "(CC(N)=O)"),
    "Asp": ResidueSpec("Asp", "C4H5NO3", 115.026943, "(CC(=O)O)"),
    "Cys": ResidueSpec("Cys", "C3H5NOS", 103.009185, "(CS)"),
    "Gln": ResidueSpec("Gln", "C5H8N2O2", 128.058578, "(CCC(N)=O)"),
    "Glu": ResidueSpec("Glu", "C5H7NO3", 129.042593, "(CCC(=O)O)"),
    "Gly": ResidueSpec("Gly", "C2H3NO", 57.021464, ""),
    "His": ResidueSpec("His", "C6H7N3O", 137.058912, "(Cc1c[nH]cn1)"),
    "Lys": ResidueSpec("Lys", "C6H12N2O", 128.094963, "(CCCCN)"),
    "Met": ResidueSpec("Met", "C5H9NOS", 131.040485, "(CCSC)"),
    "Phe": ResidueSpec("Phe", "C9H9NO", 147.068414, "(Cc1ccccc1)"),
    "Pro": ResidueSpec("Pro", "C5H7NO", 97.052764, "pro"),
    "Pyr": ResidueSpec("Pyr", "C5H5NO2", 111.032029, "pyr", first_only=True),
    "Ser": ResidueSpec("Ser", "C3H5NO2", 87.032028, "(CO)"),
    "Thr": ResidueSpec("Thr", "C4H7NO2", 101.047679, "(C(O)C)"),
    "Trp": ResidueSpec("Trp", "C11H10N2O", 186.079313, "(Cc1c[nH]c2ccccc12)"),
    "Tyr": ResidueSpec("Tyr", "C9H9NO2", 163.063320, "(Cc1ccc(O)cc1)"),
    "Val": ResidueSpec("Val", "C5H9NO", 99.068414, "(C(C)C)"),
    # Xle is used for composition scoring; generation expands it to Leu/Ile variants.
    "Xle": ResidueSpec("Xle", "C6H11NO", 113.084064, "(CC(C)C)"),
}

XLE_VARIANTS = {"Xle": ("Leu", "Ile")}
COLLAPSE_TO_RESIDUE = {"Leu": "Xle", "Ile": "Xle"}
AA_SMILES = {
    "Ala": "NC(C)C(=O)O",
    "Arg": "NC(CCCNC(=N)N)C(=O)O",
    "Asn": "NC(CC(N)=O)C(=O)O",
    "Asp": "NC(CC(=O)O)C(=O)O",
    "Cys": "NC(CS)C(=O)O",
    "Gln": "NC(CCC(N)=O)C(=O)O",
    "Glu": "NC(CCC(=O)O)C(=O)O",
    "Gly": "NCC(=O)O",
    "His": "NC(Cc1c[nH]cn1)C(=O)O",
    "Ile": "NC(C(C)CC)C(=O)O",
    "Leu": "NC(CC(C)C)C(=O)O",
    "Lys": "NC(CCCCN)C(=O)O",
    "Met": "NC(CCSC)C(=O)O",
    "Phe": "NC(Cc1ccccc1)C(=O)O",
    "Pro": "N1CCCC1C(=O)O",
    "Pyr": "N1C(=O)CCC1C(=O)O",
    "Ser": "NC(CO)C(=O)O",
    "Thr": "NC(C(O)C)C(=O)O",
    "Trp": "NC(Cc1c[nH]c2ccccc12)C(=O)O",
    "Tyr": "NC(Cc1ccc(O)cc1)C(=O)O",
    "Val": "NC(C(C)C)C(=O)O",
}


@dataclass(frozen=True)
class PeptideModification:
    name: str
    delta_formula: str
    delta_mass: float
    supports_smiles: bool
    aliases: tuple[str, ...] = ()


MODIFICATIONS: dict[str, PeptideModification] = {
    "none": PeptideModification("none", "", 0.0, True, ("peptide", "tripeptide", "dipeptide", "tetrapeptide")),
    "n_acetyl": PeptideModification("n_acetyl", "C2H2O", 42.010565, True, ("n-acetyl", "acetyl")),
    "glycoside": PeptideModification("glycoside", "C6H10O5", 162.052823, True, ("glycosylated", "glycoside", "hexose", "sugar")),
    "glucuronide": PeptideModification("glucuronide", "C6H8O6", 176.032088, False, ("glucuronide", "glucuronidated")),
    "gamma_glutamyl": PeptideModification("gamma_glutamyl", "C5H7NO3", 129.042593, False, ("gamma-glu", "gamma glutamyl", "g-glu")),
}


def formula_string(counts: dict[str, int]) -> str:
    parts: list[str] = []
    for element in ["C", "H", "N", "O", "P", "S", "F", "Cl", "Br", "I"]:
        count = int(counts.get(element, 0))
        if count > 0:
            parts.append(element if count == 1 else f"{element}{count}")
    return "".join(parts)


def add_formula_counts(*formulas: Any) -> dict[str, int]:
    out: dict[str, int] = {}
    for formula in formulas:
        for element, count in formula_counts(formula).items():
            out[element] = out.get(element, 0) + int(count)
    return {k: v for k, v in out.items() if v}


def peptide_formula(tokens: list[str], modification: str = "none") -> str:
    counts = dict(WATER_FORMULA)
    for token in tokens:
        spec = RESIDUES[residue_token(token)]
        for element, count in formula_counts(spec.formula).items():
            counts[element] = counts.get(element, 0) + int(count)
    mod = MODIFICATIONS.get(modification, MODIFICATIONS["none"])
    for element, count in formula_counts(mod.delta_formula).items():
        counts[element] = counts.get(element, 0) + int(count)
    return formula_string(counts)


def peptide_mass(tokens: list[str], modification: str = "none") -> float:
    return float(WATER_MASS + sum(RESIDUES[residue_token(token)].mass for token in tokens) + MODIFICATIONS.get(modification, MODIFICATIONS["none"]).delta_mass)


def normalized_composition(tokens: list[str]) -> str:
    return normalize_peptide_composition("-".join(tokens))


def composition_tokens(composition: Any) -> list[str]:
    return [token for token in normalize_peptide_composition(composition).split("-") if token]


def residue_token(token: Any) -> str:
    text = str(token or "").strip()
    return COLLAPSE_TO_RESIDUE.get(text, text)


def enumerate_peptide_compositions(
    neutral_mass: float,
    *,
    lengths: tuple[int, ...] = (2, 3, 4),
    modifications: tuple[str, ...] = ("none", "n_acetyl", "glycoside", "glucuronide", "gamma_glutamyl"),
    tolerance_da: float = 0.02,
) -> pd.DataFrame:
    residues = [token for token in RESIDUES if token != "Pyr"] + ["Pyr"]
    rows: list[dict[str, Any]] = []
    for length in lengths:
        for combo in itertools.combinations_with_replacement(residues, length):
            if combo.count("Pyr") > 1:
                continue
            tokens = sorted(combo, key=lambda token: token)
            comp = normalized_composition(tokens)
            if not comp:
                continue
            base_tokens = composition_tokens(comp)
            for mod_name in modifications:
                mass = peptide_mass(base_tokens, mod_name)
                error = mass - float(neutral_mass)
                if abs(error) <= float(tolerance_da):
                    rows.append(
                        {
                            "normalized_peptide_composition": comp,
                            "peptide_length": int(length),
                            "modification": mod_name,
                            "peptide_formula": peptide_formula(base_tokens, mod_name),
                            "theoretical_neutral_mass": round(float(mass), 6),
                            "mass_error_da": round(float(error), 6),
                            "mass_score": round(float(max(0.0, 1.0 - abs(error) / max(float(tolerance_da), 1e-9))), 6),
                            "smiles_supported": bool(MODIFICATIONS[mod_name].supports_smiles),
                        }
                    )
    if not rows:
        return pd.DataFrame()
    return pd.DataFrame(rows).drop_duplicates(["normalized_peptide_composition", "modification", "theoretical_neutral_mass"]).reset_index(drop=True)


def build_peptide_composition_index(
    *,
    lengths: tuple[int, ...] = (2, 3, 4),
    modifications: tuple[str, ...] = ("none", "n_acetyl", "glycoside", "glucuronide", "gamma_glutamyl"),
) -> pd.DataFrame:
    residues = [token for token in RESIDUES if token != "Pyr"] + ["Pyr"]
    rows: list[dict[str, Any]] = []
    for length in lengths:
        for combo in itertools.combinations_with_replacement(residues, int(length)):
            if combo.count("Pyr") > 1:
                continue
            comp = normalized_composition(sorted(combo, key=lambda token: token))
            if not comp:
                continue
            tokens = composition_tokens(comp)
            for mod_name in modifications:
                rows.append(
                    {
                        "normalized_peptide_composition": comp,
                        "peptide_length": int(length),
                        "modification": mod_name,
                        "peptide_formula": peptide_formula(tokens, mod_name),
                        "theoretical_neutral_mass": round(float(peptide_mass(tokens, mod_name)), 6),
                        "mass_score": 0.0,
                        "mass_error_da": 0.0,
                        "smiles_supported": bool(MODIFICATIONS[mod_name].supports_smiles),
                    }
                )
    if not rows:
        return pd.DataFrame()
    return (
        pd.DataFrame(rows)
        .drop_duplicates(["normalized_peptide_composition", "modification", "theoretical_neutral_mass"])
        .sort_values(["theoretical_neutral_mass", "peptide_length", "normalized_peptide_composition"], kind="mergesort")
        .reset_index(drop=True)
    )


def match_peptide_composition_index(index: pd.DataFrame, neutral_mass: float, tolerance_da: float) -> pd.DataFrame:
    mass = safe_float(neutral_mass, math.nan)
    if index.empty or not np.isfinite(mass):
        return pd.DataFrame(columns=list(index.columns) + ["mass_error_da", "mass_score"])
    tol = float(tolerance_da)
    scoped = index[pd.to_numeric(index["theoretical_neutral_mass"], errors="coerce").between(mass - tol, mass + tol, inclusive="both")].copy()
    if scoped.empty:
        return scoped
    scoped["mass_error_da"] = pd.to_numeric(scoped["theoretical_neutral_mass"], errors="coerce") - mass
    scoped["mass_score"] = (1.0 - scoped["mass_error_da"].abs() / max(tol, 1e-9)).clip(lower=0.0, upper=1.0)
    return scoped.sort_values(["mass_score", "peptide_length", "normalized_peptide_composition"], ascending=[False, True, True], kind="mergesort").reset_index(drop=True)


def ordered_sequences_for_composition(composition: Any, max_permutations: int = 24) -> list[str]:
    tokens = composition_tokens(composition)
    if not tokens:
        return []
    seen: set[str] = set()
    out: list[str] = []
    for perm in itertools.permutations(tokens):
        if "Pyr" in perm and perm[0] != "Pyr":
            continue
        seq = "-".join(perm)
        if seq in seen:
            continue
        seen.add(seq)
        out.append(seq)
    return sorted(out)[: int(max_permutations)]


def xle_variant_sequences(sequence: Any, max_variants: int = 8) -> list[tuple[str, str]]:
    tokens = [token for token in normalize_ordered_peptide_sequence(sequence).split("-") if token]
    if not tokens:
        return []
    options = [XLE_VARIANTS.get(token, (token,)) for token in tokens]
    rows: list[tuple[str, str]] = []
    for variant in itertools.product(*options):
        actual = "-".join(variant)
        collapsed = "-".join(["Xle" if token in {"Leu", "Ile"} else token for token in variant])
        rows.append((collapsed, actual))
    return rows[: int(max_variants)]


def build_oligo_smiles(tokens: list[str]) -> str:
    try:
        return "".join(AA_SMILES[token][:-1] for token in tokens) + "O"
    except KeyError:
        return ""


def apply_modification_smiles(base_smiles: str, modification: str) -> str:
    if not base_smiles:
        return ""
    if modification == "none":
        return base_smiles
    if modification == "n_acetyl" and base_smiles.startswith("N") and not base_smiles.startswith("N1"):
        return "CC(=O)N" + base_smiles[1:]
    if modification == "glycoside" and base_smiles.startswith("N") and not base_smiles.startswith("N1"):
        # Simple N-glycosylated proxy; used as a review candidate, not a stereochemical claim.
        return "OC1OC(CO)C(O)C(O)C1N" + base_smiles[1:]
    return ""


def peptide_smiles_from_sequence(sequence: Any, modification: str = "none") -> str:
    tokens = [token for token in str(sequence or "").split("-") if token]
    base = build_oligo_smiles(tokens)
    modified = apply_modification_smiles(base, modification)
    return canonicalize_smiles(modified) or ""


def parse_legacy_annotations(path: Path, lib_id: int = 209) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(columns=["comp_id", "lib_id", "parent_mass", "details", "unknown_name", "legacy_composition", "legacy_modification"])
    df = pd.read_csv(path, low_memory=False, encoding="latin1")
    rows: list[dict[str, Any]] = []
    for _, row in df.iterrows():
        if safe_int(row.get("lib_id"), -1) != int(lib_id):
            continue
        details = str(row.get("details", "") or "")
        text = details.lower()
        comp = normalize_peptide_composition(details)
        is_peptide = bool(comp and len(comp.split("-")) >= 2) or "peptide" in text or "leu/ile" in text
        if not is_peptide:
            continue
        mods = [name for name, spec in MODIFICATIONS.items() if name != "none" and any(alias in text for alias in spec.aliases)]
        if not mods:
            mods = ["none"]
        rows.append(
            {
                "comp_id": safe_int(row.get("comp_id"), -1),
                "lib_id": safe_int(row.get("lib_id"), int(lib_id)),
                "parent_mass": safe_float(row.get("parent_mass"), math.nan),
                "details": details,
                "unknown_name": row.get("unknown_name", ""),
                "legacy_composition": comp,
                "legacy_modification": "|".join(sorted(set(mods))),
                "legacy_track": "annotation_assisted",
            }
        )
    return pd.DataFrame(rows).drop_duplicates().reset_index(drop=True) if rows else pd.DataFrame()


def feature_rows(feature_table: pd.DataFrame, lib_id: int = 209) -> pd.DataFrame:
    df = feature_table[pd.to_numeric(feature_table.get("lib_id"), errors="coerce").fillna(-1).astype(int).eq(int(lib_id))].copy()
    if df.empty:
        return df
    rows: list[dict[str, Any]] = []
    for feature_key, group in df.groupby("feature_key", dropna=False):
        first = group.iloc[0]
        parent_mz = safe_float(first.get("parent_mz"), safe_float(first.get("mass"), math.nan))
        source_neutral_mass = safe_float(first.get("source_neutral_mass"), math.nan)
        neutral_mass = source_neutral_mass if np.isfinite(source_neutral_mass) else neutral_mass_from_precursor(parent_mz, int(lib_id))
        ions: set[int] = set()
        losses: set[int] = set()
        for _, peak in group.iterrows():
            mz = safe_float(peak.get("mass"), math.nan)
            if not np.isfinite(mz) or mz <= 0:
                continue
            ion_bin, loss_bin = ion_and_loss_bins(mz, parent_mz, int(lib_id))
            if ion_bin:
                ions.add(int(ion_bin))
            if loss_bin:
                losses.add(int(loss_bin))
        rows.append(
            {
                "feature_key": str(feature_key),
                "comp_id": safe_int(first.get("comp_id"), -1),
                "control_id": str(first.get("compound_name", "") or xid_from_feature_key(feature_key) or ""),
                "lib_id": int(lib_id),
                "ion_mode": "negative",
                "parent_mz": parent_mz,
                "neutral_mass_fallback": neutral_mass,
                "source_parent_adduct": str(first.get("source_parent_adduct", "") or ""),
                "source_neutral_mass": source_neutral_mass if np.isfinite(source_neutral_mass) else "",
                "ri": safe_float(first.get("ri"), math.nan),
                "rt": safe_float(first.get("rt"), math.nan),
                "n_ms2_peaks": int(group.shape[0]),
                "observed_ion_bins_json": json.dumps(sorted(ions)),
                "observed_loss_bins_json": json.dumps(sorted(losses)),
            }
        )
    return pd.DataFrame(rows)


def load_adduct_hypotheses(root: Path, features: pd.DataFrame, lib_id: int = 209, top_n: int = 3) -> pd.DataFrame:
    path = root / "neutral_mass_hypotheses.csv"
    if path.exists():
        raw = pd.read_csv(path, low_memory=False)
        raw = raw[pd.to_numeric(raw.get("lib_id"), errors="coerce").fillna(-1).astype(int).eq(int(lib_id))].copy()
        if not raw.empty:
            if "adduct_rank" not in raw.columns:
                raw["adduct_rank"] = 1
            if "adduct_score" not in raw.columns:
                raw["adduct_score"] = 0.5
            raw = raw.sort_values(["feature_key", "adduct_rank", "adduct_score"], ascending=[True, True, False], kind="mergesort")
            raw = raw.drop_duplicates(["feature_key", "adduct_label", "neutral_mass"]).groupby("feature_key", dropna=False).head(int(top_n))
            return raw[["feature_key", "adduct_label", "neutral_mass", "adduct_score", "adduct_rank"]].copy()
    return pd.DataFrame(
        [
            {
                "feature_key": row.feature_key,
                "adduct_label": "[M-H]-",
                "neutral_mass": row.neutral_mass_fallback,
                "adduct_score": 0.5,
                "adduct_rank": 1,
            }
            for row in features.itertuples()
        ]
    )


def ri_plausibility(ri: Any, peptide_length: int, modification: str) -> float:
    value = safe_float(ri, math.nan)
    if not np.isfinite(value):
        return 0.5
    center = {2: 1700.0, 3: 2600.0, 4: 3600.0}.get(int(peptide_length), 2600.0)
    if modification in {"glycoside", "glucuronide"}:
        center -= 350.0
    sigma = 1650.0
    return float(np.exp(-((value - center) ** 2) / (2 * sigma * sigma)))


def spectrum_score_for_sequence(sequence: str, ion_bins: set[int], loss_bins: set[int], calibration: dict[str, Any] | None) -> float:
    if not calibration:
        return 0.0
    score, _ = calibrated_spectrum_score(sequence, ion_bins, loss_bins, calibration)
    return float(score)
