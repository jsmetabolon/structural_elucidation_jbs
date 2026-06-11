#!/usr/bin/env python3
"""Shared helpers for biotransformation-guided 209 de novo generation."""

from __future__ import annotations

import itertools
import json
import math
import re
from dataclasses import dataclass
from typing import Any, Iterable

import numpy as np
import pandas as pd

from denovo_improvement_utils import (
    HAS_RDKIT,
    canonicalize_smiles,
    exact_mass_from_smiles,
    formula_counts,
    formula_from_smiles,
    inchikey_from_smiles,
)
from nce_evidence_utils import formula_mass, mass_bin, neutral_mass_from_precursor
from nce_generation_v3_utils import mass_from_bin, operator_specific_modifications
from x_domain_analog_utils import safe_float, safe_int, xid_from_feature_key

if HAS_RDKIT:
    from rdkit import Chem
else:  # pragma: no cover
    Chem = None


GOOD_TANIMOTO = 0.7
DEFAULT_LIB_ID = 209
DEFAULT_ION_MODE = "negative"
TRACK_CONTROL_BLIND = "control_blind"
TRACK_ANNOTATION_ASSISTED = "annotation_assisted"


@dataclass(frozen=True)
class TransformSpec:
    name: str
    delta_mass: float
    delta_formula: str
    aliases: tuple[str, ...]
    diagnostic_losses: tuple[float, ...]
    diagnostic_ions: tuple[float, ...] = ()
    generation_operators: tuple[str, ...] = ()


TRANSFORMS: dict[str, TransformSpec] = {
    "sulfate": TransformSpec(
        "sulfate",
        79.956815,
        "SO3",
        ("sulfate", "sulfated", "sulfo", "sulfonate"),
        (79.9568, 96.9601),
        (79.9568, 96.9601),
        ("sulfate",),
    ),
    "disulfate": TransformSpec(
        "disulfate",
        159.913630,
        "S2O6",
        ("disulfate", "di-sulfate"),
        (79.9568, 96.9601, 159.9136),
        (79.9568, 96.9601),
        ("sulfate",),
    ),
    "glucuronide": TransformSpec(
        "glucuronide",
        176.032088,
        "C6H8O6",
        ("glucuronide", "glucuronidated", "glucuronic", "uronide"),
        (176.0321,),
        (),
        ("glucuronide",),
    ),
    "glycoside": TransformSpec(
        "glycoside",
        162.052823,
        "C6H10O5",
        ("glycoside", "glycosyl", "glucoside", "hexose", "ribose", "sugar"),
        (162.0528, 132.0423, 120.0423),
        (),
        ("glycoside",),
    ),
    "phosphate": TransformSpec(
        "phosphate",
        79.966331,
        "HPO3",
        ("phosphate", "phospho", "phosphoryl"),
        (79.9663, 97.9769),
        (79.9663, 97.9769),
        ("phosphate",),
    ),
    "acetyl": TransformSpec(
        "acetyl",
        42.010565,
        "C2H2O",
        ("acetyl", "n-acetyl", "acetylated"),
        (42.0106,),
        (),
        ("acetyl",),
    ),
    "glycine": TransformSpec(
        "glycine",
        57.021464,
        "C2H3NO",
        ("glycine", "glycyl", "acylglycine", "gly-"),
        (57.0215, 75.0320),
        (),
        ("glycine",),
    ),
    "taurine": TransformSpec(
        "taurine",
        107.004099,
        "C2H5NO2S",
        ("taurine", "tauryl"),
        (107.0041, 124.0068),
        (),
        ("taurine",),
    ),
    "carnitine": TransformSpec(
        "carnitine",
        145.110279,
        "C7H15NO2",
        ("carnitine", "acylcarnitine"),
        (59.0735, 85.0284, 145.1103),
        (),
        ("carnitine",),
    ),
    "glutamine_like": TransformSpec(
        "glutamine_like",
        128.058578,
        "C5H8N2O2",
        ("glutamine", "glutathione", "gamma-glu", "gln"),
        (128.0586, 146.0691),
        (),
        ("glutamine",),
    ),
}


COMMON_NEGATIVE_ADDUCTS: dict[str, float] = {
    "[M-H]-": -1.007276466812,
    "[M+Cl]-": 34.969402,
    "[M+HCOO]-": 44.998201,
    "[M+CH3COO]-": 59.013851,
    "[M+Br]-": 78.918885,
    "[M+CF3COO]-": 112.985586,
}

COMMON_POSITIVE_ADDUCTS: dict[str, float] = {
    "[M+H]+": 1.007276466812,
    "[M+NH4]+": 18.03382,
    "[M+Na]+": 22.989218,
    "[M+CH3OH+H]+": 33.033489,
    "[M+K]+": 38.9632,
    "[M+CH3CN+H]+": 42.003823,
}

NEGATIVE_BACKGROUND_IONS: dict[str, tuple[float, ...]] = {
    "cyanide_acetonitrile": (26.0,),
    "formate": (45.0,),
    "acetate": (59.0,),
    "tfa_related": (69.0, 113.0, 127.0, 227.0, 249.0),
    "phosphate": (79.0, 97.0),
    "sulfate": (80.0, 97.0),
    "methanesulfonate": (95.0,),
}

POSITIVE_BACKGROUND_IONS: dict[str, tuple[float, ...]] = {
    "methanol": (33.0, 65.0),
    "acetonitrile": (42.0, 64.0, 83.0, 120.0, 137.0),
    "dmso": (79.0, 102.0, 157.0, 179.0),
    "phthalate_plasticizer": (149.0, 163.0, 167.0, 279.0, 391.0),
    "siloxane": (355.0, 371.0, 445.0, 519.0, 593.0),
}


AA_RESIDUES: dict[str, tuple[str, float, str]] = {
    "Ala": ("C3H5NO", 71.037114, "C"),
    "Arg": ("C6H12N4O", 156.101111, "CCCNC(=N)N"),
    "Asn": ("C4H6N2O2", 114.042927, "CC(N)=O"),
    "Asp": ("C4H5NO3", 115.026943, "CC(=O)O"),
    "Cys": ("C3H5NOS", 103.009185, "CS"),
    "Gln": ("C5H8N2O2", 128.058578, "CCC(N)=O"),
    "Glu": ("C5H7NO3", 129.042593, "CCC(=O)O"),
    "Gly": ("C2H3NO", 57.021464, ""),
    "His": ("C6H7N3O", 137.058912, "Cc1c[nH]cn1"),
    "Ile": ("C6H11NO", 113.084064, "C(C)CC"),
    "Leu": ("C6H11NO", 113.084064, "CC(C)C"),
    "Lys": ("C6H12N2O", 128.094963, "CCCCN"),
    "Met": ("C5H9NOS", 131.040485, "CCSC"),
    "Phe": ("C9H9NO", 147.068414, "Cc1ccccc1"),
    "Pro": ("C5H7NO", 97.052764, "pro"),
    "Ser": ("C3H5NO2", 87.032028, "CO"),
    "Thr": ("C4H7NO2", 101.047679, "C(O)C"),
    "Trp": ("C11H10N2O", 186.079313, "Cc1c[nH]c2ccccc12"),
    "Tyr": ("C9H9NO2", 163.063320, "Cc1ccc(O)cc1"),
    "Val": ("C5H9NO", 99.068414, "C(C)C"),
    "Pyr": ("C5H5NO2", 111.032029, "pyr"),
}

WATER_MASS = 18.010565


def parse_xid(value: Any) -> str:
    match = re.search(r"(X-\d+)", str(value or ""))
    return match.group(1) if match else ""


def transform_by_alias(text: Any) -> list[str]:
    lower = str(text or "").lower()
    names: list[str] = []
    for name, spec in TRANSFORMS.items():
        if any(alias in lower for alias in spec.aliases):
            names.append(name)
    if re.search(r"\b(di|tri|tetra)?peptide\b|pyroglu|pyr-|leu/ile|tripeptide|tetrapeptide", lower):
        names.append("peptide")
    return sorted(set(names))


def formula_subtract(formula: Any, delta_formula: Any) -> str:
    counts = formula_counts(formula)
    delta = formula_counts(delta_formula)
    if not counts or not delta:
        return ""
    out: dict[str, int] = {}
    for element, count in counts.items():
        value = int(count) - int(delta.get(element, 0))
        if value < 0:
            return ""
        if value:
            out[element] = value
    return formula_string(out)


def formula_string(counts: dict[str, int]) -> str:
    parts: list[str] = []
    for element in ["C", "H", "N", "O", "P", "S", "F", "Cl", "Br", "I"]:
        count = int(counts.get(element, 0))
        if count <= 0:
            continue
        parts.append(element if count == 1 else f"{element}{count}")
    return "".join(parts)


def parse_json_int_list(value: Any) -> list[int]:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return []
    if isinstance(value, list):
        raw = value
    else:
        try:
            raw = json.loads(str(value or "[]"))
        except Exception:
            raw = [part for part in re.split(r"[;,\s]+", str(value or "")) if part]
    out: list[int] = []
    for item in raw if isinstance(raw, list) else []:
        number = safe_int(item, 0)
        if number:
            out.append(number)
    return out


def mass_evidence_score(target_masses: Iterable[float], ion_bins: Iterable[int], loss_bins: Iterable[int], tolerance_da: float = 0.4) -> float:
    observed = [mass_from_bin(int(bin_id)) for bin_id in list(ion_bins) + list(loss_bins) if int(bin_id) != 0]
    if not observed:
        return 0.0
    best = 0.0
    for target in target_masses:
        for mass in observed:
            err = abs(float(mass) - float(target))
            best = max(best, 1.0 / (1.0 + err / max(float(tolerance_da), 1e-6)))
    return float(np.clip(best, 0.0, 1.0))


def background_ion_matches(ion_bins: Iterable[int], lib_id: int = DEFAULT_LIB_ID, tolerance_da: float = 0.5) -> dict[str, Any]:
    observed = [mass_from_bin(int(bin_id)) for bin_id in ion_bins if int(bin_id) != 0]
    table = NEGATIVE_BACKGROUND_IONS if int(lib_id) in {209, 305} else POSITIVE_BACKGROUND_IONS
    labels: list[str] = []
    masses: list[float] = []
    for label, targets in table.items():
        for target in targets:
            if any(abs(float(mass) - float(target)) <= float(tolerance_da) for mass in observed):
                labels.append(label)
                masses.append(float(target))
                break
    return {
        "background_ion_labels_json": json.dumps(sorted(set(labels))),
        "background_ion_masses_json": json.dumps(sorted(set(round(mass, 4) for mass in masses))),
        "background_ion_match_count": int(len(set(labels))),
    }


def transform_evidence_score(spec: TransformSpec, row: pd.Series) -> float:
    ion_bins = parse_json_int_list(row.get("observed_ion_bins_json", "[]")) + parse_json_int_list(row.get("ion_bins_json", "[]"))
    loss_bins = parse_json_int_list(row.get("observed_loss_bins_json", "[]")) + parse_json_int_list(row.get("loss_bins_json", "[]"))
    return max(
        mass_evidence_score(spec.diagnostic_ions, ion_bins, (), tolerance_da=0.4),
        mass_evidence_score(spec.diagnostic_losses, (), loss_bins, tolerance_da=0.4),
    )


def evidence_track_allowed(source: str, track: str) -> bool:
    if str(track) == TRACK_ANNOTATION_ASSISTED:
        return True
    return "annotation" not in str(source).lower() and "unknown_list" not in str(source).lower() and "known_x" not in str(source).lower()


def neutral_mass_for_feature(parent_mz: Any, lib_id: int = DEFAULT_LIB_ID) -> float:
    return neutral_mass_from_precursor(parent_mz, int(lib_id))


def neutral_mass_hypotheses(parent_mz: Any, lib_id: int = DEFAULT_LIB_ID, adduct_hint: Any = "") -> list[dict[str, Any]]:
    mz = safe_float(parent_mz, math.nan)
    if not np.isfinite(mz) or mz <= 0:
        return []
    adducts = COMMON_NEGATIVE_ADDUCTS if int(lib_id) in {209, 305} else COMMON_POSITIVE_ADDUCTS
    hint = str(adduct_hint or "").lower()
    rows: list[dict[str, Any]] = []
    for label, offset in adducts.items():
        neutral = mz - float(offset)
        priority = 0.50
        if label in {"[M-H]-", "[M+H]+"}:
            priority = 1.00
        if "cl" in hint and "Cl" in label:
            priority = max(priority, 0.90)
        if ("form" in hint or "cho2" in hint) and "HCOO" in label:
            priority = max(priority, 0.90)
        if ("acet" in hint or "ch3coo" in hint) and "CH3COO" in label:
            priority = max(priority, 0.90)
        if "br" in hint and "Br" in label:
            priority = max(priority, 0.90)
        if "tfa" in hint or "cf3" in hint:
            priority = max(priority, 0.90 if "CF3" in label else priority)
        if "(m-2h)/2" in hint or "2h" in hint:
            neutral = 2.0 * mz + 2.0 * COMMON_NEGATIVE_ADDUCTS["[M-H]-"] * -1.0
            label = "[M-2H]2-"
            priority = 0.95
        rows.append(
            {
                "mass_hypothesis_adduct": label,
                "neutral_mass": float(neutral),
                "adduct_priority": float(priority),
                "adduct_source": "common_adduct_reference",
            }
        )
    rows = sorted(rows, key=lambda row: (-row["adduct_priority"], abs(row["neutral_mass"] - mz), row["mass_hypothesis_adduct"]))
    dedup: list[dict[str, Any]] = []
    seen: set[str] = set()
    for row in rows:
        key = str(row["mass_hypothesis_adduct"])
        if key not in seen:
            dedup.append(row)
            seen.add(key)
    return dedup


def peptide_compositions_for_mass(
    target_mass: float,
    *,
    tolerance_da: float = 0.05,
    min_len: int = 2,
    max_len: int = 4,
    max_results: int = 40,
) -> pd.DataFrame:
    if not np.isfinite(float(target_mass)) or target_mass <= 0:
        return pd.DataFrame()
    rows: list[dict[str, Any]] = []
    names = sorted(AA_RESIDUES)
    residue_masses = {name: value[1] for name, value in AA_RESIDUES.items()}
    for length in range(int(min_len), int(max_len) + 1):
        for combo in itertools.combinations_with_replacement(names, length):
            if "Pyr" in combo and combo.count("Pyr") > 1:
                continue
            mass = WATER_MASS + sum(residue_masses[name] for name in combo)
            err = abs(float(mass) - float(target_mass))
            if err <= float(tolerance_da):
                rows.append(
                    {
                        "peptide_composition": "-".join(combo),
                        "peptide_length": int(length),
                        "peptide_mass": float(mass),
                        "peptide_mass_error_da": float(err),
                        "leu_ile_ambiguous": bool("Leu" in combo or "Ile" in combo),
                    }
                )
    if not rows:
        return pd.DataFrame()
    return pd.DataFrame(rows).sort_values(["peptide_mass_error_da", "peptide_composition"], kind="mergesort").head(int(max_results)).reset_index(drop=True)


def peptide_smiles_for_sequence(sequence: Iterable[str]) -> str:
    parts = list(sequence)
    if not parts:
        return ""
    chunks: list[str] = []
    for idx, aa in enumerate(parts):
        side = AA_RESIDUES.get(str(aa), ("", 0.0, ""))[2]
        if side == "pro":
            unit = "N1CCCC1C(=O)" if idx < len(parts) - 1 else "N1CCCC1C(=O)O"
        elif side == "pyr":
            unit = "N1C(=O)CCC1C(=O)" if idx < len(parts) - 1 else "N1C(=O)CCC1C(=O)O"
        else:
            branch = f"({side})" if side else ""
            unit = f"NC{branch}C(=O)" if idx < len(parts) - 1 else f"NC{branch}C(=O)O"
        chunks.append(unit)
    smiles = "".join(chunks).replace("(=O)N", "(=O)N")
    return canonicalize_smiles(smiles) or ""


def peptide_candidates_from_composition(composition: str, *, max_permutations: int = 24) -> list[str]:
    names = [part for part in str(composition or "").split("-") if part]
    if not names:
        return []
    out: list[str] = []
    seen_sequences: set[tuple[str, ...]] = set()
    for perm in itertools.permutations(names):
        if perm in seen_sequences:
            continue
        seen_sequences.add(perm)
        smi = peptide_smiles_for_sequence(perm)
        if smi and smi not in out:
            out.append(smi)
        if len(out) >= int(max_permutations):
            break
    return out


def reattach_transform(core_smiles: str, transform_name: str, limit: int = 80) -> list[tuple[str, str]]:
    canonical = canonicalize_smiles(core_smiles) or ""
    if not canonical:
        return []
    spec = TRANSFORMS.get(str(transform_name))
    if spec is None:
        return []
    products = operator_specific_modifications(canonical, spec.generation_operators, limit=limit)
    out: list[tuple[str, str]] = []
    for candidate, operators in products.items():
        product = canonicalize_smiles(candidate)
        if product and product != canonical:
            out.append((product, "|".join(sorted(operators))))
    return out[: int(limit)]


def structure_key(smiles: Any) -> tuple[str, str]:
    canonical = canonicalize_smiles(smiles) or ""
    return canonical, inchikey_from_smiles(canonical) if canonical else ""


def candidate_mass_score(candidate_mass: float, target_mass: float, tolerance_da: float) -> float:
    if not np.isfinite(candidate_mass) or not np.isfinite(target_mass):
        return 0.0
    err = abs(float(candidate_mass) - float(target_mass))
    return float(np.clip(1.0 / (1.0 + err / max(float(tolerance_da), 1e-6)), 0.0, 1.0))


def profile_tolerance(profile: str) -> float:
    return 1.0 if str(profile).replace("non_x_holdout_", "") == "strict_1da" else 5.0


def xid_from_any(row: pd.Series) -> str:
    for col in ["control_id", "xid", "Unknown", "feature_key"]:
        if col in row.index:
            xid = parse_xid(row.get(col))
            if xid:
                return xid
    return xid_from_feature_key(row.get("feature_key", ""))


def summarize_transform_counts(evidence: pd.DataFrame) -> pd.DataFrame:
    if evidence.empty:
        return pd.DataFrame(columns=["track", "transform_name", "n_hypotheses", "n_controls"])
    return (
        evidence.groupby(["track", "transform_name"], dropna=False)
        .agg(n_hypotheses=("transform_name", "size"), n_controls=("control_id", "nunique"))
        .reset_index()
        .sort_values(["track", "n_hypotheses", "transform_name"], ascending=[True, False, True], kind="mergesort")
    )


def canonical_smiles_ok(smiles: Any) -> bool:
    if not smiles:
        return False
    if not HAS_RDKIT or Chem is None:
        return bool(canonicalize_smiles(smiles))
    return Chem.MolFromSmiles(str(smiles)) is not None


def formula_distance_score(left_formula: Any, right_formula: Any) -> float:
    left = formula_counts(left_formula)
    right = formula_counts(right_formula)
    if not left or not right:
        return 0.5
    keys = set(left) | set(right)
    distance = sum(abs(int(left.get(k, 0)) - int(right.get(k, 0))) for k in keys)
    size = sum(max(int(left.get(k, 0)), int(right.get(k, 0))) for k in keys)
    return float(np.clip(1.0 - distance / max(size, 1), 0.0, 1.0))
