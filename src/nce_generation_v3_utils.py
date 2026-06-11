#!/usr/bin/env python3
"""Portable NCE-style evidence and generation helpers for V3.

This module keeps the old NewCompoundElucidation ideas as small deterministic
pieces that can run inside the current spec2mol workflow.  It intentionally
does not import the old app runtime.
"""

from __future__ import annotations

import json
import math
from collections import Counter
from functools import lru_cache
from typing import Any, Iterable

import numpy as np
import pandas as pd

from denovo_improvement_utils import (
    HAS_RDKIT,
    FORMULA_ELEMENTS,
    canonicalize_smiles,
    exact_mass_from_smiles,
    formula_counts,
    formula_from_smiles,
    inchikey_from_smiles,
)
from nce_evidence_utils import (
    MASS_BIN_PPM,
    MONOISOTOPIC_MASS,
    candidate_fragment_bins,
    formula_harmony_score,
    formula_mass,
    ion_and_loss_bins,
    mass_bin,
    motif_labels_for_legacy_fragment,
    motif_labels_for_smiles,
    neutral_mass_from_precursor,
    overlap_score,
    parse_json_list,
    targeted_modification_candidates,
)
from nce_workflow_leverage_v2 import LEGACY_MOTIF_KEYWORDS, legacy_motif_candidate_score
from x_domain_analog_utils import safe_float, safe_int

if HAS_RDKIT:
    from rdkit import Chem
    from rdkit.Chem import AllChem
else:  # pragma: no cover - dependency guard
    Chem = None
    AllChem = None


TARGET_MOTIFS = (
    "sulfate",
    "glucuronide_like",
    "glycoside_like",
    "phosphate",
    "peptide_glutathione",
    "taurine",
    "carnitine",
    "bile_acid_steroid",
    "glycine",
)

OPERATOR_REACTION_PREFIXES = {
    "sulfate": ("hydroxyl_sulfation",),
    "phosphate": ("hydroxyl_phosphorylation",),
    "acetyl": ("hydroxyl_acetylation", "amine_acetylation"),
    "glycine": ("amine_glycine_conjugation",),
    "glutamine": ("amine_glutamine_like_conjugation",),
    "glutathione": ("amine_glutamine_like_conjugation",),
    "taurine": ("carboxyl_taurine_amide",),
    "carnitine": ("carboxyl_carnitine_ester",),
    "carboxyethyl": ("thiol_carboxyethylation",),
}

CUSTOM_REACTION_SMARTS = {
    "hydroxyl_glucosylation": "[O;H1:1]>>[O:1]C1OC(O)C(O)C(O)C1O",
    "hydroxyl_glucuronidation": "[O;H1:1]>>[O:1]C1OC(C(=O)O)C(O)C(O)C1O",
}


def safe_json_list(values: Iterable[Any]) -> str:
    return json.dumps([value for value in values if str(value) not in {"", "nan", "None"}], sort_keys=True)


def mass_from_bin(bin_id: int, ppm: float = MASS_BIN_PPM) -> float:
    """Inverse of the legacy logarithmic mass bin.

    `nce_evidence_utils.bin_base` preserves the historical helper name but is
    not the inverse of `mass_bin`; V3 scoring needs the actual approximate mass
    represented by a bin.
    """

    return float(1000.0 / ((1.0 - float(ppm) / 1_000_000.0) ** (int(bin_id) - 299999)))


def formula_to_string(counts: dict[str, int]) -> str:
    parts: list[str] = []
    for element in FORMULA_ELEMENTS:
        count = int(counts.get(element, 0))
        if count <= 0:
            continue
        parts.append(element if count == 1 else f"{element}{count}")
    return "".join(parts)


def dbe_from_counts(counts: dict[str, int]) -> float:
    c = counts.get("C", 0)
    h = counts.get("H", 0)
    n = counts.get("N", 0)
    hal = counts.get("F", 0) + counts.get("Cl", 0) + counts.get("Br", 0) + counts.get("I", 0)
    return float(c - h / 2.0 + n / 2.0 + 1.0 - hal / 2.0)


def formula_plausible(counts: dict[str, int], neutral_mass: float | None = None) -> bool:
    if not counts or counts.get("C", 0) <= 0:
        return False
    c = counts.get("C", 0)
    h = counts.get("H", 0)
    hetero = sum(counts.get(element, 0) for element in ["N", "O", "P", "S", "F", "Cl", "Br", "I"])
    dbe = dbe_from_counts(counts)
    if h < 0 or h > 4 * c + 2 + counts.get("N", 0):
        return False
    if dbe < -0.5 or dbe > max(35.0, c + 8.0):
        return False
    if hetero > max(4 * c + 6, 24):
        return False
    if neutral_mass is not None:
        mass = formula_mass(formula_to_string(counts))
        if not np.isfinite(mass) or abs(mass - float(neutral_mass)) > 0.5:
            return False
    return True


def known_formula_catalog(candidate_source: pd.DataFrame) -> pd.DataFrame:
    if candidate_source.empty:
        return pd.DataFrame(columns=["candidate_formula", "formula_mass", "formula_frequency", "libs_json"])
    formula_col = "candidate_formula" if "candidate_formula" in candidate_source.columns else "formula"
    if formula_col not in candidate_source.columns:
        return pd.DataFrame(columns=["candidate_formula", "formula_mass", "formula_frequency", "libs_json"])
    rows: list[dict[str, Any]] = []
    work = candidate_source.copy()
    work["_formula"] = work[formula_col].fillna("").astype(str)
    work = work[work["_formula"].str.strip().ne("") & ~work["_formula"].str.lower().eq("nan")].copy()
    for formula, group in work.groupby("_formula", dropna=True):
        mass = formula_mass(formula)
        if not np.isfinite(mass):
            continue
        libs = sorted({safe_int(value, 0) for value in group.get("lib_id", pd.Series(dtype=int)).tolist() if safe_int(value, 0) > 0})
        rows.append(
            {
                "candidate_formula": str(formula),
                "formula_mass": float(mass),
                "formula_frequency": int(group.shape[0]),
                "libs_json": json.dumps(libs),
            }
        )
    return pd.DataFrame(rows).sort_values(["formula_mass", "candidate_formula"], kind="mergesort").reset_index(drop=True)


@lru_cache(maxsize=20_000)
def enumerate_formula_candidates_synthetic(
    neutral_mass_bucket: str,
    tolerance_da: float,
    max_candidates: int,
) -> tuple[tuple[str, float, float, float, str], ...]:
    """Enumerate plausible formulas by solving H after choosing heavy atoms.

    The search is intentionally bounded for metabolite-like CHNOPS/halogen
    formulas.  It is a fallback when known non-X formulas do not cover a mass.
    """

    target = float(neutral_mass_bucket)
    tol = float(tolerance_da)
    max_c = min(80, max(1, int(target / 12.0) + 2))
    max_o = min(26, max(2, int(target / 16.0)))
    max_n = min(10, max(2, int(target / 28.0) + 4))
    max_s = min(4, max(1, int(target / 96.0) + 2))
    max_p = min(2, max(1, int(target / 180.0) + 1))
    halogen_sets = [
        {},
        {"F": 1},
        {"F": 2},
        {"Cl": 1},
        {"Cl": 2},
        {"Br": 1},
        {"I": 1},
    ]
    rows: list[tuple[str, float, float, float, str]] = []
    h_mass = MONOISOTOPIC_MASS["H"]
    for c in range(1, max_c + 1):
        c_mass = c * MONOISOTOPIC_MASS["C"]
        if c_mass > target + tol:
            break
        for n in range(0, max_n + 1):
            cn_mass = c_mass + n * MONOISOTOPIC_MASS["N"]
            if cn_mass > target + tol:
                break
            for o in range(0, max_o + 1):
                cno_mass = cn_mass + o * MONOISOTOPIC_MASS["O"]
                if cno_mass > target + tol:
                    break
                for s in range(0, max_s + 1):
                    cnos_mass = cno_mass + s * MONOISOTOPIC_MASS["S"]
                    if cnos_mass > target + tol:
                        break
                    for p in range(0, max_p + 1):
                        base_mass = cnos_mass + p * MONOISOTOPIC_MASS["P"]
                        if base_mass > target + tol:
                            break
                        for halogens in halogen_sets:
                            heavy_mass = base_mass + sum(MONOISOTOPIC_MASS[element] * count for element, count in halogens.items())
                            remaining = target - heavy_mass
                            if remaining < -tol:
                                continue
                            h = int(round(remaining / h_mass))
                            counts = {"C": c, "H": h, "N": n, "O": o, "P": p, "S": s, **halogens}
                            if not formula_plausible(counts):
                                continue
                            formula = formula_to_string(counts)
                            mass = formula_mass(formula)
                            err = abs(mass - target)
                            if err <= tol:
                                rows.append(
                                    (
                                        formula,
                                        float(mass),
                                        float(err),
                                        formula_harmony_score(formula, target),
                                        "synthetic_chnops_halogen",
                                    )
                                )
    rows = sorted(rows, key=lambda item: (item[2], -item[3], item[0]))[: int(max_candidates)]
    return tuple(rows)


def formula_candidates_for_mass(
    *,
    neutral_mass: float,
    lib_id: int,
    formula_catalog: pd.DataFrame,
    tolerance_da: float = 0.05,
    max_candidates: int = 30,
    synthetic_min: int = 0,
) -> pd.DataFrame:
    if not np.isfinite(float(neutral_mass)) or neutral_mass <= 0:
        return pd.DataFrame()
    rows: list[dict[str, Any]] = []
    if not formula_catalog.empty:
        cat = formula_catalog.copy()
        masses = pd.to_numeric(cat.get("formula_mass"), errors="coerce")
        hits = cat[(masses - float(neutral_mass)).abs().le(float(tolerance_da))].copy()
        if hits.empty and tolerance_da < 0.5:
            hits = cat[(masses - float(neutral_mass)).abs().le(0.5)].copy()
        for _, row in hits.iterrows():
            formula = str(row.get("candidate_formula", "") or "")
            mass = safe_float(row.get("formula_mass"), math.nan)
            if not formula or not np.isfinite(mass):
                continue
            rows.append(
                {
                    "candidate_formula": formula,
                    "formula_mass": mass,
                    "formula_mass_error_da": abs(mass - float(neutral_mass)),
                    "formula_harmony_score": formula_harmony_score(formula, neutral_mass),
                    "formula_source": "known_non_x_formula",
                    "formula_frequency": safe_int(row.get("formula_frequency"), 1),
                }
            )
    if int(synthetic_min) > 0 and len(rows) < int(synthetic_min):
        tol = max(float(tolerance_da), 0.02)
        synthetic = enumerate_formula_candidates_synthetic(
            f"{float(neutral_mass):.4f}",
            tol,
            max(int(max_candidates), int(synthetic_min) * 4),
        )
        for formula, mass, err, harmony, source in synthetic:
            rows.append(
                {
                    "candidate_formula": formula,
                    "formula_mass": mass,
                    "formula_mass_error_da": err,
                    "formula_harmony_score": harmony,
                    "formula_source": source,
                    "formula_frequency": 0,
                }
            )
    if not rows:
        return pd.DataFrame()
    out = pd.DataFrame(rows).drop_duplicates("candidate_formula")
    out["lib_id"] = int(lib_id)
    out = out.sort_values(
        ["formula_mass_error_da", "formula_source", "formula_harmony_score", "candidate_formula"],
        ascending=[True, True, False, True],
        kind="mergesort",
    ).head(int(max_candidates))
    for rank in range(1, out.shape[0] + 1):
        out.iloc[rank - 1, out.columns.get_loc("lib_id")] = int(lib_id)
    out.insert(0, "formula_rank", list(range(1, out.shape[0] + 1)))
    return out.reset_index(drop=True)


def consensus_bins_for_feature(rows: pd.DataFrame, *, parent_mz: float, lib_id: int, max_bins: int = 80) -> dict[str, Any]:
    if rows.empty:
        return {
            "observed_ion_bins_json": "[]",
            "observed_loss_bins_json": "[]",
            "observed_major_ion_bins_json": "[]",
            "observed_intensity_by_bin_json": "{}",
            "peak_count": 0,
        }
    mass_col = "mass" if "mass" in rows.columns else "mz" if "mz" in rows.columns else ""
    if not mass_col:
        return {
            "observed_ion_bins_json": "[]",
            "observed_loss_bins_json": "[]",
            "observed_major_ion_bins_json": "[]",
            "observed_intensity_by_bin_json": "{}",
            "peak_count": int(rows.shape[0]),
        }
    intensity_col = ""
    for col in ["intensity", "rel_intensity", "abundance", "height"]:
        if col in rows.columns:
            intensity_col = col
            break
    counter: Counter[int] = Counter()
    loss_counter: Counter[int] = Counter()
    intensity: dict[str, float] = {}
    for _, row in rows.iterrows():
        mz = safe_float(row.get(mass_col), math.nan)
        if not np.isfinite(mz) or mz <= 0:
            continue
        ion_bin, loss_bin = ion_and_loss_bins(mz, parent_mz, lib_id, MASS_BIN_PPM)
        weight = safe_float(row.get(intensity_col), 1.0) if intensity_col else 1.0
        weight = max(weight, 1.0)
        if ion_bin:
            counter[int(ion_bin)] += 1
            intensity[str(int(ion_bin))] = intensity.get(str(int(ion_bin)), 0.0) + float(weight)
        if loss_bin:
            loss_counter[int(loss_bin)] += 1
    ion_bins = [item for item, _ in counter.most_common(max_bins)]
    loss_bins = [item for item, _ in loss_counter.most_common(max_bins)]
    total = sum(intensity.values()) or 1.0
    intensity_norm = {key: round(value / total, 6) for key, value in sorted(intensity.items())}
    return {
        "observed_ion_bins_json": json.dumps(sorted(ion_bins)),
        "observed_loss_bins_json": json.dumps(sorted(loss_bins)),
        "observed_major_ion_bins_json": json.dumps(sorted([item for item, _ in counter.most_common(12)])),
        "observed_intensity_by_bin_json": json.dumps(intensity_norm, sort_keys=True),
        "peak_count": int(rows.shape[0]),
    }


def formula_ion_consistency(
    *,
    candidate_formula: Any,
    observed_ion_bins: Iterable[Any],
    observed_loss_bins: Iterable[Any],
    intensity_by_bin: dict[str, float] | None = None,
) -> dict[str, Any]:
    formula = str(candidate_formula or "")
    f_mass = formula_mass(formula)
    if not np.isfinite(f_mass) or f_mass <= 0:
        return {
            "formula_ion_match_count": 0,
            "weighted_formula_ion_score": 0.0,
            "formula_loss_score": 0.0,
            "unexplained_major_ion_count": 0,
        }
    ions = [safe_int(value, 0) for value in observed_ion_bins if safe_int(value, 0)]
    losses = [safe_int(value, 0) for value in observed_loss_bins if safe_int(value, 0)]
    intensity_by_bin = intensity_by_bin or {}
    ion_matches: list[int] = []
    weighted = 0.0
    for bin_id in ions:
        mass = mass_from_bin(bin_id)
        if 10.0 <= mass <= f_mass - 1.0:
            ion_matches.append(bin_id)
            weighted += safe_float(intensity_by_bin.get(str(bin_id)), 1.0 / max(len(ions), 1))
    loss_matches = [bin_id for bin_id in losses if 5.0 <= mass_from_bin(bin_id) <= f_mass - 1.0]
    return {
        "formula_ion_match_count": int(len(set(ion_matches))),
        "weighted_formula_ion_score": float(np.clip(weighted if intensity_by_bin else len(set(ion_matches)) / max(len(set(ions)), 1), 0.0, 1.0))
        if ions
        else 0.5,
        "formula_loss_score": float(np.clip(len(set(loss_matches)) / max(len(set(losses)), 1), 0.0, 1.0)) if losses else 0.5,
        "unexplained_major_ion_count": int(max(0, len(set(ions[:12])) - len(set(ion_matches)))),
    }


def parse_int_list_json(value: Any) -> list[int]:
    out: list[int] = []
    for item in parse_json_list(value):
        number = safe_int(item, 0)
        if number:
            out.append(number)
    return out


def parse_float_dict_json(value: Any) -> dict[str, float]:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return {}
    if isinstance(value, dict):
        return {str(k): safe_float(v, 0.0) for k, v in value.items()}
    try:
        parsed = json.loads(str(value or "{}"))
    except Exception:
        return {}
    if not isinstance(parsed, dict):
        return {}
    return {str(k): safe_float(v, 0.0) for k, v in parsed.items()}


def candidate_forward_match(
    *,
    smiles: Any,
    lib_id: int,
    observed_ion_bins: Iterable[Any],
    observed_loss_bins: Iterable[Any],
    legacy_ion_bins: Iterable[Any] = (),
    legacy_loss_bins: Iterable[Any] = (),
) -> dict[str, Any]:
    frag_bins, loss_bins = candidate_fragment_bins(str(smiles or ""), int(lib_id))
    observed_ions = {safe_int(value, 0) for value in observed_ion_bins if safe_int(value, 0)}
    observed_losses = {safe_int(value, 0) for value in observed_loss_bins if safe_int(value, 0)}
    legacy_ions = {safe_int(value, 0) for value in legacy_ion_bins if safe_int(value, 0)}
    legacy_losses = {safe_int(value, 0) for value in legacy_loss_bins if safe_int(value, 0)}
    ion_target = observed_ions | legacy_ions
    loss_target = observed_losses | legacy_losses
    ion_score = overlap_score(frag_bins, ion_target)
    loss_score = overlap_score(loss_bins, loss_target)
    legacy_bonus = 0.0
    if legacy_ions or legacy_losses:
        legacy_bonus = 0.5 * overlap_score(frag_bins, legacy_ions) + 0.5 * overlap_score(loss_bins, legacy_losses)
    score = float(np.clip(0.48 * ion_score + 0.38 * loss_score + 0.14 * legacy_bonus, 0.0, 1.0))
    return {
        "forward_spectrum_score": score,
        "matched_ion_bins_json": json.dumps(sorted(set(frag_bins) & ion_target)),
        "matched_loss_bins_json": json.dumps(sorted(set(loss_bins) & loss_target)),
        "candidate_fragment_bins_json": json.dumps(sorted(frag_bins)[:80]),
        "candidate_loss_bins_json": json.dumps(sorted(loss_bins)[:80]),
    }


def motif_class_scores_from_evidence(row: pd.Series) -> dict[str, Any]:
    labels = set()
    source = "none"
    for col in ["legacy_top_motif_v2", "legacy_motifs_json_v2", "legacy_historical_notes_v2", "unknown_name", "details"]:
        if col not in row.index:
            continue
        value = row.get(col)
        if str(value or "").strip():
            labels.update(motif_labels_for_legacy_fragment(value))
            text = str(value or "").lower()
            for label, needles in LEGACY_MOTIF_KEYWORDS.items():
                if any(needle in text for needle in needles):
                    labels.add(label)
    try:
        parsed = json.loads(str(row.get("legacy_motifs_json_v2", "{}") or "{}"))
        if isinstance(parsed, dict):
            labels.update(str(k) for k, v in parsed.items() if safe_float(v, 0.0) > 0)
    except Exception:
        pass
    if labels:
        source = "legacy_ms2smiles_rule_proxy"
    scores: dict[str, float] = {motif: 0.0 for motif in TARGET_MOTIFS}
    for motif in labels:
        if motif in scores:
            scores[motif] = max(scores[motif], 0.8)
    top = sorted(scores.items(), key=lambda item: (-item[1], item[0]))
    return {
        **{f"motif_{motif}_score": float(score) for motif, score in scores.items()},
        "motif_top_label": top[0][0] if top and top[0][1] > 0 else "",
        "motif_score_source": source,
        "motif_labels_json": json.dumps(sorted(label for label in labels if label)),
    }


def select_generation_operators(evidence_row: pd.Series, worklist_operator: str = "") -> list[str]:
    text = " ".join(
        str(value or "")
        for value in [
            worklist_operator,
            evidence_row.get("legacy_top_motif_v2", ""),
            evidence_row.get("legacy_motifs_json_v2", ""),
            evidence_row.get("legacy_historical_notes_v2", ""),
            evidence_row.get("motif_top_label", ""),
            evidence_row.get("motif_labels_json", ""),
        ]
    ).lower()
    operators: list[str] = []
    mapping = [
        ("sulfate", ["sulfate", "sulfo"]),
        ("phosphate", ["phosphate", "phospho"]),
        ("glucuronide", ["glucuron", "uronic"]),
        ("glycoside", ["glycoside", "glycos", "sugar", "hexose"]),
        ("glycine", ["glycine"]),
        ("glutamine", ["glutamine", "glutathione", "peptide", "gamma-glu"]),
        ("taurine", ["taurine"]),
        ("carnitine", ["carnitine"]),
        ("acetyl", ["acetyl"]),
        ("carboxyethyl", ["carboxyethyl"]),
    ]
    for operator, needles in mapping:
        if any(needle in text for needle in needles) and operator not in operators:
            operators.append(operator)
    if not operators:
        operators.extend(["mass_delta_local", "mcs_rgroup", "brics_crossover"])
    return operators[:6]


def _run_custom_reaction(smiles: str, reaction_name: str, limit: int = 40) -> list[str]:
    if not HAS_RDKIT or Chem is None or AllChem is None:
        return []
    smarts = CUSTOM_REACTION_SMARTS.get(reaction_name)
    if not smarts:
        return []
    canonical = canonicalize_smiles(smiles)
    if not canonical:
        return []
    mol = Chem.MolFromSmiles(canonical)
    if mol is None:
        return []
    out: list[str] = []
    try:
        rxn = AllChem.ReactionFromSmarts(smarts)
        products = rxn.RunReactants((mol,))
    except Exception:
        return out
    for product_tuple in products[: int(limit)]:
        for product in product_tuple:
            try:
                Chem.SanitizeMol(product)
                candidate = canonicalize_smiles(Chem.MolToSmiles(product, isomericSmiles=True))
            except Exception:
                candidate = None
            if candidate and candidate != canonical and candidate not in out:
                out.append(candidate)
            if len(out) >= int(limit):
                return out
    return out


def operator_specific_modifications(
    anchor_smiles: str,
    operators: Iterable[str],
    *,
    limit: int = 120,
) -> dict[str, set[str]]:
    generated: dict[str, set[str]] = {}
    selected = {str(operator) for operator in operators if str(operator)}
    all_targeted = targeted_modification_candidates(anchor_smiles)
    for candidate, reaction_name in all_targeted:
        keep = False
        for operator in selected:
            prefixes = OPERATOR_REACTION_PREFIXES.get(operator, ())
            if reaction_name in prefixes:
                keep = True
                break
        if keep:
            generated.setdefault(candidate, set()).add(reaction_name)
        if len(generated) >= int(limit):
            break
    if "glycoside" in selected:
        for candidate in _run_custom_reaction(anchor_smiles, "hydroxyl_glucosylation", limit=30):
            generated.setdefault(candidate, set()).add("hydroxyl_glucosylation")
    if "glucuronide" in selected:
        for candidate in _run_custom_reaction(anchor_smiles, "hydroxyl_glucuronidation", limit=30):
            generated.setdefault(candidate, set()).add("hydroxyl_glucuronidation")
    return generated


def formulas_match_or_near(candidate_formula: str, top_formulas: Iterable[str], min_overlap: float = 0.88) -> tuple[bool, str, float]:
    candidate_counts = formula_counts(candidate_formula)
    if not candidate_counts:
        return False, "", 0.0
    best_formula = ""
    best_score = 0.0
    for formula in top_formulas:
        formula = str(formula or "")
        counts = formula_counts(formula)
        if not counts:
            continue
        keys = set(candidate_counts) | set(counts)
        distance = sum(abs(candidate_counts.get(key, 0) - counts.get(key, 0)) for key in keys)
        size = sum(max(candidate_counts.get(key, 0), counts.get(key, 0)) for key in keys)
        score = 1.0 - distance / max(size, 1)
        if score > best_score:
            best_formula = formula
            best_score = float(score)
    return bool(best_score >= float(min_overlap)), best_formula, best_score


def structure_key(smiles: Any) -> tuple[str, str]:
    canonical = canonicalize_smiles(smiles) or ""
    return canonical, inchikey_from_smiles(canonical) if canonical else ""


def candidate_motif_score(smiles: str, evidence_row: pd.Series) -> float:
    labels = motif_labels_for_smiles(smiles)
    return legacy_motif_candidate_score(labels, evidence_row.get("legacy_motifs_json_v2", "{}"), evidence_row.get("legacy_top_motif_v2", ""))
