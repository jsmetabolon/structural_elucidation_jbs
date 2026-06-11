#!/usr/bin/env python3
# pylint: disable=broad-exception-caught,no-member,unused-import,unused-argument
"""Dependency-light chemistry evidence helpers for NCE integration experiments.

The functions here intentionally mirror the useful algorithmic pieces from the
old NewCompoundElucidation workflow without importing that application.  They
are deterministic, small enough for unit tests, and safe to run inside the
current spec2mol pipeline.
"""

from __future__ import annotations

import json
import math
import re
from functools import lru_cache
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd

from structural_elucidation.utils.denovo import (
    FORMULA_ELEMENTS,
    HAS_RDKIT,
    canonicalize_smiles,
    exact_mass_from_smiles,
    formula_counts,
    formula_from_smiles,
    inchikey_from_smiles,
    pick_col,
    read_table,
)

if HAS_RDKIT:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
else:  # pragma: no cover - runtime dependency guard
    Chem = None
    AllChem = None
    Descriptors = None
    rdMolDescriptors = None


PROTON_MASS = 1.007276466812
ELECTRON_MASS = 0.000548579909065
NEGATIVE_LIBS = {209, 305}
POSITIVE_LIBS = {400, 402}
MASS_BIN_PPM = 10.0

MONOISOTOPIC_MASS = {
    "C": 12.0,
    "H": 1.00782503223,
    "N": 14.00307400443,
    "O": 15.99491461957,
    "P": 30.97376199842,
    "S": 31.9720711744,
    "F": 18.99840316273,
    "Cl": 34.968852682,
    "Br": 78.9183376,
    "I": 126.9044719,
}

MOTIF_SMARTS = {
    "sulfate": "[SX4](=O)(=O)(O)O",
    "phosphate": "[PX4](=O)(O)(O)O",
    "carboxyl": "C(=O)[O;H,-]",
    "amide": "C(=O)N",
    "glycine": "NCC(=O)[O;H,-]",
    "taurine": "NCCS(=O)(=O)O",
    "carnitine": "C[N+](C)(C)CC(O)CC(=O)[O;H,-]",
    "glucuronide_like": "C(=O)[O;H,-]",
    "glycoside_like": "COC",
    "sterol_core": "C1CCC2(C(C1)CCC1C2CCC2(C)C(CC1)CC2)C",
    "ceramide_backbone": "CCCCCCCCCCCCCC(=O)NC(CO)C(O)",
    "glucuronide": "OC1OC(C(=O)O)C(O)C(O)C1O",
    "sulfonate": "CS(=O)(=O)O",
    "isoalloxazine": "c1cc2nc3c(=O)[nH]c(=O)nc3n(c2cc1)",
}

TARGETED_REACTION_SMARTS = {
    "hydroxyl_sulfation": "[O:1]>>[O:1]S(=O)(=O)O",
    "hydroxyl_phosphorylation": "[O:1]>>[O:1]P(=O)(O)O",
    "hydroxyl_acetylation": "[O:1]>>[O:1]C(C)=O",
    "lactoyl_amine_conjugation": "[N;H1,H2:1]>>[N:1]C(=O)C(O)C",
    "lactoyl_hydroxyl_conjugation": "[O:1]>>[O:1]C(=O)C(O)C",
    "glucuronide_hydroxyl_conjugation": "[O:1]>>[O:1]C1OC(C(=O)O)C(O)C(O)C1O",
    "glucuronide_carboxyl_conjugation": "[C:1](=O)[O;H1:2]>>[C:1](=O)[O:2]C1OC(C(=O)O)C(O)C(O)C1O",
    "glucuronide_amine_conjugation": "[N;H1,H2:1]>>[N:1]C1OC(C(=O)O)C(O)C(O)C1O",
    "amine_formylation": "[N;H1,H2:1]>>[N:1]C=O",
    "amine_acetylation": "[N;H1,H2:1]>>[N:1]C(C)=O",
    "amine_glycine_conjugation": "[N;H1,H2:1]>>[N:1]CC(=O)O",
    "amine_glutamine_like_conjugation": "[N;H1,H2:1]>>[N:1]C(=O)CCC(N)C(=O)O",
    "carboxyl_taurine_amide": "[C:1](=O)[O;H1:2]>>[C:1](=O)NCCS(=O)(=O)O",
    "carboxyl_carnitine_ester": "[C:1](=O)[O;H1:2]>>[C:1](=O)OCC(O)C[N+](C)(C)C",
    "carboxyl_methylation": "[C:1](=O)[O;H1:2]>>[C:1](=O)[O:2]C",
    "thiol_carboxyethylation": "[S;H1:1]>>[S:1]CC(C(=O)O)C(=O)O",
}


def normalize_ion_mode(value: Any) -> str:
    text = str(value or "").strip().lower()
    if text in {"209", "305", "negative", "neg", "0", "0.0"}:
        return "negative"
    if text in {"400", "402", "positive", "pos", "1", "1.0"}:
        return "positive"
    return text


def lib_ion_mode(lib_id: int) -> str:
    return "negative" if int(lib_id) in NEGATIVE_LIBS else "positive"


def neutral_mass_from_precursor(parent_mz: Any, lib_id: int) -> float:
    """Convert the observed precursor m/z to a neutral mass hypothesis.

    This matches the old NCE convention: positive parent m/z is treated as
    [M+H]+, negative parent m/z as [M-H]-.
    """

    try:
        mz = float(parent_mz)
    except Exception:
        return float("nan")
    if not np.isfinite(mz):
        return float("nan")
    if int(lib_id) in POSITIVE_LIBS:
        return mz - PROTON_MASS
    if int(lib_id) in NEGATIVE_LIBS:
        return mz + PROTON_MASS
    return mz


def neutral_ion_mass(mz: Any, lib_id: int) -> float:
    """Neutralize an MS2 ion m/z with the old NCE electron-mass convention."""

    try:
        value = float(mz)
    except Exception:
        return float("nan")
    if not np.isfinite(value):
        return float("nan")
    if int(lib_id) in POSITIVE_LIBS:
        return value + ELECTRON_MASS
    if int(lib_id) in NEGATIVE_LIBS:
        return value - ELECTRON_MASS
    return value


def mass_bin(mass: Any, ppm: float = MASS_BIN_PPM) -> int:
    value = float(mass)
    if value <= 0 or not np.isfinite(value):
        return 0
    # Old CompoundDiscovery/MS2SMILES binning used a 10 ppm logarithmic grid
    # with 1000 amu assigned to bin 299999 and smaller masses decreasing from
    # there.  Keep the same convention so legacy fragment features line up.
    return int(299999 + np.int32(math.log(1000.0 / value) / math.log(1.0 - float(ppm) / 1_000_000.0)))


def bin_base(bin_id: int, ppm: float = MASS_BIN_PPM) -> float:
    return float(1000.0 / ((1.0 - float(ppm) / 1_000_000.0) ** (int(bin_id) - 299999)))


def mass_interval(mass: Any, ppm: float = MASS_BIN_PPM) -> tuple[int, int]:
    value = float(mass)
    if value <= 0 or not np.isfinite(value):
        return (0, 0)
    delta = value * (float(ppm) / 2.0) / 1_000_000.0
    return (mass_bin(value - delta, ppm), mass_bin(value + delta, ppm))


def ion_and_loss_bins(ion_mz: Any, parent_mz: Any, lib_id: int, ppm: float = MASS_BIN_PPM) -> tuple[int, int]:
    ion_mass = neutral_ion_mass(ion_mz, lib_id)
    parent_mass = neutral_mass_from_precursor(parent_mz, lib_id)
    if not np.isfinite(ion_mass) or not np.isfinite(parent_mass):
        return (0, 0)
    return (mass_bin(ion_mass, ppm), mass_bin(max(parent_mass - ion_mass, 1e-9), ppm))


def parse_json_list(value: Any) -> list[Any]:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return []
    if isinstance(value, list):
        return value
    text = str(value or "").strip()
    if not text:
        return []
    try:
        parsed = json.loads(text)
        return parsed if isinstance(parsed, list) else []
    except Exception:
        return [part.strip() for part in re.split(r"[;,\s]+", text) if part.strip()]


def formula_mass(formula: Any) -> float:
    counts = formula_counts(formula)
    if not counts:
        return float("nan")
    mass = 0.0
    for element, count in counts.items():
        if element not in MONOISOTOPIC_MASS:
            return float("nan")
        mass += MONOISOTOPIC_MASS[element] * int(count)
    return float(mass)


def formula_harmony_score(formula: Any, neutral_mass: Any = None) -> float:
    """Simple deterministic formula plausibility proxy.

    The original NCE formula scorer used trained LOF models.  Those private
    caches are not portable, so this proxy captures the strongest cheap checks:
    mass consistency, sane H/C and heteroatom ratios, halogen balance, and DBE.
    """

    counts = formula_counts(formula)
    if not counts:
        return 0.0
    c = counts.get("C", 0)
    h = counts.get("H", 0)
    n = counts.get("N", 0)
    o = counts.get("O", 0)
    p = counts.get("P", 0)
    s = counts.get("S", 0)
    hal = counts.get("F", 0) + counts.get("Cl", 0) + counts.get("Br", 0) + counts.get("I", 0)
    score = 1.0
    if c <= 0:
        score *= 0.45
    else:
        h_c = h / max(c, 1)
        hetero_c = (n + o + p + s + hal) / max(c, 1)
        if not (0.2 <= h_c <= 3.4):
            score *= 0.55
        if hetero_c > 2.6:
            score *= 0.65
        if o > max(4 * c + 4, 12):
            score *= 0.70
    dbe = c - h / 2.0 + n / 2.0 + 1.0 - hal / 2.0
    if dbe < -0.5:
        score *= 0.35
    elif dbe < 0.0:
        score *= 0.70
    if dbe > max(30.0, c + 6.0):
        score *= 0.75
    if neutral_mass is not None:
        try:
            target = float(neutral_mass)
        except Exception:
            target = float("nan")
        f_mass = formula_mass(formula)
        if np.isfinite(target) and np.isfinite(f_mass):
            err = abs(f_mass - target)
            if err <= 0.01:
                score *= 1.0
            elif err <= 0.1:
                score *= 0.92
            elif err <= 1.0:
                score *= 0.72
            else:
                score *= 0.35
    return float(np.clip(score, 0.0, 1.0))


@lru_cache(maxsize=200_000)
def cached_formula_harmony(formula: str, neutral_mass_bucket: str = "") -> float:
    neutral_mass = None
    if neutral_mass_bucket:
        try:
            neutral_mass = float(neutral_mass_bucket)
        except Exception:
            neutral_mass = None
    return formula_harmony_score(formula, neutral_mass)


@lru_cache(maxsize=200_000)
def cached_smiles_formula(smiles: str) -> str:
    return formula_from_smiles(smiles)


@lru_cache(maxsize=200_000)
def cached_smiles_mass(smiles: str) -> float:
    return exact_mass_from_smiles(smiles)


@lru_cache(maxsize=200_000)
def cached_smiles_inchikey(smiles: str) -> str:
    return inchikey_from_smiles(smiles)


@lru_cache(maxsize=100_000)
def ri_proxy_from_smiles(smiles: str) -> float:
    if not HAS_RDKIT or Chem is None or Descriptors is None:
        return float("nan")
    mol = Chem.MolFromSmiles(str(smiles or ""))
    if mol is None:
        return float("nan")
    try:
        mass = float(Descriptors.ExactMolWt(mol))
        logp = float(Descriptors.MolLogP(mol))
        tpsa = float(Descriptors.TPSA(mol))
        rings = float(rdMolDescriptors.CalcNumRings(mol)) if rdMolDescriptors is not None else 0.0
        return float(250.0 + 1.7 * mass + 55.0 * logp + 8.0 * rings - 1.4 * tpsa)
    except Exception:
        return float("nan")


def ri_similarity_score(candidate_smiles: Any, observed_ri: Any, candidate_ri: Any = None) -> float:
    try:
        observed = float(observed_ri)
    except Exception:
        return 0.5
    if not np.isfinite(observed) or observed <= 0:
        return 0.5
    try:
        predicted = float(candidate_ri)
    except Exception:
        predicted = ri_proxy_from_smiles(str(candidate_smiles or ""))
    if not np.isfinite(predicted):
        return 0.5
    scale = 300.0 if candidate_ri is not None and str(candidate_ri) != "" else 800.0
    return float(np.clip(math.exp(-abs(observed - predicted) / scale), 0.0, 1.0))


@lru_cache(maxsize=100_000)
def motif_labels_for_smiles(smiles: str) -> tuple[str, ...]:
    labels: set[str] = set()
    if not HAS_RDKIT or Chem is None:
        return tuple()
    mol = Chem.MolFromSmiles(str(smiles or ""))
    if mol is None:
        return tuple()
    for label, smarts in MOTIF_SMARTS.items():
        patt = Chem.MolFromSmarts(smarts)
        if patt is not None and mol.HasSubstructMatch(patt):
            labels.add(label)
    formula = formula_from_smiles(smiles)
    counts = formula_counts(formula)
    if counts.get("S", 0) > 0 and counts.get("O", 0) >= 4:
        labels.add("sulfur_oxygen_rich")
    if counts.get("P", 0) > 0 and counts.get("O", 0) >= 3:
        labels.add("phosphate_like")
    if counts.get("O", 0) >= 6 and counts.get("C", 0) >= 6:
        labels.add("glycoside_like")
    return tuple(sorted(labels))


def motif_labels_for_legacy_fragment(fragment: Any) -> tuple[str, ...]:
    text = str(fragment or "").strip().lower()
    labels: set[str] = set()
    keyword_map = {
        "sulfate": ["sulfate", "s(=o)(=o)", "[s"],
        "glucuronide_like": ["glucuron", "c(=o)o", "o)co"],
        "glycoside_like": ["glucose", "ribose", "sugar", "oc(o)"],
        "glycine": ["glycine", "ncc(=o)"],
        "taurine": ["taurine", "nccs"],
        "carnitine": ["carnitine", "[n+](c)(c)c"],
        "amide": ["c(=o)n"],
        "carboxyl": ["c(=o)o"],
        "phosphate": ["phosphate", "p(=o)"],
    }
    for label, keywords in keyword_map.items():
        if any(keyword in text for keyword in keywords):
            labels.add(label)
    return tuple(sorted(labels))


@lru_cache(maxsize=100_000)
def candidate_fragment_bins(smiles: str, lib_id: int, ppm: float = MASS_BIN_PPM) -> tuple[tuple[int, ...], tuple[int, ...]]:
    if not HAS_RDKIT or Chem is None:
        return (tuple(), tuple())
    mol = Chem.MolFromSmiles(str(smiles or ""))
    if mol is None:
        return (tuple(), tuple())
    parent_mass = exact_mass_from_smiles(smiles)
    if not np.isfinite(parent_mass):
        return (tuple(), tuple())
    bins: set[int] = set()
    losses: set[int] = set()
    for bond in mol.GetBonds():
        if bond.IsInRing() or bond.GetBondTypeAsDouble() != 1.0:
            continue
        try:
            fragged = Chem.FragmentOnBonds(mol, [bond.GetIdx()], addDummies=False)
            for frag in Chem.GetMolFrags(fragged, asMols=True, sanitizeFrags=True):
                mass = float(Descriptors.ExactMolWt(frag))
                if mass <= 12.0 or mass >= parent_mass - 1.0:
                    continue
                bins.add(mass_bin(mass, ppm))
                losses.add(mass_bin(max(parent_mass - mass, 1e-9), ppm))
        except Exception:
            continue
    return (tuple(sorted(bins)), tuple(sorted(losses)))


def overlap_score(left: Iterable[Any], right: Iterable[Any]) -> float:
    left_set = {str(item) for item in left if str(item) not in {"", "0", "nan"}}
    right_set = {str(item) for item in right if str(item) not in {"", "0", "nan"}}
    if not left_set or not right_set:
        return 0.5
    overlap = len(left_set & right_set)
    denom = math.sqrt(len(left_set) * len(right_set))
    return float(np.clip(overlap / max(denom, 1.0), 0.0, 1.0))


def forward_spectrum_score(candidate_smiles: Any, observed_ion_bins: Iterable[Any], observed_loss_bins: Iterable[Any], lib_id: int) -> float:
    frag_bins, loss_bins = candidate_fragment_bins(str(candidate_smiles or ""), int(lib_id))
    ion_score = overlap_score(frag_bins, observed_ion_bins)
    loss_score = overlap_score(loss_bins, observed_loss_bins)
    return float(np.clip(0.55 * ion_score + 0.45 * loss_score, 0.0, 1.0))


def legacy_ms2smiles_summary(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(
            columns=["comp_id", "legacy_motif_count", "legacy_top_motif", "legacy_motifs_json", "legacy_mean_probability"]
        )
    df = read_table(path)
    comp_col = pick_col(df, ["comp_id", "CompID", "compound_id"])
    frag_col = pick_col(df, ["fragment", "frag_smiles", "frag_chem_name"])
    prob_col = pick_col(df, ["probability", "score", "prediction_probability"])
    lib_col = pick_col(df, ["library", "lib_id"])
    if comp_col is None or frag_col is None:
        return pd.DataFrame()
    work = df.copy()
    if lib_col is not None:
        work = work[pd.to_numeric(work[lib_col], errors="coerce").eq(209)].copy()
    work["_comp_id"] = pd.to_numeric(work[comp_col], errors="coerce")
    work = work[work["_comp_id"].notna()].copy()
    work["_motifs"] = work[frag_col].map(lambda value: list(motif_labels_for_legacy_fragment(value)))
    work["_prob"] = pd.to_numeric(work[prob_col], errors="coerce") if prob_col is not None else 1.0
    rows: list[dict[str, Any]] = []
    for comp_id, group in work.groupby("_comp_id", dropna=True):
        motifs: dict[str, int] = {}
        for values in group["_motifs"]:
            for motif in values:
                motifs[motif] = motifs.get(motif, 0) + 1
        top = sorted(motifs.items(), key=lambda item: (-item[1], item[0]))
        rows.append(
            {
                "comp_id": int(comp_id),
                "legacy_motif_count": int(sum(motifs.values())),
                "legacy_top_motif": top[0][0] if top else "",
                "legacy_motifs_json": json.dumps(motifs, sort_keys=True),
                "legacy_mean_probability": float(pd.to_numeric(group["_prob"], errors="coerce").fillna(0.0).mean()),
            }
        )
    return pd.DataFrame(rows).sort_values("comp_id").reset_index(drop=True)


def clean_structure_key(smiles: Any) -> tuple[str, str]:
    canonical = canonicalize_smiles(smiles) or ""
    return canonical, inchikey_from_smiles(canonical) if canonical else ""


def load_excluded_structure_keys(paths: Iterable[Path]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for path in paths:
        if not path.exists():
            continue
        try:
            df = read_table(path)
        except Exception:
            continue
        smiles_col = pick_col(
            df,
            [
                "selected_canonical_smiles",
                "canonical_smiles",
                "smiles",
                "SMILES",
                "structure_smiles",
            ],
        )
        if smiles_col is None:
            continue
        id_col = pick_col(df, ["XID", "xid", "control_id", "comp_id", "CompID"])
        tier_col = pick_col(df, ["evidence_tier", "truth_tier", "tier"])
        for _, row in df.iterrows():
            canonical, inchikey = clean_structure_key(row.get(smiles_col))
            if not canonical:
                continue
            rows.append(
                {
                    "source_path": str(path),
                    "source_id": "" if id_col is None else str(row.get(id_col, "")),
                    "evidence_tier": "" if tier_col is None else str(row.get(tier_col, "")),
                    "canonical_smiles": canonical,
                    "inchikey": inchikey,
                }
            )
    out = pd.DataFrame(rows)
    if out.empty:
        return pd.DataFrame(columns=["source_path", "source_id", "evidence_tier", "canonical_smiles", "inchikey"])
    return out.drop_duplicates(["canonical_smiles", "inchikey"]).reset_index(drop=True)


@lru_cache(maxsize=10_000)
def targeted_modification_candidates(smiles: str) -> tuple[tuple[str, str], ...]:
    """Generate a small set of valid targeted local edits from one seed."""

    if not HAS_RDKIT or Chem is None or AllChem is None:
        return tuple()
    canonical = canonicalize_smiles(smiles)
    if not canonical:
        return tuple()
    mol = Chem.MolFromSmiles(canonical)
    if mol is None:
        return tuple()
    out: set[tuple[str, str]] = set()
    for name, smarts in TARGETED_REACTION_SMARTS.items():
        try:
            rxn = AllChem.ReactionFromSmarts(smarts)
            products = rxn.RunReactants((mol,))
        except Exception:
            continue
        for product_tuple in products[:50]:
            for product in product_tuple:
                try:
                    Chem.SanitizeMol(product)
                    product_smiles = Chem.MolToSmiles(product, isomericSmiles=True)
                    canonical_product = canonicalize_smiles(product_smiles)
                    if canonical_product and canonical_product != canonical:
                        out.add((canonical_product, name))
                except Exception:
                    continue
    return tuple(sorted(out))
