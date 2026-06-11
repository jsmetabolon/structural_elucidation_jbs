#!/usr/bin/env python3
"""Shared helpers for lactoyl-amino-acid heuristic scoring."""

from __future__ import annotations

import json
import math
import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from x_domain_analog_utils import safe_float


PROTON_MASS = 1.007276466812
LACTOYL_DELTA_MASS = 72.021129368  # C3H4O2 = lactic acid - H2O
GOOD_TANIMOTO = 0.7


AMINO_ACIDS: dict[str, dict[str, Any]] = {
    "Ala": {"names": ["alanine", "ala"], "mass": 89.047678, "formula": "C3H7NO2", "classes": ["small"]},
    "Arg": {"names": ["arginine", "arg"], "mass": 174.111676, "formula": "C6H14N4O2", "classes": ["basic"]},
    "Asn": {"names": ["asparagine", "asn"], "mass": 132.053493, "formula": "C4H8N2O3", "classes": ["amide"]},
    "Asp": {"names": ["aspartic", "aspartate", "asp"], "mass": 133.037508, "formula": "C4H7NO4", "classes": ["acidic"]},
    "Cys": {"names": ["cysteine", "cys"], "mass": 121.019749, "formula": "C3H7NO2S", "classes": ["sulfur"]},
    "Gln": {"names": ["glutamine", "gln"], "mass": 146.069143, "formula": "C5H10N2O3", "classes": ["amide"]},
    "Glu": {"names": ["glutamic", "glutamate", "glu"], "mass": 147.053158, "formula": "C5H9NO4", "classes": ["acidic"]},
    "Gly": {"names": ["glycine", "gly"], "mass": 75.032028, "formula": "C2H5NO2", "classes": ["small"]},
    "His": {"names": ["histidine", "his"], "mass": 155.069477, "formula": "C6H9N3O2", "classes": ["basic", "aromatic"]},
    "Ile": {"names": ["isoleucine", "ile"], "mass": 131.094629, "formula": "C6H13NO2", "classes": ["branched"]},
    "Leu": {"names": ["leucine", "leu"], "mass": 131.094629, "formula": "C6H13NO2", "classes": ["branched"]},
    "Lys": {"names": ["lysine", "lys"], "mass": 146.105528, "formula": "C6H14N2O2", "classes": ["basic"]},
    "Met": {"names": ["methionine", "met"], "mass": 149.051049, "formula": "C5H11NO2S", "classes": ["sulfur"]},
    "Phe": {"names": ["phenylalanine", "phe"], "mass": 165.078979, "formula": "C9H11NO2", "classes": ["aromatic"]},
    "Pro": {"names": ["proline", "pro"], "mass": 115.063329, "formula": "C5H9NO2", "classes": ["cyclic"]},
    "Ser": {"names": ["serine", "ser"], "mass": 105.042593, "formula": "C3H7NO3", "classes": ["hydroxy"]},
    "Thr": {"names": ["threonine", "thr"], "mass": 119.058243, "formula": "C4H9NO3", "classes": ["hydroxy"]},
    "Trp": {"names": ["tryptophan", "trp"], "mass": 204.089878, "formula": "C11H12N2O2", "classes": ["aromatic"]},
    "Tyr": {"names": ["tyrosine", "tyr"], "mass": 181.073893, "formula": "C9H11NO3", "classes": ["aromatic", "hydroxy"]},
    "Val": {"names": ["valine", "val"], "mass": 117.078979, "formula": "C5H11NO2", "classes": ["branched"]},
}

THEORETICAL_N_LACTOYL_SMILES: dict[str, str] = {
    "Ala": "CC(O)C(=O)NC(C)C(=O)O",
    "Arg": "CC(O)C(=O)NC(CCCN=C(N)N)C(=O)O",
    "Asn": "CC(O)C(=O)NC(CC(=O)N)C(=O)O",
    "Asp": "CC(O)C(=O)NC(CC(=O)O)C(=O)O",
    "Cys": "CC(O)C(=O)NC(CS)C(=O)O",
    "Gln": "CC(O)C(=O)NC(CCC(=O)N)C(=O)O",
    "Glu": "CC(O)C(=O)NC(CCC(=O)O)C(=O)O",
    "Gly": "CC(O)C(=O)NCC(=O)O",
    "His": "CC(O)C(=O)NC(CC1=CN=CN1)C(=O)O",
    "Ile": "CC(O)C(=O)NC(C(C)CC)C(=O)O",
    "Leu": "CC(O)C(=O)NC(CC(C)C)C(=O)O",
    "Lys": "CC(O)C(=O)NC(CCCCN)C(=O)O",
    "Met": "CC(O)C(=O)NC(CCSC)C(=O)O",
    "Phe": "CC(O)C(=O)NC(CC1=CC=CC=C1)C(=O)O",
    "Pro": "CC(O)C(=O)N1CCCC1C(=O)O",
    "Ser": "CC(O)C(=O)NC(CO)C(=O)O",
    "Thr": "CC(O)C(=O)NC(C(O)C)C(=O)O",
    "Trp": "CC(O)C(=O)NC(CC1=CNC2=CC=CC=C12)C(=O)O",
    "Tyr": "CC(O)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)O",
    "Val": "CC(O)C(=O)NC(C(C)C)C(=O)O",
}


AA_BY_NAME: dict[str, str] = {}
for aa, meta in AMINO_ACIDS.items():
    for name in meta["names"]:
        AA_BY_NAME[name.lower()] = aa


def parse_xid(value: Any) -> str:
    match = re.search(r"(X-\d+)", str(value or ""))
    return match.group(1) if match else ""


def normalize_name(value: Any) -> str:
    return re.sub(r"[^a-z0-9]+", " ", str(value or "").lower()).strip()


def infer_lactoyl_amino_acid(value: Any) -> str:
    text = normalize_name(value)
    if "lactoyl" not in text and "lactyl" not in text and "n lac" not in text:
        return ""
    for name, aa in sorted(AA_BY_NAME.items(), key=lambda item: len(item[0]), reverse=True):
        compact = text.replace(" ", "")
        if re.search(rf"\b{name}\b", text) or f"lactoyl{name}" in compact or f"lactyl{name}" in compact:
            return aa
    return ""


def formula_add(formula: str, delta: str = "C3H4O2") -> str:
    counts: dict[str, int] = {}
    for source in [formula, delta]:
        for elem, number in re.findall(r"([A-Z][a-z]?)(\d*)", source):
            counts[elem] = counts.get(elem, 0) + int(number or 1)
    order = ["C", "H", "N", "O", "S", "P", "Cl", "Br", "F", "I"]
    return "".join(f"{elem}{counts[elem] if counts[elem] != 1 else ''}" for elem in order if counts.get(elem, 0))


def lactoyl_mass(aa: str) -> float:
    return float(AMINO_ACIDS[aa]["mass"]) + LACTOYL_DELTA_MASS


def lactoyl_formula(aa: str) -> str:
    return formula_add(str(AMINO_ACIDS[aa]["formula"]))


def neutral_mass_from_parent(parent_mz: Any, ion_mode: str = "negative", adduct: str | None = None) -> float:
    mz = safe_float(parent_mz, math.nan)
    if math.isnan(mz):
        return math.nan
    text = str(adduct or "").upper()
    mode = str(ion_mode or "").lower()
    if text in {"[M-H]-", "M-H", "M-H-"} or mode.startswith("neg"):
        return mz + PROTON_MASS
    if text in {"[M+H]+", "M+H", "M+H+"} or mode.startswith("pos"):
        return mz - PROTON_MASS
    return mz + PROTON_MASS


def normalize_peak_table(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, low_memory=False)
    if df.empty:
        return pd.DataFrame(columns=["feature_key", "control_id", "comp_id", "parent_mz", "rt", "fragment_mz", "intensity", "source_path"])
    colmap = {col.lower(): col for col in df.columns}
    def get(*names: str) -> pd.Series:
        for name in names:
            if name.lower() in colmap:
                return df[colmap[name.lower()]]
        return pd.Series(np.nan, index=df.index)
    out = pd.DataFrame(
        {
            "feature_key": get("feature_key", "FEATURE_KEY").astype(str),
            "control_id": get("control_id", "compound_name", "COMPOUND_NAME").map(parse_xid),
            "comp_id": pd.to_numeric(get("comp_id", "COMP_ID"), errors="coerce"),
            "parent_mz": pd.to_numeric(get("parent_mz", "PARENT_MASS", "parent_mass"), errors="coerce"),
            "rt": pd.to_numeric(get("rt", "RT"), errors="coerce"),
            "fragment_mz": pd.to_numeric(get("fragment_mz", "mass", "MASS"), errors="coerce"),
            "intensity": pd.to_numeric(get("intensity", "MASS_RATIO"), errors="coerce").fillna(0.0),
            "source_path": str(path),
        }
    )
    out = out[out["fragment_mz"].notna()].copy()
    return out


def load_peak_tables(paths: list[Path]) -> pd.DataFrame:
    frames = [normalize_peak_table(path) for path in paths if Path(path).exists()]
    if not frames:
        return pd.DataFrame(columns=["feature_key", "control_id", "comp_id", "parent_mz", "rt", "fragment_mz", "intensity", "source_path"])
    out = pd.concat(frames, ignore_index=True, sort=False)
    out = out.drop_duplicates(["feature_key", "fragment_mz", "intensity", "source_path"]).reset_index(drop=True)
    return out


def peak_bins(peaks: pd.DataFrame, bin_width: float = 0.01) -> pd.DataFrame:
    if peaks.empty:
        return peaks.assign(fragment_bin=pd.Series(dtype=float), rel_intensity=pd.Series(dtype=float))
    out = peaks.copy()
    out["fragment_bin"] = (pd.to_numeric(out["fragment_mz"], errors="coerce") / float(bin_width)).round().astype("Int64").astype(float) * float(bin_width)
    max_int = out.groupby("feature_key")["intensity"].transform("max").replace(0, np.nan)
    out["rel_intensity"] = (out["intensity"] / max_int).fillna(0.0).clip(0.0, 1.0)
    return out


def bin_support_for_feature(peaks: pd.DataFrame, feature_key: str, bins: list[float], tolerance: float = 0.02) -> tuple[int, float, str]:
    scoped = peaks[peaks["feature_key"].astype(str).eq(str(feature_key))]
    if scoped.empty or not bins:
        return 0, 0.0, ""
    values = pd.to_numeric(scoped["fragment_mz"], errors="coerce").dropna().to_numpy()
    intensities = pd.to_numeric(scoped["rel_intensity"] if "rel_intensity" in scoped.columns else scoped["intensity"], errors="coerce").fillna(0.0).to_numpy()
    matched: list[float] = []
    score = 0.0
    for target in bins:
        diffs = np.abs(values - float(target))
        if len(diffs) and float(diffs.min()) <= float(tolerance):
            idx = int(diffs.argmin())
            matched.append(float(target))
            score += float(intensities[idx])
    return len(matched), min(score / max(len(bins), 1), 1.0), json.dumps(matched)


def dataframe_to_md(df: pd.DataFrame, max_rows: int = 40) -> str:
    if df.empty:
        return "_No rows._"
    from x_domain_analog_utils import dataframe_to_markdown

    return dataframe_to_markdown(df.head(int(max_rows)))
