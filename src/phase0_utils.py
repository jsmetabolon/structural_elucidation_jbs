from __future__ import annotations

import ast
import math
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

try:
    from rdkit import Chem
    from rdkit import RDLogger
    from rdkit.Chem import Descriptors, rdMolDescriptors

    HAS_RDKIT = True
    RDLogger.DisableLog("rdApp.warning")
    RDLogger.DisableLog("rdApp.error")
except ImportError:
    HAS_RDKIT = False


PROTON_MASS = 1.007276
AMMONIUM_MASS = 18.03382
SODIUM_MINUS_PROTON_MASS = 22.989218
POTASSIUM_MINUS_PROTON_MASS = 38.9632
C13_SHIFT_MASS = 1.003355
VALID_ION_MODES = {"negative", "positive"}


def clean_col_name(value: object) -> str:
    return str(value).strip().strip('"').lstrip("\ufeff")


def clean_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out.columns = [clean_col_name(c) for c in out.columns]
    return out


def read_csv_clean(path: Path, **kwargs) -> pd.DataFrame:
    try:
        df = pd.read_csv(path, **kwargs)
    except Exception:
        fallback_kwargs = dict(kwargs)
        fallback_kwargs.setdefault("engine", "python")
        df = pd.read_csv(path, **fallback_kwargs)
    return clean_columns(df)


def first_nonempty(values: Iterable[object]) -> object:
    for value in values:
        if pd.isna(value):
            continue
        if isinstance(value, str) and not value.strip():
            continue
        return value
    return np.nan


def canonicalize_smiles(smiles: str) -> str | None:
    if not HAS_RDKIT or not smiles or not str(smiles).strip():
        return None
    mol = Chem.MolFromSmiles(str(smiles).strip())
    if mol is None:
        return None
    try:
        Chem.SanitizeMol(mol)
        return Chem.MolToSmiles(mol, isomericSmiles=True)
    except Exception:
        return None


def exact_mol_weight(smiles: str) -> float | None:
    if not HAS_RDKIT:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        return float(Descriptors.ExactMolWt(mol))
    except Exception:
        return None


def molecular_formula(smiles: str) -> str | None:
    if not HAS_RDKIT:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        return str(rdMolDescriptors.CalcMolFormula(mol))
    except Exception:
        return None


def normalize_ion_mode(ion_mode: str) -> str:
    mode = str(ion_mode).strip().lower()
    if mode not in VALID_ION_MODES:
        raise ValueError(f"Unsupported ion mode {ion_mode!r}. Expected one of {sorted(VALID_ION_MODES)}.")
    return mode


def neutral_mass_from_precursor(parent_mz: float, ion_mode: str = "negative") -> float:
    mode = normalize_ion_mode(ion_mode)
    if mode == "negative":
        return float(parent_mz) + PROTON_MASS  # [M-H]-
    return float(parent_mz) - PROTON_MASS  # [M+H]+


def adduct_neutral_mass_error(predicted_mass: float, parent_mz: float, ion_mode: str = "negative") -> float:
    neutral_mass = neutral_mass_from_precursor(parent_mz, ion_mode=ion_mode)
    return abs(float(predicted_mass) - neutral_mass)


def closest_neutral_mass(parent_mz: float, predicted_mass: float, ion_mode: str = "negative") -> float:
    _ = float(predicted_mass)
    return neutral_mass_from_precursor(parent_mz, ion_mode=ion_mode)


def adduct_neutral_mass_error_ppm(predicted_mass: float, parent_mz: float, ion_mode: str = "negative") -> float:
    reference_mass = closest_neutral_mass(parent_mz, predicted_mass, ion_mode=ion_mode)
    denom = max(abs(reference_mass), 1e-12)
    return abs(float(predicted_mass) - reference_mass) / denom * 1e6


def _adduct_neutral_mass_hypotheses(parent_mz: float, ion_mode: str) -> list[tuple[str, float]]:
    mode = normalize_ion_mode(ion_mode)
    mz = float(parent_mz)
    if mode == "negative":
        return [
            ("[M-H]-", mz + PROTON_MASS),
            ("[M+Cl]-", mz - 34.969402),
            ("[M+HCOO]-", mz - 44.998201),
            ("[M+CH3COO]-", mz - 59.013851),
            ("[M-2H]2-", 2.0 * mz + 2.0 * PROTON_MASS),
            ("[2M-H]-", (mz + PROTON_MASS) / 2.0),
        ]
    return [
        ("[M+H]+", mz - PROTON_MASS),
        ("[M+NH4]+", mz - AMMONIUM_MASS),
        ("[M+Na]+", mz - SODIUM_MINUS_PROTON_MASS),
        ("[M+K]+", mz - POTASSIUM_MINUS_PROTON_MASS),
        ("[M+2H]2+", 2.0 * mz - 2.0 * PROTON_MASS),
        ("[2M+H]+", (mz - PROTON_MASS) / 2.0),
    ]


def best_adduct_neutral_mass_match(
    predicted_mass: float,
    parent_mz: float,
    ion_mode: str = "negative",
    include_isotopes: bool = True,
) -> dict[str, object]:
    """Return the closest neutral-mass/adduct hypothesis for a predicted molecule."""
    pred = float(predicted_mass)
    best: dict[str, object] | None = None
    for adduct, neutral_mass in _adduct_neutral_mass_hypotheses(float(parent_mz), ion_mode):
        if not np.isfinite(neutral_mass) or neutral_mass <= 0:
            continue
        shifts = [0, 1, 2] if include_isotopes else [0]
        for isotope_shift_count in shifts:
            adjusted_neutral = float(neutral_mass) - float(isotope_shift_count) * C13_SHIFT_MASS
            abs_err = abs(pred - adjusted_neutral)
            row = {
                "abs_mass_error_da": float(abs_err),
                "hypothesis": adduct,
                "adduct": adduct,
                "neutral_mass": float(adjusted_neutral),
                "unshifted_neutral_mass": float(neutral_mass),
                "isotope_shift_count": int(isotope_shift_count),
            }
            if best is None or float(row["abs_mass_error_da"]) < float(best["abs_mass_error_da"]):
                best = row
    if best is None:
        neutral_mass = neutral_mass_from_precursor(parent_mz, ion_mode=ion_mode)
        best = {
            "abs_mass_error_da": abs(pred - float(neutral_mass)),
            "hypothesis": "[M-H]-" if normalize_ion_mode(ion_mode) == "negative" else "[M+H]+",
            "adduct": "[M-H]-" if normalize_ion_mode(ion_mode) == "negative" else "[M+H]+",
            "neutral_mass": float(neutral_mass),
            "unshifted_neutral_mass": float(neutral_mass),
            "isotope_shift_count": 0,
        }
    return best


def strict_tolerance_da(
    parent_mz: float,
    ppm: float = 10.0,
    da_floor: float = 0.005,
    ion_mode: str = "negative",
) -> float:
    neutral_approx = max(neutral_mass_from_precursor(parent_mz, ion_mode=ion_mode), 0.0)
    return max(neutral_approx * float(ppm) * 1e-6, float(da_floor))


def passes_mass_filter(
    predicted_mass: float,
    parent_mz: float,
    ppm: float = 10.0,
    da_floor: float = 0.005,
    ion_mode: str = "negative",
) -> bool:
    err = adduct_neutral_mass_error(predicted_mass, parent_mz, ion_mode=ion_mode)
    return err <= strict_tolerance_da(parent_mz, ppm=ppm, da_floor=da_floor, ion_mode=ion_mode)


def soft_mass_score(
    predicted_mass: float,
    parent_mz: float,
    window_da: float = 0.5,
    ion_mode: str = "negative",
) -> float:
    err = adduct_neutral_mass_error(predicted_mass, parent_mz, ion_mode=ion_mode)
    if err <= 0.001:
        return 1.0
    if err >= float(window_da):
        return 0.0
    return 1.0 - (err / float(window_da))


def pre_strict_shortlist_score(
    mass_score: float,
    mean_log_prob: float,
    support_frac: float,
    latent_cosine: float = float("nan"),
) -> float:
    norm_lp = 1.0 / (1.0 + np.exp(-float(mean_log_prob) - 1.0))
    cos_term = max(0.0, float(latent_cosine)) if not np.isnan(latent_cosine) else 0.0
    return 0.40 * float(mass_score) + 0.25 * norm_lp + 0.20 * float(support_frac) + 0.15 * cos_term


def extract_xid(text: str) -> str | None:
    match = re.search(r"(X-\d+)", str(text))
    return match.group(1) if match else None


def safe_literal_list(value: object) -> list[str]:
    if pd.isna(value):
        return []
    text = str(value).strip()
    if not text:
        return []
    try:
        parsed = ast.literal_eval(text)
    except Exception:
        return []
    if isinstance(parsed, list):
        return [str(v) for v in parsed]
    return []


def load_jj_mm_merged_info(path: Path, sheet_name: str = "JJ + MM Merged Info") -> pd.DataFrame:
    raw = pd.read_excel(path, sheet_name=sheet_name, header=None)
    if raw.empty:
        return pd.DataFrame()

    header_row_idx = 0
    expected_header_tokens = {"LIB_ID", "CompID", "Unknown"}
    for idx in range(min(len(raw), 10)):
        values = {
            clean_col_name(value)
            for value in raw.iloc[idx].tolist()
            if not pd.isna(value) and str(value).strip()
        }
        if expected_header_tokens.issubset(values):
            header_row_idx = idx
            break

    header = [
        clean_col_name(value) if not pd.isna(value) and str(value).strip() else f"col_{idx}"
        for idx, value in enumerate(raw.iloc[header_row_idx].tolist())
    ]
    df = raw.iloc[header_row_idx + 1 :].copy()
    df.columns = header
    df = clean_columns(df)
    return df.dropna(how="all").reset_index(drop=True)


def load_comp_id_list(path: Path) -> set[int]:
    df = read_csv_clean(path)
    if "comp_id" in df.columns:
        values = df["comp_id"]
    else:
        values = df.iloc[:, 0]
    return {
        int(v)
        for v in pd.to_numeric(values, errors="coerce").dropna().astype(int).tolist()
    }


def load_embedding_id_set(path: Path) -> set[int]:
    if not path.exists():
        return set()
    ids: set[int] = set()
    for emb_path in path.glob("*.pt"):
        try:
            ids.add(int(emb_path.stem))
        except ValueError:
            continue
    return ids


def reconstruct_shortlisted_candidates(
    raw_samples_csv: Path,
    feature_meta_csv: Path,
    soft_window_da: float = 0.5,
    pre_strict_topn: int = 50,
) -> pd.DataFrame:
    raw_df = read_csv_clean(raw_samples_csv)
    feature_meta_df = read_csv_clean(feature_meta_csv)

    raw_df["feature_key"] = raw_df["feature_key"].astype(str)
    raw_df["smiles"] = raw_df["smiles"].astype(str)
    raw_df["mean_log_prob"] = pd.to_numeric(raw_df["mean_log_prob"], errors="coerce")

    canonical_counts: dict[str, Counter] = defaultdict(Counter)
    canonical_log_probs: dict[str, dict[str, float]] = defaultdict(dict)
    draw_counts = raw_df.groupby("feature_key").size().to_dict()
    canonical_cache: dict[str, str | None] = {}

    for row in raw_df.itertuples(index=False):
        raw_smi = str(row.smiles)
        feature_key = str(row.feature_key)
        mean_log_prob = float(row.mean_log_prob) if not pd.isna(row.mean_log_prob) else float("nan")

        if raw_smi not in canonical_cache:
            canonical_cache[raw_smi] = canonicalize_smiles(raw_smi)
        canonical_smiles = canonical_cache[raw_smi]
        if canonical_smiles is None:
            continue

        canonical_counts[feature_key][canonical_smiles] += 1
        prev = canonical_log_probs[feature_key].get(canonical_smiles, float("-inf"))
        if not np.isnan(mean_log_prob) and mean_log_prob > prev:
            canonical_log_probs[feature_key][canonical_smiles] = mean_log_prob

    feature_meta_df["feature_key"] = feature_meta_df["feature_key"].astype(str)
    rows: list[dict[str, object]] = []

    for meta in feature_meta_df.to_dict(orient="records"):
        feature_key = str(meta["feature_key"])
        parent_mz = float(meta["parent_mz"])
        feature_label = str(meta.get("feature_label", extract_xid(feature_key) or ""))
        comp_id = int(meta["comp_id"])
        draws = int(draw_counts.get(feature_key, 0))

        scored: list[dict[str, object]] = []
        for canonical_smiles, support_count in canonical_counts.get(feature_key, Counter()).items():
            pred_mass = exact_mol_weight(canonical_smiles)
            if pred_mass is None:
                continue

            support_frac = float(support_count) / max(draws, 1)
            mean_log_prob = canonical_log_probs.get(feature_key, {}).get(canonical_smiles, float("nan"))
            mass_score = soft_mass_score(pred_mass, parent_mz, window_da=soft_window_da)
            shortlist_score = pre_strict_shortlist_score(
                mass_score=mass_score,
                mean_log_prob=mean_log_prob,
                support_frac=support_frac,
                latent_cosine=float("nan"),
            )

            scored.append(
                {
                    "feature_key": feature_key,
                    "feature_label": feature_label,
                    "comp_id": comp_id,
                    "parent_mz": parent_mz,
                    "predicted_smiles": canonical_smiles,
                    "predicted_formula": molecular_formula(canonical_smiles),
                    "predicted_exact_mass": pred_mass,
                    "mean_log_prob": mean_log_prob,
                    "generation_support_count": int(support_count),
                    "generation_support_fraction": support_frac,
                    "soft_mass_score": mass_score,
                    "shortlist_score": shortlist_score,
                    "strict_mass_pass_10ppm_0p005da": passes_mass_filter(
                        pred_mass, parent_mz, ppm=10.0, da_floor=0.005
                    ),
                    "samples_drawn": draws,
                }
            )

        scored.sort(key=lambda record: float(record["shortlist_score"]), reverse=True)
        for shortlist_rank, record in enumerate(scored, start=1):
            record["shortlist_rank"] = shortlist_rank
            record["is_shortlisted"] = shortlist_rank <= int(pre_strict_topn)
            record["n_unique_valid_canonical"] = len(scored)
            rows.append(record)

    if not rows:
        return pd.DataFrame()

    return pd.DataFrame(rows)


def finite_float_or_none(value: object) -> float | None:
    try:
        parsed = float(value)
    except Exception:
        return None
    if math.isfinite(parsed):
        return parsed
    return None
