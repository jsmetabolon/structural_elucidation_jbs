#!/usr/bin/env python3
"""Utilities for recursive substructure discovery and conditioning.

The v1 model is intentionally small and deterministic around feature creation:
it learns multi-label substructure/motif evidence from known non-X compounds,
then emits fixed-width embeddings that can be used for frozen-pool reranking or
SELFIES decoder conditioning.
"""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
import torch
import torch.nn as nn

from denovo_improvement_utils import (
    HAS_RDKIT,
    canonicalize_smiles,
    exact_mass_from_smiles,
    formula_from_smiles,
    inchikey_from_smiles,
    pick_col,
    read_table,
)
from nce_evidence_utils import (
    candidate_fragment_bins,
    ion_and_loss_bins,
    lib_ion_mode,
    mass_bin,
    motif_labels_for_smiles,
    neutral_mass_from_precursor,
)
from x_domain_analog_utils import safe_float, safe_int

if HAS_RDKIT:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
else:  # pragma: no cover - optional runtime dependency
    Chem = None
    Descriptors = None


DEFAULT_OUTPUT_ROOT = Path("models/recursive_substructure_decoder_v1")
DEFAULT_ACTIVE_LIBS = (209, 305, 400)
DEFAULT_DIAGNOSTIC_LIBS = (402,)
LIB_ORDER = [209, 305, 400, 402]
MOTIF_PREFIX = "motif:"
FRAGMENT_PREFIX = "frag:"
FORMULA_PREFIX = "formula:"
LABEL_GROUPS = ("motif", "fragment", "formula", "other")


def parse_csv(text: Any) -> list[str]:
    return [part.strip() for part in str(text or "").split(",") if part.strip()]


def parse_libs(text: Any) -> list[int]:
    return [int(part) for part in parse_csv(text)]


def resolve(path: Path, repo_root: Path) -> Path:
    return path if path.is_absolute() else repo_root / path


def stable_bucket(value: Any, modulo: int) -> int:
    digest = hashlib.sha1(str(value).encode("utf-8")).hexdigest()
    return int(digest[:12], 16) % max(int(modulo), 1)


def normalize_text(value: Any) -> str:
    text = str(value or "").strip()
    return "" if text.lower() == "nan" else text


def parse_json_list(value: Any) -> list[Any]:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return []
    if isinstance(value, list):
        return value
    text = normalize_text(value)
    if not text:
        return []
    try:
        parsed = json.loads(text)
        return parsed if isinstance(parsed, list) else []
    except Exception:
        return [part for part in text.replace(";", ",").split(",") if part.strip()]


def json_dumps(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def structure_keys(smiles: Any) -> tuple[str, str]:
    canonical = canonicalize_smiles(smiles) or ""
    return canonical, inchikey_from_smiles(canonical) if canonical else ""


def fragment_smiles_for_labels(smiles: Any, *, max_fragments: int = 48) -> list[str]:
    if not HAS_RDKIT or Chem is None:
        return []
    canonical = canonicalize_smiles(smiles)
    if not canonical:
        return []
    mol = Chem.MolFromSmiles(canonical)
    if mol is None:
        return []
    fragments: set[str] = set()
    for bond in mol.GetBonds():
        if bond.IsInRing() or bond.GetBondTypeAsDouble() != 1.0:
            continue
        try:
            fragged = Chem.FragmentOnBonds(mol, [bond.GetIdx()], addDummies=False)
            for frag in Chem.GetMolFrags(fragged, asMols=True, sanitizeFrags=True):
                if frag.GetNumHeavyAtoms() < 2:
                    continue
                smi = Chem.MolToSmiles(frag, isomericSmiles=True)
                can = canonicalize_smiles(smi)
                if can:
                    fragments.add(can)
        except Exception:
            continue
    return sorted(fragments)[: int(max_fragments)]


@lru_cache(maxsize=300_000)
def labels_from_smiles(smiles: Any) -> list[str]:
    canonical = canonicalize_smiles(smiles)
    if not canonical:
        return []
    labels: set[str] = set()
    for motif in motif_labels_for_smiles(canonical):
        labels.add(f"{MOTIF_PREFIX}{motif}")
    for fragment in fragment_smiles_for_labels(canonical):
        labels.add(f"{FRAGMENT_PREFIX}{fragment}")
        formula = formula_from_smiles(fragment)
        if formula:
            labels.add(f"{FORMULA_PREFIX}{formula}")
    parent_formula = formula_from_smiles(canonical)
    if parent_formula:
        labels.add(f"{FORMULA_PREFIX}{parent_formula}")
    return sorted(labels)


def motif_label_set(labels: Iterable[str]) -> set[str]:
    return {str(label)[len(MOTIF_PREFIX) :] for label in labels if str(label).startswith(MOTIF_PREFIX)}


def label_group(label: Any) -> str:
    text = str(label)
    if text.startswith(MOTIF_PREFIX):
        return "motif"
    if text.startswith(FRAGMENT_PREFIX):
        return "fragment"
    if text.startswith(FORMULA_PREFIX):
        return "formula"
    return "other"


def label_groups_for_vocab(vocab: Iterable[Any]) -> np.ndarray:
    return np.asarray([label_group(label) for label in vocab], dtype=object)


def thresholds_by_group(
    thresholds: dict[str, Any] | None = None,
    *,
    default: float = 0.35,
) -> dict[str, float]:
    out = {group: float(default) for group in LABEL_GROUPS}
    out["global"] = float(default)
    if thresholds:
        source = thresholds.get("thresholds_by_group", thresholds) if isinstance(thresholds, dict) else {}
        for key, value in source.items():
            try:
                out[str(key)] = float(value)
            except Exception:
                continue
    return out


def load_excluded_truth(
    *,
    known_x: Path | None = None,
    known_x_all: Path | None = None,
    x_known_catalog: Path | None = None,
    extra_paths: Iterable[Path] = (),
) -> tuple[set[int], set[str], set[str]]:
    """Return excluded comp_ids, canonical SMILES, and InChIKeys.

    This intentionally treats broad inferred X structures as truth for leakage
    avoidance when the curation catalog is present.
    """

    comp_ids: set[int] = set()
    canonicals: set[str] = set()
    inchikeys: set[str] = set()

    paths: list[Path] = []
    if known_x is not None:
        paths.append(Path(known_x))
    if known_x_all is not None:
        paths.append(Path(known_x_all))
    if x_known_catalog is not None:
        paths.append(Path(x_known_catalog))
    paths.extend(Path(path) for path in extra_paths)

    for path in paths:
        if not path.exists():
            continue
        try:
            df = read_table(path)
        except Exception:
            continue
        comp_col = pick_col(df, ["CompID", "comp_id", "COMP_ID"])
        if comp_col is not None:
            for value in pd.to_numeric(df[comp_col], errors="coerce").dropna().astype(int).tolist():
                comp_ids.add(int(value))
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
        for value in df[smiles_col].tolist():
            canonical, inchikey = structure_keys(value)
            if canonical:
                canonicals.add(canonical)
            if inchikey:
                inchikeys.add(inchikey)
    return comp_ids, canonicals, inchikeys


def load_structure_rows(paths: Iterable[Path], libs: Iterable[int]) -> pd.DataFrame:
    wanted_libs = {int(lib) for lib in libs}
    rows: list[dict[str, Any]] = []
    for path in paths:
        if not Path(path).exists():
            continue
        df = read_table(Path(path))
        comp_col = pick_col(df, ["comp_id", "CompID", "COMP_ID", "id"])
        smiles_col = pick_col(df, ["canonical_smiles", "smiles", "SMILES"])
        if comp_col is None or smiles_col is None:
            continue
        lib_col = pick_col(df, ["lib_id", "LIB_ID", "parent_mz_source_lib"])
        lib_ids_col = pick_col(df, ["lib_ids", "library_ids"])
        parent_col = pick_col(df, ["parent_mz", "precursor_mz", "mz"])
        ri_col = pick_col(df, ["ri", "RI", "rt", "RT"])
        name_col = pick_col(df, ["compound_name", "name", "known_name"])
        for _, row in df.iterrows():
            try:
                comp_id = int(row[comp_col])
            except Exception:
                continue
            smiles = normalize_text(row.get(smiles_col))
            canonical = canonicalize_smiles(smiles)
            if not canonical:
                continue
            row_libs: set[int] = set()
            if lib_ids_col is not None:
                for token in str(row.get(lib_ids_col, "")).replace(";", ",").split(","):
                    try:
                        row_libs.add(int(float(token)))
                    except Exception:
                        pass
            if lib_col is not None:
                try:
                    row_libs.add(int(float(row.get(lib_col))))
                except Exception:
                    pass
            if not row_libs:
                row_libs = set(wanted_libs)
            for lib_id in sorted(row_libs & wanted_libs):
                rows.append(
                    {
                        "comp_id": int(comp_id),
                        "lib_id": int(lib_id),
                        "ion_mode": lib_ion_mode(int(lib_id)),
                        "canonical_smiles": canonical,
                        "inchikey": inchikey_from_smiles(canonical),
                        "parent_mz": safe_float(row.get(parent_col), math.nan) if parent_col else math.nan,
                        "ri": safe_float(row.get(ri_col), math.nan) if ri_col else math.nan,
                        "compound_name": normalize_text(row.get(name_col)) if name_col else "",
                        "source_csv": str(path),
                    }
                )
    if not rows:
        return pd.DataFrame()
    out = pd.DataFrame(rows).drop_duplicates(["comp_id", "lib_id", "canonical_smiles"]).reset_index(drop=True)
    return out


def feature_summary_from_table(feature_table: Path, libs: Iterable[int]) -> pd.DataFrame:
    if not Path(feature_table).exists():
        return pd.DataFrame()
    df = read_table(Path(feature_table))
    feature_col = pick_col(df, ["feature_key", "FeatureKey", "feature_id"])
    comp_col = pick_col(df, ["comp_id", "CompID", "COMP_ID"])
    lib_col = pick_col(df, ["lib_id", "LIB_ID", "library", "parent_mz_source_lib"])
    parent_col = pick_col(df, ["parent_mz", "precursor_mz", "parent_mass", "mz"])
    mz_col = pick_col(df, ["mass", "mz", "m/z", "ion_mz"])
    intensity_col = pick_col(df, ["intensity", "rel_intensity", "rel_abundance"])
    ri_col = pick_col(df, ["ri", "RI", "rt", "RT", "retention_index", "retention_time"])
    name_col = pick_col(df, ["compound_name", "name", "known_name"])
    if comp_col is None or parent_col is None:
        return pd.DataFrame()
    wanted = {int(lib) for lib in libs}
    work = df.copy()
    if lib_col is not None:
        work["_lib_id"] = pd.to_numeric(work[lib_col], errors="coerce")
        work = work[work["_lib_id"].isin(wanted)].copy()
    else:
        work["_lib_id"] = next(iter(sorted(wanted))) if wanted else 0
    if feature_col is None:
        work["_feature_key"] = [
            f"known_{safe_int(lib, 0)}_{safe_int(comp, -1)}"
            for lib, comp in zip(work["_lib_id"], work[comp_col])
        ]
        feature_col = "_feature_key"
    rows: list[dict[str, Any]] = []
    for feature_key, group in work.groupby(feature_col, dropna=False, sort=True):
        lib_id = safe_int(group["_lib_id"].dropna().iloc[0], 0) if group["_lib_id"].notna().any() else 0
        if lib_id not in wanted:
            continue
        parent_mz = safe_float(group[parent_col].dropna().iloc[0], math.nan) if group[parent_col].notna().any() else math.nan
        ion_bins: list[int] = []
        loss_bins: list[int] = []
        intensity_by_bin: dict[int, float] = {}
        if mz_col is not None and np.isfinite(parent_mz):
            peaks = group.copy()
            peaks["_mz"] = pd.to_numeric(peaks[mz_col], errors="coerce")
            if intensity_col is not None:
                peaks["_intensity"] = pd.to_numeric(peaks[intensity_col], errors="coerce").fillna(0.0)
                peaks = peaks.sort_values("_intensity", ascending=False, kind="mergesort")
            else:
                peaks["_intensity"] = 1.0
            for _, peak in peaks.head(80).iterrows():
                mz = safe_float(peak.get("_mz"), math.nan)
                if not np.isfinite(mz):
                    continue
                ion_bin, loss_bin = ion_and_loss_bins(mz, parent_mz, int(lib_id))
                if ion_bin:
                    ion_bins.append(int(ion_bin))
                    intensity_by_bin[int(ion_bin)] = max(float(intensity_by_bin.get(int(ion_bin), 0.0)), safe_float(peak.get("_intensity"), 1.0))
                if loss_bin:
                    loss_bins.append(int(loss_bin))
        rows.append(
            {
                "feature_key": normalize_text(feature_key),
                "comp_id": safe_int(group[comp_col].dropna().iloc[0], -1) if group[comp_col].notna().any() else -1,
                "lib_id": int(lib_id),
                "ion_mode": lib_ion_mode(int(lib_id)),
                "parent_mz": float(parent_mz),
                "ri": safe_float(group[ri_col].dropna().iloc[0], math.nan) if ri_col and group[ri_col].notna().any() else math.nan,
                "feature_label": normalize_text(group[name_col].dropna().iloc[0]) if name_col and group[name_col].notna().any() else "",
                "peak_count": int(group.shape[0]),
                "observed_ion_bins_json": json_dumps(sorted(set(ion_bins))[:120]),
                "observed_loss_bins_json": json_dumps(sorted(set(loss_bins))[:120]),
                "observed_intensity_by_bin_json": json_dumps({str(k): round(v, 6) for k, v in intensity_by_bin.items()}),
            }
        )
    return pd.DataFrame(rows)


def fallback_feature_summary_from_structures(structures: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for _, row in structures.iterrows():
        lib_id = safe_int(row.get("lib_id"), 0)
        comp_id = safe_int(row.get("comp_id"), -1)
        parent_mz = safe_float(row.get("parent_mz"), math.nan)
        if not np.isfinite(parent_mz):
            mass = exact_mass_from_smiles(row.get("canonical_smiles"))
            if np.isfinite(mass):
                parent_mz = mass - 1.007276466812 if lib_id in {209, 305} else mass + 1.007276466812
        rows.append(
            {
                "feature_key": f"known_{lib_id}_{comp_id}",
                "comp_id": comp_id,
                "lib_id": lib_id,
                "ion_mode": lib_ion_mode(lib_id),
                "parent_mz": parent_mz,
                "ri": safe_float(row.get("ri"), math.nan),
                "feature_label": normalize_text(row.get("compound_name")),
                "peak_count": 0,
                "observed_ion_bins_json": "[]",
                "observed_loss_bins_json": "[]",
                "observed_intensity_by_bin_json": "{}",
            }
        )
    return pd.DataFrame(rows)


def load_external_embedding_map(index_csv: Path | None, matrix_path: Path | None = None) -> tuple[dict[str, np.ndarray], int]:
    if index_csv is None or not Path(index_csv).exists():
        return {}, 0
    index = read_table(Path(index_csv))
    vector_cols = [
        col
        for col in index.columns
        if str(col).startswith("substructure_emb_")
        or (
            str(col).startswith("embedding_")
            and str(col) not in {"embedding_row", "embedding_dim"}
        )
    ]
    if vector_cols:
        out: dict[str, np.ndarray] = {}
        for _, row in index.iterrows():
            keys = [normalize_text(row.get("feature_key")), str(safe_int(row.get("comp_id"), -1))]
            vec = pd.to_numeric(row[vector_cols], errors="coerce").fillna(0.0).to_numpy(dtype=np.float32)
            for key in keys:
                if key and key != "-1":
                    out[key] = vec
        dim = len(vector_cols)
        return out, dim
    if matrix_path is None or not Path(matrix_path).exists():
        for suffix in ["feature_ms2smiles_embedding.npy", "substructure_embeddings.npy", "x_substructure_embeddings.npy"]:
            candidate = Path(index_csv).with_name(suffix)
            if candidate.exists():
                matrix_path = candidate
                break
    if matrix_path is None or not Path(matrix_path).exists():
        return {}, 0
    if str(matrix_path).endswith(".pt"):
        matrix_obj = torch.load(Path(matrix_path), map_location="cpu")
        matrix = matrix_obj.detach().cpu().numpy() if torch.is_tensor(matrix_obj) else np.asarray(matrix_obj)
    else:
        matrix = np.load(Path(matrix_path))
    if matrix.ndim != 2:
        return {}, 0
    row_col = pick_col(index, ["embedding_row", "row", "index"])
    out: dict[str, np.ndarray] = {}
    for idx, row in index.iterrows():
        row_idx = safe_int(row.get(row_col), idx) if row_col is not None else int(idx)
        if row_idx < 0 or row_idx >= matrix.shape[0]:
            continue
        vec = np.asarray(matrix[row_idx], dtype=np.float32)
        for key in [normalize_text(row.get("feature_key")), str(safe_int(row.get("comp_id"), -1))]:
            if key and key != "-1":
                out[key] = vec
    return out, int(matrix.shape[1])


def build_feature_vector(
    row: pd.Series,
    *,
    spectrum_bins: int,
    external_embedding: np.ndarray | None = None,
    external_embedding_dim: int = 0,
    feature_latent: np.ndarray | None = None,
    feature_latent_dim: int = 0,
) -> np.ndarray:
    ion_vec = np.zeros(int(spectrum_bins), dtype=np.float32)
    loss_vec = np.zeros(int(spectrum_bins), dtype=np.float32)
    for bin_id in parse_json_list(row.get("observed_ion_bins_json")):
        value = safe_int(bin_id, 0)
        if value:
            ion_vec[value % int(spectrum_bins)] = 1.0
    for bin_id in parse_json_list(row.get("observed_loss_bins_json")):
        value = safe_int(bin_id, 0)
        if value:
            loss_vec[value % int(spectrum_bins)] = 1.0
    parent_mz = safe_float(row.get("parent_mz"), math.nan)
    lib_id = safe_int(row.get("lib_id"), 0)
    neutral = neutral_mass_from_precursor(parent_mz, lib_id) if np.isfinite(parent_mz) else math.nan
    ri = safe_float(row.get("ri"), math.nan)
    peak_count = safe_float(row.get("peak_count"), 0.0)
    meta = np.array(
        [
            np.clip((parent_mz if np.isfinite(parent_mz) else 0.0) / 1000.0, 0.0, 2.0),
            np.clip((neutral if np.isfinite(neutral) else 0.0) / 1000.0, 0.0, 2.0),
            np.clip((ri if np.isfinite(ri) else 0.0) / 5000.0, 0.0, 2.0),
            np.clip(peak_count / 200.0, 0.0, 2.0),
        ],
        dtype=np.float32,
    )
    lib_one_hot = np.zeros(len(LIB_ORDER), dtype=np.float32)
    if lib_id in LIB_ORDER:
        lib_one_hot[LIB_ORDER.index(lib_id)] = 1.0
    ion_mode = lib_ion_mode(lib_id)
    mode_vec = np.array([1.0, 0.0] if ion_mode == "negative" else [0.0, 1.0], dtype=np.float32)
    ext = np.zeros(int(external_embedding_dim), dtype=np.float32)
    if external_embedding is not None and int(external_embedding_dim) > 0:
        ext[: min(int(external_embedding_dim), len(external_embedding))] = np.asarray(external_embedding, dtype=np.float32)[: int(external_embedding_dim)]
    latent = np.zeros(int(feature_latent_dim), dtype=np.float32)
    if feature_latent is not None and int(feature_latent_dim) > 0:
        latent[: min(int(feature_latent_dim), len(feature_latent))] = np.asarray(feature_latent, dtype=np.float32)[: int(feature_latent_dim)]
    return np.concatenate([ion_vec, loss_vec, meta, lib_one_hot, mode_vec, ext, latent]).astype(np.float32)


def build_input_matrix(
    features: pd.DataFrame,
    *,
    spectrum_bins: int,
    external_embedding_by_key: dict[str, np.ndarray] | None = None,
    external_embedding_dim: int = 0,
    feature_latent_by_key: dict[str, np.ndarray] | None = None,
    feature_latent_dim: int = 0,
) -> np.ndarray:
    vectors: list[np.ndarray] = []
    external_embedding_by_key = external_embedding_by_key or {}
    feature_latent_by_key = feature_latent_by_key or {}
    for _, row in features.iterrows():
        keys = [normalize_text(row.get("feature_key")), str(safe_int(row.get("comp_id"), -1))]
        ext = next((external_embedding_by_key[key] for key in keys if key in external_embedding_by_key), None)
        latent = next((feature_latent_by_key[key] for key in keys if key in feature_latent_by_key), None)
        vectors.append(
            build_feature_vector(
                row,
                spectrum_bins=int(spectrum_bins),
                external_embedding=ext,
                external_embedding_dim=int(external_embedding_dim),
                feature_latent=latent,
                feature_latent_dim=int(feature_latent_dim),
            )
        )
    if not vectors:
        return np.zeros((0, int(2 * spectrum_bins + 10 + external_embedding_dim + feature_latent_dim)), dtype=np.float32)
    return np.vstack(vectors).astype(np.float32)


@dataclass
class RecursiveModelOutput:
    logits_pass1: torch.Tensor
    logits_pass2: torch.Tensor
    embedding_pass1: torch.Tensor
    embedding_pass2: torch.Tensor


class RecursiveSubstructureNet(nn.Module):
    def __init__(
        self,
        input_dim: int,
        target_dim: int,
        hidden_dim: int = 512,
        embedding_dim: int = 128,
        dropout: float = 0.1,
    ) -> None:
        super().__init__()
        self.input_dim = int(input_dim)
        self.target_dim = int(target_dim)
        self.hidden_dim = int(hidden_dim)
        self.embedding_dim = int(embedding_dim)
        self.encoder = nn.Sequential(
            nn.Linear(self.input_dim, self.hidden_dim),
            nn.ReLU(),
            nn.Dropout(float(dropout)),
            nn.Linear(self.hidden_dim, self.hidden_dim),
            nn.ReLU(),
        )
        self.embedding_head1 = nn.Linear(self.hidden_dim, self.embedding_dim)
        self.logit_head1 = nn.Linear(self.embedding_dim, self.target_dim)
        self.recursive_encoder = nn.Sequential(
            nn.Linear(self.input_dim + self.target_dim + self.embedding_dim, self.hidden_dim),
            nn.ReLU(),
            nn.Dropout(float(dropout)),
            nn.Linear(self.hidden_dim, self.hidden_dim),
            nn.ReLU(),
        )
        self.embedding_head2 = nn.Linear(self.hidden_dim, self.embedding_dim)
        self.logit_head2 = nn.Linear(self.embedding_dim, self.target_dim)

    def forward(self, x: torch.Tensor) -> RecursiveModelOutput:
        x = x.float()
        h1 = self.encoder(x)
        emb1 = torch.tanh(self.embedding_head1(h1))
        logits1 = self.logit_head1(emb1)
        pass1_probs = torch.sigmoid(logits1).detach()
        h2_input = torch.cat([x, pass1_probs, emb1], dim=1)
        h2 = self.recursive_encoder(h2_input)
        emb2 = torch.tanh(self.embedding_head2(h2))
        logits2 = self.logit_head2(emb2)
        return RecursiveModelOutput(logits1, logits2, emb1, emb2)


def multilabel_metrics(y_true: np.ndarray, y_prob: np.ndarray, threshold: float = 0.35) -> dict[str, float]:
    if y_true.size == 0 or y_prob.size == 0:
        return {"micro_precision": 0.0, "micro_recall": 0.0, "micro_f1": 0.0, "mean_labels_per_row": 0.0}
    y_pred = (y_prob >= float(threshold)).astype(np.int32)
    y_true_i = (y_true > 0).astype(np.int32)
    tp = int(((y_pred == 1) & (y_true_i == 1)).sum())
    fp = int(((y_pred == 1) & (y_true_i == 0)).sum())
    fn = int(((y_pred == 0) & (y_true_i == 1)).sum())
    precision = tp / max(tp + fp, 1)
    recall = tp / max(tp + fn, 1)
    f1 = (2 * precision * recall / max(precision + recall, 1e-12)) if (precision + recall) > 0 else 0.0
    return {
        "micro_precision": float(precision),
        "micro_recall": float(recall),
        "micro_f1": float(f1),
        "mean_labels_per_row": float(y_true_i.sum(axis=1).mean()) if y_true_i.ndim == 2 and y_true_i.shape[0] else 0.0,
    }


def top_predictions_rows(
    *,
    feature_index: pd.DataFrame,
    vocab: list[str],
    probs: np.ndarray,
    embeddings: np.ndarray,
    top_k: int = 25,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for row_idx, row in feature_index.reset_index(drop=True).iterrows():
        values = probs[row_idx] if row_idx < probs.shape[0] else np.zeros(len(vocab), dtype=np.float32)
        order = np.argsort(-values)[: int(top_k)]
        top_labels = [
            {"label": str(vocab[idx]), "probability": round(float(values[idx]), 6)}
            for idx in order
            if float(values[idx]) > 0
        ]
        motifs = {
            str(vocab[idx])[len(MOTIF_PREFIX) :]: round(float(values[idx]), 6)
            for idx in order
            if str(vocab[idx]).startswith(MOTIF_PREFIX)
        }
        rows.append(
            {
                **row.to_dict(),
                "embedding_row": int(row_idx),
                "embedding_dim": int(embeddings.shape[1]) if embeddings.ndim == 2 else 0,
                "top_substructure_labels_json": json_dumps(top_labels),
                "top_motif_probs_json": json_dumps(motifs),
                "max_substructure_probability": float(np.max(values)) if values.size else 0.0,
                "mean_substructure_probability": float(np.mean(values)) if values.size else 0.0,
                "positive_label_count_0p35": int((values >= 0.35).sum()) if values.size else 0,
            }
        )
    return pd.DataFrame(rows)


def candidate_label_scores(candidate_smiles: Any, vocab: list[str], probs: np.ndarray) -> dict[str, Any]:
    labels = set(labels_from_smiles(candidate_smiles))
    candidate_vec = np.zeros(len(vocab), dtype=np.float32)
    label_to_idx = {label: idx for idx, label in enumerate(vocab)}
    for label in labels:
        idx = label_to_idx.get(label)
        if idx is not None:
            candidate_vec[idx] = 1.0
    probs_vec = np.asarray(probs, dtype=np.float32)
    if probs_vec.shape[0] != len(vocab):
        probs_vec = np.zeros(len(vocab), dtype=np.float32)
    hit_probs = [float(probs_vec[label_to_idx[label]]) for label in labels if label in label_to_idx]
    weighted = float(np.mean(hit_probs)) if hit_probs else 0.0
    denom = float(np.linalg.norm(candidate_vec) * np.linalg.norm(probs_vec))
    cosine = float(np.dot(candidate_vec, probs_vec) / denom) if denom > 0 else 0.0
    motif_labels = motif_label_set(labels)
    predicted_motifs = {
        vocab[idx][len(MOTIF_PREFIX) :]
        for idx, value in enumerate(probs_vec)
        if vocab[idx].startswith(MOTIF_PREFIX) and float(value) >= 0.35
    }
    motif_overlap = len(motif_labels & predicted_motifs) / max(math.sqrt(len(motif_labels) * len(predicted_motifs)), 1.0) if motif_labels and predicted_motifs else 0.5
    return {
        "recursive_candidate_label_count": int(len(labels)),
        "recursive_substructure_weighted_score": float(np.clip(weighted, 0.0, 1.0)),
        "recursive_label_cosine": float(np.clip(cosine, 0.0, 1.0)),
        "recursive_motif_score": float(np.clip(motif_overlap, 0.0, 1.0)),
        "recursive_candidate_labels_json": json_dumps(sorted(labels)[:80]),
        "recursive_candidate_motifs_json": json_dumps(sorted(motif_labels)),
        "recursive_predicted_motifs_json": json_dumps(sorted(predicted_motifs)),
    }


def grouped_candidate_label_scores(
    candidate_smiles: Any,
    vocab: list[str],
    probs: np.ndarray,
    thresholds: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Score candidate labels against calibrated recursive probabilities.

    The older V1 score treats every substructure label as one pool. V1.1 keeps
    that score available, but also reports motif, fragment-SMILES, and
    formula-fragment agreement separately so rare high-precision evidence does
    not get washed out by the larger fragment vocabulary.
    """

    base = candidate_label_scores(candidate_smiles, vocab, probs)
    labels = set(labels_from_smiles(candidate_smiles))
    label_to_idx = {label: idx for idx, label in enumerate(vocab)}
    probs_vec = np.asarray(probs, dtype=np.float32)
    if probs_vec.shape[0] != len(vocab):
        probs_vec = np.zeros(len(vocab), dtype=np.float32)
    threshold_map = thresholds_by_group(thresholds)

    group_rows: dict[str, dict[str, Any]] = {}
    for group in LABEL_GROUPS:
        group_labels = sorted(label for label in labels if label_group(label) == group and label in label_to_idx)
        hit_probs = [float(probs_vec[label_to_idx[label]]) for label in group_labels]
        threshold = float(threshold_map.get(group, threshold_map.get("global", 0.35)))
        if hit_probs:
            mean_prob = float(np.mean(hit_probs))
            max_prob = float(np.max(hit_probs))
            calibrated_fraction = float(np.mean(np.asarray(hit_probs, dtype=float) >= threshold))
        else:
            mean_prob = 0.5 if group in {"motif", "formula"} else 0.0
            max_prob = 0.0
            calibrated_fraction = 0.0
        group_rows[group] = {
            "count": int(len(group_labels)),
            "mean_prob": float(np.clip(mean_prob, 0.0, 1.0)),
            "max_prob": float(np.clip(max_prob, 0.0, 1.0)),
            "calibrated_fraction": float(np.clip(calibrated_fraction, 0.0, 1.0)),
        }

    motif_component = 0.65 * group_rows["motif"]["mean_prob"] + 0.35 * group_rows["motif"]["calibrated_fraction"]
    fragment_component = 0.65 * group_rows["fragment"]["mean_prob"] + 0.35 * group_rows["fragment"]["calibrated_fraction"]
    formula_component = 0.65 * group_rows["formula"]["mean_prob"] + 0.35 * group_rows["formula"]["calibrated_fraction"]
    calibrated_score = 0.35 * motif_component + 0.45 * fragment_component + 0.20 * formula_component
    base.update(
        {
            "recursive_motif_group_score": float(np.clip(motif_component, 0.0, 1.0)),
            "recursive_fragment_group_score": float(np.clip(fragment_component, 0.0, 1.0)),
            "recursive_formula_group_score": float(np.clip(formula_component, 0.0, 1.0)),
            "recursive_calibrated_group_score": float(np.clip(calibrated_score, 0.0, 1.0)),
            "recursive_motif_candidate_label_count": int(group_rows["motif"]["count"]),
            "recursive_fragment_candidate_label_count": int(group_rows["fragment"]["count"]),
            "recursive_formula_candidate_label_count": int(group_rows["formula"]["count"]),
            "recursive_group_scores_json": json_dumps(group_rows),
        }
    )
    return base
