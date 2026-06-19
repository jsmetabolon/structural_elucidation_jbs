#!/usr/bin/env python3
"""Batch SMILES prediction for merged-library spectral features.

Architecture overview
---------------------
The pipeline encodes MS2 feature tables with a mass-conditioned ``Net1D``
encoder (optionally calibrated), then runs two complementary tracks:

**Track A — mass-gated retrieval**
    Cosine KNN over a comp_id-keyed embedding bank within progressive
    precursor-mass gates (0.05 → 0.5 → 1.0 Da). Optional reranking blends
    latent cosine, fragment explanation, HTML formula bonus, prevalence, and
    structural-neighbor scores into ``retrieval_score``.

**Track B — constrained de novo generation**
    Mass- (and optionally RI-) conditioned decoder sampling from query
    latents, with soft-mass shortlisting and strict adduct-aware gating.
    Optional retrieval-guided latent proposals interpolate toward top analogs.

**Merged output**
    Track A rank-1 rows and Track B survivors are combined via
    ``--merge_policy`` (default ``pin_high_rank1_joint_rest``) into
    ``{prefix}_merged_predictions.csv``. Policies may pull SELFIES-scored
    candidates from ``candidate_metrics_all`` and joint-sort retrieval +
    de novo rows per feature. ``--skip_constrained_denovo`` emits
    retrieval-only merged outputs.

Primary artifacts (under ``--output_dir``):
    ``{prefix}_retrieval_report.csv``, ``{prefix}_denovo_report.csv``,
    ``{prefix}_merged_predictions.csv``, ``{prefix}_feature_latents.pt``,
    ``{prefix}_run_summary.json``.
"""

from __future__ import annotations

from structural_elucidation.paths import REPO_ROOT

import argparse
import importlib
import json
import re
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    from rdkit.DataStructs import TanimotoSimilarity
    from rdkit.Chem import AllChem

    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False

from structural_elucidation.utils.device import DEVICE_HELP_TEXT, resolve_device
from structural_elucidation.inference.model_conv import (
    Net1D,
    infer_mass_conditioned_compatible_length,
    load_compatible_state_dict,
)
from structural_elucidation.inference.calibration import apply_latent_calibrator, load_latent_calibrator
from structural_elucidation.utils.html import (
    audit_candidate_formula_against_html_signals,
    parse_html_signals,
)
from structural_elucidation.utils.spectrum import (
    SpectrumPreprocessingConfig,
    preprocess_and_vectorize,
)
from structural_elucidation.utils.substructure import compute_structural_neighbor_scores
from structural_elucidation.utils.phase0 import (
    PROTON_MASS,
    adduct_neutral_mass_error as shared_adduct_neutral_mass_error,
    molecular_formula,
    neutral_mass_from_precursor,
    normalize_ion_mode,
    strict_tolerance_da as shared_strict_tolerance_da,
)
from structural_elucidation.inference.retrieval_scoring import (
    RetrievalScoringConfig,
    compute_retrieval_score,
    resolve_scoring_config,
)
from structural_elucidation.inference.ri_audit import load_feature_source_audit_summary
from structural_elucidation.utils.merge import (
    SELFIES_SCORE_FORMULA_VERSION,
    build_policy_high_protected_medium_low_joint,
    build_policy_pin_high_rank1_joint_rest,
    build_policy_retrieval_rank1_plus_best_selfies,
    build_policy_selfies_only,
    build_policy_top5_selfies_then_joint_rest,
    finalize_policy_rows,
    map_retrieval_rows_for_merge,
    map_selfies_rows_for_merge,
    rank_selfies_pool,
    select_selfies_rows_per_feature,
    to_merged_predictions_csv_schema,
)


# ── torch / checkpoint helpers ────────────────────────────────────────────────

def robust_torch_load(path: Path, map_location: str = "cpu", weights_only: bool = False):
    """Load with compatibility across torch versions."""
    original_rebuild = getattr(torch._utils, "_rebuild_device_tensor_from_numpy", None)

    def load() -> object:
        try:
            return torch.load(path, map_location=map_location, weights_only=weights_only)
        except TypeError:
            return torch.load(path, map_location=map_location)

    if original_rebuild is None or str(map_location) != "cpu":
        return load()

    def rebuild_on_cpu(data, dtype, device, requires_grad):
        if str(device).startswith("mps"):
            tensor = torch.from_numpy(data).to(dtype=dtype, device="cpu")
            tensor.requires_grad = requires_grad
            return tensor
        return original_rebuild(data, dtype, device, requires_grad)

    torch._utils._rebuild_device_tensor_from_numpy = rebuild_on_cpu
    try:
        return load()
    finally:
        torch._utils._rebuild_device_tensor_from_numpy = original_rebuild


def register_moses_aliases() -> None:
    """
    Older decoder artifacts reference `moses.*` modules during unpickling.
    Register aliases to the local `decoder.moses.*` package before torch.load.
    """
    vendor_root = REPO_ROOT / "vendor"
    vendor_root_str = str(vendor_root)
    if vendor_root_str not in sys.path:
        sys.path.insert(0, vendor_root_str)
    alias_map = {
        "moses": "decoder.moses",
        "moses.utils": "decoder.moses.utils",
        "moses.trans": "decoder.moses.trans",
        "moses.models_storage": "decoder.moses.models_storage",
    }
    for alias, real_name in alias_map.items():
        if alias in sys.modules:
            continue
        sys.modules[alias] = importlib.import_module(real_name)


# ── feature table I/O ─────────────────────────────────────────────────────────

def read_table(path: Path) -> pd.DataFrame:
    """Load a CSV with tolerant delimiter detection and BOM stripping."""
    try:
        df = pd.read_csv(path, sep=None, engine="python")
    except Exception:
        df = pd.read_csv(path)
    df.columns = [str(c).strip().strip('"').lstrip("\ufeff") for c in df.columns]
    return df


def load_feature_table(path: Path) -> pd.DataFrame:
    """Load encoder-ready feature table and enforce required columns."""
    if path.suffix.lower() in {".pkl", ".pickle"}:
        df = pd.read_pickle(path)
    else:
        df = read_table(path)

    expected = {
        "feature_key",
        "comp_id",
        "compound_name",
        "rt",
        "parent_mz",
        "mass",
        "intensity",
        "source_html",
    }
    missing = [c for c in expected if c not in df.columns]
    if missing:
        raise KeyError(f"Feature table missing columns: {missing}. Found: {df.columns.tolist()}")

    out = df.copy()
    out["feature_key"] = out["feature_key"].astype(str)
    out["comp_id"] = out["comp_id"].astype(int)
    out["compound_name"] = out["compound_name"].astype(str)
    out["rt"] = out["rt"].astype(float)
    out["parent_mz"] = out["parent_mz"].astype(float)
    out["mass"] = out["mass"].astype(float)
    out["intensity"] = out["intensity"].astype(float)
    out["source_html"] = out["source_html"].astype(str)
    return out


def vectorize_peaks(
    peaks: Iterable[tuple[float, float]],
    minmass: float,
    maxmass: float,
    resolution: int,
    preprocessing_config: SpectrumPreprocessingConfig,
) -> np.ndarray:
    """Vectorize peak list via shared spectrum preprocessing pipeline."""
    return preprocess_and_vectorize(
        peaks=peaks,
        config=preprocessing_config,
        minmass=minmass,
        maxmass=maxmass,
        resolution=resolution,
    )["vector"]


def normalize_scalar(value: float, min_value: float, max_value: float) -> float:
    """Linearly scale a scalar into [0, 1] using inclusive bounds."""
    if max_value <= min_value:
        return 0.0
    scaled = (float(value) - float(min_value)) / (float(max_value) - float(min_value))
    return float(np.clip(scaled, 0.0, 1.0))


# ── chemistry helpers ─────────────────────────────────────────────────────────

def validate_smiles(smiles: str) -> bool:
    """Return True when SMILES parses and sanitizes under RDKit."""
    if not HAS_RDKIT:
        return True
    if not smiles or not smiles.strip():
        return False
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    try:
        Chem.SanitizeMol(mol)
        return True
    except Exception:
        return False


def exact_mol_weight(smiles: str) -> float | None:
    """Exact molecular weight from SMILES, or None when invalid."""
    if not HAS_RDKIT:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        return Descriptors.ExactMolWt(mol)
    except Exception:
        return None


def adduct_neutral_mass_error(predicted_mass: float, parent_mz: float, ion_mode: str = "negative") -> float:
    """Delegate to shared adduct-aware neutral-mass error in Da."""
    return shared_adduct_neutral_mass_error(predicted_mass, parent_mz, ion_mode=ion_mode)


def passes_strict_mass_filter(predicted_mass: float, parent_mz: float, ion_mode: str = "negative") -> bool:
    """max(10 ppm, 0.005 Da) mode-aware tolerance."""
    err = adduct_neutral_mass_error(predicted_mass, parent_mz, ion_mode=ion_mode)
    return err <= shared_strict_tolerance_da(parent_mz, ion_mode=ion_mode)


_MASS_SCORE_PLATEAU_DA: float = 0.001  # full mass score below this Da error


def mass_match_score(
    predicted_mass: float | None,
    parent_mz: float,
    tolerance_da: float = 0.005,
    ion_mode: str = "negative",
) -> float:
    """Score [0, 1] from mass error vs adduct-corrected neutral mass.

    Uses ppm relative to the reference neutral mass with the same Da windows as
    before: full score below ``_MASS_SCORE_PLATEAU_DA``, zero at or above
    ``tolerance_da``, and ``1 - err_ppm / tol_ppm`` in between (equivalent to
    ``1 - err_da / tolerance_da`` since ``err_ppm / tol_ppm == err_da / tolerance_da``).
    """
    if predicted_mass is None:
        return 0.0
    err_da = adduct_neutral_mass_error(predicted_mass, parent_mz, ion_mode=ion_mode)
    ref = max(
        abs(float(neutral_mass_from_precursor(parent_mz, ion_mode=ion_mode))),
        1e-12,
    )
    err_ppm = err_da / ref * 1e6
    plateau_ppm = _MASS_SCORE_PLATEAU_DA / ref * 1e6
    tol_ppm = float(tolerance_da) / ref * 1e6
    if err_ppm <= plateau_ppm:
        return 1.0
    if err_ppm >= tol_ppm:
        return 0.0
    return max(0.0, 1.0 - err_ppm / tol_ppm)


def tanimoto_similarity(smi_a: str, smi_b: str, radius: int = 2) -> float:
    """Morgan fingerprint Tanimoto between two SMILES. Returns 0 on any failure."""
    if not HAS_RDKIT:
        return 0.0
    try:
        mol_a = Chem.MolFromSmiles(smi_a)
        mol_b = Chem.MolFromSmiles(smi_b)
        if mol_a is None or mol_b is None:
            return 0.0
        fp_a = AllChem.GetMorganFingerprintAsBitVect(mol_a, radius, nBits=2048)
        fp_b = AllChem.GetMorganFingerprintAsBitVect(mol_b, radius, nBits=2048)
        return float(TanimotoSimilarity(fp_a, fp_b))
    except Exception:
        return 0.0


def composite_confidence(
    generation_support: float,
    is_valid: bool,
    mass_score: float,
    mean_log_prob: float,
    latent_cosine: float,
    retrieval_tanimoto: float = 0.0,
) -> float:
    """Combine signals into a single [0, 1] confidence.
    Weights are heuristic and can be tuned."""
    validity = 1.0 if is_valid else 0.0
    norm_lp = 1.0 / (1.0 + np.exp(-mean_log_prob - 1.0))  # sigmoid centering around -1
    cos_term = max(0.0, latent_cosine) if not np.isnan(latent_cosine) else 0.0
    retrieval_bonus = min(1.0, retrieval_tanimoto / 0.7) if retrieval_tanimoto > 0.3 else 0.0

    score = (
        0.15 * generation_support
        + 0.10 * validity
        + 0.25 * mass_score
        + 0.20 * norm_lp
        + 0.15 * cos_term
        + 0.15 * retrieval_bonus
    )
    return float(np.clip(score, 0.0, 1.0))


# ── encoder ─────────────────────────────────────────────────────────────────────

def load_encoder_checkpoint(
    ckpt_path: Path,
    n_bins: int,
    minmass: float,
    maxmass: float,
    parent_mz_min: float,
    parent_mz_max: float,
    device: torch.device,
    apply_tanh: bool = True,
) -> tuple[Net1D, int, dict[str, object]]:
    """Load mass-conditioned ``Net1D`` encoder with shape-compatible bin width."""
    state = robust_torch_load(ckpt_path, map_location="cpu", weights_only=False)
    if not isinstance(state, dict):
        raise RuntimeError(f"Expected state_dict in {ckpt_path}, got {type(state)}")
    model_n_bins = int(infer_mass_conditioned_compatible_length(state, preferred_length_in=n_bins) or n_bins)
    model = Net1D(
        length_in=model_n_bins,
        minmass=minmass,
        maxmass=maxmass,
        parent_mz_min=parent_mz_min,
        parent_mz_max=parent_mz_max,
        apply_tanh=apply_tanh,
    ).to(device)
    load_info = load_compatible_state_dict(model, state)
    model.eval()
    return model, model_n_bins, load_info


def align_vector_length(vec: np.ndarray, target_n_bins: int) -> np.ndarray:
    """Truncate or zero-pad a binned spectrum to ``target_n_bins``."""
    current_n_bins = int(vec.shape[0])
    if current_n_bins == int(target_n_bins):
        return vec
    if current_n_bins > int(target_n_bins):
        return vec[: int(target_n_bins)]
    out = np.zeros(int(target_n_bins), dtype=vec.dtype)
    out[:current_n_bins] = vec
    return out


def encode_features(
    features_df: pd.DataFrame,
    encoder: Net1D,
    minmass: float,
    maxmass: float,
    resolution: int,
    device: torch.device,
    preprocessing_config: SpectrumPreprocessingConfig,
    encoder_n_bins: int | None = None,
) -> tuple[list[dict], torch.Tensor, dict[str, torch.Tensor], dict[str, list[tuple[float, float]]]]:
    """Encode grouped feature peaks to latents and per-feature metadata rows."""
    feature_rows = []
    latents = []
    latent_map: dict[str, torch.Tensor] = {}
    processed_peak_map: dict[str, list[tuple[float, float]]] = {}

    grouped = features_df.groupby("feature_key", sort=True)
    for feature_key, sub in grouped:
        peaks = list(zip(sub["mass"].tolist(), sub["intensity"].tolist()))
        parent_mz = float(sub["parent_mz"].iloc[0])
        comp_id = int(sub["comp_id"].iloc[0])
        compound_name = str(sub["compound_name"].iloc[0])
        rt = float(sub["rt"].iloc[0])
        source_html = str(sub["source_html"].iloc[0])

        processed = preprocess_and_vectorize(
            peaks=peaks,
            config=preprocessing_config,
            minmass=minmass,
            maxmass=maxmass,
            resolution=resolution,
        )
        vec = processed["vector"]
        if encoder_n_bins is not None:
            vec = align_vector_length(vec, encoder_n_bins)
        x = torch.from_numpy(vec).unsqueeze(0).unsqueeze(0).to(device)
        pmz = torch.tensor([parent_mz], dtype=torch.float32, device=device)

        with torch.no_grad():
            z = encoder(x, parent_mz=pmz).squeeze(0).detach().cpu()

        feature_rows.append(
            {
                "feature_key": feature_key,
                "comp_id": comp_id,
                "feature_label": compound_name,
                "rt_or_ri": rt,
                "parent_mz": parent_mz,
                "source_html": source_html,
                "n_peaks": int(processed["n_raw_peaks"]),
                "n_raw_peaks": int(processed["n_raw_peaks"]),
                "n_retained_peaks": int(processed["n_retained_peaks"]),
            }
        )
        latents.append(z)
        latent_map[feature_key] = z
        processed_peak_map[feature_key] = list(processed["peaks"])

    latent_tensor = torch.stack(latents, dim=0) if latents else torch.empty((0, 512))
    return feature_rows, latent_tensor, latent_map, processed_peak_map


# ── decoder artifacts ─────────────────────────────────────────────────────────

def derive_model_specific_artifact_name(model_name: str, prefix: str, extension: str) -> str:
    """Build ``{prefix}_{model_stem}{extension}`` from a decoder checkpoint name."""
    stem = Path(model_name).stem
    suffix = stem[len("model_") :] if stem.startswith("model_") else stem
    return f"{prefix}_{suffix}{extension}"


def resolve_decoder_artifact_path(
    decoder_dir: Path,
    explicit_name: str,
    generic_name: str,
    model_name: str,
    prefix: str,
) -> Path:
    """Pick model-specific or generic decoder artifact path in ``decoder_dir``."""
    if explicit_name:
        return (decoder_dir / explicit_name).resolve()
    specific = decoder_dir / derive_model_specific_artifact_name(model_name, prefix, ".nb")
    if specific.exists():
        return specific.resolve()
    return (decoder_dir / generic_name).resolve()


def derive_decoder_summary_path(decoder_dir: Path, model_name: str, explicit_path: str) -> Path:
    """Resolve fine-tune summary JSON beside decoder checkpoints."""
    if explicit_path:
        return Path(explicit_path).expanduser().resolve()
    return (decoder_dir / derive_model_specific_artifact_name(model_name, "finetune", "_summary.json")).resolve()


def maybe_load_decoder_summary_payload(summary_path: Path | None) -> dict[str, object] | None:
    """Load decoder summary JSON when the path exists and is a dict."""
    if summary_path is None or not summary_path.exists():
        return None
    loaded = json.loads(summary_path.read_text())
    if not isinstance(loaded, dict):
        return None
    return loaded


def resolve_decoder_mass_normalization(
    decoder_summary_payload: dict[str, object] | None,
) -> tuple[float | None, float | None, str]:
    """Read mass norm bounds from decoder summary JSON when available."""
    if not isinstance(decoder_summary_payload, dict):
        return None, None, "runtime_args"
    mass_norm = decoder_summary_payload.get("mass_normalization")
    if not isinstance(mass_norm, dict):
        return None, None, "runtime_args"
    min_value = mass_norm.get("min")
    max_value = mass_norm.get("max")
    if min_value is None or max_value is None:
        return None, None, "runtime_args"
    return float(min_value), float(max_value), "decoder_summary"


def derive_decoder_diagnosis_output_dir(output_dir: Path, explicit_path: str) -> Path:
    """Default sibling directory for post-run decoder mass diagnosis."""
    if explicit_path:
        return Path(explicit_path).expanduser().resolve()
    return (output_dir.parent / f"{output_dir.name}_decoder_mass_diagnosis").resolve()


def load_decoder(
    decoder_dir: Path,
    device: torch.device,
    model_name: str = "model_mass_209.pt",
    config_name: str = "",
    vocab_name: str = "",
):
    """Load constrained de novo decoder (vocab, config, weights) onto device."""
    register_moses_aliases()
    from decoder.moses.models_storage import ModelsStorage

    vocab_path = resolve_decoder_artifact_path(
        decoder_dir=decoder_dir,
        explicit_name=vocab_name,
        generic_name="vocab.nb",
        model_name=model_name,
        prefix="vocab",
    )
    config_path = resolve_decoder_artifact_path(
        decoder_dir=decoder_dir,
        explicit_name=config_name,
        generic_name="config.nb",
        model_name=model_name,
        prefix="config",
    )
    vocab = robust_torch_load(vocab_path, map_location="cpu", weights_only=False)
    config = robust_torch_load(config_path, map_location="cpu", weights_only=False)
    preferred = decoder_dir / model_name
    fallback = decoder_dir / "model.pt"
    state_path = preferred if preferred.exists() else fallback
    if not state_path.exists():
        raise FileNotFoundError(f"Decoder weights not found. Tried: {preferred} and {fallback}")
    print(f"Loading decoder vocab: {vocab_path}")
    print(f"Loading decoder config: {config_path}")
    print(f"Loading decoder weights: {state_path}")
    state = robust_torch_load(state_path, map_location="cpu", weights_only=False)

    models = ModelsStorage()
    decoder = models.get_model_class("translation")(vocab, config)
    # Some decoder checkpoints were trained with a smaller conditioning MLP
    # hidden size than the current model definition. Filter shape-mismatched
    # tensors so core encoder/decoder weights still load.
    model_state = decoder.state_dict()
    filtered_state = {}
    skipped_shape: list[tuple[str, tuple[int, ...], tuple[int, ...]]] = []
    for key, value in state.items():
        target = model_state.get(key)
        if target is None:
            filtered_state[key] = value
            continue
        if torch.is_tensor(value) and torch.is_tensor(target):
            if tuple(value.shape) != tuple(target.shape):
                skipped_shape.append((key, tuple(value.shape), tuple(target.shape)))
                continue
        filtered_state[key] = value

    incompatible = decoder.load_state_dict(filtered_state, strict=False)
    missing = list(getattr(incompatible, "missing_keys", []))
    unexpected = list(getattr(incompatible, "unexpected_keys", []))
    if skipped_shape:
        preview = ", ".join(
            f"{name} ckpt={src_shape} model={dst_shape}" for name, src_shape, dst_shape in skipped_shape[:8]
        )
        print(f"Decoder load skipped {len(skipped_shape)} shape-mismatched tensor(s): {preview}")
    if missing:
        print(f"Decoder load missing keys: {missing}")
    if unexpected:
        print(f"Decoder load unexpected keys: {unexpected}")
    decoder.to(device)
    decoder.eval()
    return decoder, vocab, state_path.resolve(), config_path, vocab_path


def run_decoder_mass_diagnosis(
    python_executable: str,
    repo_root: Path,
    prediction_dir: Path,
    decoder_dir: Path,
    decoder_model_name: str,
    decoder_summary_json: Path,
    diagnosis_output_dir: Path,
    enforce_gate: bool,
    gate_allowed_root_causes: str,
    gate_min_bank_targets_with_strict_pass: int,
    gate_min_raw_features_within_1da: int,
    gate_fail_on_ri_mismatch: bool,
    ion_mode: str,
) -> dict:
    """Invoke ``analyze_decoder_mass_conditioning.py`` and return gate metadata."""
    diagnosis_script = repo_root / "analyze_decoder_mass_conditioning.py"
    cmd = [
        python_executable,
        str(diagnosis_script),
        "--prediction-dir",
        str(prediction_dir),
        "--decoder-dir",
        str(decoder_dir),
        "--decoder-model-name",
        decoder_model_name,
        "--decoder-summary-json",
        str(decoder_summary_json),
        "--output-dir",
        str(diagnosis_output_dir),
        "--gate-allowed-root-causes",
        str(gate_allowed_root_causes),
        "--gate-min-bank-targets-with-strict-pass",
        str(gate_min_bank_targets_with_strict_pass),
        "--gate-min-raw-features-within-1da",
        str(gate_min_raw_features_within_1da),
        "--ion-mode",
        str(ion_mode),
    ]
    if enforce_gate:
        cmd.append("--enforce-gate")
    if gate_fail_on_ri_mismatch:
        cmd.append("--gate-fail-on-ri-mismatch")

    print("\n=== Post-run decoder mass diagnosis ===", flush=True)
    print("Running:", " ".join(cmd), flush=True)
    completed = subprocess.run(cmd, check=False)

    gate_status_path = diagnosis_output_dir / "gate_status.json"
    gate_status = None
    if gate_status_path.exists():
        gate_status = json.loads(gate_status_path.read_text())

    return {
        "requested": True,
        "enforced": bool(enforce_gate),
        "return_code": int(completed.returncode),
        "decoder_summary_json": str(decoder_summary_json),
        "diagnosis_output_dir": str(diagnosis_output_dir),
        "gate_status_json": str(gate_status_path),
        "gate_pass": gate_status.get("pass") if isinstance(gate_status, dict) else None,
        "ion_mode": str(ion_mode),
    }


# ── Track B sampling ────────────────────────────────────────────────────────────

def sample_generation_support(
    decoder,
    feature_keys: list[str],
    latent_tensor: torch.Tensor,
    mass_conditions: list[float] | None,
    ri_conditions: list[float] | None,
    num_samples: int,
    temp: float,
    batch_size: int,
    max_len: int,
    top_p: float = 0.0,
    temps: list[float] | None = None,
) -> tuple[dict[str, Counter], dict[str, int], pd.DataFrame, dict[str, dict[str, float]]]:
    """Batch decoder sampling; aggregate per-SMILES counts and log-probs."""
    counts: dict[str, Counter] = {k: Counter() for k in feature_keys}
    draws = {k: 0 for k in feature_keys}
    log_prob_accum: dict[str, dict[str, list[float]]] = defaultdict(lambda: defaultdict(list))
    raw_rows = []
    n = len(feature_keys)

    if n == 0:
        return (
            counts,
            draws,
            pd.DataFrame(
                columns=[
                    "sample_idx",
                    "feature_key",
                    "smiles",
                    "mean_log_prob",
                    "temp",
                    "mass_condition",
                    "ri_condition",
                ],
            ),
            {},
        )

    z_all = latent_tensor.to(decoder.device)
    mass_all = None
    ri_all = None
    if mass_conditions is not None:
        if len(mass_conditions) != n:
            raise ValueError(f"mass_conditions length {len(mass_conditions)} does not match n_features={n}")
        mass_all = torch.tensor(mass_conditions, dtype=torch.float32, device=decoder.device).unsqueeze(1)
    if ri_conditions is not None:
        if len(ri_conditions) != n:
            raise ValueError(f"ri_conditions length {len(ri_conditions)} does not match n_features={n}")
        ri_all = torch.tensor(ri_conditions, dtype=torch.float32, device=decoder.device).unsqueeze(1)

    if temps and len(temps) > 1:
        per_temp = max(1, num_samples // len(temps))
        schedule = []
        for t in temps:
            schedule.extend([t] * per_temp)
        while len(schedule) < num_samples:
            schedule.append(temps[-1])
        schedule = schedule[:num_samples]
    else:
        schedule = [temp] * num_samples

    for sample_idx, sample_temp in enumerate(schedule):
        for start in range(0, n, batch_size):
            end = min(start + batch_size, n)
            z_batch = z_all[start:end]
            mass_batch = mass_all[start:end] if mass_all is not None else None
            ri_batch = ri_all[start:end] if ri_all is not None else None
            result = decoder.sample(
                end - start,
                z=z_batch,
                temp=sample_temp,
                max_len=max_len,
                top_p=top_p,
                return_log_probs=True,
                mass=mass_batch,
                ri=ri_batch,
            )
            decoded, mean_lps = result
            for local_i, smiles in enumerate(decoded):
                feature_key = feature_keys[start + local_i]
                counts[feature_key][smiles] += 1
                draws[feature_key] += 1
                log_prob_accum[feature_key][smiles].append(mean_lps[local_i])
                raw_rows.append(
                    {
                        "sample_idx": sample_idx,
                        "feature_key": feature_key,
                        "smiles": smiles,
                        "mean_log_prob": mean_lps[local_i],
                        "temp": sample_temp,
                        "mass_condition": float(mass_batch[local_i].item()) if mass_batch is not None else float("nan"),
                        "ri_condition": float(ri_batch[local_i].item()) if ri_batch is not None else float("nan"),
                    }
                )

    raw_df = pd.DataFrame(raw_rows)
    best_log_probs: dict[str, dict[str, float]] = {}
    for fk, smi_dict in log_prob_accum.items():
        best_log_probs[fk] = {smi: max(lps) for smi, lps in smi_dict.items()}

    return counts, draws, raw_df, best_log_probs


# ── embedding bank / retrieval tables ───────────────────────────────────────────

def _extract_embedding_tensor(obj, source: Path) -> torch.Tensor:
    """Normalize per-comp_id ``.pt`` payloads to a 1-D float tensor."""
    if torch.is_tensor(obj):
        z = obj
    elif isinstance(obj, (tuple, list)):
        z = next((v for v in obj if torch.is_tensor(v)), None)
        if z is None:
            raise ValueError(f"No tensor found in tuple/list embedding: {source}")
    elif isinstance(obj, dict):
        z = obj.get("z", None) or obj.get("latent", None) or obj.get("embedding", None)
        if z is None or not torch.is_tensor(z):
            raise ValueError(f"No tensor z/latent/embedding in dict embedding: {source}")
    else:
        raise ValueError(f"Unsupported embedding object type {type(obj)} in {source}")

    z = z.detach().cpu().float()
    if z.ndim == 2 and z.shape[0] == 1:
        z = z.squeeze(0)
    if z.ndim != 1:
        raise ValueError(f"Expected 1D embedding for {source}, got shape={tuple(z.shape)}")
    return z


def load_embedding_bank(emb_dir: Path) -> tuple[list[int], torch.Tensor]:
    """Load comp_id-keyed latent bank from ``emb_dir/*.pt`` files."""
    bank_ids = []
    bank_rows = []
    for emb_file in sorted(emb_dir.glob("*.pt")):
        try:
            emb_id = int(emb_file.stem)
        except ValueError:
            continue
        obj = robust_torch_load(emb_file, map_location="cpu", weights_only=False)
        z = _extract_embedding_tensor(obj, emb_file)
        bank_ids.append(emb_id)
        bank_rows.append(z)
    if not bank_rows:
        raise RuntimeError(f"No usable .pt embeddings found in {emb_dir}")
    return bank_ids, torch.stack(bank_rows, dim=0)


def _pick_col(df: pd.DataFrame, options: list[str]) -> str | None:
    """Case-insensitive column lookup from an ordered candidate list."""
    lookup = {c.lower(): c for c in df.columns}
    for candidate in options:
        col = lookup.get(candidate.lower())
        if col is not None:
            return col
    return None


def load_retrieval_smiles(path: Path) -> dict[int, str]:
    """Load comp_id -> SMILES from the canonical retrieval table (comp_id_to_smiles_209.csv)."""
    df = read_table(path)
    id_col = _pick_col(df, ["comp_id", "COMP_ID", "id"])
    smiles_col = _pick_col(df, ["smiles", "SMILES", "canonical_smiles"])
    if id_col is None:
        raise KeyError(f"No comp_id column in {path}; found: {df.columns.tolist()}")
    if smiles_col is None:
        raise KeyError(f"No SMILES column in {path}; found: {df.columns.tolist()}")
    out: dict[int, str] = {}
    for _, row in df.iterrows():
        try:
            cid = int(row[id_col])
        except Exception:
            continue
        smi = str(row[smiles_col]).strip()
        if smi and smi.lower() != "nan" and cid not in out:
            out[cid] = smi
    return out


def load_retrieval_inchikey14(path: Path) -> dict[int, str]:
    """Load comp_id -> InChIKey first block from retrieval SMILES table."""
    if not path.exists():
        return {}
    df = read_table(path)
    id_col = _pick_col(df, ["comp_id", "COMP_ID", "id"])
    ik14_col = _pick_col(df, ["inchikey14", "InChIKey_smiles_firstBlock", "InChIKey14"])
    smiles_col = _pick_col(df, ["smiles", "SMILES", "canonical_smiles"])
    if id_col is None:
        return {}
    out: dict[int, str] = {}
    for _, row in df.iterrows():
        try:
            cid = int(row[id_col])
        except Exception:
            continue
        ik14 = str(row[ik14_col]).strip() if ik14_col is not None else ""
        if not ik14 and HAS_RDKIT and smiles_col is not None:
            try:
                mol = Chem.MolFromSmiles(str(row[smiles_col]).strip())
                if mol is not None:
                    ik14 = str(Chem.MolToInchiKey(mol))[:14]
            except Exception:
                ik14 = ""
        if ik14 and cid not in out:
            out[cid] = ik14
    return out


def load_prevalence_map(path: Path) -> dict[str, float]:
    """Load inchikey14 -> prevalence score from FASSTrecords prevalence CSV."""
    if not path.exists():
        return {}
    df = read_table(path)
    key_col = _pick_col(df, ["inchikey14", "InChIKey_smiles_firstBlock"])
    score_col = _pick_col(df, ["prevalence_score"])
    if key_col is None or score_col is None:
        return {}
    out: dict[str, float] = {}
    for _, row in df.iterrows():
        key = str(row[key_col]).strip()
        if not key:
            continue
        try:
            score = float(row[score_col])
        except Exception:
            score = 0.0
        if not np.isfinite(score):
            score = 0.0
        out[key] = float(np.clip(score, 0.0, 1.0))
    return out


def load_retrieval_mass_hints(path: Path, ion_mode: str = "negative") -> dict[int, float]:
    """Load optional comp_id -> parent_mz hints from retrieval SMILES CSV."""
    if not path.exists():
        return {}
    df = read_table(path)
    id_col = _pick_col(df, ["comp_id", "COMP_ID", "id"])
    if str(ion_mode).strip().lower() == "positive":
        mz_col = _pick_col(
            df,
            ["theoretical_parent_mz_pos", "parent_mz", "theoretical_parent_mz"],
        )
    else:
        mz_col = _pick_col(
            df,
            ["theoretical_parent_mz_neg", "parent_mz", "theoretical_parent_mz"],
        )
    exact_mass_col = _pick_col(df, ["exact_mass", "exact_molecular_weight", "exact_mol_weight"])
    if id_col is None:
        return {}

    out: dict[int, float] = {}
    for _, row in df.iterrows():
        try:
            cid = int(row[id_col])
        except Exception:
            continue
        if mz_col is not None:
            try:
                mz_val = float(row[mz_col])
                if np.isfinite(mz_val):
                    out[cid] = mz_val
                    continue
            except Exception:
                pass
        if exact_mass_col is not None:
            try:
                exact_mass = float(row[exact_mass_col])
            except Exception:
                continue
            if not np.isfinite(exact_mass):
                continue
            if str(ion_mode).strip().lower() == "negative":
                out[cid] = float(exact_mass - PROTON_MASS)
            else:
                out[cid] = float(exact_mass + PROTON_MASS)
    return out


def load_retrieval_names(path: Path) -> dict[int, str]:
    """Load comp_id -> display name from name_dict.csv."""
    df = read_table(path)
    id_col = _pick_col(df, ["comp_id", "COMP_ID", "id"])
    name_col = _pick_col(df, ["name", "SHORT_NAME", "CHEMICAL_NAME", "COMPOUND_NAME"])
    if id_col is None or name_col is None:
        return {}
    out: dict[int, str] = {}
    for _, row in df.iterrows():
        try:
            cid = int(row[id_col])
        except Exception:
            continue
        name = str(row[name_col]).strip()
        if name and name.lower() != "nan" and cid not in out:
            out[cid] = name
    return out


def load_retrieval_mz(path: Path, fallback_hints: dict[int, float] | None = None) -> dict[int, float]:
    """Load comp_id -> precursor m/z from pickle, with optional fallback hints."""
    import pickle

    fallback_hints = fallback_hints or {}
    if not path.exists():
        return {int(k): float(v) for k, v in fallback_hints.items()}

    with open(path, "rb") as fh:
        d = pickle.load(fh)
    out = {int(k): float(v) for k, v in d.items()}
    for cid, mz in fallback_hints.items():
        out.setdefault(int(cid), float(mz))
    return out


def load_chemical_map(path: Path) -> tuple[dict[int, str], dict[int, str]]:
    """Legacy loader kept for optional enrichment from MASTER_METABOLITE_LIST.csv."""
    df = read_table(path)
    id_col = _pick_col(df, ["CHEMICAL_ID", "chemical_id", "id"])
    smiles_col = _pick_col(df, ["SMILES", "smiles", "canonical_smiles"])
    name_col = _pick_col(df, ["SHORT_NAME", "CHEMICAL_NAME", "COMPOUND_NAME", "name"])

    if id_col is None or smiles_col is None:
        return {}, {}

    id_to_smiles: dict[int, str] = {}
    id_to_name: dict[int, str] = {}
    for _, row in df.iterrows():
        try:
            chem_id = int(row[id_col])
        except Exception:
            continue
        smiles = str(row[smiles_col]).strip()
        if smiles and smiles.lower() != "nan" and chem_id not in id_to_smiles:
            id_to_smiles[chem_id] = smiles

        if name_col is not None:
            name = str(row[name_col]).strip()
            if name and name.lower() != "nan" and chem_id not in id_to_name:
                id_to_name[chem_id] = name

    return id_to_smiles, id_to_name


# ── HTML / isotope audit ────────────────────────────────────────────────────────

def serialize_list_columns_for_csv(df: pd.DataFrame, columns: Sequence[str]) -> pd.DataFrame:
    """JSON-encode list/dict columns so CSV exports round-trip cleanly."""
    out = df.copy()
    for column in columns:
        if column in out.columns:
            out[column] = out[column].apply(
                lambda value: json.dumps(value) if isinstance(value, (list, dict, set)) else value
            )
    return out


def build_html_signal_feature_tables(
    feature_meta_df: pd.DataFrame,
    candidate_formula_max: int,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, dict[str, object]]]:
    """Parse HTML signals into formula and isotope feature tables per feature_key."""
    formula_rows: list[dict[str, object]] = []
    isotope_rows: list[dict[str, object]] = []
    signal_lookup: dict[str, dict[str, object]] = {}
    keep_cols = ["feature_key", "feature_label", "source_html"]
    unique_features = feature_meta_df[keep_cols].drop_duplicates(subset=["feature_key"]).copy()
    for row in unique_features.itertuples(index=False):
        try:
            record = parse_html_signals(
                source_html=row.source_html,
                feature_key=str(row.feature_key),
                feature_label=str(row.feature_label),
                candidate_formula_max=int(candidate_formula_max),
            )
        except Exception:
            record = {
                "feature_key": str(row.feature_key),
                "feature_label": str(row.feature_label),
                "source_html": str(row.source_html),
                "usable_formula_source": "error",
                "usable_formula_set": [],
                "feature_level_formula_set": [],
                "candidate_level_formula_set": [],
                "candidate_formula_unique_count": 0,
                "feature_is_likely_isotopologue": False,
                "feature_isotopologue_reason": None,
                "isotope_element_hints": {},
                "isotope_signal_count": 0,
                "observed_isotope_rows_json": [],
                "c_predict": None,
                "parse_status": "error",
                "observed_isotope_rows": [],
            }
        signal_lookup[str(record["feature_key"])] = dict(record)
        formula_rows.append(
            {
                "feature_key": record["feature_key"],
                "feature_label": record["feature_label"],
                "source_html": record["source_html"],
                "usable_formula_source": record["usable_formula_source"],
                "usable_formula_set": record["usable_formula_set"],
                "feature_level_formula_set": record.get("feature_level_formula_set", []),
                "candidate_level_formula_set": record.get("candidate_level_formula_set", []),
                "candidate_formula_unique_count": int(record.get("candidate_formula_unique_count", 0)),
            }
        )
        isotope_rows.append(
            {
                "feature_key": record["feature_key"],
                "source_html": record["source_html"],
                "feature_is_likely_isotopologue": bool(record.get("feature_is_likely_isotopologue", False)),
                "feature_isotopologue_reason": record.get("feature_isotopologue_reason"),
                "isotope_element_hints": dict(record.get("isotope_element_hints") or {}),
                "isotope_signal_count": int(record.get("isotope_signal_count", 0) or 0),
                "observed_isotope_rows_json": record.get("observed_isotope_rows_json", []),
                "c_predict": record.get("c_predict"),
                "parse_status": record.get("parse_status", "unavailable"),
            }
        )
    html_formula_df = pd.DataFrame(formula_rows)
    if not html_formula_df.empty:
        html_formula_df["usable_formula_count"] = html_formula_df["usable_formula_set"].apply(len)
        html_formula_df = html_formula_df.sort_values("feature_key").reset_index(drop=True)
    html_isotope_df = pd.DataFrame(isotope_rows)
    if not html_isotope_df.empty:
        html_isotope_df = html_isotope_df.sort_values("feature_key").reset_index(drop=True)
    return html_formula_df, html_isotope_df, signal_lookup


def annotate_dataframe_with_isotope_audit(
    df: pd.DataFrame,
    *,
    feature_signal_lookup: dict[str, dict[str, object]],
    formula_column: str,
    smiles_column: str | None = None,
) -> pd.DataFrame:
    """Attach isotope audit columns from HTML signal lookup to a dataframe."""
    if df.empty:
        return df.copy()

    out = df.copy()
    audit_rows: list[dict[str, object]] = []
    for row in out.itertuples(index=False):
        row_dict = row._asdict()
        feature_signal = feature_signal_lookup.get(str(row_dict.get("feature_key", "")), {})
        formula_value = row_dict.get(formula_column)
        if (formula_value is None or str(formula_value).strip() in {"", "nan", "None"}) and smiles_column is not None:
            smiles_value = row_dict.get(smiles_column)
            if isinstance(smiles_value, str) and smiles_value.strip():
                formula_value = molecular_formula(smiles_value)
        audit_rows.append(
            audit_candidate_formula_against_html_signals(
                formula=formula_value,
                feature_signal=feature_signal,
            )
        )
    audit_df = pd.DataFrame(audit_rows, index=out.index)
    for column in audit_df.columns:
        out[column] = audit_df[column]
    return out


def build_isotope_audit_summary(
    *,
    artifact_prefix: str,
    html_isotope_df: pd.DataFrame,
    retrieval_df: pd.DataFrame,
    denovo_df: pd.DataFrame,
    merged_df: pd.DataFrame,
    html_isotope_csv_path: Path,
) -> dict[str, object]:
    """Summarize isotope audit status counts across retrieval/denovo/merged tables."""
    parse_status_counts = (
        html_isotope_df.get("parse_status", pd.Series(dtype=object))
        .fillna("unavailable")
        .astype(str)
        .value_counts()
        .to_dict()
        if not html_isotope_df.empty
        else {}
    )
    feature_status_counts = (
        html_isotope_df.get("feature_is_likely_isotopologue", pd.Series(dtype=bool))
        .fillna(False)
        .astype(bool)
        .value_counts()
        .to_dict()
        if not html_isotope_df.empty
        else {}
    )
    retrieval_status_counts = (
        retrieval_df.get("isotope_audit_status", pd.Series(dtype=object))
        .fillna("unavailable")
        .astype(str)
        .value_counts()
        .to_dict()
    )
    denovo_status_counts = (
        denovo_df.get("isotope_audit_status", pd.Series(dtype=object))
        .fillna("unavailable")
        .astype(str)
        .value_counts()
        .to_dict()
    )
    merged_status_counts = (
        merged_df.get("isotope_audit_status", pd.Series(dtype=object))
        .fillna("unavailable")
        .astype(str)
        .value_counts()
        .to_dict()
    )
    return {
        "artifact_prefix": artifact_prefix,
        "feature_rows": int(len(html_isotope_df)),
        "feature_parse_status_counts": {str(key): int(value) for key, value in parse_status_counts.items()},
        "feature_isotopologue_counts": {str(key): int(value) for key, value in feature_status_counts.items()},
        "feature_rows_with_signals": (
            int(pd.to_numeric(html_isotope_df.get("isotope_signal_count"), errors="coerce").fillna(0).gt(0).sum())
            if not html_isotope_df.empty
            else 0
        ),
        "retrieval_audit_status_counts": {str(key): int(value) for key, value in retrieval_status_counts.items()},
        "denovo_audit_status_counts": {str(key): int(value) for key, value in denovo_status_counts.items()},
        "merged_audit_status_counts": {str(key): int(value) for key, value in merged_status_counts.items()},
        "html_isotope_signals_csv": str(html_isotope_csv_path),
    }


# ── fragment DB parsing ─────────────────────────────────────────────────────

SPECTRUM_PAIR_RE: re.Pattern[str] = re.compile(r"\(\s*([-+0-9.eE]+)\s*,\s*([-+0-9.eE]+)\s*\)")
FRAGMENT_BLOCK_RE: re.Pattern[str] = re.compile(r'^\s*"(?P<fragments>.*)"\s*,\s*"(?P<spectrum>\[.*\])"\s*$')


def parse_spectrum_literal(text: str) -> list[tuple[float, float]]:
    """Parse ``(mz, intensity)`` pairs from a spectrum literal string."""
    peaks: list[tuple[float, float]] = []
    for mz_text, intensity_text in SPECTRUM_PAIR_RE.findall(str(text)):
        try:
            peaks.append((float(mz_text), float(intensity_text)))
        except Exception:
            continue
    return peaks


def load_fragment_database(path: Path) -> dict[int, dict[str, object]]:
    """Load comp_id-keyed fragment spectra from library CSV exports."""
    if not path.exists():
        return {}

    records: dict[int, dict[str, object]] = {}
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        header = next(handle, None)
        if header is None:
            return records
        for line in handle:
            raw = line.rstrip("\n")
            if not raw.strip():
                continue
            prefix, sep, tail = raw.rpartition(', "')
            if not sep:
                continue
            prefix_parts = prefix.split(",", 8)
            if len(prefix_parts) < 9:
                continue
            spectrum_text = f'"{tail.strip()}' if not tail.lstrip().startswith('"') else tail.strip()
            fragments_and_meta = prefix_parts[8]
            match = FRAGMENT_BLOCK_RE.match(f"{fragments_and_meta},{spectrum_text}")
            if match is None:
                continue
            try:
                cid = int(prefix_parts[0].strip())
            except Exception:
                continue
            formula = str(prefix_parts[5]).strip()
            canonical_smiles = str(prefix_parts[7]).strip()
            spectrum = parse_spectrum_literal(match.group("spectrum"))
            records[cid] = {
                "candidate_formula": formula if formula and formula.lower() != "nan" else None,
                "canonical_smiles": (
                    canonical_smiles if canonical_smiles and canonical_smiles.lower() != "nan" else None
                ),
                "spectrum": spectrum,
            }
    return records


def compute_fragment_match_metrics(
    observed_peaks: list[tuple[float, float]],
    candidate_peaks: list[tuple[float, float]],
    tolerance_da: float = 0.02,
) -> dict[str, float | int]:
    """Score observed vs candidate fragment spectra within Da tolerance."""
    if not observed_peaks:
        return {
            "fragment_match_count": 0,
            "fragment_matched_peak_fraction": 0.0,
            "fragment_explained_intensity_fraction": 0.0,
        }
    if not candidate_peaks:
        return {
            "fragment_match_count": 0,
            "fragment_matched_peak_fraction": 0.0,
            "fragment_explained_intensity_fraction": 0.0,
        }

    candidate_mz = sorted(float(mz) for mz, _ in candidate_peaks)
    total_intensity = sum(float(intensity) for _, intensity in observed_peaks)
    matched_count = 0
    matched_intensity = 0.0
    for mz, intensity in observed_peaks:
        mz = float(mz)
        lo = mz - float(tolerance_da)
        hi = mz + float(tolerance_da)
        hit = False
        for candidate in candidate_mz:
            if candidate < lo:
                continue
            if candidate > hi:
                break
            hit = True
            break
        if hit:
            matched_count += 1
            matched_intensity += float(intensity)

    return {
        "fragment_match_count": int(matched_count),
        "fragment_matched_peak_fraction": float(matched_count / max(len(observed_peaks), 1)),
        "fragment_explained_intensity_fraction": float(matched_intensity / max(total_intensity, 1e-12)),
    }


MASS_GATES_DA: list[float] = [0.05, 0.5, 1.0]  # progressive retrieval mass gates

RETRIEVAL_SCORE_TIERS: list[tuple[float, str]] = [  # min score -> tier label
    (0.70, "high"),
    (0.55, "medium"),
    (0.40, "low"),
]


def precursor_mass_score(abs_mass_diff: float, gate_da: float) -> float:
    """Triangular [0, 1] score from absolute mass error within a gate window."""
    if not np.isfinite(abs_mass_diff) or gate_da <= 0:
        return 0.0
    if abs_mass_diff <= 0.001:
        return 1.0
    if abs_mass_diff >= gate_da:
        return 0.0
    return float(max(0.0, 1.0 - (abs_mass_diff / gate_da)))


def retrieval_confidence_tier(retrieval_score: float) -> str:
    """Map retrieval score to high/medium/low/reject tier labels."""
    for min_score, label in RETRIEVAL_SCORE_TIERS:
        if float(retrieval_score) >= min_score:
            return label
    return "reject"


# ── Track A retrieval ───────────────────────────────────────────────────────────

def compute_mass_gated_retrieval(
    feature_keys: list[str],
    feature_parent_mzs: list[float],
    feature_latents: torch.Tensor,
    feature_peak_map: dict[str, list[tuple[float, float]]],
    bank_ids: list[int],
    bank_latents: torch.Tensor,
    bank_mz: dict[int, float],
    id_to_smiles: dict[int, str],
    id_to_name: dict[int, str],
    id_to_inchikey14: dict[int, str] | None,
    fragment_db: dict[int, dict[str, object]] | None,
    html_formula_lookup: dict[str, set[str]] | None,
    retrieval_scoring_config: RetrievalScoringConfig | None,
    prevalence_map: dict[str, float] | None,
    topk: int = 10,
    fragment_match_tolerance_da: float = 0.02,
    fragment_match_target: float = 5.0,
    structural_neighbor_top_n: int = 5,
    structural_min_mcs_atoms: int = 6,
    structural_mcs_timeout_sec: int = 1,
) -> pd.DataFrame:
    """Mass-gated KNN retrieval with fragment/intensity/formula-aware reranking."""
    zq = F.normalize(feature_latents.float(), dim=1)
    zb = F.normalize(bank_latents.float(), dim=1)
    cos_matrix = zq @ zb.T  # (N_features, N_bank)
    scoring_cfg = retrieval_scoring_config or resolve_scoring_config(None)

    bank_mz_arr = np.array([bank_mz.get(bid, float("nan")) for bid in bank_ids])
    rows = []

    for i, feature_key in enumerate(feature_keys):
        parent_mz = feature_parent_mzs[i]
        cos_scores = cos_matrix[i].cpu().numpy()  # (N_bank,)
        abs_diffs = np.abs(bank_mz_arr - parent_mz)

        best_gate = None
        gate_mask = None
        for gate in MASS_GATES_DA:
            mask = abs_diffs <= gate
            if mask.any():
                best_gate = gate
                gate_mask = mask
                break

        if best_gate is None:
            rows.append(
                {
                    "feature_key": feature_key,
                    "rank": 1,
                    "comp_id": "",
                    "smiles": "",
                    "name": "",
                    "latent_cosine": float("nan"),
                    "gate_used": float("nan"),
                    "abs_mass_diff": float("nan"),
                    "precursor_mass_score": 0.0,
                    "candidate_formula": None,
                    "usable_formula_count": 0,
                    "html_formula_bonus": 0.0,
                    "fragment_match_count": 0,
                    "fragment_matched_peak_fraction": 0.0,
                    "fragment_explained_intensity_fraction": 0.0,
                    "fragment_match_count_score": 0.0,
                    "isotope_formula_score": 0.0,
                    "substructure_consensus_score": 0.0,
                    "tanimoto_neighbor_score": 0.0,
                    "prevalence_score": 0.0,
                    "retrieval_score": 0.0,
                    "confidence_tier": "no_gated_analog",
                    "recommended": False,
                    "status": "no_gated_analog",
                }
            )
            continue

        candidate_indices = np.where(gate_mask)[0]
        candidate_cos = cos_scores[candidate_indices]
        top_order = np.argsort(-candidate_cos)[: max(int(topk), 1)]

        scored_rows: list[dict[str, object]] = []
        rescoring_active = (fragment_db is not None) or bool(html_formula_lookup)
        usable_formula_set = set(html_formula_lookup.get(feature_key, set())) if html_formula_lookup else set()
        observed_peaks = feature_peak_map.get(feature_key, [])
        for oi in top_order:
            bi = int(candidate_indices[oi])
            cid = int(bank_ids[bi])
            cos_val = float(candidate_cos[oi])
            diff = float(abs_diffs[bi])
            candidate_record = (fragment_db or {}).get(cid, {})
            candidate_formula = candidate_record.get("candidate_formula") or molecular_formula(
                id_to_smiles.get(cid, "")
            )
            prevalence_score = 0.0
            if prevalence_map:
                ik14 = str((id_to_inchikey14 or {}).get(cid, "")).strip()
                prevalence_score = float(prevalence_map.get(ik14, 0.0))
            fragment_metrics = {
                "fragment_match_count": 0,
                "fragment_matched_peak_fraction": 0.0,
                "fragment_explained_intensity_fraction": 0.0,
            }
            precursor_score = 0.0
            html_formula_bonus = 0.0
            fragment_match_count_score = 0.0
            isotope_formula_score = 0.0
            substructure_consensus_score = 0.0
            tanimoto_neighbor_score = 0.0
            if rescoring_active:
                fragment_metrics = compute_fragment_match_metrics(
                    observed_peaks=observed_peaks,
                    candidate_peaks=candidate_record.get("spectrum", []),
                    tolerance_da=float(fragment_match_tolerance_da),
                )
                precursor_score = precursor_mass_score(diff, best_gate)
                html_formula_bonus = float(
                    1.0 if usable_formula_set and candidate_formula in usable_formula_set else 0.0
                )
                fragment_match_count_score = float(
                    np.clip(
                        float(fragment_metrics["fragment_match_count"]) / max(float(fragment_match_target), 1.0),
                        0.0,
                        1.0,
                    )
                )
                retrieval_score = compute_retrieval_score(
                    latent_cosine=cos_val,
                    precursor_mass_score=precursor_score,
                    fragment_explained_intensity_fraction=float(
                        fragment_metrics["fragment_explained_intensity_fraction"]
                    ),
                    fragment_matched_peak_fraction=float(fragment_metrics["fragment_matched_peak_fraction"]),
                    html_formula_bonus=html_formula_bonus,
                    isotope_formula_score=isotope_formula_score,
                    fragment_match_count_score=fragment_match_count_score,
                    substructure_consensus_score=substructure_consensus_score,
                    tanimoto_neighbor_score=tanimoto_neighbor_score,
                    prevalence_score=prevalence_score,
                    config=scoring_cfg,
                )
            else:
                retrieval_score = float(cos_val)
            scored_rows.append(
                {
                    "feature_key": feature_key,
                    "comp_id": cid,
                    "smiles": id_to_smiles.get(cid, ""),
                    "name": id_to_name.get(cid, ""),
                    "latent_cosine": cos_val,
                    "gate_used": best_gate,
                    "abs_mass_diff": diff,
                    "precursor_mass_score": precursor_score,
                    "candidate_formula": candidate_formula,
                    "usable_formula_count": int(len(usable_formula_set)),
                    "html_formula_bonus": html_formula_bonus,
                    "isotope_formula_score": isotope_formula_score,
                    "fragment_match_count_score": fragment_match_count_score,
                    "substructure_consensus_score": substructure_consensus_score,
                    "tanimoto_neighbor_score": tanimoto_neighbor_score,
                    "prevalence_score": prevalence_score,
                    "retrieval_score": float(retrieval_score),
                    **fragment_metrics,
                    "status": "ok",
                }
            )

        if rescoring_active and scored_rows:
            ordered_rows = sorted(
                scored_rows,
                key=lambda row: (
                    -float(row["retrieval_score"]),
                    -float(row["latent_cosine"]),
                    float(row["abs_mass_diff"]),
                    str(row["smiles"]),
                ),
            )
            structural_scores = compute_structural_neighbor_scores(
                [str(row["smiles"]) for row in ordered_rows],
                top_n=int(structural_neighbor_top_n),
                min_mcs_atoms=int(structural_min_mcs_atoms),
                mcs_timeout_sec=int(structural_mcs_timeout_sec),
            )
            for row, structural in zip(ordered_rows, structural_scores):
                row["substructure_consensus_score"] = float(structural["substructure_consensus_score"])
                row["tanimoto_neighbor_score"] = float(structural["tanimoto_neighbor_score"])
                row["retrieval_score"] = compute_retrieval_score(
                    latent_cosine=float(row.get("latent_cosine", 0.0)),
                    precursor_mass_score=float(row.get("precursor_mass_score", 0.0)),
                    fragment_explained_intensity_fraction=float(row.get("fragment_explained_intensity_fraction", 0.0)),
                    fragment_matched_peak_fraction=float(row.get("fragment_matched_peak_fraction", 0.0)),
                    html_formula_bonus=float(row.get("html_formula_bonus", 0.0)),
                    isotope_formula_score=float(row.get("isotope_formula_score", 0.0)),
                    fragment_match_count_score=float(row.get("fragment_match_count_score", 0.0)),
                    substructure_consensus_score=float(row.get("substructure_consensus_score", 0.0)),
                    tanimoto_neighbor_score=float(row.get("tanimoto_neighbor_score", 0.0)),
                    prevalence_score=float(row.get("prevalence_score", 0.0)),
                    config=scoring_cfg,
                )
            scored_rows = ordered_rows

        scored_rows.sort(
            key=lambda row: (
                -float(row["retrieval_score"]),
                -float(row["latent_cosine"]),
                float(row["abs_mass_diff"]),
                str(row["smiles"]),
            )
        )
        for rank, row in enumerate(scored_rows[:topk], start=1):
            row["rank"] = rank
            row["confidence_tier"] = retrieval_confidence_tier(float(row["retrieval_score"]))
            row["recommended"] = rank == 1
            rows.append(row)

    return pd.DataFrame(rows)


# ── scoring / artifact writers ──────────────────────────────────────────────────

def cosine_value(a: torch.Tensor, b: torch.Tensor) -> float:
    """Cosine similarity between two 1-D tensors."""
    a = a.float()
    b = b.float()
    denom = float(a.norm().item() * b.norm().item())
    if denom <= 0:
        return float("nan")
    return float(torch.dot(a, b).item() / denom)


def build_smiles_mu_cache(decoder, vocab, smiles_set: Iterable[str]) -> dict[str, torch.Tensor]:
    """Cache decoder encoder means ``mu`` for generated SMILES strings."""
    cache: dict[str, torch.Tensor] = {}
    for smiles in sorted(set(smiles_set)):
        if not isinstance(smiles, str) or not smiles.strip():
            continue
        try:
            token_ids = vocab.string2ids(smiles)
            if not token_ids:
                continue
            toks = torch.LongTensor([token_ids]).to(decoder.device)
            with torch.no_grad():
                out = decoder.forward_encoder(toks)
                mu = out[0] if isinstance(out, (tuple, list)) else out
            cache[smiles] = mu.squeeze(0).detach().cpu().float()
        except Exception:
            continue
    return cache


def build_alignment_report(
    parsed_feature_comp_ids: Iterable[int],
    embedding_ids: Iterable[int],
    retrieval_smiles_ids: Iterable[int],
    retrieval_mz_ids: Iterable[int],
) -> dict:
    """Measure overlap against comp_id-keyed bank and canonical SMILES table."""
    emb_set = set(int(i) for i in embedding_ids)
    parsed_set = set(int(i) for i in parsed_feature_comp_ids)
    smiles_set = set(int(i) for i in retrieval_smiles_ids)
    mz_set = set(int(i) for i in retrieval_mz_ids)

    return {
        "embedding_bank_count": len(emb_set),
        "retrieval_smiles_table_count": len(smiles_set),
        "retrieval_mz_table_count": len(mz_set),
        "feature_comp_ids_count": len(parsed_set),
        "feature_overlap_with_embedding_bank": len(parsed_set & emb_set),
        "smiles_table_overlap_with_embedding_bank": len(smiles_set & emb_set),
        "mz_table_overlap_with_embedding_bank": len(mz_set & emb_set),
        "embedding_bank_missing_smiles": len(emb_set - smiles_set),
        "embedding_bank_missing_mz": len(emb_set - mz_set),
    }


def summarize_numeric_array(values: np.ndarray) -> dict[str, float | int]:
    """Return count/min/median/mean/max summary for a numeric array."""
    arr = np.asarray(values, dtype=np.float32)
    if arr.size == 0:
        return {"count": 0, "min": 0.0, "median": 0.0, "mean": 0.0, "max": 0.0}
    return {
        "count": int(arr.size),
        "min": float(np.min(arr)),
        "median": float(np.median(arr)),
        "mean": float(np.mean(arr)),
        "max": float(np.max(arr)),
    }


def artifact_filename(prefix: str, stem: str, suffix: str) -> str:
    """Build ``{prefix}_{stem}{suffix}`` when prefix is non-empty."""
    normalized = str(prefix).strip()
    if not normalized:
        return f"{stem}{suffix}"
    return f"{normalized}_{stem}{suffix}"


def write_dataframe_with_aliases(df: pd.DataFrame, primary_path: Path, alias_paths: Iterable[Path]) -> None:
    """Write CSV to primary path and legacy alias filenames."""
    df.to_csv(primary_path, index=False)
    for alias_path in alias_paths:
        if alias_path == primary_path:
            continue
        df.to_csv(alias_path, index=False)


def write_json_with_aliases(payload: dict[str, object], primary_path: Path, alias_paths: Iterable[Path]) -> None:
    """Write JSON to primary path and legacy alias filenames."""
    text = json.dumps(payload, indent=2)
    primary_path.write_text(text)
    for alias_path in alias_paths:
        if alias_path == primary_path:
            continue
        alias_path.write_text(text)


def write_tensor_with_aliases(payload: object, primary_path: Path, alias_paths: Iterable[Path]) -> None:
    """Write ``torch.save`` payload to primary path and aliases."""
    torch.save(payload, primary_path)
    for alias_path in alias_paths:
        if alias_path == primary_path:
            continue
        torch.save(payload, alias_path)


def choose_encoder_checkpoint(path: Path) -> Path:
    """Resolve encoder checkpoint with historical fallback search paths."""
    if path.exists():
        return path

    candidates = [
        Path("models/spectra_encoder_209_trueZ_mass_filtered_ft.pt"),
        Path("models/spectra_encoder_209_trueZ_infonce_ft_conservative.pt"),
        Path("models/spectra_encoder_209_trueZ.pt"),
        Path("models/spectra_encoder_209.pt"),
        Path("models/spectra_encoder_new.pt"),
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()
    raise FileNotFoundError(f"Encoder checkpoint not found at {path} and no fallback checkpoint exists.")


# ── CLI ─────────────────────────────────────────────────────────────────────────

def get_parser() -> argparse.ArgumentParser:
    """Build CLI parser for Track A/B prediction and merge policies."""
    parser = argparse.ArgumentParser(
        description="Predict and score SMILES for all X_209 features using encoder+decoder."
    )
    parser.add_argument(
        "--feature_table",
        type=str,
        default="models/x209_features_encoder.pkl",
        help="Path to encoder-ready feature table (from parse_x209_html.py).",
    )
    parser.add_argument(
        "--encoder_ckpt",
        type=str,
        default="models/spectra_encoder_209_trueZ_mass_filtered_ft.pt",
        help="Path to Net1D checkpoint. If missing, fallback checkpoints are auto-tried.",
    )
    parser.add_argument("--decoder_dir", type=str, default="decoder/models")
    parser.add_argument(
        "--decoder_model_name",
        type=str,
        default="model_mass_209.pt",
        help="Preferred decoder checkpoint in decoder_dir; falls back to model.pt if missing.",
    )
    parser.add_argument(
        "--decoder_config_name",
        type=str,
        default="",
        help="Optional decoder config file in decoder_dir. Defaults to model-specific config if present, else config.nb.",
    )
    parser.add_argument(
        "--decoder_vocab_name",
        type=str,
        default="",
        help="Optional decoder vocab file in decoder_dir. Defaults to model-specific vocab if present, else vocab.nb.",
    )
    parser.add_argument("--embedding_bank_dir", type=str, default="embeddings_latent_from_smiles_expanded")
    parser.add_argument(
        "--retrieval_smiles_csv",
        type=str,
        default="SMILES_DB/comp_id_to_smiles_expanded.csv",
        help="comp_id-keyed SMILES table for retrieval identity.",
    )
    parser.add_argument(
        "--retrieval_name_csv", type=str, default="SMILES_DB/name_dict.csv", help="comp_id-keyed display-name table."
    )
    parser.add_argument(
        "--retrieval_mz_pkl",
        type=str,
        default="SMILES_DB/cid_to_parent_mz_expanded.pkl",
        help="comp_id-keyed precursor m/z pickle for mass gating.",
    )
    parser.add_argument(
        "--chemical_map_csv", type=str, default="", help="Optional enrichment from MASTER_METABOLITE_LIST.csv (legacy)."
    )
    parser.add_argument("--train_csv", type=str, default="train.csv")
    parser.add_argument("--output_dir", type=str, default="models/x209_smiles_predictions")
    parser.add_argument(
        "--artifact_prefix",
        type=str,
        default="x209",
        help="Filename prefix for retrieval/denovo/summary artifacts written to output_dir.",
    )
    parser.add_argument(
        "--run_decoder_mass_diagnosis",
        action="store_true",
        help="Run analyze_decoder_mass_conditioning.py after prediction outputs are written.",
    )
    parser.add_argument(
        "--enforce_decoder_mass_gate",
        action="store_true",
        help="When post-run diagnosis is enabled, exit non-zero if the decoder gate fails.",
    )
    parser.add_argument(
        "--decoder_diagnosis_output_dir",
        type=str,
        default="",
        help="Optional output directory for post-run decoder diagnosis. Defaults to a sibling of output_dir.",
    )
    parser.add_argument(
        "--decoder_summary_json",
        type=str,
        default="",
        help="Optional decoder fine-tune summary JSON used by the post-run diagnosis.",
    )
    parser.add_argument(
        "--decoder_gate_allowed_root_causes",
        type=str,
        default="spectrum_latent_off_manifold,conditioning_metadata_mismatch",
        help="Comma-separated root causes that count as passing states for the enforced decoder gate.",
    )
    parser.add_argument(
        "--decoder_gate_min_bank_targets_with_strict_pass",
        type=int,
        default=1,
        help="Minimum bank-side probe targets with any strict mass pass for the enforced decoder gate.",
    )
    parser.add_argument(
        "--decoder_gate_min_raw_features_within_1da",
        type=int,
        default=1,
        help="Minimum features with any valid raw draw within 1 Da for the enforced decoder gate.",
    )
    parser.add_argument(
        "--decoder_gate_fail_on_ri_mismatch",
        action="store_true",
        help="Treat RI training/inference mismatch as a hard failure in the enforced decoder gate.",
    )

    parser.add_argument(
        "--device",
        type=str,
        default="auto",
        help=DEVICE_HELP_TEXT,
    )
    parser.add_argument(
        "--ion_mode",
        type=str,
        choices=("negative", "positive"),
        default="negative",
        help="Ionization mode used by strict mass matching logic.",
    )
    parser.add_argument("--minmass", type=float, default=50.0)
    parser.add_argument("--maxmass", type=float, default=499.97)
    parser.add_argument(
        "--parent_mz_min",
        type=float,
        default=None,
        help="Parent-m/z normalization min. Defaults to --minmass.",
    )
    parser.add_argument(
        "--parent_mz_max",
        type=float,
        default=None,
        help="Parent-m/z normalization max. Defaults to --maxmass.",
    )
    parser.add_argument("--resolution", type=int, default=2)
    parser.add_argument(
        "--min_rel_intensity_pct",
        type=float,
        default=2.0,
        help="Drop peaks below this percentage of the spectrum base peak before vectorization.",
    )
    parser.add_argument(
        "--top_peaks",
        type=int,
        default=25,
        help="Retain at most this many peaks per spectrum after thresholding.",
    )
    parser.add_argument(
        "--disable_basepeak_norm",
        action="store_true",
        help="Bypass basepeak normalization after thresholding. Intended for the gated no-norm branch only.",
    )
    parser.add_argument("--decode_batch_size", type=int, default=128)
    parser.add_argument(
        "--num_samples",
        type=int,
        default=60,
        help="Total samples per feature across all schedule modes (increased for mass-first).",
    )
    parser.add_argument(
        "--skip_constrained_denovo",
        action="store_true",
        help="Skip Track B constrained de novo generation and emit retrieval-only merged outputs.",
    )
    parser.add_argument(
        "--decode_temp", type=float, default=1.0, help="Decoding temperature (default lowered from 2.0 to 1.0)."
    )
    parser.add_argument(
        "--decode_temps",
        type=str,
        default="0.3,0.5,0.8,1.0,1.2",
        help="Comma-separated temperature sweep. Samples split evenly across temps.",
    )
    parser.add_argument(
        "--top_p", type=float, default=0.0, help="Nucleus sampling threshold (0 disables, 0.9 recommended)."
    )
    parser.add_argument(
        "--top_p_schedule",
        type=str,
        default="0.0,0.85,0.95",
        help="Comma-separated top-p values to interleave with temp sweep.",
    )
    parser.add_argument(
        "--soft_mass_window_da", type=float, default=0.5, help="Soft mass window in Da for pre-strict shortlisting."
    )
    parser.add_argument(
        "--pre_strict_topn", type=int, default=50, help="Keep top-N candidates by soft score before final strict gate."
    )
    parser.add_argument("--decode_max_len", type=int, default=100)
    parser.add_argument(
        "--disable_ri_conditioning",
        action="store_true",
        help="Disable RT/RI conditioning for decoder sampling.",
    )
    parser.add_argument(
        "--ri_min",
        type=float,
        default=None,
        help="Manual RI min for normalization. Defaults to feature table min.",
    )
    parser.add_argument(
        "--ri_max",
        type=float,
        default=None,
        help="Manual RI max for normalization. Defaults to feature table max.",
    )
    parser.add_argument("--retrieval_topk", type=int, default=10)
    parser.add_argument(
        "--merge_policy",
        type=str,
        choices=(
            "retrieval_default",
            "pin_high_rank1_joint_rest",
            "high_protected_medium_low_joint",
            "retrieval_rank1_plus_best_selfies",
            "selfies_only_5da",
            "top5_selfies_then_joint_rest",
        ),
        default="pin_high_rank1_joint_rest",
        help=(
            "Merge policy for x209_merged_predictions.csv. "
            "retrieval_default preserves the historical retrieval-first ordering."
        ),
    )
    parser.add_argument(
        "--selfies_policy_mass_window_da",
        type=float,
        default=5.0,
        help="Mass window for policy-level SELFIES pool construction.",
    )
    parser.add_argument(
        "--selfies_policy_cap_per_feature",
        type=int,
        default=3,
        help="Maximum SELFIES rows per feature contributed to joint-sort policies.",
    )
    parser.add_argument(
        "--retrieval_rescoring",
        type=str,
        choices=("off", "on"),
        default="on",
        help="Enable fragment/intensity/html-formula retrieval rescoring. Set to off for baseline-refresh isolation.",
    )
    parser.add_argument(
        "--retrieval_weights_json",
        type=str,
        default="",
        help="Optional JSON with retrieval reranking weights and latent_cos_transform for rescoring runs.",
    )
    parser.add_argument(
        "--fragment_db_csv",
        type=str,
        default="SMILES_DB/library_209_fragment_database_element0.csv",
        help="Local fragment/intensity database used for retrieval rescoring.",
    )
    parser.add_argument(
        "--fragment_match_tolerance_da",
        type=float,
        default=0.02,
        help="Product-ion tolerance in Da for retrieval fragment/intensity matching.",
    )
    parser.add_argument(
        "--fragment_match_target",
        type=float,
        default=5.0,
        help="Target fragment match count used to normalize fragment_match_count_score.",
    )
    parser.add_argument(
        "--html_formula_mode",
        type=str,
        choices=("off", "soft"),
        default="soft",
        help="Use existing HTML-derived formulas as a soft reranking bonus when retrieval rescoring is on.",
    )
    parser.add_argument(
        "--prevalence_csv",
        type=str,
        default="SMILES_DB/fasstrecords_prevalence.csv",
        help="Optional inchikey14 prevalence table with prevalence_score from FASSTrecords.",
    )
    parser.add_argument(
        "--structural_neighbor_top_n",
        type=int,
        default=5,
        help="Top-N retrieval candidates used for substructure/tanimoto neighborhood scoring.",
    )
    parser.add_argument(
        "--structural_min_mcs_atoms",
        type=int,
        default=6,
        help="Minimum MCS atom count for substructure consensus scoring.",
    )
    parser.add_argument(
        "--structural_mcs_timeout_sec",
        type=int,
        default=1,
        help="Pairwise/top-N MCS timeout in seconds for structural scoring.",
    )
    parser.add_argument(
        "--candidate_formula_max",
        type=int,
        default=5,
        help="Maximum candidate-level HTML formula-set size treated as directly usable.",
    )
    parser.add_argument(
        "--encoder_calibrator",
        type=str,
        default="",
        help="Optional affine latent calibrator checkpoint applied after encoder inference.",
    )
    parser.add_argument("--max_candidates_per_feature", type=int, default=20)
    parser.add_argument(
        "--latent_proposal_alphas",
        type=str,
        default="0.1,0.2,0.3",
        help="Interpolation alphas for retrieval-guided latent proposals.",
    )
    parser.add_argument("--proposal_samples_per", type=int, default=3, help="Decode samples per latent proposal.")
    parser.add_argument(
        "--enable_latent_proposals", action="store_true", help="Enable retrieval-guided latent proposals for de novo."
    )
    parser.add_argument(
        "--compute_candidate_latent_cosine",
        action="store_true",
        help="Optional: re-encode generated SMILES and compute candidate latent cosine (slower).",
    )
    parser.add_argument(
        "--mass_filter_tolerance",
        type=float,
        default=5.0,
        help="Mass tolerance in Da for filtering candidates (0 disables).",
    )
    parser.add_argument(
        "--mass_filter_hard", action="store_true", help="If set, discard candidates outside mass tolerance entirely."
    )
    parser.add_argument(
        "--apply_tanh",
        action="store_true",
        default=True,
        help="Apply tanh to encoder output to match decoder's [-1,1] input range (default: True).",
    )
    parser.add_argument(
        "--no_tanh",
        action="store_true",
        help="Disable tanh on encoder output (use only for checkpoints trained without tanh).",
    )
    return parser


# ── Track B ranking helpers ─────────────────────────────────────────────────────

def canonicalize_smiles(smiles: str) -> str | None:
    """Return canonical SMILES or None if invalid."""
    if not HAS_RDKIT or not smiles or not smiles.strip():
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        Chem.SanitizeMol(mol)
        return Chem.MolToSmiles(mol, isomericSmiles=True)
    except Exception:
        return None


def soft_mass_score(
    predicted_mass: float,
    parent_mz: float,
    window_da: float = 0.5,
    ion_mode: str = "negative",
) -> float:
    """Continuous [0, 1] score: 1.0 at perfect match, decays over window_da."""
    err = adduct_neutral_mass_error(predicted_mass, parent_mz, ion_mode=ion_mode)
    if err <= 0.001:
        return 1.0
    if err >= window_da:
        return 0.0
    return 1.0 - (err / window_da)


def pre_strict_shortlist_score(
    mass_score: float,
    mean_log_prob: float,
    support_frac: float,
    latent_cosine: float = float("nan"),
) -> float:
    """Composite score for pre-strict shortlisting. Higher is better."""
    norm_lp = 1.0 / (1.0 + np.exp(-mean_log_prob - 1.0))
    cos_term = max(0.0, latent_cosine) if not np.isnan(latent_cosine) else 0.0
    return 0.40 * mass_score + 0.25 * norm_lp + 0.20 * support_frac + 0.15 * cos_term


def build_latent_proposals(
    z_query: torch.Tensor,
    retrieval_latents: list[torch.Tensor],
    alphas: list[float] = (0.1, 0.2, 0.3),
) -> list[torch.Tensor]:
    """Create interpolated latent proposals from retrieval analogs.

    Given z_query and retrieved z vectors, produce:
      - z_query (always included)
      - z_query + alpha * (z_retrieved - z_query)  for each alpha and retrieved z
    """
    proposals = [z_query]
    for z_ret in retrieval_latents:
        diff = z_ret - z_query
        for alpha in alphas:
            proposals.append(z_query + alpha * diff)
    return proposals


def decode_from_proposals(
    decoder,
    proposals: list[torch.Tensor],
    mass_cond: float,
    ri_cond: float | None,
    num_samples_per_proposal: int,
    temp: float,
    max_len: int,
    top_p: float = 0.0,
) -> tuple[Counter, dict[str, float]]:
    """Decode from multiple latent proposals, returning merged counts and best log-probs."""
    merged_counts: Counter = Counter()
    best_lps: dict[str, float] = {}

    device = decoder.device
    for z in proposals:
        z_batch = z.unsqueeze(0).to(device)
        mass_t = torch.tensor([[mass_cond]], dtype=torch.float32, device=device)
        ri_t = torch.tensor([[ri_cond or 0.0]], dtype=torch.float32, device=device)

        for _ in range(num_samples_per_proposal):
            result = decoder.sample(
                1,
                z=z_batch,
                temp=temp,
                max_len=max_len,
                top_p=top_p,
                return_log_probs=True,
                mass=mass_t,
                ri=ri_t,
            )
            decoded, mean_lps = result
            smi = decoded[0]
            lp = mean_lps[0]
            merged_counts[smi] += 1
            prev = best_lps.get(smi, float("-inf"))
            if lp > prev:
                best_lps[smi] = lp

    return merged_counts, best_lps


def denovo_rerank_score(
    mass_err_da: float,
    latent_cosine: float,
    log_prob: float,
    support_frac: float,
    retrieval_tanimoto: float = 0.0,
) -> float:
    """De novo reranker: composite score optimized for mass error + latent consistency + confidence."""
    mass_term = max(0.0, 1.0 - mass_err_da / 0.5)
    cos_term = max(0.0, latent_cosine) if not np.isnan(latent_cosine) else 0.0
    lp_term = 1.0 / (1.0 + np.exp(-log_prob - 1.0))
    ret_bonus = min(1.0, retrieval_tanimoto / 0.7) if retrieval_tanimoto > 0.3 else 0.0

    return 0.30 * mass_term + 0.25 * cos_term + 0.20 * lp_term + 0.10 * support_frac + 0.15 * ret_bonus


def rank_denovo_candidates(
    candidates: list[dict],
    max_candidates: int = 20,
) -> list[dict]:
    """Rank surviving de novo candidates by mass consistency, log-prob, cosine, support."""
    candidates.sort(
        key=lambda r: (
            -r["mass_match_score"],
            -np.nan_to_num(r["mean_log_prob"], nan=-99.0),
            -np.nan_to_num(r["candidate_latent_cosine"], nan=-1.0),
            -r["generation_support_fraction"],
            r["predicted_smiles"],
        )
    )
    return candidates[:max_candidates]


def main() -> None:
    """Run encoder, Track A retrieval, Track B de novo, and merged output."""
    args = get_parser().parse_args()
    ion_mode = normalize_ion_mode(str(args.ion_mode))

    feature_table_path = Path(args.feature_table).expanduser().resolve()
    decoder_dir = Path(args.decoder_dir).expanduser().resolve()
    embedding_bank_dir = Path(args.embedding_bank_dir).expanduser().resolve()
    retrieval_smiles_csv = Path(args.retrieval_smiles_csv).expanduser().resolve()
    retrieval_name_csv = Path(args.retrieval_name_csv).expanduser().resolve()
    retrieval_mz_pkl = Path(args.retrieval_mz_pkl).expanduser().resolve()
    fragment_db_csv = Path(args.fragment_db_csv).expanduser().resolve()
    prevalence_csv = Path(args.prevalence_csv).expanduser().resolve()
    output_dir = Path(args.output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    artifact_prefix = str(args.artifact_prefix).strip() or "x209"
    retrieval_rescoring_enabled = str(args.retrieval_rescoring).strip().lower() == "on"
    effective_html_formula_mode = str(args.html_formula_mode) if retrieval_rescoring_enabled else "off"
    retrieval_scoring_config = (
        resolve_scoring_config(args.retrieval_weights_json) if retrieval_rescoring_enabled else None
    )
    preprocessing_config = SpectrumPreprocessingConfig(
        min_rel_intensity_pct=float(args.min_rel_intensity_pct),
        top_peaks=int(args.top_peaks),
        disable_basepeak_norm=bool(args.disable_basepeak_norm),
    )

    device = resolve_device(args.device)
    n_bins = int((float(args.maxmass) - float(args.minmass)) * int(args.resolution)) + 1

    # --- Load retrieval identity tables (comp_id-keyed) ---
    print("Loading retrieval identity tables...")
    id_to_smiles = load_retrieval_smiles(retrieval_smiles_csv)
    id_to_name = load_retrieval_names(retrieval_name_csv) if retrieval_name_csv.exists() else {}
    id_to_inchikey14 = load_retrieval_inchikey14(retrieval_smiles_csv)
    mz_hints = load_retrieval_mass_hints(retrieval_smiles_csv, ion_mode=ion_mode)
    bank_mz = load_retrieval_mz(retrieval_mz_pkl, fallback_hints=mz_hints)
    prevalence_map = (
        load_prevalence_map(prevalence_csv) if retrieval_rescoring_enabled and prevalence_csv.exists() else {}
    )
    print(f"  SMILES table: {len(id_to_smiles)} comp_ids")
    print(f"  Name table: {len(id_to_name)} comp_ids")
    print(f"  InChIKey14 table: {len(id_to_inchikey14)} comp_ids")
    print(f"  Precursor m/z table: {len(bank_mz)} comp_ids")
    if mz_hints:
        print(f"  Precursor m/z hints from SMILES table: {len(mz_hints)} comp_ids")
    if retrieval_rescoring_enabled:
        print(f"  Prevalence table rows: {len(prevalence_map)}")

    if args.chemical_map_csv:
        chem_path = Path(args.chemical_map_csv).expanduser().resolve()
        if chem_path.exists():
            legacy_smiles, legacy_names = load_chemical_map(chem_path)
            for cid, smi in legacy_smiles.items():
                id_to_smiles.setdefault(cid, smi)
            for cid, nm in legacy_names.items():
                id_to_name.setdefault(cid, nm)
            print(f"  Enriched from legacy map: +{len(legacy_smiles)} SMILES, +{len(legacy_names)} names")

    # --- Encoder ---
    encoder_ckpt = choose_encoder_checkpoint(Path(args.encoder_ckpt).expanduser())
    print(f"Using encoder checkpoint: {encoder_ckpt}")
    print(f"Using decoder directory: {decoder_dir}")
    print(f"Using embedding bank: {embedding_bank_dir}")

    features_df = load_feature_table(feature_table_path)
    use_tanh = bool(args.apply_tanh) and not bool(args.no_tanh)
    parent_mz_min = float(args.parent_mz_min) if args.parent_mz_min is not None else float(args.minmass)
    parent_mz_max = float(args.parent_mz_max) if args.parent_mz_max is not None else float(args.maxmass)
    print(f"Encoder apply_tanh: {use_tanh}")
    print(f"Spectrum preprocessing: {preprocessing_config.to_dict()}")
    print(f"Retrieval rescoring: {'enabled' if retrieval_rescoring_enabled else 'disabled'}")
    if retrieval_rescoring_enabled:
        print(f"Retrieval scoring config: {retrieval_scoring_config.to_dict()}")
    print(f"Ion mode: {ion_mode}")

    encoder, encoder_n_bins, encoder_load_info = load_encoder_checkpoint(
        ckpt_path=encoder_ckpt,
        n_bins=n_bins,
        minmass=float(args.minmass),
        maxmass=float(args.maxmass),
        parent_mz_min=parent_mz_min,
        parent_mz_max=parent_mz_max,
        device=device,
        apply_tanh=use_tanh,
    )
    if encoder_n_bins != n_bins:
        print(
            "Encoder input width override: "
            f"requested_n_bins={n_bins}, checkpoint_compatible_n_bins={encoder_n_bins}. "
            "Spectra will be zero-padded or truncated after preprocessing."
        )
    skipped_mismatched = dict(encoder_load_info.get("skipped_mismatched") or {})
    if skipped_mismatched:
        print(
            "Encoder checkpoint load skipped mismatched tensors: "
            + ", ".join(
                f"{name}({payload['checkpoint']}->{payload['model']})"
                for name, payload in sorted(skipped_mismatched.items())
            )
        )
    calibrator = None
    calibrator_path = None
    if str(args.encoder_calibrator).strip():
        calibrator_path = Path(args.encoder_calibrator).expanduser().resolve()
        if not calibrator_path.exists():
            raise FileNotFoundError(f"Encoder calibrator not found: {calibrator_path}")
        calibrator = load_latent_calibrator(calibrator_path, device=device)
        print(f"Using encoder calibrator: {calibrator_path}")

    feature_rows, feature_latents_raw, feature_latent_map_raw, feature_peak_map = encode_features(
        features_df=features_df,
        encoder=encoder,
        minmass=float(args.minmass),
        maxmass=float(args.maxmass),
        resolution=int(args.resolution),
        device=device,
        preprocessing_config=preprocessing_config,
        encoder_n_bins=encoder_n_bins,
    )
    feature_latents = apply_latent_calibrator(feature_latents_raw.to(device), calibrator).detach().cpu()
    feature_latent_map = {
        feature_key: feature_latents[idx].detach().cpu()
        for idx, feature_key in enumerate(feature_latent_map_raw.keys())
    }
    feature_meta_df = pd.DataFrame(feature_rows)
    feature_keys = feature_meta_df["feature_key"].tolist()
    feature_parent_mzs = feature_meta_df["parent_mz"].tolist()
    print(f"Encoded features: {len(feature_keys)}")

    latent_norms_pre = feature_latents_raw.norm(dim=1).cpu().numpy() if len(feature_latents_raw) else np.array([])
    latent_norms_post = feature_latents.norm(dim=1).cpu().numpy() if len(feature_latents) else np.array([])

    feature_latents_pt = output_dir / artifact_filename(artifact_prefix, "feature_latents", ".pt")
    write_tensor_with_aliases(
        feature_latent_map,
        feature_latents_pt,
        [output_dir / "feature_latents.pt"],
    )
    feature_meta_csv = output_dir / artifact_filename(artifact_prefix, "feature_meta", ".csv")
    write_dataframe_with_aliases(
        feature_meta_df,
        feature_meta_csv,
        [output_dir / "feature_meta.csv"],
    )

    print("Parsing HTML-derived isotope audit signals...")
    html_feature_df, html_isotope_df, html_signal_lookup = build_html_signal_feature_tables(
        feature_meta_df=feature_meta_df,
        candidate_formula_max=int(args.candidate_formula_max),
    )
    html_isotope_csv = output_dir / artifact_filename(artifact_prefix, "html_isotope_signals", ".csv")
    html_isotope_serialized = serialize_list_columns_for_csv(
        html_isotope_df,
        ["isotope_element_hints", "observed_isotope_rows_json"],
    )
    write_dataframe_with_aliases(
        html_isotope_serialized,
        html_isotope_csv,
        [output_dir / "html_isotope_signals.csv"],
    )

    html_formula_lookup: dict[str, set[str]] | None = None
    if retrieval_rescoring_enabled and effective_html_formula_mode == "soft":
        print("Parsing HTML-derived formula evidence for soft retrieval bonus...")
        if not html_feature_df.empty:
            html_formula_lookup = {
                str(row.feature_key): set(row.usable_formula_set if isinstance(row.usable_formula_set, list) else [])
                for row in html_feature_df.itertuples(index=False)
            }
            html_feature_csv = output_dir / artifact_filename(artifact_prefix, "html_formula_signals", ".csv")
            html_feature_serialized = serialize_list_columns_for_csv(
                html_feature_df,
                ["usable_formula_set", "feature_level_formula_set", "candidate_level_formula_set"],
            )
            write_dataframe_with_aliases(
                html_feature_serialized,
                html_feature_csv,
                [output_dir / "html_formula_signals.csv"],
            )
        else:
            html_feature_csv = None
    else:
        if not retrieval_rescoring_enabled and str(args.html_formula_mode) != "off":
            print("Skipping HTML formula parsing because retrieval rescoring is disabled.")
        html_feature_csv = None

    # --- Embedding bank ---
    bank_ids, bank_latents = load_embedding_bank(embedding_bank_dir)
    bank_latents = apply_latent_calibrator(bank_latents.to(device), calibrator).detach().cpu()
    fragment_db = load_fragment_database(fragment_db_csv) if retrieval_rescoring_enabled else None

    # ====================================================================
    # TRACK A: Mass-gated retrieval
    # ====================================================================
    print("\n=== Track A: Mass-gated retrieval ===")
    retrieval_df = compute_mass_gated_retrieval(
        feature_keys=feature_keys,
        feature_parent_mzs=feature_parent_mzs,
        feature_latents=feature_latents,
        feature_peak_map=feature_peak_map,
        bank_ids=bank_ids,
        bank_latents=bank_latents,
        bank_mz=bank_mz,
        id_to_smiles=id_to_smiles,
        id_to_name=id_to_name,
        id_to_inchikey14=id_to_inchikey14,
        fragment_db=fragment_db,
        html_formula_lookup=html_formula_lookup,
        retrieval_scoring_config=retrieval_scoring_config,
        prevalence_map=prevalence_map,
        topk=int(args.retrieval_topk),
        fragment_match_tolerance_da=float(args.fragment_match_tolerance_da),
        fragment_match_target=float(args.fragment_match_target),
        structural_neighbor_top_n=int(args.structural_neighbor_top_n),
        structural_min_mcs_atoms=int(args.structural_min_mcs_atoms),
        structural_mcs_timeout_sec=int(args.structural_mcs_timeout_sec),
    )
    retrieval_df = annotate_dataframe_with_isotope_audit(
        retrieval_df,
        feature_signal_lookup=html_signal_lookup,
        formula_column="candidate_formula",
        smiles_column="smiles",
    )
    retrieval_csv_path = output_dir / artifact_filename(artifact_prefix, "retrieval_report", ".csv")
    retrieval_serialized = serialize_list_columns_for_csv(
        retrieval_df,
        ["isotope_element_hints"],
    )
    write_dataframe_with_aliases(
        retrieval_serialized,
        retrieval_csv_path,
        [output_dir / "retrieval_report.csv"],
    )

    tier_counts = retrieval_df[retrieval_df["rank"] == 1]["confidence_tier"].value_counts()
    print(f"  Retrieval tier distribution (rank-1): {tier_counts.to_dict()}")
    rec_count = int(retrieval_df["recommended"].sum())
    print(f"  Recommended analogs: {rec_count}")
    retrieval_rank1 = retrieval_df[retrieval_df["rank"] == 1].copy()
    retrieval_rank1_cos = pd.to_numeric(retrieval_rank1.get("latent_cosine"), errors="coerce").dropna()
    retrieval_rank1_score = pd.to_numeric(retrieval_rank1.get("retrieval_score"), errors="coerce").dropna()
    skip_constrained_denovo = bool(args.skip_constrained_denovo) or int(args.num_samples) <= 0

    # ====================================================================
    # TRACK B: Constrained de novo generation
    # ====================================================================
    print("\n=== Track B: Constrained de novo generation ===")
    decoder = None
    vocab = None
    decoder_model_path = None
    decoder_config_path = None
    decoder_vocab_path = None
    candidate_smiles_set = set()
    candidate_metrics_df = pd.DataFrame()
    denovo_df = pd.DataFrame()
    features_with_denovo = 0
    raw_samples_csv = output_dir / artifact_filename(artifact_prefix, "generation_samples_raw", ".csv")
    denovo_csv_path = output_dir / artifact_filename(artifact_prefix, "denovo_report", ".csv")
    candidate_metrics_csv_path = output_dir / artifact_filename(artifact_prefix, "candidate_metrics_all", ".csv")
    RI_COVERAGE_THRESHOLD = 0.05  # require >5% RI coverage to use conditioning
    ri_vals = pd.to_numeric(feature_meta_df["rt_or_ri"], errors="coerce").tolist()
    ri_finite = [float(v) for v in ri_vals if np.isfinite(v) and float(v) > 0]
    ri_coverage = len(ri_finite) / max(len(ri_vals), 1)
    ri_auto_disabled = False
    ri_min_used = None
    ri_max_used = None
    feature_ri_norm = None

    if skip_constrained_denovo:
        if bool(args.skip_constrained_denovo):
            print("  Skipping constrained de novo generation by user flag.")
        else:
            print("  Skipping constrained de novo generation because --num_samples <= 0.")
        raw_samples_df = pd.DataFrame()
        write_dataframe_with_aliases(
            raw_samples_df,
            raw_samples_csv,
            [output_dir / "generation_samples_raw.csv"],
        )
        write_dataframe_with_aliases(
            denovo_df,
            denovo_csv_path,
            [output_dir / "denovo_report.csv"],
        )
        write_dataframe_with_aliases(
            candidate_metrics_df,
            candidate_metrics_csv_path,
            [output_dir / "candidate_metrics_all.csv"],
        )
    else:
        decoder, vocab, decoder_model_path, decoder_config_path, decoder_vocab_path = load_decoder(
            decoder_dir=decoder_dir,
            device=device,
            model_name=str(args.decoder_model_name),
            config_name=str(args.decoder_config_name),
            vocab_name=str(args.decoder_vocab_name),
        )
        decoder_summary_path = derive_decoder_summary_path(
            decoder_dir=decoder_dir,
            model_name=str(args.decoder_model_name),
            explicit_path=str(args.decoder_summary_json),
        )
        decoder_summary_payload = maybe_load_decoder_summary_payload(decoder_summary_path)
        decoder_mass_min, decoder_mass_max, decoder_mass_source = resolve_decoder_mass_normalization(
            decoder_summary_payload
        )
        mass_condition_min = float(decoder_mass_min) if decoder_mass_min is not None else float(parent_mz_min)
        mass_condition_max = float(decoder_mass_max) if decoder_mass_max is not None else float(parent_mz_max)
        print(
            "  Mass conditioning normalization:",
            f"min={mass_condition_min:.4f}",
            f"max={mass_condition_max:.4f}",
            f"source={decoder_mass_source}",
        )

        feature_mass_norm = [
            normalize_scalar(pmz, mass_condition_min, mass_condition_max) for pmz in feature_parent_mzs
        ]
        if not bool(args.disable_ri_conditioning):
            if ri_coverage < RI_COVERAGE_THRESHOLD:
                ri_auto_disabled = True
                print(
                    f"  RI conditioning AUTO-DISABLED: coverage={ri_coverage:.2%} "
                    f"({len(ri_finite)}/{len(ri_vals)}) < threshold={RI_COVERAGE_THRESHOLD:.0%}"
                )
                print("  Mass conditioning remains active.")
            else:
                ri_min_used = float(args.ri_min) if args.ri_min is not None else float(np.min(ri_finite))
                ri_max_used = float(args.ri_max) if args.ri_max is not None else float(np.max(ri_finite))
                feature_ri_norm = [
                    normalize_scalar(v, ri_min_used, ri_max_used) if np.isfinite(v) and float(v) > 0 else 0.0
                    for v in ri_vals
                ]
                print(
                    f"  RI conditioning enabled: coverage={ri_coverage:.2%}, "
                    f"normalized with min={ri_min_used:.4f}, max={ri_max_used:.4f}"
                )
        else:
            print("  RI conditioning disabled (user flag).")

        temps = None
        if args.decode_temps:
            temps = [float(t.strip()) for t in args.decode_temps.split(",") if t.strip()]
            print(f"  Temperature sweep: {temps}")

        gen_counts, gen_draws, raw_samples_df, best_log_probs = sample_generation_support(
            decoder=decoder,
            feature_keys=feature_keys,
            latent_tensor=feature_latents,
            mass_conditions=feature_mass_norm,
            ri_conditions=feature_ri_norm,
            num_samples=int(args.num_samples),
            temp=float(args.decode_temp),
            batch_size=int(args.decode_batch_size),
            max_len=int(args.decode_max_len),
            top_p=float(args.top_p),
            temps=temps,
        )
        write_dataframe_with_aliases(
            raw_samples_df,
            raw_samples_csv,
            [output_dir / "generation_samples_raw.csv"],
        )

        if bool(args.enable_latent_proposals):
            proposal_alphas = [float(a.strip()) for a in args.latent_proposal_alphas.split(",") if a.strip()]
            samples_per = int(args.proposal_samples_per)
            print(f"  Generating latent proposals: alphas={proposal_alphas}, samples_per={samples_per}")

            retrieval_latent_cache: dict[int, torch.Tensor] = {}
            for bid, bz in zip(bank_ids, bank_latents):
                retrieval_latent_cache[bid] = bz

            n_proposal_generated = 0
            for i, feature_key in enumerate(feature_keys):
                fk_retrieval = retrieval_df[retrieval_df["feature_key"] == feature_key]
                fk_retrieval = fk_retrieval[fk_retrieval["status"] == "ok"].head(3)

                if fk_retrieval.empty:
                    continue

                z_query = feature_latents[i]
                ret_latents = []
                for _, rrow in fk_retrieval.iterrows():
                    cid = int(rrow["comp_id"]) if rrow["comp_id"] != "" else None
                    if cid is not None and cid in retrieval_latent_cache:
                        ret_latents.append(retrieval_latent_cache[cid])

                if not ret_latents:
                    continue

                proposals = build_latent_proposals(z_query, ret_latents, alphas=proposal_alphas)
                mass_cond_val = feature_mass_norm[i]
                ri_cond_val = feature_ri_norm[i] if feature_ri_norm is not None else None

                proposal_counts, proposal_lps = decode_from_proposals(
                    decoder=decoder,
                    proposals=proposals[1:],
                    mass_cond=mass_cond_val,
                    ri_cond=ri_cond_val,
                    num_samples_per_proposal=samples_per,
                    temp=float(args.decode_temp),
                    max_len=int(args.decode_max_len),
                    top_p=0.85,
                )

                for smi, cnt in proposal_counts.items():
                    gen_counts[feature_key][smi] += cnt
                    gen_draws[feature_key] = gen_draws.get(feature_key, 0) + cnt
                    n_proposal_generated += cnt
                for smi, lp in proposal_lps.items():
                    prev = best_log_probs.get(feature_key, {}).get(smi, float("-inf"))
                    if lp > prev:
                        best_log_probs.setdefault(feature_key, {})[smi] = lp

            print(f"  Latent proposals generated: {n_proposal_generated} additional samples")

        print("  Canonicalizing and deduplicating raw generations...")
        canonical_counts: dict[str, Counter] = {k: Counter() for k in feature_keys}
        canonical_log_probs: dict[str, dict[str, float]] = defaultdict(dict)
        raw_to_canonical: dict[str, str | None] = {}
        n_raw_total = 0
        n_valid_unique = 0

        for fk in feature_keys:
            fk_lp = best_log_probs.get(fk, {})
            for raw_smi, count in gen_counts[fk].items():
                n_raw_total += count
                if raw_smi not in raw_to_canonical:
                    raw_to_canonical[raw_smi] = canonicalize_smiles(raw_smi)
                can = raw_to_canonical[raw_smi]
                if can is None:
                    continue
                canonical_counts[fk][can] += count
                lp = fk_lp.get(raw_smi, float("nan"))
                prev = canonical_log_probs[fk].get(can, float("-inf"))
                if not np.isnan(lp) and lp > prev:
                    canonical_log_probs[fk][can] = lp

        for fk in feature_keys:
            n_valid_unique += len(canonical_counts[fk])

        print(f"  Raw samples: {n_raw_total}, valid unique canonical: {n_valid_unique}")

        for counter in canonical_counts.values():
            candidate_smiles_set.update(counter.keys())
        smiles_mu_cache: dict[str, torch.Tensor] = {}
        if bool(args.compute_candidate_latent_cosine):
            smiles_mu_cache = build_smiles_mu_cache(
                decoder=decoder,
                vocab=vocab,
                smiles_set=candidate_smiles_set,
            )

        if HAS_RDKIT:
            print("  RDKit available -- mass-first scoring enabled.")
        else:
            print("  WARNING: RDKit not installed. Validity/mass checks disabled.")

        feature_latent_lookup = {feature_keys[i]: feature_latents[i].detach().cpu() for i in range(len(feature_keys))}
        feature_meta_lookup = {row["feature_key"]: row for row in feature_rows}
        max_candidates = int(args.max_candidates_per_feature)
        soft_window = float(args.soft_mass_window_da)
        pre_strict_n = int(args.pre_strict_topn)

        n_soft_pass = 0
        n_strict_pass = 0
        denovo_rows = []
        candidate_metrics_rows = []
        for feature_key in feature_keys:
            counter = canonical_counts[feature_key]
            draws = max(gen_draws.get(feature_key, 0), 1)
            meta = feature_meta_lookup[feature_key]
            parent_mz = float(meta["parent_mz"])
            feature_latent = feature_latent_lookup[feature_key]
            fk_log_probs_can = canonical_log_probs.get(feature_key, {})

            scored = []
            for can_smi, count in counter.items():
                pred_mass = exact_mol_weight(can_smi)
                if pred_mass is None:
                    continue
                abs_mass_error_da = adduct_neutral_mass_error(pred_mass, parent_mz, ion_mode=ion_mode)

                s_mass = soft_mass_score(pred_mass, parent_mz, window_da=soft_window, ion_mode=ion_mode)
                support = float(count) / float(draws)
                mean_lp = fk_log_probs_can.get(can_smi, float("nan"))

                candidate_cos = float("nan")
                if bool(args.compute_candidate_latent_cosine):
                    mu = smiles_mu_cache.get(can_smi)
                    if mu is not None:
                        candidate_cos = cosine_value(feature_latent, mu)

                shortlist_score = pre_strict_shortlist_score(
                    s_mass,
                    mean_lp,
                    support,
                    candidate_cos,
                )

                scored.append(
                    {
                        "predicted_smiles": can_smi,
                        "is_valid": True,
                        "predicted_exact_mass": pred_mass,
                        "predicted_formula": molecular_formula(can_smi),
                        "abs_mass_error_da": float(abs_mass_error_da),
                        "soft_mass_score": s_mass,
                        "mass_match_score": mass_match_score(pred_mass, parent_mz, ion_mode=ion_mode),
                        "strict_mass_pass": passes_strict_mass_filter(pred_mass, parent_mz, ion_mode=ion_mode),
                        "mean_log_prob": mean_lp,
                        "generation_support_count": int(count),
                        "generation_support_fraction": support,
                        "candidate_latent_cosine": candidate_cos,
                        "shortlist_score": shortlist_score,
                    }
                )

            for row in scored:
                candidate_metrics_rows.append(
                    {
                        "feature_key": feature_key,
                        "comp_id": meta["comp_id"],
                        "parent_mz": parent_mz,
                        **{
                            key: value
                            for key, value in row.items()
                            if key not in ("strict_mass_pass", "shortlist_score", "soft_mass_score")
                        },
                    }
                )

            scored.sort(key=lambda r: -r["shortlist_score"])
            shortlisted = scored[:pre_strict_n]
            n_soft_pass += len(shortlisted)

            strict_survivors = [c for c in shortlisted if c["strict_mass_pass"]]
            n_strict_pass += len(strict_survivors)

            ranked = rank_denovo_candidates(strict_survivors, max_candidates=max_candidates)
            for rank, row in enumerate(ranked, start=1):
                denovo_rows.append(
                    {
                        "feature_key": feature_key,
                        "comp_id": meta["comp_id"],
                        "feature_label": meta["feature_label"],
                        "parent_mz": parent_mz,
                        "rank": rank,
                        **{
                            k: v
                            for k, v in row.items()
                            if k not in ("strict_mass_pass", "shortlist_score", "soft_mass_score")
                        },
                        "samples_drawn": draws,
                    }
                )

        print(f"  Mass-first pipeline: {n_soft_pass} soft-pass -> {n_strict_pass} strict-pass")

        denovo_df = pd.DataFrame(denovo_rows)
        if not denovo_df.empty:
            denovo_df = denovo_df.sort_values(["feature_key", "rank"]).reset_index(drop=True)
        denovo_df = annotate_dataframe_with_isotope_audit(
            denovo_df,
            feature_signal_lookup=html_signal_lookup,
            formula_column="predicted_formula",
            smiles_column="predicted_smiles",
        )
        denovo_serialized = serialize_list_columns_for_csv(
            denovo_df,
            ["isotope_element_hints"],
        )
        write_dataframe_with_aliases(
            denovo_serialized,
            denovo_csv_path,
            [output_dir / "denovo_report.csv"],
        )
        candidate_metrics_df = pd.DataFrame(candidate_metrics_rows)
        if not candidate_metrics_df.empty:
            candidate_metrics_df = candidate_metrics_df.sort_values(
                ["feature_key", "mass_match_score", "abs_mass_error_da", "predicted_smiles"],
                ascending=[True, False, True, True],
                kind="mergesort",
            ).reset_index(drop=True)
        candidate_metrics_df = annotate_dataframe_with_isotope_audit(
            candidate_metrics_df,
            feature_signal_lookup=html_signal_lookup,
            formula_column="predicted_formula",
            smiles_column="predicted_smiles",
        )
        candidate_metrics_serialized = serialize_list_columns_for_csv(
            candidate_metrics_df,
            ["isotope_element_hints"],
        )
        write_dataframe_with_aliases(
            candidate_metrics_serialized,
            candidate_metrics_csv_path,
            [output_dir / "candidate_metrics_all.csv"],
        )

        features_with_denovo = denovo_df["feature_key"].nunique() if not denovo_df.empty else 0
        print(f"  De novo: {len(denovo_df)} candidates across {features_with_denovo}/{len(feature_keys)} features")

    # ====================================================================
    # Merged output
    # ====================================================================
    print("\n=== Merged output ===")
    authoritative_track_default = "retrieval"
    merge_policy_name = "retrieval_default" if skip_constrained_denovo else str(args.merge_policy)
    merge_policy_joint_sort_used = False
    merge_policy_selfies_cap = 0
    merge_policy_selfies_pool_rows = 0

    if merge_policy_name == "retrieval_default":
        merged_rows = []
        for _, rrow in retrieval_df.iterrows():
            merged_rows.append(
                {
                    "feature_key": rrow["feature_key"],
                    "track": "retrieval",
                    "rank": rrow["rank"],
                    "authoritative_track_default": authoritative_track_default,
                    "comp_id": rrow.get("comp_id", ""),
                    "smiles": rrow.get("smiles", ""),
                    "name": rrow.get("name", ""),
                    "latent_cosine": rrow.get("latent_cosine", float("nan")),
                    "gate_used": rrow.get("gate_used", float("nan")),
                    "abs_mass_diff": rrow.get("abs_mass_diff", float("nan")),
                    "confidence_tier": rrow.get("confidence_tier", ""),
                    "recommended": rrow.get("recommended", False),
                    "candidate_formula": rrow.get("candidate_formula", None),
                    "retrieval_score": rrow.get("retrieval_score", None),
                    "fragment_match_count": rrow.get("fragment_match_count", None),
                    "fragment_matched_peak_fraction": rrow.get("fragment_matched_peak_fraction", None),
                    "fragment_explained_intensity_fraction": rrow.get("fragment_explained_intensity_fraction", None),
                    "predicted_exact_mass": None,
                    "mass_match_score": None,
                    "mean_log_prob": None,
                    "generation_support_fraction": None,
                    "candidate_latent_cosine": None,
                }
            )

        for _, drow in denovo_df.iterrows():
            merged_rows.append(
                {
                    "feature_key": drow["feature_key"],
                    "track": "de_novo",
                    "rank": drow["rank"],
                    "authoritative_track_default": authoritative_track_default,
                    "comp_id": drow.get("comp_id", ""),
                    "smiles": drow.get("predicted_smiles", ""),
                    "name": "",
                    "latent_cosine": None,
                    "gate_used": None,
                    "abs_mass_diff": None,
                    "confidence_tier": None,
                    "recommended": None,
                    "candidate_formula": drow.get("predicted_formula", None),
                    "retrieval_score": None,
                    "fragment_match_count": None,
                    "fragment_matched_peak_fraction": None,
                    "fragment_explained_intensity_fraction": None,
                    "predicted_exact_mass": drow.get("predicted_exact_mass", None),
                    "mass_match_score": drow.get("mass_match_score", None),
                    "mean_log_prob": drow.get("mean_log_prob", None),
                    "generation_support_fraction": drow.get("generation_support_fraction", None),
                    "candidate_latent_cosine": drow.get("candidate_latent_cosine", None),
                }
            )

        merged_df = pd.DataFrame(merged_rows)
        if not merged_df.empty:
            merged_priority = {"retrieval": 0, "de_novo": 1}
            merged_df["merged_rank"] = (
                merged_df.assign(_track_priority=merged_df["track"].map(merged_priority).fillna(99).astype(int))
                .sort_values(["feature_key", "_track_priority", "rank", "smiles"], kind="mergesort")
                .groupby("feature_key")
                .cumcount()
                .add(1)
            )
            merged_df["is_authoritative_track"] = merged_df["track"].eq(authoritative_track_default)
            merged_df = merged_df.drop(columns=["_track_priority"], errors="ignore")
    else:
        retrieval_topk_df = retrieval_df.copy()
        retrieval_topk_df["feature_key"] = retrieval_topk_df["feature_key"].astype(str)
        retrieval_topk_df["rank"] = pd.to_numeric(retrieval_topk_df["rank"], errors="coerce").fillna(9999).astype(int)
        retrieval_topk_df = retrieval_topk_df[retrieval_topk_df["rank"] <= int(args.retrieval_topk)].copy()

        policy_candidates_df = candidate_metrics_df.copy()
        if not policy_candidates_df.empty:
            policy_candidates_df["feature_key"] = policy_candidates_df["feature_key"].astype(str)
            policy_candidates_df["predicted_smiles"] = policy_candidates_df["predicted_smiles"].astype(str)

        selfies_pool = rank_selfies_pool(
            policy_candidates_df,
            tolerance_da=float(args.selfies_policy_mass_window_da),
        )
        merge_policy_selfies_pool_rows = int(len(selfies_pool))

        feature_keys_for_policy = sorted(set(retrieval_topk_df["feature_key"].tolist()))
        rank1_df = retrieval_topk_df[retrieval_topk_df["rank"] == 1].copy()
        rank1_tier_map = {
            str(row.feature_key): str(row.confidence_tier)
            for row in rank1_df[["feature_key", "confidence_tier"]].itertuples(index=False)
        }

        retrieval_merge_df = map_retrieval_rows_for_merge(retrieval_topk_df)
        selfies_merge_df = map_selfies_rows_for_merge(selfies_pool, feature_meta_df=feature_meta_df)
        selfies_cap_df = select_selfies_rows_per_feature(
            selfies_merge_df, max_per_feature=int(args.selfies_policy_cap_per_feature)
        )

        if merge_policy_name == "pin_high_rank1_joint_rest":
            raw_ranked = build_policy_pin_high_rank1_joint_rest(
                feature_keys_for_policy,
                retrieval_merge_df,
                selfies_cap_df,
                rank1_tier_map,
            )
            merge_policy_joint_sort_used = True
            merge_policy_selfies_cap = int(args.selfies_policy_cap_per_feature)
        elif merge_policy_name == "high_protected_medium_low_joint":
            raw_ranked = build_policy_high_protected_medium_low_joint(
                feature_keys_for_policy,
                retrieval_merge_df,
                selfies_cap_df,
                rank1_tier_map,
            )
            merge_policy_joint_sort_used = True
            merge_policy_selfies_cap = int(args.selfies_policy_cap_per_feature)
        elif merge_policy_name == "retrieval_rank1_plus_best_selfies":
            raw_ranked = build_policy_retrieval_rank1_plus_best_selfies(
                feature_keys_for_policy,
                retrieval_merge_df,
                selfies_merge_df,
            )
        elif merge_policy_name == "selfies_only_5da":
            raw_ranked = build_policy_selfies_only(
                feature_keys_for_policy,
                selfies_merge_df,
            )
        elif merge_policy_name == "top5_selfies_then_joint_rest":
            raw_ranked = build_policy_top5_selfies_then_joint_rest(
                feature_keys_for_policy,
                retrieval_merge_df,
                selfies_cap_df,
            )
            merge_policy_joint_sort_used = True
            merge_policy_selfies_cap = int(args.selfies_policy_cap_per_feature)
        else:
            raise ValueError(f"Unsupported merge policy: {merge_policy_name}")

        ranked_with_marks = finalize_policy_rows(raw_ranked, topk=int(args.retrieval_topk))
        merged_df = to_merged_predictions_csv_schema(ranked_with_marks, topk=int(args.retrieval_topk))

    merged_df = annotate_dataframe_with_isotope_audit(
        merged_df,
        feature_signal_lookup=html_signal_lookup,
        formula_column="candidate_formula",
        smiles_column="smiles",
    )
    merged_csv_path = output_dir / artifact_filename(artifact_prefix, "merged_predictions", ".csv")
    merged_serialized = serialize_list_columns_for_csv(
        merged_df,
        ["isotope_element_hints"],
    )
    write_dataframe_with_aliases(
        merged_serialized,
        merged_csv_path,
        [output_dir / "merged_predictions.csv"],
    )
    print(
        "  Merged:",
        f"{len(merged_df)} rows",
        f"(policy={merge_policy_name}, retrieval_rows={len(retrieval_df)}, denovo_rows={len(denovo_df)})",
    )

    # ====================================================================
    # Alignment report
    # ====================================================================
    alignment_report = build_alignment_report(
        parsed_feature_comp_ids=feature_meta_df["comp_id"].tolist(),
        embedding_ids=bank_ids,
        retrieval_smiles_ids=id_to_smiles.keys(),
        retrieval_mz_ids=bank_mz.keys(),
    )
    alignment_report_path = output_dir / artifact_filename(artifact_prefix, "id_alignment_report", ".json")
    write_json_with_aliases(
        alignment_report,
        alignment_report_path,
        [output_dir / "id_alignment_report.json"],
    )

    # ====================================================================
    # Run summary
    # ====================================================================
    diagnosis_result = {
        "requested": False,
        "enforced": False,
        "return_code": None,
        "decoder_summary_json": None,
        "diagnosis_output_dir": None,
        "gate_status_json": None,
        "gate_pass": None,
        "ion_mode": ion_mode,
    }
    diagnosis_output_dir = None
    if bool(args.run_decoder_mass_diagnosis):
        decoder_summary_json = derive_decoder_summary_path(
            decoder_dir=decoder_dir,
            model_name=str(args.decoder_model_name),
            explicit_path=str(args.decoder_summary_json),
        )
        if not decoder_summary_json.exists():
            raise FileNotFoundError(f"Decoder summary JSON not found for post-run diagnosis: {decoder_summary_json}")
        diagnosis_output_dir = derive_decoder_diagnosis_output_dir(
            output_dir=output_dir,
            explicit_path=str(args.decoder_diagnosis_output_dir),
        )
        diagnosis_output_dir.mkdir(parents=True, exist_ok=True)
        diagnosis_result = run_decoder_mass_diagnosis(
            python_executable=sys.executable,
            repo_root=REPO_ROOT,
            prediction_dir=output_dir,
            decoder_dir=decoder_dir,
            decoder_model_name=str(args.decoder_model_name),
            decoder_summary_json=decoder_summary_json,
            diagnosis_output_dir=diagnosis_output_dir,
            enforce_gate=bool(args.enforce_decoder_mass_gate),
            gate_allowed_root_causes=str(args.decoder_gate_allowed_root_causes),
            gate_min_bank_targets_with_strict_pass=int(args.decoder_gate_min_bank_targets_with_strict_pass),
            gate_min_raw_features_within_1da=int(args.decoder_gate_min_raw_features_within_1da),
            gate_fail_on_ri_mismatch=bool(args.decoder_gate_fail_on_ri_mismatch),
            ion_mode=ion_mode,
        )

    feature_source_audit_path, feature_source_audit_summary = load_feature_source_audit_summary(feature_table_path)
    ri_audit = {
        "artifact_prefix": artifact_prefix,
        "feature_table": str(feature_table_path),
        "feature_count": int(len(ri_vals)),
        "finite_positive_ri_count": int(len(ri_finite)),
        "ri_coverage_fraction": float(ri_coverage),
        "ri_conditioning_requested": not bool(args.disable_ri_conditioning),
        "ri_conditioning_auto_disabled": bool(ri_auto_disabled),
        "ri_conditioning_enabled": feature_ri_norm is not None,
        "ri_coverage_threshold": float(RI_COVERAGE_THRESHOLD),
        "ri_normalization_min": float(ri_min_used) if ri_min_used is not None else None,
        "ri_normalization_max": float(ri_max_used) if ri_max_used is not None else None,
        "ri_value_summary": summarize_numeric_array(np.asarray(ri_finite, dtype=np.float32)),
        "feature_source_audit_json": str(feature_source_audit_path) if feature_source_audit_path is not None else None,
        "feature_source_provenance_summary": feature_source_audit_summary,
        "decoder_mass_conditioning_source": decoder_mass_source if not skip_constrained_denovo else None,
        "decoder_summary_json": (
            str(decoder_summary_path) if not skip_constrained_denovo and decoder_summary_payload is not None else None
        ),
    }
    ri_audit_path = output_dir / artifact_filename(artifact_prefix, "ri_conditioning_audit", ".json")
    write_json_with_aliases(
        ri_audit,
        ri_audit_path,
        [output_dir / "ri_conditioning_audit.json"],
    )
    isotope_audit_summary = build_isotope_audit_summary(
        artifact_prefix=artifact_prefix,
        html_isotope_df=html_isotope_df,
        retrieval_df=retrieval_df,
        denovo_df=denovo_df,
        merged_df=merged_df,
        html_isotope_csv_path=html_isotope_csv,
    )
    isotope_audit_summary_path = output_dir / artifact_filename(artifact_prefix, "isotope_audit_summary", ".json")
    write_json_with_aliases(
        isotope_audit_summary,
        isotope_audit_summary_path,
        [output_dir / "isotope_audit_summary.json"],
    )

    summary = {
        "artifact_prefix": artifact_prefix,
        "feature_table": str(feature_table_path),
        "encoder_checkpoint": str(encoder_ckpt),
        "encoder_calibrator": str(calibrator_path) if calibrator_path is not None else None,
        "decoder_dir": str(decoder_dir),
        "decoder_model_file": str(decoder_model_path) if decoder_model_path is not None else None,
        "decoder_config_file": str(decoder_config_path) if decoder_config_path is not None else None,
        "decoder_vocab_file": str(decoder_vocab_path) if decoder_vocab_path is not None else None,
        "skip_constrained_denovo": bool(skip_constrained_denovo),
        "embedding_bank_dir": str(embedding_bank_dir),
        "retrieval_smiles_csv": str(retrieval_smiles_csv),
        "retrieval_mz_pkl": str(retrieval_mz_pkl),
        "fragment_db_csv": str(fragment_db_csv),
        "prevalence_csv": str(prevalence_csv) if prevalence_csv.exists() else None,
        "prevalence_rows_loaded": int(len(prevalence_map)),
        "retrieval_rescoring": str(args.retrieval_rescoring),
        "retrieval_weights_json": (
            str(Path(args.retrieval_weights_json).expanduser().resolve())
            if str(args.retrieval_weights_json).strip()
            else None
        ),
        "retrieval_scoring_config": (
            retrieval_scoring_config.to_dict() if retrieval_scoring_config is not None else None
        ),
        "merge_policy": merge_policy_name,
        "merge_policy_joint_sort_used": bool(merge_policy_joint_sort_used),
        "merge_policy_selfies_cap_per_feature": int(merge_policy_selfies_cap),
        "merge_policy_selfies_pool_rows": int(merge_policy_selfies_pool_rows),
        "merge_policy_selfies_pool_mass_window_da": float(args.selfies_policy_mass_window_da),
        "merge_policy_selfies_score_formula": SELFIES_SCORE_FORMULA_VERSION,
        "ion_mode": ion_mode,
        "spectrum_bin_minmass": float(args.minmass),
        "spectrum_bin_maxmass": float(args.maxmass),
        "parent_mz_min": parent_mz_min,
        "parent_mz_max": parent_mz_max,
        "mass_conditioning": {
            "normalization_min": float(mass_condition_min) if not skip_constrained_denovo else None,
            "normalization_max": float(mass_condition_max) if not skip_constrained_denovo else None,
            "normalization_source": decoder_mass_source if not skip_constrained_denovo else None,
            "decoder_summary_json": (
                str(decoder_summary_path)
                if not skip_constrained_denovo and decoder_summary_payload is not None
                else None
            ),
        },
        "preprocessing": preprocessing_config.to_dict(),
        "html_formula_mode": effective_html_formula_mode,
        "candidate_formula_max": int(args.candidate_formula_max),
        "fragment_match_tolerance_da": float(args.fragment_match_tolerance_da),
        "fragment_match_target": float(args.fragment_match_target),
        "structural_neighbor_top_n": int(args.structural_neighbor_top_n),
        "structural_min_mcs_atoms": int(args.structural_min_mcs_atoms),
        "structural_mcs_timeout_sec": int(args.structural_mcs_timeout_sec),
        "num_features_encoded": int(feature_meta_df.shape[0]),
        "num_unique_raw_candidates": int(len(candidate_smiles_set)),
        "retrieval_report_rows": int(len(retrieval_df)),
        "denovo_report_rows": int(len(denovo_df)),
        "merged_report_rows": int(len(merged_df)),
        "features_with_recommended_analog": rec_count,
        "features_with_denovo_candidate": features_with_denovo,
        "retrieval_tier_distribution": tier_counts.to_dict(),
        "feature_latent_norm_pre_calibration": summarize_numeric_array(latent_norms_pre),
        "feature_latent_norm_post_calibration": summarize_numeric_array(latent_norms_post),
        "retrieval_rank1_latent_cosine_summary": summarize_numeric_array(retrieval_rank1_cos.to_numpy()),
        "retrieval_rank1_score_summary": summarize_numeric_array(retrieval_rank1_score.to_numpy()),
        "isotope_audit": isotope_audit_summary,
        "ri_conditioning": {
            "enabled": feature_ri_norm is not None,
            "requested": not bool(args.disable_ri_conditioning),
            "auto_disabled": bool(ri_auto_disabled),
            "coverage_fraction": float(ri_coverage),
            "finite_positive_count": int(len(ri_finite)),
            "feature_count": int(len(ri_vals)),
            "threshold": float(RI_COVERAGE_THRESHOLD),
            "normalization_min": float(ri_min_used) if ri_min_used is not None else None,
            "normalization_max": float(ri_max_used) if ri_max_used is not None else None,
            "audit_json": str(ri_audit_path),
        },
        "authoritative_track_default": authoritative_track_default,
        "outputs": {
            "retrieval_report": str(retrieval_csv_path),
            "denovo_report": str(denovo_csv_path),
            "candidate_metrics_all": str(candidate_metrics_csv_path),
            "merged_predictions": str(merged_csv_path),
            "generation_samples_raw": str(raw_samples_csv),
            "feature_latents": str(feature_latents_pt),
            "feature_meta": str(feature_meta_csv),
            "alignment_report": str(alignment_report_path),
            "html_formula_signals": str(html_feature_csv) if html_feature_csv is not None else None,
            "html_isotope_signals": str(html_isotope_csv),
            "isotope_audit_summary": str(isotope_audit_summary_path),
            "ri_conditioning_audit": str(ri_audit_path),
        },
        "postrun_decoder_mass_diagnosis": diagnosis_result,
    }
    summary_path = output_dir / artifact_filename(artifact_prefix, "run_summary", ".json")
    write_json_with_aliases(
        summary,
        summary_path,
        [output_dir / "run_summary.json"],
    )

    print(f"\nWrote: {feature_latents_pt}")
    print(f"Wrote: {feature_meta_csv}")
    print(f"Wrote: {raw_samples_csv}")
    print(f"Wrote: {retrieval_csv_path}")
    print(f"Wrote: {denovo_csv_path}")
    print(f"Wrote: {candidate_metrics_csv_path}")
    print(f"Wrote: {merged_csv_path}")
    print(f"Wrote: {alignment_report_path}")
    print(f"Wrote: {html_isotope_csv}")
    print(f"Wrote: {isotope_audit_summary_path}")
    print(f"Wrote: {summary_path}")
    if diagnosis_output_dir is not None:
        print(f"Wrote: {diagnosis_output_dir}")

    if (
        bool(args.run_decoder_mass_diagnosis)
        and bool(args.enforce_decoder_mass_gate)
        and int(diagnosis_result["return_code"] or 0) != 0
    ):
        raise SystemExit(int(diagnosis_result["return_code"] or 1))


if __name__ == "__main__":
    main()
