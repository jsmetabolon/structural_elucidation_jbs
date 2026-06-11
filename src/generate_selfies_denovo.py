#!/usr/bin/env python3
"""Generate de novo candidates from trained SELFIES decoder on X209 features."""

from __future__ import annotations

import argparse
import json
import random
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import torch

from denovo_improvement_utils import (
    FORMULA_ELEMENTS,
    adduct_vector,
    load_formula_feature_map,
    zero_formula_vector,
)
from device_utils import DEVICE_HELP_TEXT, resolve_device
from latent_calibration import apply_latent_calibrator, load_latent_calibrator
from recursive_substructure_utils import load_external_embedding_map
from model1D2Conv import (
    Net1D,
    infer_mass_conditioned_compatible_length,
    load_compatible_state_dict,
)
from selfies_decoder_model import SelfiesDecoder, SelfiesVocab
from substructure_search import (
    build_feature_consensus_smarts_from_retrieval,
    smiles_matches_substructure,
    tanimoto_similarity_smiles,
)
from spectrum_preprocessing import (
    SpectrumPreprocessingConfig,
    preprocess_and_vectorize,
)
from phase0_utils import (
    adduct_neutral_mass_error,
    best_adduct_neutral_mass_match,
    exact_mol_weight,
    molecular_formula,
    neutral_mass_from_precursor,
    normalize_ion_mode,
    strict_tolerance_da,
)

try:
    from merged_model import (
        MergedNet1D,
        ion_mode_label_to_scalar,
        lib_id_to_index,
        load_compatible_state_dict_for_merged,
    )
except ImportError:
    MergedNet1D = None
    ion_mode_label_to_scalar = None
    lib_id_to_index = None
    load_compatible_state_dict_for_merged = None


def robust_torch_load(path: Path, map_location: str = "cpu"):
    original_rebuild = getattr(torch._utils, "_rebuild_device_tensor_from_numpy", None)

    def load() -> object:
        try:
            return torch.load(path, map_location=map_location, weights_only=False)
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


def resolve_generation_device(raw_device: str) -> torch.device:
    device = resolve_device(raw_device)
    if device.type == "mps":
        print(
            "MPS selected for SELFIES generation, but torch.multinomial is not supported "
            "reliably here; falling back to CPU for decoder sampling."
        )
        return torch.device("cpu")
    return device


def read_table(path: Path) -> pd.DataFrame:
    if path.suffix.lower() in {".pkl", ".pickle"}:
        df = pd.read_pickle(path)
    else:
        try:
            df = pd.read_csv(path, sep=None, engine="python")
        except Exception:
            df = pd.read_csv(path)
    df.columns = [
        str(col).strip().strip('"').lstrip("\ufeff") for col in df.columns
    ]
    return df


def choose_encoder_checkpoint(path: Path) -> Path:
    if path.exists():
        return path.resolve()
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
    raise FileNotFoundError(
        f"Encoder checkpoint not found at {path} and no fallback checkpoint exists."
    )


def detect_encoder_kind(
    state: dict[str, object],
    requested_kind: str,
) -> str:
    requested = str(requested_kind).strip().lower()
    if requested in {"simple", "merged"}:
        return requested

    if any(
        str(key).startswith("mode_embedding.") or str(key).startswith("lib_embedding.")
        for key in state.keys()
    ):
        return "merged"
    return "simple"


def load_encoder_checkpoint(
    ckpt_path: Path,
    n_bins: int,
    minmass: float,
    maxmass: float,
    parent_mz_min: float,
    parent_mz_max: float,
    device: torch.device,
    requested_kind: str = "auto",
) -> tuple[Net1D | MergedNet1D, int, dict[str, object], str]:
    state = robust_torch_load(ckpt_path, map_location="cpu")
    if not isinstance(state, dict):
        raise RuntimeError(
            f"Expected state_dict in {ckpt_path}, got {type(state)}"
        )
    encoder_kind = detect_encoder_kind(state, requested_kind)
    model_n_bins = int(
        infer_mass_conditioned_compatible_length(state, preferred_length_in=n_bins) or n_bins
    )
    if encoder_kind == "merged":
        if MergedNet1D is None or load_compatible_state_dict_for_merged is None:
            raise RuntimeError("Merged encoder support is unavailable in this runtime.")
        model = MergedNet1D(
            length_in=model_n_bins,
            minmass=minmass,
            maxmass=maxmass,
            apply_tanh=True,
            parent_mz_min=parent_mz_min,
            parent_mz_max=parent_mz_max,
        ).to(device)
        load_info = load_compatible_state_dict_for_merged(
            model,
            state,
            ignore_new_conditioning_params=False,
        )
    else:
        model = Net1D(
            length_in=model_n_bins,
            minmass=minmass,
            maxmass=maxmass,
            apply_tanh=True,
            parent_mz_min=parent_mz_min,
            parent_mz_max=parent_mz_max,
        ).to(device)
        load_info = load_compatible_state_dict(model, state)
    model.eval()
    return model, model_n_bins, load_info, encoder_kind


def align_vector_length(vector: np.ndarray, target_n_bins: int) -> np.ndarray:
    current_n_bins = int(vector.shape[0])
    if current_n_bins == int(target_n_bins):
        return vector
    if current_n_bins > int(target_n_bins):
        return vector[: int(target_n_bins)]
    out = np.zeros(int(target_n_bins), dtype=vector.dtype)
    out[:current_n_bins] = vector
    return out


def load_feature_table(path: Path) -> pd.DataFrame:
    df = read_table(path)
    expected = {
        "feature_key",
        "comp_id",
        "compound_name",
        "rt",
        "parent_mz",
        "mass",
        "intensity",
    }
    missing = sorted(expected - set(df.columns))
    if missing:
        raise KeyError(f"Feature table missing columns: {missing}")
    out = df.copy()
    if "source_html" not in out.columns:
        out["source_html"] = ""
    out["feature_key"] = out["feature_key"].astype(str)
    out["comp_id"] = (
        pd.to_numeric(out["comp_id"], errors="coerce").fillna(-1).astype(int)
    )
    out["compound_name"] = out["compound_name"].astype(str)
    out["source_html"] = out["source_html"].fillna("").astype(str)
    out["parent_mz"] = pd.to_numeric(out["parent_mz"], errors="coerce")
    out["mass"] = pd.to_numeric(out["mass"], errors="coerce")
    out["intensity"] = pd.to_numeric(out["intensity"], errors="coerce")
    out = out.dropna(subset=["parent_mz", "mass", "intensity"])
    return out


def _normalize_mode_token(value: object) -> str:
    if value is None:
        text = ""
    elif isinstance(value, float) and np.isnan(value):
        text = ""
    else:
        text = str(value).strip().lower()
    if text in {"negative", "neg", "-", "0", "0.0"}:
        return "negative"
    if text in {"positive", "pos", "+", "1", "1.0"}:
        return "positive"
    return text


def filter_feature_table_scope(
    features_df: pd.DataFrame,
    *,
    lib_id: int | None,
    ion_mode: str,
    require_scope: bool = True,
) -> tuple[pd.DataFrame, dict[str, object]]:
    original_rows = int(features_df.shape[0])
    original_features = int(features_df["feature_key"].nunique())
    work = features_df.copy()
    filters_applied: list[str] = []
    filters_missing: list[str] = []

    if lib_id is not None:
        if "lib_id" in work.columns:
            lib_values = pd.to_numeric(work["lib_id"], errors="coerce")
            work = work[lib_values.eq(int(lib_id))].copy()
            filters_applied.append("lib_id")
        else:
            filters_missing.append("lib_id")

    mode_col = next(
        (
            col
            for col in ("ion_mode", "ion_mode_label", "mode")
            if col in work.columns
        ),
        None,
    )
    if mode_col is not None:
        requested_mode = _normalize_mode_token(ion_mode)
        mode_values = work[mode_col].map(_normalize_mode_token)
        work = work[mode_values.eq(requested_mode)].copy()
        filters_applied.append(str(mode_col))
    else:
        filters_missing.append("ion_mode")

    if require_scope and filters_missing:
        raise KeyError(
            "Feature table is missing required per-lib/per-mode scope columns: "
            f"missing={filters_missing}, lib_id={lib_id}, ion_mode={ion_mode}, "
            f"available_columns={list(features_df.columns)}"
        )

    if work.empty:
        raise ValueError(
            "Feature table scope filter removed all rows: "
            f"lib_id={lib_id}, ion_mode={ion_mode}, "
            f"original_rows={original_rows}, original_features={original_features}, "
            f"filters_applied={filters_applied}, filters_missing={filters_missing}"
        )

    summary = {
        "requested_lib_id": int(lib_id) if lib_id is not None else None,
        "requested_ion_mode": str(ion_mode),
        "input_rows": original_rows,
        "input_features": original_features,
        "output_rows": int(work.shape[0]),
        "output_features": int(work["feature_key"].nunique()),
        "filters_applied": filters_applied,
        "filters_missing": filters_missing,
    }
    return work, summary


def encode_feature_latents(
    features_df: pd.DataFrame,
    encoder: Net1D | MergedNet1D,
    minmass: float,
    maxmass: float,
    resolution: int,
    device: torch.device,
    preprocessing_config: SpectrumPreprocessingConfig,
    encoder_n_bins: int | None = None,
    encoder_kind: str = "simple",
    encoder_ion_mode: float | None = None,
    encoder_lib_id: int | None = None,
) -> tuple[pd.DataFrame, torch.Tensor]:
    rows: list[dict[str, object]] = []
    latents: list[torch.Tensor] = []

    grouped = features_df.groupby("feature_key", sort=True)
    for feature_key, sub in grouped:
        peaks = list(zip(sub["mass"].tolist(), sub["intensity"].tolist()))
        parent_mz = float(sub["parent_mz"].iloc[0])
        comp_id = int(sub["comp_id"].iloc[0])
        feature_label = str(sub["compound_name"].iloc[0])

        processed = preprocess_and_vectorize(
            peaks=peaks,
            config=preprocessing_config,
            minmass=minmass,
            maxmass=maxmass,
            resolution=resolution,
        )
        vector = processed["vector"]
        if encoder_n_bins is not None:
            vector = align_vector_length(vector, encoder_n_bins)
        x = torch.from_numpy(vector).unsqueeze(0).unsqueeze(0).to(device)
        pmz = torch.tensor([parent_mz], dtype=torch.float32, device=device)

        with torch.no_grad():
            if str(encoder_kind) == "merged":
                z = encoder(
                    x,
                    parent_mz=pmz,
                    ion_mode=encoder_ion_mode,
                    lib_id=encoder_lib_id,
                ).squeeze(0).detach().cpu().float()
            else:
                z = encoder(x, parent_mz=pmz).squeeze(0).detach().cpu().float()

        rows.append(
            {
                "feature_key": str(feature_key),
                "comp_id": int(comp_id),
                "feature_label": feature_label,
                "rt": float(sub["rt"].iloc[0]) if "rt" in sub.columns else float("nan"),
                "parent_mz": float(parent_mz),
                "n_raw_peaks": int(processed["n_raw_peaks"]),
                "n_retained_peaks": int(processed["n_retained_peaks"]),
                "n_vectorized_peaks": int(processed["n_vectorized_peaks"]),
                "n_out_of_range_retained_peaks": int(
                    processed["n_out_of_range_retained_peaks"]
                ),
                "dropped_out_of_range_peak_fraction": float(
                    processed["dropped_out_of_range_peak_fraction"]
                ),
            }
        )
        latents.append(z)

    meta_df = (
        pd.DataFrame(rows).sort_values("feature_key").reset_index(drop=True)
    )
    latent_tensor = (
        torch.stack(latents, dim=0) if latents else torch.empty((0, 512))
    )
    return meta_df, latent_tensor


def build_temperature_schedule(
    num_samples: int, temps: list[float]
) -> list[float]:
    if not temps:
        return [1.0] * int(num_samples)
    per_temp = max(1, int(num_samples) // len(temps))
    schedule: list[float] = []
    for temp in temps:
        schedule.extend([float(temp)] * per_temp)
    while len(schedule) < int(num_samples):
        schedule.append(float(temps[-1]))
    return schedule[: int(num_samples)]


def mass_match_score(abs_err_da: float, tolerance_da: float = 5.0) -> float:
    if not np.isfinite(abs_err_da):
        return 0.0
    if abs_err_da <= 0.001:
        return 1.0
    if abs_err_da >= float(tolerance_da):
        return 0.0
    return float(max(0.0, 1.0 - abs_err_da / float(tolerance_da)))


def mass_error_bin(abs_err_da: float, strict_pass: bool) -> int:
    if bool(strict_pass):
        return 0
    if float(abs_err_da) <= 0.5:
        return 1
    if float(abs_err_da) <= 1.0:
        return 2
    if float(abs_err_da) <= 2.0:
        return 3
    return 4


def artifact_filename(prefix: str, stem: str, suffix: str) -> str:
    normalized = str(prefix).strip()
    if not normalized:
        return f"{stem}{suffix}"
    return f"{normalized}_{stem}{suffix}"


def write_dataframe_with_aliases(
    df: pd.DataFrame, primary_path: Path, alias_paths: Iterable[Path]
) -> None:
    df.to_csv(primary_path, index=False)
    for alias_path in alias_paths:
        if alias_path == primary_path:
            continue
        df.to_csv(alias_path, index=False)


def write_json_with_aliases(
    payload: dict[str, object], primary_path: Path, alias_paths: Iterable[Path]
) -> None:
    text = json.dumps(payload, indent=2)
    primary_path.write_text(text)
    for alias_path in alias_paths:
        if alias_path == primary_path:
            continue
        alias_path.write_text(text)


def write_tensor_with_aliases(
    payload: object, primary_path: Path, alias_paths: Iterable[Path]
) -> None:
    torch.save(payload, primary_path)
    for alias_path in alias_paths:
        if alias_path == primary_path:
            continue
        torch.save(payload, alias_path)


def summarize_numeric_array(values: np.ndarray) -> dict[str, float | int]:
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


def load_selfies_decoder(
    model_ckpt: Path,
    vocab_override: Path | None,
    device: torch.device,
) -> tuple[SelfiesDecoder, SelfiesVocab, dict[str, object]]:
    payload = robust_torch_load(model_ckpt, map_location="cpu")
    if not isinstance(payload, dict):
        raise RuntimeError(
            f"Unexpected checkpoint payload type: {type(payload)}"
        )

    model_config = dict(payload.get("model_config", {}))
    if not model_config:
        raise RuntimeError("Checkpoint missing model_config.")
    state_dict = payload.get("model_state_dict")
    if not isinstance(state_dict, dict):
        raise RuntimeError("Checkpoint missing model_state_dict.")

    vocab_path = (
        vocab_override
        if vocab_override is not None
        else Path(str(payload.get("vocab_json", "")))
    )
    if vocab_path is None or not str(vocab_path):
        raise RuntimeError(
            "No vocabulary path provided and checkpoint has no vocab_json."
        )
    vocab = SelfiesVocab.load_json(Path(vocab_path).expanduser().resolve())

    decoder_model_config = {
        key: value for key, value in model_config.items() if key != "conditioning_version"
    }
    model = SelfiesDecoder(vocab=vocab, **decoder_model_config).to(device)
    incompatible = model.load_state_dict(state_dict, strict=False)
    missing = list(getattr(incompatible, "missing_keys", []))
    unexpected = list(getattr(incompatible, "unexpected_keys", []))
    if missing:
        print(f"Model load missing keys: {missing}")
    if unexpected:
        print(f"Model load unexpected keys: {unexpected}")
    model.eval()
    return model, vocab, payload


def run_control_eval(
    repo_root: Path,
    output_dir: Path,
    denovo_csv: Path,
    control_lib_id: int,
) -> dict[str, object]:
    script = repo_root / "evaluate_denovo_controls.py"
    cmd = [
        sys.executable,
        str(script),
        "--prediction_dir",
        str(output_dir),
        "--denovo_csv",
        str(denovo_csv),
        "--output_dir",
        str(output_dir),
        "--control_lib_id",
        str(int(control_lib_id)),
    ]
    print("Running control evaluation:", " ".join(cmd))
    completed = subprocess.run(cmd, check=False)
    summary_path = output_dir / "control_eval_summary.json"
    return {
        "command": cmd,
        "return_code": int(completed.returncode),
        "summary_json": str(summary_path),
        "summary_exists": bool(summary_path.exists()),
    }


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Generate X209 de novo candidates from a trained SELFIES decoder."
    )
    parser.add_argument(
        "--model_ckpt",
        type=str,
        default="models/selfies_decoder_209/selfies_decoder_209.pt",
    )
    parser.add_argument(
        "--vocab_json",
        type=str,
        default="",
        help="Optional override vocabulary JSON.",
    )
    parser.add_argument(
        "--feature_table", type=str, default="models/x209_features_encoder.pkl"
    )
    parser.add_argument(
        "--encoder_ckpt",
        type=str,
        default="models/spectra_encoder_209_trueZ_mass_filtered_ft.pt",
    )
    parser.add_argument(
        "--encoder_calibrator",
        type=str,
        default="",
        help="Optional affine latent calibrator checkpoint applied after encoder inference.",
    )
    parser.add_argument(
        "--encoder_kind",
        type=str,
        choices=("auto", "simple", "merged"),
        default="auto",
        help="Encoder architecture. 'auto' detects merged checkpoints by conditioning layers.",
    )
    parser.add_argument(
        "--lib_id",
        type=int,
        default=None,
        help="Library id conditioning for merged encoders (209, 305, 400, or 402).",
    )
    parser.add_argument(
        "--encoder_ion_mode",
        type=float,
        default=None,
        help="Optional numeric ion-mode conditioning override for merged encoders.",
    )
    parser.add_argument(
        "--output_dir", type=str, default="models/x209_selfies_denovo"
    )
    parser.add_argument(
        "--artifact_prefix",
        type=str,
        default="x209",
        help="Filename prefix for denovo feature/report artifacts written to output_dir.",
    )
    parser.add_argument("--num_samples", type=int, default=60)
    parser.add_argument(
        "--target_mass_matched_per_feature",
        type=int,
        default=0,
        help=(
            "Optional adaptive decode target. When >0, keep sampling features "
            "until this many unique candidates pass --mass_tolerance_da or "
            "--max_samples_per_feature is reached."
        ),
    )
    parser.add_argument(
        "--max_samples_per_feature",
        type=int,
        default=0,
        help=(
            "Maximum decode draws per feature when adaptive mass sampling is enabled. "
            "0 uses 5x --num_samples."
        ),
    )
    parser.add_argument(
        "--ranking_mode",
        type=str,
        choices=("mass_first", "mass_bin_support"),
        default="mass_first",
        help=(
            "Candidate ranking. mass_first preserves the legacy exact-mass-first "
            "sort; mass_bin_support ranks by strict/0.5/1/2/5 Da bins, then "
            "generation support and log probability."
        ),
    )
    parser.add_argument(
        "--decode_temps", type=str, default="0.3,0.5,0.8,1.0,1.2"
    )
    parser.add_argument(
        "--top_p",
        type=float,
        default=0.0,
        help="Optional nucleus sampling threshold for standalone SELFIES decoding (0 disables).",
    )
    parser.add_argument("--decode_batch_size", type=int, default=128)
    parser.add_argument("--decode_max_len", type=int, default=120)
    parser.add_argument("--mass_tolerance_da", type=float, default=5.0)
    parser.add_argument(
        "--disable_adduct_hypotheses",
        action="store_true",
        help=(
            "Use only the legacy primary-adduct neutral mass. By default, "
            "generation ranking uses common adduct/isotope neutral-mass hypotheses."
        ),
    )
    parser.add_argument(
        "--allow_missing_feature_scope",
        action="store_true",
        help="Allow missing lib_id/mode columns in feature_table. Intended only for legacy debugging.",
    )
    parser.add_argument(
        "--substructure_filter_mode",
        type=str,
        choices=("off", "consensus_substructure", "retrieval_tanimoto"),
        default="off",
        help="Optional post-generation filter driven by retrieval structure priors.",
    )
    parser.add_argument(
        "--retrieval_report_csv",
        type=str,
        default="",
        help="Optional retrieval report used to derive consensus substructure constraints per feature.",
    )
    parser.add_argument("--substructure_consensus_top_n", type=int, default=5)
    parser.add_argument("--substructure_min_mcs_atoms", type=int, default=6)
    parser.add_argument("--substructure_mcs_timeout_sec", type=int, default=3)
    parser.add_argument("--retrieval_tanimoto_threshold", type=float, default=0.5)
    parser.add_argument("--control_lib_id", type=int, default=209)
    parser.add_argument(
        "--seed",
        type=int,
        default=-1,
        help="Optional RNG seed for deterministic decoding. Negative values leave sampling unseeded.",
    )
    parser.add_argument("--skip_control_eval", action="store_true")
    parser.add_argument(
        "--max_features",
        type=int,
        default=0,
        help="Optional limit for smoke runs.",
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
    )
    parser.add_argument("--minmass", type=float, default=50.0)
    parser.add_argument("--maxmass", type=float, default=499.97)
    parser.add_argument("--parent_mz_min", type=float, default=50.0)
    parser.add_argument("--parent_mz_max", type=float, default=499.97)
    parser.add_argument("--resolution", type=int, default=2)
    parser.add_argument("--min_rel_intensity_pct", type=float, default=2.0)
    parser.add_argument("--top_peaks", type=int, default=25)
    parser.add_argument("--disable_basepeak_norm", action="store_true")
    parser.add_argument(
        "--disable_ri_conditioning",
        action="store_true",
        help="Compatibility flag for sweep runners. The standalone SELFIES decoder does not support RI conditioning.",
    )
    parser.add_argument(
        "--conditioning_version",
        type=str,
        choices=("auto", "mass_v0", "formula_v1", "substructure_v1", "biotransform_v1"),
        default="auto",
        help="Expected decoder conditioning version. auto uses checkpoint metadata.",
    )
    parser.add_argument(
        "--formula_csv",
        type=str,
        default="",
        help="Optional feature_key or comp_id keyed formula metadata for formula_v1 generation.",
    )
    parser.add_argument(
        "--adduct_conditioning",
        action="store_true",
        help="Compatibility flag; adduct conditioning is used when the checkpoint declares support.",
    )
    parser.add_argument(
        "--substructure_embedding_csv",
        type=str,
        default="",
        help="Optional feature_key/comp_id keyed external embedding index CSV for substructure_v1 or biotransform_v1.",
    )
    parser.add_argument(
        "--substructure_embedding_pt",
        type=str,
        default="",
        help="Optional .npy/.pt matrix paired with --substructure_embedding_csv.",
    )
    return parser


def main() -> None:
    args = get_parser().parse_args()
    device = resolve_generation_device(args.device)
    ion_mode = normalize_ion_mode(str(args.ion_mode))
    if not bool(args.disable_ri_conditioning):
        print("RI conditioning requested when supported by checkpoint config.")

    if int(args.seed) >= 0:
        seed = int(args.seed)
        random.seed(seed)
        np.random.seed(seed)
        torch.manual_seed(seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(seed)
        print(f"Using decode seed: {seed}")

    repo_root = Path(__file__).resolve().parent
    model_ckpt = Path(args.model_ckpt).expanduser().resolve()
    vocab_override = (
        Path(args.vocab_json).expanduser().resolve()
        if str(args.vocab_json).strip()
        else None
    )
    feature_table = Path(args.feature_table).expanduser().resolve()
    output_dir = Path(args.output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    artifact_prefix = str(args.artifact_prefix).strip() or "x209"

    model, _, checkpoint_payload = load_selfies_decoder(
        model_ckpt=model_ckpt,
        vocab_override=vocab_override,
        device=device,
    )
    model_config_payload = dict(checkpoint_payload.get("model_config") or {})
    checkpoint_ri_supported = bool(model_config_payload.get("use_ri_conditioning", False))
    ri_conditioning_effective = bool(checkpoint_ri_supported and not bool(args.disable_ri_conditioning))
    checkpoint_conditioning_version = str(
        model_config_payload.get("conditioning_version", "mass_v0")
    )
    requested_conditioning_version = str(args.conditioning_version)
    if requested_conditioning_version != "auto" and requested_conditioning_version != checkpoint_conditioning_version:
        raise RuntimeError(
            "Requested conditioning version does not match checkpoint metadata: "
            f"requested={requested_conditioning_version}, checkpoint={checkpoint_conditioning_version}"
        )
    formula_conditioning_effective = bool(
        model_config_payload.get("use_formula_conditioning", False)
    )
    adduct_conditioning_effective = bool(
        model_config_payload.get("use_adduct_conditioning", False)
    )
    substructure_conditioning_effective = bool(
        model_config_payload.get("use_substructure_conditioning", False)
    )
    if not bool(args.disable_ri_conditioning) and not checkpoint_ri_supported:
        print("Checkpoint config does not declare RI support; continuing without RI conditioning.")
    mass_norm_payload = checkpoint_payload.get("mass_normalization", {})
    train_mass_min = float(mass_norm_payload.get("min", 50.0))
    train_mass_max = float(mass_norm_payload.get("max", 500.0))
    print(f"Loaded decoder checkpoint: {model_ckpt}")
    print(
        f"Mass normalization from checkpoint: min={train_mass_min:.6f}, max={train_mass_max:.6f}"
    )
    formula_csv = (
        Path(args.formula_csv).expanduser().resolve()
        if str(args.formula_csv).strip()
        else None
    )
    substructure_embedding_csv = (
        Path(args.substructure_embedding_csv).expanduser().resolve()
        if str(args.substructure_embedding_csv).strip()
        else None
    )
    substructure_embedding_pt = (
        Path(args.substructure_embedding_pt).expanduser().resolve()
        if str(args.substructure_embedding_pt).strip()
        else None
    )
    formula_feature_map: dict[str, dict[str, object]] = {}
    if formula_conditioning_effective and formula_csv is not None:
        formula_feature_map.update(
            load_formula_feature_map(formula_csv, ["feature_key", "comp_id", "CompID"])
        )
    if formula_conditioning_effective:
        print(
            "Formula conditioning:",
            f"version={checkpoint_conditioning_version}",
            f"formula_csv={formula_csv if formula_csv is not None else 'none'}",
            f"metadata_rows={len(formula_feature_map)}",
        )
    substructure_embedding_map: dict[str, np.ndarray] = {}
    substructure_input_dim = int(model_config_payload.get("substructure_input_dim", 0) or 0)
    if substructure_conditioning_effective:
        if substructure_embedding_csv is None:
            raise FileNotFoundError(
                "Checkpoint requires substructure conditioning; pass --substructure_embedding_csv."
            )
        substructure_embedding_map, loaded_substructure_dim = load_external_embedding_map(
            substructure_embedding_csv, substructure_embedding_pt
        )
        if int(loaded_substructure_dim) != int(substructure_input_dim):
            raise RuntimeError(
                "Substructure embedding dimension mismatch: "
                f"checkpoint={substructure_input_dim}, supplied={loaded_substructure_dim}"
            )
        print(
            "External vector conditioning:",
            f"embedding_csv={substructure_embedding_csv}",
            f"embedding_matrix={substructure_embedding_pt if substructure_embedding_pt is not None else 'auto'}",
            f"metadata_rows={len(substructure_embedding_map)}",
            f"dim={substructure_input_dim}",
        )

    features_df = load_feature_table(feature_table)
    features_df, feature_scope_summary = filter_feature_table_scope(
        features_df,
        lib_id=args.lib_id,
        ion_mode=ion_mode,
        require_scope=not bool(args.allow_missing_feature_scope),
    )
    print(
        "Feature scope:",
        f"input_features={feature_scope_summary['input_features']}",
        f"output_features={feature_scope_summary['output_features']}",
        f"filters_applied={feature_scope_summary['filters_applied']}",
    )
    if int(args.max_features) > 0:
        keep_keys = sorted(features_df["feature_key"].unique())[
            : int(args.max_features)
        ]
        features_df = features_df[
            features_df["feature_key"].isin(keep_keys)
        ].copy()
        print(f"Limiting to first {len(keep_keys)} features for smoke run.")

    n_bins = (
        int((float(args.maxmass) - float(args.minmass)) * int(args.resolution))
        + 1
    )
    encoder_ckpt = choose_encoder_checkpoint(
        Path(args.encoder_ckpt).expanduser()
    )
    encoder, encoder_n_bins, encoder_load_info, encoder_kind = load_encoder_checkpoint(
        ckpt_path=encoder_ckpt,
        n_bins=n_bins,
        minmass=float(args.minmass),
        maxmass=float(args.maxmass),
        parent_mz_min=float(args.parent_mz_min),
        parent_mz_max=float(args.parent_mz_max),
        device=device,
        requested_kind=str(args.encoder_kind),
    )
    encoder_ion_mode = None
    encoder_lib_id = None
    if str(encoder_kind) == "merged":
        raw_encoder_mode = args.encoder_ion_mode
        if raw_encoder_mode is None:
            if ion_mode_label_to_scalar is None:
                raise RuntimeError("Merged encoder ion-mode utilities are unavailable.")
            raw_encoder_mode = ion_mode_label_to_scalar(ion_mode)
        encoder_ion_mode = float(raw_encoder_mode)

        if args.lib_id is None:
            raise ValueError("--lib_id is required when using a merged encoder checkpoint.")
        if lib_id_to_index is None:
            raise RuntimeError("Merged encoder lib-id utilities are unavailable.")
        encoder_lib_id = int(args.lib_id)
        _ = lib_id_to_index(encoder_lib_id)

        print(
            "Using merged encoder conditioning:",
            f"lib_id={encoder_lib_id}",
            f"encoder_ion_mode={encoder_ion_mode:.1f}",
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
        calibrator_path = (
            Path(args.encoder_calibrator).expanduser().resolve()
        )
        if not calibrator_path.exists():
            raise FileNotFoundError(
                f"Encoder calibrator not found: {calibrator_path}"
            )
        calibrator = load_latent_calibrator(
            calibrator_path, device=device
        )
        print(f"Using encoder calibrator: {calibrator_path}")

    preprocessing_config = SpectrumPreprocessingConfig(
        min_rel_intensity_pct=float(args.min_rel_intensity_pct),
        top_peaks=int(args.top_peaks),
        disable_basepeak_norm=bool(args.disable_basepeak_norm),
    )

    feature_meta_df, feature_latents_raw = encode_feature_latents(
        features_df=features_df,
        encoder=encoder,
        minmass=float(args.minmass),
        maxmass=float(args.maxmass),
        resolution=int(args.resolution),
        device=device,
        preprocessing_config=preprocessing_config,
        encoder_n_bins=encoder_n_bins,
        encoder_kind=encoder_kind,
        encoder_ion_mode=encoder_ion_mode,
        encoder_lib_id=encoder_lib_id,
    )
    feature_latents = (
        apply_latent_calibrator(
            feature_latents_raw.to(device), calibrator
        )
        .detach()
        .cpu()
    )
    latent_norms_pre = (
        feature_latents_raw.norm(dim=1).cpu().numpy()
        if len(feature_latents_raw)
        else np.array([])
    )
    latent_norms_post = (
        feature_latents.norm(dim=1).cpu().numpy()
        if len(feature_latents)
        else np.array([])
    )
    feature_latents = feature_latents.to(device)
    print(f"Encoded features: {len(feature_meta_df)}")

    feature_keys = feature_meta_df["feature_key"].tolist()
    feature_index_by_key = {str(key): idx for idx, key in enumerate(feature_keys)}
    parent_mzs = feature_meta_df["parent_mz"].astype(float).tolist()
    requested_filter_mode = str(args.substructure_filter_mode).strip().lower()
    retrieval_report_csv = (
        Path(args.retrieval_report_csv).expanduser().resolve()
        if str(args.retrieval_report_csv).strip()
        else None
    )
    consensus_smarts_by_feature: dict[str, str] = {}
    top_retrieval_smiles_by_feature: dict[str, str] = {}
    effective_filter_mode = requested_filter_mode
    if requested_filter_mode != "off":
        if retrieval_report_csv is None or not retrieval_report_csv.exists():
            print(
                "Substructure filtering requested but retrieval report is missing; "
                "disabling substructure filtering."
            )
            effective_filter_mode = "off"
        else:
            retrieval_df = read_table(retrieval_report_csv)
            consensus_smarts_by_feature, top_retrieval_smiles_by_feature = (
                build_feature_consensus_smarts_from_retrieval(
                    retrieval_df,
                    top_n=int(args.substructure_consensus_top_n),
                    min_mcs_atoms=int(args.substructure_min_mcs_atoms),
                    timeout_sec=int(args.substructure_mcs_timeout_sec),
                )
            )
            print(
                "Loaded retrieval constraints: "
                f"consensus_features={len(consensus_smarts_by_feature)}, "
                f"top_retrieval_features={len(top_retrieval_smiles_by_feature)}"
            )
    mass_conditions = [
        float(
            np.clip(
                (pmz - train_mass_min)
                / max(train_mass_max - train_mass_min, 1e-12),
                0.0,
                1.0,
            )
        )
        for pmz in parent_mzs
    ]
    mass_tensor = torch.tensor(
        mass_conditions, dtype=torch.float32, device=device
    ).unsqueeze(1)
    ri_tensor = None
    ri_norm_payload = dict(checkpoint_payload.get("ri_normalization", {}) or {})
    if ri_conditioning_effective:
        ri_min = float(ri_norm_payload.get("min", 0.0))
        ri_max = float(ri_norm_payload.get("max", 1.0))
        raw_ri_values = []
        for _, row in feature_meta_df.iterrows():
            value = row.get("ri", row.get("rt", np.nan))
            try:
                numeric = float(value)
            except Exception:
                numeric = float("nan")
            if not np.isfinite(numeric):
                numeric = ri_min
            raw_ri_values.append(
                float(
                    np.clip(
                        (numeric - ri_min) / max(ri_max - ri_min, 1e-12),
                        0.0,
                        1.0,
                    )
                )
            )
        ri_tensor = torch.tensor(
            raw_ri_values, dtype=torch.float32, device=device
        ).unsqueeze(1)
    formula_tensor = None
    formula_condition_sources: list[str] = ["disabled" for _ in feature_keys]
    formula_condition_formulas: list[str] = ["" for _ in feature_keys]
    formula_condition_nonzero: list[bool] = [False for _ in feature_keys]
    if formula_conditioning_effective:
        vectors: list[list[float]] = []
        formula_condition_sources = []
        formula_condition_formulas = []
        formula_condition_nonzero = []
        for _, row in feature_meta_df.iterrows():
            feature_key = str(row.get("feature_key", ""))
            comp_id_key = str(int(row.get("comp_id", -1)))
            payload = formula_feature_map.get(feature_key) or formula_feature_map.get(comp_id_key)
            if payload:
                vector = [float(value) for value in payload.get("vector", zero_formula_vector())]
                formula_text = str(payload.get("formula", ""))
                source = "formula_csv"
            else:
                vector = zero_formula_vector()
                formula_text = ""
                source = "missing_formula"
            vectors.append(vector)
            formula_condition_sources.append(source)
            formula_condition_formulas.append(formula_text)
            formula_condition_nonzero.append(any(abs(float(value)) > 0 for value in vector))
        formula_tensor = torch.tensor(vectors, dtype=torch.float32, device=device)
    adduct_tensor = None
    if adduct_conditioning_effective:
        adduct_values = [adduct_vector(ion_mode) for _ in feature_keys]
        adduct_tensor = torch.tensor(adduct_values, dtype=torch.float32, device=device)
    substructure_tensor = None
    substructure_condition_sources: list[str] = ["disabled" for _ in feature_keys]
    substructure_condition_nonzero: list[bool] = [False for _ in feature_keys]
    if substructure_conditioning_effective:
        vectors: list[np.ndarray] = []
        substructure_condition_sources = []
        substructure_condition_nonzero = []
        for _, row in feature_meta_df.iterrows():
            feature_key = str(row.get("feature_key", ""))
            comp_id_key = str(int(row.get("comp_id", -1)))
            vector = substructure_embedding_map.get(feature_key)
            if vector is None:
                vector = substructure_embedding_map.get(comp_id_key)
            if vector is None:
                arr = np.zeros(int(substructure_input_dim), dtype=np.float32)
                source = "missing_substructure_embedding"
            else:
                arr = np.asarray(vector, dtype=np.float32).reshape(-1)
                if arr.shape[0] != int(substructure_input_dim):
                    raise RuntimeError(
                        f"Substructure embedding dim mismatch for {feature_key}: "
                        f"got {arr.shape[0]}, expected {substructure_input_dim}"
                    )
                source = "substructure_embedding_csv"
            vectors.append(arr)
            substructure_condition_sources.append(source)
            substructure_condition_nonzero.append(bool(np.any(np.abs(arr) > 0)))
        substructure_tensor = torch.tensor(np.vstack(vectors), dtype=torch.float32, device=device)

    schedule = build_temperature_schedule(
        num_samples=int(args.num_samples),
        temps=[
            float(token.strip())
            for token in str(args.decode_temps).split(",")
            if token.strip()
        ],
    )
    print(f"Decode schedule size: {len(schedule)} samples per feature")

    generation_counts: dict[str, Counter] = {
        key: Counter() for key in feature_keys
    }
    generation_draws: dict[str, int] = {key: 0 for key in feature_keys}
    best_log_probs: dict[str, dict[str, float]] = defaultdict(dict)
    raw_rows: list[dict[str, object]] = []

    def record_decoded_batch(
        *,
        sample_idx: int,
        temp: float,
        batch_positions: list[int],
    ) -> None:
        if not batch_positions:
            return
        position_tensor = torch.tensor(
            batch_positions, dtype=torch.long, device=device
        )
        z_batch = feature_latents.index_select(0, position_tensor)
        mass_batch = mass_tensor.index_select(0, position_tensor)
        ri_batch = (
            ri_tensor.index_select(0, position_tensor)
            if ri_tensor is not None
            else None
        )
        formula_batch = (
            formula_tensor.index_select(0, position_tensor)
            if formula_tensor is not None
            else None
        )
        adduct_batch = (
            adduct_tensor.index_select(0, position_tensor)
            if adduct_tensor is not None
            else None
        )
        substructure_batch = (
            substructure_tensor.index_select(0, position_tensor)
            if substructure_tensor is not None
            else None
        )
        selfies_list, mean_lps = model.generate(
            z=z_batch,
            mass=mass_batch,
            ri=ri_batch,
            formula=formula_batch,
            adduct=adduct_batch,
            substructure=substructure_batch,
            temp=float(temp),
            max_len=int(args.decode_max_len),
            top_p=float(args.top_p),
            return_log_probs=True,
        )
        for local_idx, selfies_text in enumerate(selfies_list):
            feature_pos = int(batch_positions[local_idx])
            feature_key = feature_keys[feature_pos]
            mean_lp = float(mean_lps[local_idx])
            mass_condition = float(mass_batch[local_idx].item())
            formula_source = str(formula_condition_sources[feature_pos])
            formula_text = str(formula_condition_formulas[feature_pos])
            formula_nonzero = bool(formula_condition_nonzero[feature_pos])
            substructure_source = str(substructure_condition_sources[feature_pos])
            substructure_nonzero = bool(substructure_condition_nonzero[feature_pos])
            smiles = model.selfies_to_smiles(selfies_text)
            if smiles is None:
                raw_rows.append(
                    {
                        "sample_idx": int(sample_idx),
                        "feature_key": feature_key,
                        "selfies": selfies_text,
                        "smiles": "",
                        "mean_log_prob": mean_lp,
                        "temp": float(temp),
                        "mass_condition": mass_condition,
                        "formula_condition_source": formula_source,
                        "formula_condition_formula": formula_text,
                        "formula_condition_nonzero": formula_nonzero,
                        "substructure_condition_source": substructure_source,
                        "substructure_condition_nonzero": substructure_nonzero,
                    }
                )
                generation_draws[feature_key] += 1
                continue

            generation_counts[feature_key][smiles] += 1
            generation_draws[feature_key] += 1
            prev = best_log_probs[feature_key].get(smiles, float("-inf"))
            if mean_lp > prev:
                best_log_probs[feature_key][smiles] = mean_lp
            raw_rows.append(
                {
                    "sample_idx": int(sample_idx),
                    "feature_key": feature_key,
                    "selfies": selfies_text,
                    "smiles": smiles,
                    "mean_log_prob": mean_lp,
                    "temp": float(temp),
                    "mass_condition": mass_condition,
                    "formula_condition_source": formula_source,
                    "formula_condition_formula": formula_text,
                    "formula_condition_nonzero": formula_nonzero,
                    "substructure_condition_source": substructure_source,
                    "substructure_condition_nonzero": substructure_nonzero,
                }
            )

    for sample_idx, temp in enumerate(schedule):
        if sample_idx % 5 == 0 or sample_idx == (len(schedule) - 1):
            print(
                f"[decode sample {sample_idx + 1:02d}/{len(schedule):02d}] "
                f"temp={float(temp):.2f}"
            )
        for start in range(0, len(feature_keys), int(args.decode_batch_size)):
            end = min(start + int(args.decode_batch_size), len(feature_keys))
            record_decoded_batch(
                sample_idx=int(sample_idx),
                temp=float(temp),
                batch_positions=list(range(start, end)),
            )

    adaptive_summary: dict[str, object] = {
        "enabled": False,
        "target_mass_matched_per_feature": int(args.target_mass_matched_per_feature),
        "max_samples_per_feature": int(args.max_samples_per_feature),
        "additional_rounds": 0,
        "additional_draws": 0,
    }
    target_mass_matched = max(int(args.target_mass_matched_per_feature), 0)
    if target_mass_matched > 0 and feature_keys:
        parent_mz_by_feature = {
            str(key): float(parent_mz)
            for key, parent_mz in zip(feature_keys, parent_mzs)
        }

        def count_mass_matched_candidates(feature_key: str) -> int:
            parent_mz = float(parent_mz_by_feature[feature_key])
            matched = 0
            for smiles in generation_counts[feature_key].keys():
                pred_mass = exact_mol_weight(smiles)
                if pred_mass is None:
                    continue
                abs_err = adduct_neutral_mass_error(
                    pred_mass, parent_mz, ion_mode=ion_mode
                )
                if not bool(args.disable_adduct_hypotheses):
                    mass_match = best_adduct_neutral_mass_match(
                        pred_mass,
                        parent_mz,
                        ion_mode=ion_mode,
                        include_isotopes=True,
                    )
                    abs_err = float(mass_match["abs_mass_error_da"])
                if abs_err <= float(args.mass_tolerance_da):
                    matched += 1
            return matched

        max_samples_per_feature = int(args.max_samples_per_feature)
        if max_samples_per_feature <= 0:
            max_samples_per_feature = max(int(args.num_samples) * 5, int(args.num_samples))
        if max_samples_per_feature < int(args.num_samples):
            raise ValueError(
                "--max_samples_per_feature must be >= --num_samples when "
                "adaptive mass sampling is enabled."
            )

        next_sample_idx = len(schedule)
        additional_draws = 0
        additional_rounds = 0
        while True:
            active_positions = [
                idx
                for idx, feature_key in enumerate(feature_keys)
                if int(generation_draws.get(feature_key, 0)) < max_samples_per_feature
                and count_mass_matched_candidates(feature_key) < target_mass_matched
            ]
            if not active_positions:
                break
            temp = float(schedule[next_sample_idx % len(schedule)] if schedule else 1.0)
            additional_rounds += 1
            if additional_rounds == 1 or additional_rounds % 10 == 0:
                print(
                    "[adaptive mass decode "
                    f"round {additional_rounds}] active_features={len(active_positions)} "
                    f"temp={temp:.2f}"
                )
            for start in range(0, len(active_positions), int(args.decode_batch_size)):
                batch_positions = active_positions[start : start + int(args.decode_batch_size)]
                record_decoded_batch(
                    sample_idx=int(next_sample_idx),
                    temp=temp,
                    batch_positions=batch_positions,
                )
                additional_draws += len(batch_positions)
            next_sample_idx += 1

        adaptive_summary = {
            "enabled": True,
            "target_mass_matched_per_feature": int(target_mass_matched),
            "max_samples_per_feature": int(max_samples_per_feature),
            "additional_rounds": int(additional_rounds),
            "additional_draws": int(additional_draws),
            "features_satisfying_target": int(
                sum(
                    count_mass_matched_candidates(feature_key) >= target_mass_matched
                    for feature_key in feature_keys
                )
            ),
        }

    raw_df = pd.DataFrame(raw_rows)
    raw_csv = output_dir / artifact_filename(
        artifact_prefix, "generation_samples_raw", ".csv"
    )
    write_dataframe_with_aliases(
        raw_df,
        raw_csv,
        [
            output_dir / "selfies_generation_samples_raw.csv",
            output_dir / "generation_samples_raw.csv",
        ],
    )

    feature_meta_lookup = {
        str(row["feature_key"]): row
        for row in feature_meta_df.to_dict(orient="records")
    }

    denovo_rows: list[dict[str, object]] = []
    candidate_metrics_rows: list[dict[str, object]] = []
    features_with_any_valid = 0
    features_with_any_within_1da = 0
    substructure_filter_features_applied = 0
    substructure_filter_candidates_removed = 0
    substructure_filter_missing_constraint_features = 0

    for feature_key in feature_keys:
        parent_mz = float(feature_meta_lookup[feature_key]["parent_mz"])
        feature_pos = int(feature_index_by_key[str(feature_key)])
        feature_formula_source = str(formula_condition_sources[feature_pos])
        feature_formula_text = str(formula_condition_formulas[feature_pos])
        feature_formula_nonzero = bool(formula_condition_nonzero[feature_pos])
        feature_substructure_source = str(substructure_condition_sources[feature_pos])
        feature_substructure_nonzero = bool(substructure_condition_nonzero[feature_pos])
        draws = max(int(generation_draws.get(feature_key, 0)), 1)
        candidates: list[dict[str, object]] = []

        for smiles, support_count in generation_counts[feature_key].items():
            pred_mass = exact_mol_weight(smiles)
            if pred_mass is None:
                continue
            primary_abs_err = adduct_neutral_mass_error(
                pred_mass, parent_mz, ion_mode=ion_mode
            )
            mass_match = {
                "abs_mass_error_da": float(primary_abs_err),
                "hypothesis": "[M-H]-" if ion_mode == "negative" else "[M+H]+",
                "adduct": "[M-H]-" if ion_mode == "negative" else "[M+H]+",
                "neutral_mass": float(
                    neutral_mass_from_precursor(parent_mz, ion_mode=ion_mode)
                ),
                "isotope_shift_count": 0,
            }
            if not bool(args.disable_adduct_hypotheses):
                mass_match = best_adduct_neutral_mass_match(
                    pred_mass,
                    parent_mz,
                    ion_mode=ion_mode,
                    include_isotopes=True,
                )
            abs_err = float(mass_match["abs_mass_error_da"])
            strict_pass = bool(
                abs_err <= strict_tolerance_da(parent_mz, ion_mode=ion_mode)
            )
            err_bin = int(mass_error_bin(abs_err, strict_pass))
            candidates.append(
                {
                    "comp_id": int(feature_meta_lookup[feature_key]["comp_id"]),
                    "parent_mz": float(parent_mz),
                    "predicted_smiles": smiles,
                    "predicted_exact_mass": float(pred_mass),
                    "predicted_formula": molecular_formula(smiles),
                    "formula_condition_source": feature_formula_source,
                    "formula_condition_formula": feature_formula_text,
                    "formula_condition_nonzero": feature_formula_nonzero,
                    "substructure_condition_source": feature_substructure_source,
                    "substructure_condition_nonzero": feature_substructure_nonzero,
                    "abs_mass_error_da": float(abs_err),
                    "primary_abs_mass_error_da": float(primary_abs_err),
                    "mass_error_bin": int(err_bin),
                    "mass_error_hypothesis": str(mass_match.get("hypothesis", "")),
                    "mass_hypothesis_adduct": str(mass_match.get("adduct", "")),
                    "mass_hypothesis_neutral_mass": float(mass_match.get("neutral_mass", float("nan"))),
                    "mass_hypothesis_isotope_shift_count": int(mass_match.get("isotope_shift_count", 0) or 0),
                    "strict_mass_pass": strict_pass,
                    "mass_match_score": mass_match_score(
                        abs_err, tolerance_da=float(args.mass_tolerance_da)
                    ),
                    "mean_log_prob": float(
                        best_log_probs.get(feature_key, {}).get(
                            smiles, float("nan")
                        )
                    ),
                    "retrieval_tanimoto": float("nan"),
                    "substructure_filter_pass": True,
                    "generation_support_count": int(support_count),
                    "generation_support_fraction": float(
                        support_count / draws
                    ),
                }
            )

        if candidates:
            features_with_any_valid += 1
        if any(float(row["abs_mass_error_da"]) <= 1.0 for row in candidates):
            features_with_any_within_1da += 1

        for row in candidates:
            candidate_metrics_rows.append({"feature_key": feature_key, **row})

        kept = [
            row
            for row in candidates
            if float(row["abs_mass_error_da"]) <= float(args.mass_tolerance_da)
        ]

        if effective_filter_mode == "consensus_substructure":
            consensus_smarts = consensus_smarts_by_feature.get(str(feature_key), "")
            if consensus_smarts:
                substructure_filter_features_applied += 1
                before_count = len(kept)
                constrained = []
                for row in kept:
                    passed = bool(
                        smiles_matches_substructure(
                            str(row.get("predicted_smiles", "")),
                            consensus_smarts,
                        )
                    )
                    row["substructure_filter_pass"] = passed
                    if passed:
                        constrained.append(row)
                kept = constrained
                substructure_filter_candidates_removed += max(before_count - len(kept), 0)
            else:
                substructure_filter_missing_constraint_features += 1
        elif effective_filter_mode == "retrieval_tanimoto":
            retrieval_anchor_smiles = top_retrieval_smiles_by_feature.get(str(feature_key), "")
            if retrieval_anchor_smiles:
                substructure_filter_features_applied += 1
                before_count = len(kept)
                constrained = []
                for row in kept:
                    tanimoto_val = tanimoto_similarity_smiles(
                        str(row.get("predicted_smiles", "")),
                        retrieval_anchor_smiles,
                    )
                    row["retrieval_tanimoto"] = float(tanimoto_val)
                    passed = bool(tanimoto_val >= float(args.retrieval_tanimoto_threshold))
                    row["substructure_filter_pass"] = passed
                    if passed:
                        constrained.append(row)
                kept = constrained
                substructure_filter_candidates_removed += max(before_count - len(kept), 0)
            else:
                substructure_filter_missing_constraint_features += 1

        if str(args.ranking_mode) == "mass_bin_support":
            kept.sort(
                key=lambda row: (
                    mass_error_bin(
                        float(row["abs_mass_error_da"]),
                        bool(row["strict_mass_pass"]),
                    ),
                    -float(row["generation_support_fraction"]),
                    -np.nan_to_num(float(row["mean_log_prob"]), nan=-99.0),
                    float(row["abs_mass_error_da"]),
                    str(row["predicted_smiles"]),
                )
            )
        else:
            kept.sort(
                key=lambda row: (
                    float(row["abs_mass_error_da"]),
                    -np.nan_to_num(float(row["mean_log_prob"]), nan=-99.0),
                    -float(row["generation_support_fraction"]),
                    str(row["predicted_smiles"]),
                )
            )

        for rank, row in enumerate(kept, start=1):
            denovo_rows.append(
                {
                    "feature_key": feature_key,
                    "comp_id": int(
                        feature_meta_lookup[feature_key]["comp_id"]
                    ),
                    "feature_label": str(
                        feature_meta_lookup[feature_key]["feature_label"]
                    ),
                    "parent_mz": parent_mz,
                    "rank": int(rank),
                    "predicted_smiles": str(row["predicted_smiles"]),
                    "is_valid": True,
                    "predicted_exact_mass": float(row["predicted_exact_mass"]),
                    "predicted_formula": row["predicted_formula"],
                    "formula_condition_source": str(row["formula_condition_source"]),
                    "formula_condition_formula": str(row["formula_condition_formula"]),
                    "formula_condition_nonzero": bool(row["formula_condition_nonzero"]),
                    "mass_match_score": float(row["mass_match_score"]),
                    "abs_mass_error_da": float(row["abs_mass_error_da"]),
                    "primary_abs_mass_error_da": float(row["primary_abs_mass_error_da"]),
                    "mass_error_bin": int(row["mass_error_bin"]),
                    "mass_error_hypothesis": str(row["mass_error_hypothesis"]),
                    "mass_hypothesis_adduct": str(row["mass_hypothesis_adduct"]),
                    "mass_hypothesis_neutral_mass": float(row["mass_hypothesis_neutral_mass"]),
                    "mass_hypothesis_isotope_shift_count": int(row["mass_hypothesis_isotope_shift_count"]),
                    "strict_mass_pass": bool(row["strict_mass_pass"]),
                    "mean_log_prob": float(row["mean_log_prob"]),
                    "generation_support_count": int(
                        row["generation_support_count"]
                    ),
                    "generation_support_fraction": float(
                        row["generation_support_fraction"]
                    ),
                    "candidate_latent_cosine": float("nan"),
                    "samples_drawn": int(draws),
                }
            )

    denovo_df = pd.DataFrame(denovo_rows)
    if not denovo_df.empty:
        denovo_df = denovo_df.sort_values(["feature_key", "rank"]).reset_index(
            drop=True
        )
    denovo_csv = output_dir / artifact_filename(
        artifact_prefix, "denovo_report", ".csv"
    )
    selfies_denovo_csv = output_dir / "selfies_denovo_report.csv"
    write_dataframe_with_aliases(
        denovo_df,
        denovo_csv,
        [
            selfies_denovo_csv,
            output_dir / "denovo_report.csv",
        ],
    )

    candidate_metrics_df = pd.DataFrame(candidate_metrics_rows)
    candidate_metrics_csv = output_dir / artifact_filename(
        artifact_prefix, "candidate_metrics_all", ".csv"
    )
    write_dataframe_with_aliases(
        candidate_metrics_df,
        candidate_metrics_csv,
        [
            output_dir / "selfies_candidate_metrics_all.csv",
            output_dir / "candidate_metrics_all.csv",
        ],
    )

    mass_accuracy = {
        "n_candidates": int(candidate_metrics_df.shape[0]),
        "within_0p5_da": (
            int((candidate_metrics_df["abs_mass_error_da"] <= 0.5).sum())
            if not candidate_metrics_df.empty
            else 0
        ),
        "within_1p0_da": (
            int((candidate_metrics_df["abs_mass_error_da"] <= 1.0).sum())
            if not candidate_metrics_df.empty
            else 0
        ),
        "within_5p0_da": (
            int((candidate_metrics_df["abs_mass_error_da"] <= 5.0).sum())
            if not candidate_metrics_df.empty
            else 0
        ),
    }
    n_candidates = max(int(mass_accuracy["n_candidates"]), 1)
    mass_accuracy.update(
        {
            "within_0p5_da_rate": float(
                mass_accuracy["within_0p5_da"] / n_candidates
            ),
            "within_1p0_da_rate": float(
                mass_accuracy["within_1p0_da"] / n_candidates
            ),
            "within_5p0_da_rate": float(
                mass_accuracy["within_5p0_da"] / n_candidates
            ),
        }
    )

    control_eval = {"requested": not bool(args.skip_control_eval)}
    if not bool(args.skip_control_eval):
        control_eval = run_control_eval(
            repo_root=repo_root,
            output_dir=output_dir,
            denovo_csv=denovo_csv,
            control_lib_id=int(args.control_lib_id),
        )

    summary = {
        "model_ckpt": str(model_ckpt),
        "artifact_prefix": artifact_prefix,
        "encoder_ckpt": str(encoder_ckpt),
        "encoder_kind": str(encoder_kind),
        "encoder_calibrator": (
            str(calibrator_path) if calibrator_path is not None else None
        ),
        "encoder_conditioning": {
            "lib_id": int(encoder_lib_id) if encoder_lib_id is not None else None,
            "encoder_ion_mode": float(encoder_ion_mode) if encoder_ion_mode is not None else None,
        },
        "seed": int(args.seed),
        "feature_table": str(feature_table),
        "feature_scope": feature_scope_summary,
        "device": str(device),
        "ion_mode": ion_mode,
        "parent_mz_normalization": {
            "min": float(args.parent_mz_min),
            "max": float(args.parent_mz_max),
        },
        "decoder_mass_normalization": {
            "min": float(train_mass_min),
            "max": float(train_mass_max),
        },
        "feature_latent_norm_pre_calibration": summarize_numeric_array(
            latent_norms_pre
        ),
        "feature_latent_norm_post_calibration": summarize_numeric_array(
            latent_norms_post
        ),
        "num_features": int(len(feature_meta_df)),
        "num_samples_per_feature": int(args.num_samples),
        "adaptive_mass_sampling": adaptive_summary,
        "decode_temps": [float(x) for x in schedule],
        "top_p": float(args.top_p),
        "ranking_mode": str(args.ranking_mode),
        "mass_tolerance_da": float(args.mass_tolerance_da),
        "adduct_hypotheses": {
            "enabled": not bool(args.disable_adduct_hypotheses),
            "primary_only_when_disabled": "[M-H]-/[M+H]+",
        },
        "substructure_filter": {
            "requested_mode": requested_filter_mode,
            "effective_mode": effective_filter_mode,
            "retrieval_report_csv": str(retrieval_report_csv) if retrieval_report_csv is not None else None,
            "consensus_features_available": int(len(consensus_smarts_by_feature)),
            "top_retrieval_features_available": int(len(top_retrieval_smiles_by_feature)),
            "features_with_constraints_applied": int(substructure_filter_features_applied),
            "features_missing_constraint": int(substructure_filter_missing_constraint_features),
            "candidates_removed": int(substructure_filter_candidates_removed),
            "consensus_top_n": int(args.substructure_consensus_top_n),
            "consensus_min_mcs_atoms": int(args.substructure_min_mcs_atoms),
            "consensus_mcs_timeout_sec": int(args.substructure_mcs_timeout_sec),
            "retrieval_tanimoto_threshold": float(args.retrieval_tanimoto_threshold),
        },
        "ri_conditioning": {
            "requested": not bool(args.disable_ri_conditioning),
            "supported": bool(checkpoint_ri_supported),
            "effective": bool(ri_conditioning_effective),
            "normalization": ri_norm_payload,
        },
        "formula_conditioning": {
            "requested_version": requested_conditioning_version,
            "checkpoint_version": checkpoint_conditioning_version,
            "effective": bool(formula_conditioning_effective),
            "formula_csv": str(formula_csv) if formula_csv is not None else None,
            "metadata_rows": int(len(formula_feature_map)),
            "features_with_formula_condition": int(sum(formula_condition_nonzero)),
            "formula_elements": list(FORMULA_ELEMENTS),
        },
        "adduct_conditioning": {
            "requested": bool(args.adduct_conditioning),
            "effective": bool(adduct_conditioning_effective),
        },
        "substructure_conditioning": {
            "effective": bool(substructure_conditioning_effective),
            "embedding_csv": str(substructure_embedding_csv) if substructure_embedding_csv is not None else None,
            "embedding_matrix": str(substructure_embedding_pt) if substructure_embedding_pt is not None else None,
            "metadata_rows": int(len(substructure_embedding_map)),
            "input_dim": int(substructure_input_dim),
            "features_with_substructure_condition": int(sum(substructure_condition_nonzero)),
        },
        "features_with_any_valid_candidate": int(features_with_any_valid),
        "features_with_candidate_within_1da": int(
            features_with_any_within_1da
        ),
        "features_with_output_rows": (
            int(denovo_df["feature_key"].nunique())
            if not denovo_df.empty
            else 0
        ),
        "denovo_rows": int(len(denovo_df)),
        "mass_accuracy": mass_accuracy,
        "control_eval": control_eval,
        "outputs": {
            "feature_meta_csv": str(
                output_dir / artifact_filename(artifact_prefix, "feature_meta", ".csv")
            ),
            "feature_latents_pt": str(
                output_dir / artifact_filename(artifact_prefix, "feature_latents", ".pt")
            ),
            "raw_samples_csv": str(raw_csv),
            "candidate_metrics_csv": str(candidate_metrics_csv),
            "denovo_csv": str(denovo_csv),
            "selfies_denovo_csv": str(selfies_denovo_csv),
        },
    }

    feature_meta_out = output_dir / artifact_filename(
        artifact_prefix, "feature_meta", ".csv"
    )
    write_dataframe_with_aliases(
        feature_meta_df,
        feature_meta_out,
        [output_dir / "feature_meta.csv"],
    )
    feature_latents_out = output_dir / artifact_filename(
        artifact_prefix, "feature_latents", ".pt"
    )
    feature_latent_payload = {
        str(feature_meta_df.iloc[i]["feature_key"]): feature_latents[i]
        .detach()
        .cpu()
        for i in range(len(feature_meta_df))
    }
    write_tensor_with_aliases(
        feature_latent_payload,
        feature_latents_out,
        [output_dir / "feature_latents.pt"],
    )

    summary_out = output_dir / artifact_filename(
        artifact_prefix, "selfies_run_summary", ".json"
    )
    write_json_with_aliases(
        summary,
        summary_out,
        [output_dir / "selfies_run_summary.json"],
    )

    print(f"Wrote: {feature_meta_out}")
    print(f"Wrote: {feature_latents_out}")
    print(f"Wrote: {raw_csv}")
    print(f"Wrote: {candidate_metrics_csv}")
    print(f"Wrote: {denovo_csv}")
    print(f"Wrote: {selfies_denovo_csv}")
    print(f"Wrote: {summary_out}")

    if (
        isinstance(control_eval, dict)
        and int(control_eval.get("return_code", 0)) != 0
        and not bool(args.skip_control_eval)
    ):
        raise SystemExit(int(control_eval["return_code"]))


if __name__ == "__main__":
    main()
