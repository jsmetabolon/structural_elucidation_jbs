#!/usr/bin/env python3
"""Merged-ionization inference wrapper over predict_smiles.py."""

from __future__ import annotations

import argparse
import tempfile
from pathlib import Path
import sys
from typing import Iterable

import pandas as pd

from structural_elucidation.paths import REPO_ROOT
import torch

import structural_elucidation.inference.predict as shared_predict
from structural_elucidation.inference.model_conv import infer_mass_conditioned_compatible_length

try:
    from structural_elucidation.inference.merged_model import (
        MergedNet1D,
        ion_mode_label_to_scalar,
        lib_id_to_index,
        load_compatible_state_dict_for_merged,
    )
except ImportError:
    from structural_elucidation.inference.merged_model import (
        MergedNet1D,
        ion_mode_label_to_scalar,
        lib_id_to_index,
        load_compatible_state_dict_for_merged,
    )


_ORIGINAL_GET_PARSER = shared_predict.get_parser
_ENCODER_DEFAULT_ION_MODE = 0.0
_ENCODER_DEFAULT_LIB_IDX = 0
_LIB_TO_MODE = {
    209: "negative",
    305: "negative",
    400: "positive",
    402: "positive",
}


def _path_exists(path: Path) -> bool:
    return path.exists() and path.is_file()


def _first_existing(candidates: Iterable[Path]) -> Path | None:
    for candidate in candidates:
        path = Path(candidate)
        if path.exists():
            return path.resolve()
    return None


def _has_flag(argv: list[str], flag: str) -> bool:
    return any(token == flag or token.startswith(f"{flag}=") for token in argv)


def _extract_arg_value(argv: list[str], flag: str) -> str | None:
    for idx, token in enumerate(argv):
        if token == flag and idx + 1 < len(argv):
            return argv[idx + 1]
        if token.startswith(f"{flag}="):
            return token.split("=", 1)[1]
    return None


def _remove_arg(argv: list[str], flag: str) -> tuple[list[str], str | None]:
    out: list[str] = []
    value: str | None = None
    skip_next = False
    for idx, token in enumerate(argv):
        if skip_next:
            skip_next = False
            continue

        if token == flag:
            if idx + 1 < len(argv):
                value = argv[idx + 1]
                skip_next = True
            continue

        if token.startswith(f"{flag}="):
            value = token.split("=", 1)[1]
            continue

        out.append(token)

    return out, value


def _parse_ion_mode(argv: list[str]) -> str:
    mode = _extract_arg_value(argv, "--ion_mode")
    if mode is None:
        raise ValueError(
            "predict_merged.py requires --ion_mode (negative or positive)."
        )
    text = str(mode).strip().lower()
    if text in {"negative", "positive"}:
        return text
    raise ValueError(
        f"Unsupported --ion_mode value {mode!r}. Expected 'negative' or 'positive'."
    )


def _load_table(path: Path) -> pd.DataFrame:
    if path.suffix.lower() == ".csv":
        return pd.read_csv(path)
    if path.suffix.lower() == ".tsv":
        return pd.read_csv(path, sep="\t")
    return pd.read_pickle(path)


def _write_temp_filtered_feature_table(
    *,
    feature_table: Path,
    ion_mode: str,
    lib_id: int | None,
) -> Path:
    df = _load_table(feature_table)
    filtered = df.copy()

    if "ion_mode" in filtered.columns:
        mode_scalar = 1.0 if ion_mode == "positive" else 0.0
        filtered["ion_mode"] = pd.to_numeric(filtered["ion_mode"], errors="coerce")
        filtered = filtered[filtered["ion_mode"].eq(mode_scalar)].copy()

    if lib_id is not None:
        if "lib_id" not in filtered.columns:
            raise KeyError(
                "--lib_id was provided but feature_table does not contain lib_id column."
            )
        filtered["lib_id"] = pd.to_numeric(filtered["lib_id"], errors="coerce")
        filtered = filtered[filtered["lib_id"].eq(int(lib_id))].copy()

    if filtered.empty:
        raise RuntimeError(
            f"No rows remain after filtering feature_table by ion_mode={ion_mode}"
            + (f" and lib_id={lib_id}" if lib_id is not None else "")
        )

    tmp_dir = REPO_ROOT / "merged_encoder" / "data" / "tmp_predict_inputs"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    fd, tmp_name = tempfile.mkstemp(
        prefix="feature_table_",
        suffix=".csv",
        dir=str(tmp_dir),
    )
    Path(tmp_name).unlink(missing_ok=True)
    filtered.to_csv(tmp_name, index=False)

    # close OS fd from mkstemp
    try:
        import os

        os.close(fd)
    except Exception:
        pass

    return Path(tmp_name).resolve()


def _patched_get_parser():
    parser = _ORIGINAL_GET_PARSER()
    parser.add_argument(
        "--encoder_ion_mode",
        type=float,
        default=None,
        help=(
            "Optional scalar conditioning value passed to merged mode_embedding. "
            "Defaults to 0.0 for --ion_mode negative and 1.0 for positive."
        ),
    )
    parser.add_argument(
        "--lib_id",
        type=int,
        default=None,
        choices=(209, 305, 400, 402),
        help="Required library id conditioning value (209, 305, 400, or 402).",
    )
    parser.add_argument(
        "--decoder_model_name_neg",
        type=str,
        default="",
        help=(
            "Wrapper-only override for negative-mode constrained de novo decoder "
            "checkpoint name inside --decoder_dir."
        ),
    )
    parser.add_argument(
        "--decoder_model_name_pos",
        type=str,
        default="",
        help=(
            "Wrapper-only override for positive-mode constrained de novo decoder "
            "checkpoint name inside --decoder_dir."
        ),
    )
    return parser


def _patched_load_encoder_checkpoint(
    ckpt_path: Path,
    n_bins: int,
    minmass: float,
    maxmass: float,
    parent_mz_min: float,
    parent_mz_max: float,
    device: torch.device,
    apply_tanh: bool = True,
):
    state = shared_predict.robust_torch_load(ckpt_path, map_location="cpu", weights_only=False)
    if not isinstance(state, dict):
        raise RuntimeError(f"Expected state_dict in {ckpt_path}, got {type(state)}")

    model_n_bins = int(
        infer_mass_conditioned_compatible_length(state, preferred_length_in=n_bins)
        or n_bins
    )

    model = MergedNet1D(
        length_in=model_n_bins,
        minmass=minmass,
        maxmass=maxmass,
        parent_mz_min=parent_mz_min,
        parent_mz_max=parent_mz_max,
        apply_tanh=apply_tanh,
    ).to(device)
    model.default_ion_mode_idx = int(round(float(_ENCODER_DEFAULT_ION_MODE)))
    model.default_lib_idx = int(_ENCODER_DEFAULT_LIB_IDX)

    load_info = load_compatible_state_dict_for_merged(
        model,
        state,
        ignore_new_conditioning_params=False,
    )
    load_info["encoder_default_ion_mode"] = float(_ENCODER_DEFAULT_ION_MODE)
    load_info["encoder_default_lib_idx"] = int(_ENCODER_DEFAULT_LIB_IDX)

    model.eval()
    return model, model_n_bins, load_info


def _prepare_argv(argv: list[str]) -> list[str]:
    ion_mode = _parse_ion_mode(argv)

    # Wrapper-only args removed before shared parser sees them.
    argv, lib_id_raw = _remove_arg(argv, "--lib_id")
    argv, encoder_ion_mode_raw = _remove_arg(argv, "--encoder_ion_mode")
    argv, decoder_model_name_neg = _remove_arg(argv, "--decoder_model_name_neg")
    argv, decoder_model_name_pos = _remove_arg(argv, "--decoder_model_name_pos")

    if lib_id_raw is None:
        raise ValueError(
            "predict_merged.py requires --lib_id (one of 209, 305, 400, 402)."
        )
    lib_id = int(lib_id_raw)
    expected_mode = _LIB_TO_MODE.get(int(lib_id))
    if expected_mode is None:
        raise ValueError("Unsupported lib_id. Expected one of 209, 305, 400, 402.")
    if str(ion_mode) != str(expected_mode):
        raise ValueError(
            f"lib_id={lib_id} is mapped to ion_mode={expected_mode!r}, "
            f"but --ion_mode was {ion_mode!r}."
        )

    global _ENCODER_DEFAULT_ION_MODE, _ENCODER_DEFAULT_LIB_IDX
    if encoder_ion_mode_raw is not None:
        _ENCODER_DEFAULT_ION_MODE = float(encoder_ion_mode_raw)
    else:
        _ENCODER_DEFAULT_ION_MODE = float(1.0 if ion_mode == "positive" else 0.0)
    _ENCODER_DEFAULT_LIB_IDX = int(lib_id_to_index(lib_id))

    selected_decoder_model_name = (
        str(decoder_model_name_pos).strip()
        if ion_mode == "positive"
        else str(decoder_model_name_neg).strip()
    )
    explicit_decoder_model_name = _extract_arg_value(argv, "--decoder_model_name")
    if selected_decoder_model_name:
        if (
            explicit_decoder_model_name is not None
            and str(explicit_decoder_model_name).strip()
            and str(explicit_decoder_model_name).strip() != selected_decoder_model_name
        ):
            raise ValueError(
                "Conflicting decoder checkpoint selection: --decoder_model_name "
                f"{explicit_decoder_model_name!r} != mode-routed {selected_decoder_model_name!r} "
                f"for ion_mode={ion_mode!r}."
            )
        if explicit_decoder_model_name is None:
            argv.extend(["--decoder_model_name", selected_decoder_model_name])

    neg_mz = _first_existing(
        [
            REPO_ROOT / "SMILES_DB" / "cid_to_parent_mz_runtime_combined_neg.pkl",
            REPO_ROOT / "SMILES_DB" / "cid_to_parent_mz_expanded_neg.pkl",
            REPO_ROOT / "SMILES_DB" / "cid_to_parent_mz_expanded.pkl",
        ]
    )
    pos_mz = _first_existing(
        [
            REPO_ROOT / "SMILES_DB" / "cid_to_parent_mz_runtime_combined_pos.pkl",
            REPO_ROOT / "SMILES_DB" / "cid_to_parent_mz_expanded_pos.pkl",
            REPO_ROOT / "SMILES_DB" / "cid_to_parent_mz_x400_allowed_pos.pkl",
            REPO_ROOT / "SMILES_DB" / "cid_to_parent_mz_expanded.pkl",
        ]
    )

    default_bank = _first_existing(
        [
            REPO_ROOT / "embeddings_latent_from_smiles_runtime_combined",
            REPO_ROOT / "embeddings_latent_from_smiles_expanded_950",
            REPO_ROOT / "embeddings_latent_from_smiles_expanded",
        ]
    )
    default_smiles = _first_existing(
        [
            REPO_ROOT / "SMILES_DB" / "comp_id_to_smiles_runtime_combined.csv",
            REPO_ROOT / "SMILES_DB" / "comp_id_to_smiles_expanded_950.csv",
            REPO_ROOT / "SMILES_DB" / "comp_id_to_smiles_expanded.csv",
        ]
    )
    default_name_csv = _first_existing(
        [
            REPO_ROOT / "SMILES_DB" / "name_dict_runtime_combined.csv",
            REPO_ROOT / "SMILES_DB" / "name_dict.csv",
        ]
    )
    default_feature_table = _first_existing(
        [
            REPO_ROOT / "merged_encoder" / "data" / "ions_unknown_merged.csv",
            REPO_ROOT / "merged_encoder" / "data" / "ions_unknown_merged.pkl",
        ]
    )
    default_fragment_db = _first_existing(
        [
            REPO_ROOT / "SMILES_DB" / f"library_{int(lib_id)}_fragment_database_element0.csv",
            REPO_ROOT / "SMILES_DB" / f"library_{int(lib_id)}_fragment_database.csv",
            REPO_ROOT / "SMILES_DB" / "library_209_fragment_database_element0.csv",
        ]
    )
    default_prevalence = _first_existing(
        [
            REPO_ROOT / "SMILES_DB" / (
                "fasstrecords_prevalence_pos.csv" if ion_mode == "positive" else "fasstrecords_prevalence_neg.csv"
            ),
            REPO_ROOT / "SMILES_DB" / "fasstrecords_prevalence.csv",
        ]
    )

    if not _has_flag(argv, "--retrieval_mz_pkl"):
        selected_mz = pos_mz if ion_mode == "positive" else neg_mz
        if selected_mz is not None:
            argv.extend(["--retrieval_mz_pkl", str(selected_mz)])

    if not _has_flag(argv, "--embedding_bank_dir") and default_bank is not None:
        argv.extend(["--embedding_bank_dir", str(default_bank)])

    if not _has_flag(argv, "--retrieval_smiles_csv") and default_smiles is not None:
        argv.extend(["--retrieval_smiles_csv", str(default_smiles)])

    if not _has_flag(argv, "--retrieval_name_csv") and default_name_csv is not None:
        argv.extend(["--retrieval_name_csv", str(default_name_csv)])

    if not _has_flag(argv, "--retrieval_rescoring"):
        argv.extend(["--retrieval_rescoring", "on"])

    if not _has_flag(argv, "--fragment_db_csv") and default_fragment_db is not None:
        argv.extend(["--fragment_db_csv", str(default_fragment_db)])

    if not _has_flag(argv, "--prevalence_csv") and default_prevalence is not None:
        argv.extend(["--prevalence_csv", str(default_prevalence)])

    feature_table_raw = _extract_arg_value(argv, "--feature_table")
    feature_table_path = (
        Path(feature_table_raw).expanduser().resolve()
        if feature_table_raw is not None
        else default_feature_table
    )

    should_filter_features = feature_table_path is not None and feature_table_path.exists()
    if should_filter_features:
        filtered_path = _write_temp_filtered_feature_table(
            feature_table=feature_table_path,
            ion_mode=ion_mode,
            lib_id=int(lib_id),
        )

        if _has_flag(argv, "--feature_table"):
            argv, _ = _remove_arg(argv, "--feature_table")
        argv.extend(["--feature_table", str(filtered_path)])
    elif feature_table_path is not None and not _has_flag(argv, "--feature_table"):
        argv.extend(["--feature_table", str(feature_table_path)])

    return argv


def main() -> None:
    shared_predict.get_parser = _patched_get_parser
    shared_predict.load_encoder_checkpoint = _patched_load_encoder_checkpoint

    raw_args = sys.argv[1:]
    if any(token in {"-h", "--help"} for token in raw_args):
        shared_predict.main()
        return

    prepared = _prepare_argv(raw_args)
    sys.argv = [sys.argv[0]] + prepared
    shared_predict.main()


if __name__ == "__main__":
    main()
