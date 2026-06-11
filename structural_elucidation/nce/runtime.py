#!/usr/bin/env python3
# pylint: disable=import-error,broad-exception-caught
"""Local runtime shims for using NewCompoundElucidation assets without DB access."""

from __future__ import annotations

import os
import pickle
import re
import sys
import types
from functools import lru_cache
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from structural_elucidation.paths import REPO_ROOT

DEFAULT_NCE_ROOT = REPO_ROOT / "vendor" / "nce"
DEFAULT_REFERENCE_DIR = REPO_ROOT / "vendor" / "modelfiles"
DEFAULT_MLIMS_TABLE_CANDIDATES = (
    REPO_ROOT / "vendor" / "modelfiles" / "mlims_tables.pkl",
    REPO_ROOT / "vendor" / "modelfiles" / "mlims_tables.pkl",
)


def install_local_nce_shims() -> None:
    """Install tiny shims for optional NCE deps not needed in cache-only mode."""
    sys.modules.setdefault("regex", re)

    if "dotenv" not in sys.modules:
        dotenv = types.ModuleType("dotenv")
        dotenv.load_dotenv = lambda *args, **kwargs: None
        sys.modules["dotenv"] = dotenv

    if "mdsutils.easyQuery" not in sys.modules:
        mdsutils = types.ModuleType("mdsutils")
        easy_query = types.ModuleType("mdsutils.easyQuery")

        class EasyQuery:  # pragma: no cover - only used to guard accidental DB paths
            def __init__(self, *args: Any, **kwargs: Any) -> None:
                raise RuntimeError("EasyQuery is unavailable in local cache-only NCE mode")

        easy_query.EasyQuery = EasyQuery
        sys.modules.setdefault("mdsutils", mdsutils)
        sys.modules["mdsutils.easyQuery"] = easy_query


def prepare_nce_runtime(
    *,
    nce_root: Path | str = DEFAULT_NCE_ROOT,
    reference_dir: Path | str = DEFAULT_REFERENCE_DIR,
) -> None:
    """Point NCE imports at the local checkout and model/cache directory."""
    root = Path(str(nce_root)).expanduser()
    ref = Path(str(reference_dir)).expanduser()
    os.environ["REFERENCE_DIR"] = str(ref)
    if str(root) not in sys.path:
        sys.path.insert(0, str(root))
    install_local_nce_shims()


@lru_cache(maxsize=4)
def get_predicted_spectrum_matcher(
    lib_id: int = 209,
    *,
    nce_root: str = str(DEFAULT_NCE_ROOT),
    reference_dir: str = str(DEFAULT_REFERENCE_DIR),
):
    prepare_nce_runtime(nce_root=nce_root, reference_dir=reference_dir)
    from elucidation.predictedSpectrumMatcher import PredictedSpectrumMatcher

    return PredictedSpectrumMatcher.get_instance(int(lib_id))


@lru_cache(maxsize=1)
def get_cached_mlims_service(
    *,
    nce_root: str = str(DEFAULT_NCE_ROOT),
    reference_dir: str = str(DEFAULT_REFERENCE_DIR),
):
    prepare_nce_runtime(nce_root=nce_root, reference_dir=reference_dir)
    try:
        from elucidation.library_services.cachedMlimsLibraryService import CachedMlimsLibraryService

        return CachedMlimsLibraryService.get_instance()
    except ModuleNotFoundError:
        return get_cached_mlims_service_from_pickle(reference_dir=reference_dir)


class PickleCachedMlimsService:
    def __init__(self, tables: dict[str, pd.DataFrame]) -> None:
        self.chem_df = tables.get("chemicals", pd.DataFrame()).copy()
        self.cle_df = tables.get("chro_lib_entry", pd.DataFrame()).copy()
        self.comp_df = tables.get("chro_comp_entry", pd.DataFrame()).copy()

    def get_name_dict(self) -> dict[Any, str]:
        if self.comp_df.empty or "comp_id" not in self.comp_df.columns:
            return {}
        name_col = "report_name" if "report_name" in self.comp_df.columns else "name"
        names = self.comp_df.get(name_col, pd.Series("", index=self.comp_df.index)).fillna("")
        fallback = self.comp_df.get("name", pd.Series("", index=self.comp_df.index)).fillna("")
        out: dict[Any, str] = {}
        for comp_id, name, alt in zip(self.comp_df["comp_id"], names, fallback):
            value = str(name or alt or "").strip()
            if value:
                out[comp_id] = value
        return out


@lru_cache(maxsize=4)
def get_cached_mlims_service_from_pickle(
    *, reference_dir: str = str(DEFAULT_REFERENCE_DIR)
) -> PickleCachedMlimsService:
    candidates = [Path(str(reference_dir)).expanduser() / "mlims_tables.pkl", *DEFAULT_MLIMS_TABLE_CANDIDATES]
    for path in candidates:
        if not path.exists():
            continue
        with path.open("rb") as handle:
            tables = pickle.load(handle)
        if isinstance(tables, dict):
            return PickleCachedMlimsService(tables)
    raise ModuleNotFoundError("No NCE elucidation package or mlims_tables.pkl cache is available")


def nce_bin_name(bin_id: Any) -> str:
    value = int(float(bin_id))
    return f"bin_{value}" if value >= 0 else f"bin_m{-value}"


def observed_bin_spectrum(ion_bins: list[Any], loss_bins: list[Any]) -> dict[str, float]:
    bins = []
    for value in list(ion_bins or []) + list(loss_bins or []):
        try:
            intval = int(float(value))
        except Exception:
            continue
        if intval:
            bins.append(intval)
    if not bins:
        return {}
    weight = 1.0 / float(len(bins))
    return {nce_bin_name(value): weight for value in bins}


def score_smiles_against_observed_bins(
    smiles: Any,
    ion_bins: list[Any],
    loss_bins: list[Any],
    *,
    lib_id: int = 209,
    nce_root: Path | str = DEFAULT_NCE_ROOT,
    reference_dir: Path | str = DEFAULT_REFERENCE_DIR,
) -> float:
    if not str(smiles or "").strip():
        return 0.0
    spectrum = observed_bin_spectrum(ion_bins, loss_bins)
    if not spectrum:
        return 0.0
    matcher = get_predicted_spectrum_matcher(
        int(lib_id),
        nce_root=str(Path(str(nce_root)).expanduser()),
        reference_dir=str(Path(str(reference_dir)).expanduser()),
    )
    model_bins = set(spectrum) & set(matcher.model_dict)
    try:
        score = matcher.score_a_smiles_for_spectrum(str(smiles), spectrum, model_bins)
    except Exception:
        return 0.0
    if not np.isfinite(float(score)):
        return 0.0
    return float(max(score, 0.0))
