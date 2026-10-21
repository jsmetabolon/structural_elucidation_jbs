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


def resolve_reference_dir(reference_dir: Path | str | None = None) -> Path:
    """Explicit models path, then NCE_REFERENCE_DIR, then the bundled default."""
    value = reference_dir if reference_dir is not None else (
        os.environ.get("NCE_REFERENCE_DIR", "").strip() or DEFAULT_REFERENCE_DIR
    )
    return Path(str(value)).expanduser()


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
            """Guard class that blocks accidental database access in cache-only mode."""

            def __init__(self, *args: Any, **kwargs: Any) -> None:
                """Reject construction so cache-only runs never open a live DB connection."""
                raise RuntimeError("EasyQuery is unavailable in local cache-only NCE mode")

        easy_query.EasyQuery = EasyQuery
        sys.modules.setdefault("mdsutils", mdsutils)
        sys.modules["mdsutils.easyQuery"] = easy_query


def prepare_nce_runtime(
    *,
    nce_root: Path | str = DEFAULT_NCE_ROOT,
    reference_dir: Path | str | None = None,
) -> None:
    """Point NCE imports at the local checkout and model/cache directory."""
    root = Path(str(nce_root)).expanduser()
    ref = resolve_reference_dir(reference_dir)
    os.environ["REFERENCE_DIR"] = str(ref)
    if str(root) not in sys.path:
        sys.path.insert(0, str(root))
    install_local_nce_shims()


def get_predicted_spectrum_matcher(
    lib_id: int = 209,
    *,
    nce_root: str = str(DEFAULT_NCE_ROOT),
    reference_dir: Path | str | None = None,
):
    """Return the matcher using the selected data directory and repository code."""
    return _get_predicted_spectrum_matcher(
        int(lib_id), str(Path(str(nce_root)).expanduser()), str(resolve_reference_dir(reference_dir))
    )


@lru_cache(maxsize=4)
def _get_predicted_spectrum_matcher(lib_id: int, nce_root: str, reference_dir: str):
    prepare_nce_runtime(nce_root=nce_root, reference_dir=reference_dir)
    from elucidation.predictedSpectrumMatcher import PredictedSpectrumMatcher

    return PredictedSpectrumMatcher.get_instance(int(lib_id))


get_predicted_spectrum_matcher.cache_clear = _get_predicted_spectrum_matcher.cache_clear


def get_cached_mlims_service(
    *,
    nce_root: str = str(DEFAULT_NCE_ROOT),
    reference_dir: Path | str | None = None,
):
    """Return MLIMS tables using an explicit or environment-selected directory."""
    return _get_cached_mlims_service(
        str(Path(str(nce_root)).expanduser()), str(resolve_reference_dir(reference_dir))
    )


@lru_cache(maxsize=1)
def _get_cached_mlims_service(nce_root: str, reference_dir: str):
    prepare_nce_runtime(nce_root=nce_root, reference_dir=reference_dir)
    try:
        from elucidation.library_services.cachedMlimsLibraryService import CachedMlimsLibraryService

        return CachedMlimsLibraryService.get_instance()
    except ModuleNotFoundError:
        return get_cached_mlims_service_from_pickle(reference_dir=reference_dir)


get_cached_mlims_service.cache_clear = _get_cached_mlims_service.cache_clear


class PickleCachedMlimsService:
    """Minimal MLIMS name lookup backed by pickled reference tables."""

    def __init__(self, tables: dict[str, pd.DataFrame]) -> None:
        """Load chemical, library-entry, and compound tables from a pickle dict."""
        self.chem_df = tables.get("chemicals", pd.DataFrame()).copy()
        self.cle_df = tables.get("chro_lib_entry", pd.DataFrame()).copy()
        self.comp_df = tables.get("chro_comp_entry", pd.DataFrame()).copy()

    def get_name_dict(self) -> dict[Any, str]:
        """Map ``comp_id`` values to report or fallback compound names."""
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


def get_cached_mlims_service_from_pickle(
    *, reference_dir: Path | str | None = None
) -> PickleCachedMlimsService:
    """Load MLIMS tables from the first available ``mlims_tables.pkl`` cache."""
    return _get_cached_mlims_service_from_pickle(str(resolve_reference_dir(reference_dir)))


@lru_cache(maxsize=4)
def _get_cached_mlims_service_from_pickle(reference_dir: str) -> PickleCachedMlimsService:
    candidates = [Path(reference_dir) / "mlims_tables.pkl", *DEFAULT_MLIMS_TABLE_CANDIDATES]
    for path in candidates:
        if not path.exists():
            continue
        with path.open("rb") as handle:
            tables = pickle.load(handle)
        if isinstance(tables, dict):
            return PickleCachedMlimsService(tables)
    raise ModuleNotFoundError("No NCE elucidation package or mlims_tables.pkl cache is available")


get_cached_mlims_service_from_pickle.cache_clear = _get_cached_mlims_service_from_pickle.cache_clear


def nce_bin_name(bin_id: Any) -> str:
    """Format a bin id the way NCE spectrum dictionaries name keys."""
    value = int(float(bin_id))
    return f"bin_{value}" if value >= 0 else f"bin_m{-value}"


def observed_bin_spectrum(ion_bins: list[Any], loss_bins: list[Any]) -> dict[str, float]:
    """Build a uniform-weight NCE spectrum dict from observed fragment bins."""
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
    reference_dir: Path | str | None = None,
) -> float:
    """Score one SMILES against observed bins via the real NCE spectrum matcher."""
    if not str(smiles or "").strip():
        return 0.0
    spectrum = observed_bin_spectrum(ion_bins, loss_bins)
    if not spectrum:
        return 0.0
    matcher = get_predicted_spectrum_matcher(
        int(lib_id),
        nce_root=str(Path(str(nce_root)).expanduser()),
        reference_dir=str(resolve_reference_dir(reference_dir)),
    )
    model_bins = set(spectrum) & set(matcher.model_dict)
    try:
        score = matcher.score_a_smiles_for_spectrum(str(smiles), spectrum, model_bins)
    except Exception:
        return 0.0
    if not np.isfinite(float(score)):
        return 0.0
    return float(max(score, 0.0))
