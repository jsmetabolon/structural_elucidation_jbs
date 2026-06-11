from __future__ import annotations

import pickle
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any

import pandas as pd


REQUIRED_IONS_PKL_COLUMNS = (
    "comp_id",
    "mass",
    "intensity",
    "ri",
    "rt",
    "compound_name",
    "parent_mz",
    "feature_key",
)

PEAK_SIGNATURE_MASS_DP = 5
PEAK_SIGNATURE_INTENSITY_DP = 3
PARENT_MZ_MATCH_ATOL_DA = 1e-4


def robust_pickle_load(path: str | Path) -> Any:
    try:
        return pd.read_pickle(path)
    except Exception:
        with Path(path).expanduser().resolve().open("rb") as handle:
            return pickle.load(handle)


def _bridge_python_candidates() -> list[str]:
    candidates: list[str] = []
    for raw in (
        os.environ.get("SPEC2MOL_PICKLE_BRIDGE_PYTHON"),
        "/opt/homebrew/bin/python3",
        shutil.which("python3"),
        shutil.which("python"),
        sys.executable,
    ):
        if not raw:
            continue
        value = str(raw).strip()
        if not value or value in candidates:
            continue
        candidates.append(value)
    return candidates


def _load_pickle_via_bridge_python(resolved: Path) -> pd.DataFrame:
    with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as handle:
        temp_json = Path(handle.name)
    try:
        last_error = "no bridge python candidates attempted"
        for python_bin in _bridge_python_candidates():
            try:
                subprocess.run(
                    [
                        python_bin,
                        "-c",
                        (
                            "import pandas as pd, sys; "
                            "obj = pd.read_pickle(sys.argv[1]); "
                            "assert hasattr(obj, 'to_json'), type(obj); "
                            "obj.to_json(sys.argv[2], orient='table')"
                        ),
                        str(resolved),
                        str(temp_json),
                    ],
                    check=True,
                    capture_output=True,
                    text=True,
                )
                return pd.read_json(temp_json, orient="table")
            except Exception as exc:
                last_error = f"{python_bin}: {exc}"
        raise RuntimeError(
            f"Unable to bridge-load pickle via alternate Python interpreters for {resolved}. Last error: {last_error}"
        )
    finally:
        temp_json.unlink(missing_ok=True)


def load_pickle_dataframe_compat(path: str | Path) -> pd.DataFrame:
    resolved = Path(path).expanduser().resolve()
    try:
        obj = robust_pickle_load(resolved)
    except Exception:
        obj = _load_pickle_via_bridge_python(resolved)
    if not hasattr(obj, "columns"):
        raise TypeError(f"Expected DataFrame pickle, got {type(obj)} for {resolved}")
    return clean_columns(obj)


def clean_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out.columns = [str(col).strip().strip('"').lstrip("\ufeff") for col in out.columns]
    return out


def validate_ions_pkl_columns(df: pd.DataFrame, path: str | Path) -> None:
    missing = [col for col in REQUIRED_IONS_PKL_COLUMNS if col not in df.columns]
    if missing:
        raise KeyError(
            f"Ions PKL missing required columns {missing}: {Path(path).expanduser().resolve()} -> {df.columns.tolist()}"
        )


def load_validated_ions_pkl(path: str | Path) -> pd.DataFrame:
    resolved = Path(path).expanduser().resolve()
    df = load_pickle_dataframe_compat(resolved)
    validate_ions_pkl_columns(df, resolved)
    out = df.copy()
    out["comp_id"] = pd.to_numeric(out["comp_id"], errors="coerce")
    out["mass"] = pd.to_numeric(out["mass"], errors="coerce")
    out["intensity"] = pd.to_numeric(out["intensity"], errors="coerce")
    out["ri"] = pd.to_numeric(out["ri"], errors="coerce")
    out["rt"] = pd.to_numeric(out["rt"], errors="coerce")
    out["parent_mz"] = pd.to_numeric(out["parent_mz"], errors="coerce")
    out["compound_name"] = out["compound_name"].where(out["compound_name"].notna(), "")
    out["compound_name"] = out["compound_name"].astype(str).str.strip()
    out["feature_key"] = out["feature_key"].where(out["feature_key"].notna(), "")
    out["feature_key"] = out["feature_key"].astype(str).str.strip()
    out = out.dropna(subset=["comp_id", "mass", "intensity", "ri", "rt", "parent_mz"]).copy()
    out["comp_id"] = out["comp_id"].astype(int)
    return out


def load_library_comp_ids(path: str | Path) -> set[int]:
    obj = robust_pickle_load(path)
    if isinstance(obj, tuple):
        if not obj or not hasattr(obj[0], "columns"):
            raise TypeError(f"Unexpected library pickle structure: {type(obj)}")
        df = clean_columns(obj[0])
    elif hasattr(obj, "columns"):
        df = clean_columns(obj)
    else:
        raise TypeError(f"Unexpected library pickle payload: {type(obj)}")
    if "comp_id" not in df.columns:
        raise KeyError(f"Library pickle missing comp_id column: {Path(path).expanduser().resolve()}")
    series = pd.to_numeric(df["comp_id"], errors="coerce").dropna().astype(int)
    return set(series.tolist())


def load_parent_mz_map(path: str | Path) -> dict[int, float]:
    obj = robust_pickle_load(path)
    if not isinstance(obj, dict):
        raise TypeError(f"Expected dict from parent_mz pickle, got {type(obj)}")
    return {int(cid): float(value) for cid, value in obj.items()}


def list_embedding_comp_ids(embeddings_dir: str | Path) -> set[int]:
    emb_dir = Path(embeddings_dir).expanduser().resolve()
    out: set[int] = set()
    for path in emb_dir.glob("*.pt"):
        try:
            out.add(int(path.stem))
        except ValueError:
            continue
    return out


def peaks_signature_by_comp_id(df: pd.DataFrame) -> dict[int, tuple[tuple[float, float], ...]]:
    out: dict[int, tuple[tuple[float, float], ...]] = {}
    for comp_id, group in df.groupby("comp_id", sort=True):
        pairs = [
            (
                round(float(mass), PEAK_SIGNATURE_MASS_DP),
                round(float(intensity), PEAK_SIGNATURE_INTENSITY_DP),
            )
            for mass, intensity in zip(group["mass"], group["intensity"])
        ]
        out[int(comp_id)] = tuple(sorted(pairs))
    return out


def _sorted_int_list(values: set[int]) -> list[int]:
    return sorted(int(value) for value in values)


def compute_supported_asset_comp_ids(
    *,
    library_pkl: str | Path,
    parent_mz_pkl: str | Path,
    embeddings_dir: str | Path,
) -> set[int]:
    library_ids = load_library_comp_ids(library_pkl)
    parent_ids = set(load_parent_mz_map(parent_mz_pkl))
    embedding_ids = list_embedding_comp_ids(embeddings_dir)
    return library_ids & parent_ids & embedding_ids


def build_supported_known_ions_artifact(
    *,
    known_ions_pkl: str | Path,
    library_pkl: str | Path,
    parent_mz_pkl: str | Path,
    embeddings_dir: str | Path,
    out_pkl: str | Path | None = None,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    known_path = Path(known_ions_pkl).expanduser().resolve()
    supported_ids = compute_supported_asset_comp_ids(
        library_pkl=library_pkl,
        parent_mz_pkl=parent_mz_pkl,
        embeddings_dir=embeddings_dir,
    )
    known_df = load_validated_ions_pkl(known_path)
    filtered_df = known_df.loc[known_df["comp_id"].isin(supported_ids)].copy().reset_index(drop=True)

    out_path = None if out_pkl is None else Path(out_pkl).expanduser().resolve()
    if out_path is not None:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        filtered_df.to_pickle(out_path)

    summary: dict[str, Any] = {
        "known_ions_pkl": str(known_path),
        "library_pkl": str(Path(library_pkl).expanduser().resolve()),
        "parent_mz_pkl": str(Path(parent_mz_pkl).expanduser().resolve()),
        "embeddings_dir": str(Path(embeddings_dir).expanduser().resolve()),
        "supported_asset_comp_ids": int(len(supported_ids)),
        "input_rows": int(len(known_df)),
        "input_unique_comp_ids": int(known_df["comp_id"].nunique()),
        "output_rows": int(len(filtered_df)),
        "output_unique_comp_ids": int(filtered_df["comp_id"].nunique()),
        "output_pkl": None if out_path is None else str(out_path),
    }
    return filtered_df, summary


def build_known_ions_refresh_assessment(
    *,
    known_ions_pkl: str | Path,
    legacy_known_ions_pkl: str | Path,
    library_pkl: str | Path,
    parent_mz_pkl: str | Path,
    embeddings_dir: str | Path,
) -> dict[str, Any]:
    known_path = Path(known_ions_pkl).expanduser().resolve()
    legacy_path = Path(legacy_known_ions_pkl).expanduser().resolve()
    known_df = load_validated_ions_pkl(known_path)

    legacy_obj = robust_pickle_load(legacy_path)
    if not hasattr(legacy_obj, "columns"):
        raise TypeError(f"Expected DataFrame from legacy ions PKL, got {type(legacy_obj)} for {legacy_path}")
    legacy_df = clean_columns(legacy_obj)
    required_legacy = {"comp_id", "mass", "intensity"}
    missing_legacy = sorted(required_legacy - set(legacy_df.columns))
    if missing_legacy:
        raise KeyError(
            f"Legacy known ions PKL missing required columns {missing_legacy}: {legacy_path} -> {legacy_df.columns.tolist()}"
        )
    legacy_df = legacy_df.copy()
    legacy_df["comp_id"] = pd.to_numeric(legacy_df["comp_id"], errors="coerce")
    legacy_df["mass"] = pd.to_numeric(legacy_df["mass"], errors="coerce")
    legacy_df["intensity"] = pd.to_numeric(legacy_df["intensity"], errors="coerce")
    legacy_df = legacy_df.dropna(subset=["comp_id", "mass", "intensity"]).copy()
    legacy_df["comp_id"] = legacy_df["comp_id"].astype(int)

    library_ids = load_library_comp_ids(library_pkl)
    parent_map = load_parent_mz_map(parent_mz_pkl)
    parent_ids = set(parent_map)
    embedding_ids = list_embedding_comp_ids(embeddings_dir)

    known_ids = set(known_df["comp_id"].unique().tolist())
    legacy_ids = set(legacy_df["comp_id"].unique().tolist())
    supported_assets = compute_supported_asset_comp_ids(
        library_pkl=library_pkl,
        parent_mz_pkl=parent_mz_pkl,
        embeddings_dir=embeddings_dir,
    )
    known_supported_ids = known_ids & supported_assets
    legacy_supported_ids = legacy_ids & supported_assets

    known_supported_df = known_df.loc[known_df["comp_id"].isin(known_supported_ids), ["comp_id", "mass", "intensity"]].copy()
    legacy_supported_df = legacy_df.loc[legacy_df["comp_id"].isin(legacy_supported_ids), ["comp_id", "mass", "intensity"]].copy()
    known_signatures = peaks_signature_by_comp_id(known_supported_df)
    legacy_signatures = peaks_signature_by_comp_id(legacy_supported_df)
    changed_supported_ids = sorted(
        cid
        for cid in sorted(known_supported_ids & legacy_supported_ids)
        if known_signatures.get(cid) != legacy_signatures.get(cid)
    )

    per_cid_parent = (
        known_df.groupby("comp_id", sort=True)["parent_mz"]
        .agg(["nunique", "first"])
        .reset_index()
        .rename(columns={"nunique": "parent_mz_nunique", "first": "known_parent_mz"})
    )
    per_cid_parent["comp_id"] = per_cid_parent["comp_id"].astype(int)
    parent_overlap_df = per_cid_parent.loc[per_cid_parent["comp_id"].isin(parent_ids)].copy()
    parent_overlap_df["map_parent_mz"] = parent_overlap_df["comp_id"].map(parent_map)
    parent_overlap_df["parent_mz_matches_existing_map"] = (
        parent_overlap_df["known_parent_mz"] - parent_overlap_df["map_parent_mz"]
    ).abs().le(PARENT_MZ_MATCH_ATOL_DA)
    parent_mismatch_df = parent_overlap_df.loc[~parent_overlap_df["parent_mz_matches_existing_map"]].copy()

    supported_ids_changed = known_supported_ids != legacy_supported_ids
    supported_spectra_changed = bool(changed_supported_ids)
    retrain_required = bool(supported_ids_changed or supported_spectra_changed)

    return {
        "known_ions_pkl": str(known_path),
        "legacy_known_ions_pkl": str(legacy_path),
        "library_pkl": str(Path(library_pkl).expanduser().resolve()),
        "parent_mz_pkl": str(Path(parent_mz_pkl).expanduser().resolve()),
        "embeddings_dir": str(Path(embeddings_dir).expanduser().resolve()),
        "known_unique_comp_ids": int(len(known_ids)),
        "legacy_unique_comp_ids": int(len(legacy_ids)),
        "library_unique_comp_ids": int(len(library_ids)),
        "parent_mz_unique_comp_ids": int(len(parent_ids)),
        "embedding_unique_comp_ids": int(len(embedding_ids)),
        "supported_asset_comp_ids": int(len(supported_assets)),
        "known_supported_comp_ids": int(len(known_supported_ids)),
        "legacy_supported_comp_ids": int(len(legacy_supported_ids)),
        "supported_comp_ids_added_vs_legacy": _sorted_int_list(known_supported_ids - legacy_supported_ids),
        "supported_comp_ids_removed_vs_legacy": _sorted_int_list(legacy_supported_ids - known_supported_ids),
        "supported_spectra_changed_vs_legacy": changed_supported_ids,
        "known_only_comp_ids_outside_supported_assets": int(len(known_ids - supported_assets)),
        "legacy_only_comp_ids_outside_supported_assets": int(len(legacy_ids - supported_assets)),
        "parent_mz_singleton_comp_ids": int((per_cid_parent["parent_mz_nunique"] <= 1).sum()),
        "parent_mz_multi_value_comp_ids": int((per_cid_parent["parent_mz_nunique"] > 1).sum()),
        "parent_mz_overlap_with_existing_map": int(len(parent_overlap_df)),
        "parent_mz_map_exact_matches": int(parent_overlap_df["parent_mz_matches_existing_map"].sum()),
        "parent_mz_map_mismatches": int(len(parent_mismatch_df)),
        "parent_mz_map_mismatch_samples": [
            {
                "comp_id": int(row.comp_id),
                "known_parent_mz": float(row.known_parent_mz),
                "existing_parent_mz": float(row.map_parent_mz),
            }
            for row in parent_mismatch_df.head(10).itertuples(index=False)
        ],
        "supported_ids_changed_vs_legacy": bool(supported_ids_changed),
        "supported_spectra_changed_vs_legacy": bool(supported_spectra_changed),
        "retrain_required": retrain_required,
        "recommended_action": (
            "stop_before_training_and_use_existing_encoder"
            if not retrain_required
            else "rebuild_training_inputs_and_retrain_encoder"
        ),
    }
