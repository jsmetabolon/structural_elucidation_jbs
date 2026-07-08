#!/usr/bin/env python3
"""Build PCM per-lib pipeline inputs from the local/global CHRO-LIB reference spectra.

The PCM diagnosis seed (``docs/pcm/pcm_diagnosis_seed.csv``) keys controls by
Metabolon ``chemical_id``. The elucidation pipeline (candidates -> verify ->
pool -> rerank) is keyed by CHRO-LIB ``comp_id`` + ``lib_id`` and consumes an
observed-MS2 lookup, feature rows, and neutral-mass hypotheses. This tool
bridges the two using the MLIMS reference tables:

    chemical_id --chro_lib_entry--> (comp_id, lib_id, chro_lib_entry_id)
    (comp_id, lib_id, chro_lib_entry_id) --chro_lib_masses--> reference MS2 peaks

For each PCM entry it resolves the library it was profiled in (209/305 negative,
400/402 positive), reads the reference fragment spectrum, and emits, per
lib/ion-mode scope, exactly the artifacts ``run_inference.sh`` expects:

    data/observed_ms2_lookups/pcm_lib{L}_{ion}_observed_ms2_lookup.csv
    data/feature_rows/pcm_lib{L}_{ion}_feature_rows.csv
    data/feature_rows/pcm_lib{L}_{ion}_neutral_mass_hypotheses.csv
    docs/pcm/pcm_diagnosis_lib{L}_{ion}.csv

plus a single per-entry mapping table ``docs/pcm/pcm_entry_lib_ion_mode.csv``.

Bins are computed with ``structural_elucidation.nce.evidence`` so they are
identical to what ``pipeline.verify`` scores against. Read-only against the
MLIMS pickle; safe to run locally or on EC2.

Usage:
    PYTHONPATH=. python3 tools/build_pcm_ms2_inputs.py \
        [--mlims-pkl PATH] [--diagnosis docs/pcm/pcm_diagnosis_seed.csv] \
        [--tasks docs/pcm/pcm_tasks.csv]
"""
from __future__ import annotations

import argparse
import json
import os
import pickle
import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from structural_elucidation.nce.evidence import (
    NEGATIVE_LIBS,
    POSITIVE_LIBS,
    ion_and_loss_bins,
    neutral_mass_from_precursor,
)
from structural_elucidation.paths import REPO_ROOT

# CHRO-LIB records the molecular ion of a library entry as an "m-H"/"m+H" peak.
_BASE_ION_RE = re.compile(r"^\s*m\s*[+-]\s*h\s*$", re.IGNORECASE)
# Isotopologue confirmation peaks are ~1 Da above the parent; exclude when
# inferring the precursor by max mass so they cannot masquerade as the parent.
_ISOTOPE_RE = re.compile(r"\[|13C|C13|18O|15N|34S|37Cl|2H", re.IGNORECASE)


def ion_mode_for_lib(lib_id: int) -> str:
    """Map a CHRO-LIB library id to its ion mode (209/305 neg, 400/402 pos)."""
    if int(lib_id) in POSITIVE_LIBS:
        return "positive"
    if int(lib_id) in NEGATIVE_LIBS:
        return "negative"
    return "negative"


def resolve_mlims_pkl(explicit: str | None) -> Path:
    """Find ``mlims_tables.pkl`` from an explicit path, REFERENCE_DIR, or fallbacks."""
    candidates: list[Path] = []
    if explicit:
        candidates.append(Path(explicit).expanduser())
    ref = os.environ.get("REFERENCE_DIR")
    if ref:
        candidates.append(Path(ref).expanduser() / "mlims_tables.pkl")
    candidates += [
        REPO_ROOT / "vendor" / "modelfiles" / "mlims_tables.pkl",
        Path("~/dev/modelfiles/mlims_tables.pkl").expanduser(),
    ]
    for path in candidates:
        if path.is_file():
            return path
    raise SystemExit(
        "mlims_tables.pkl not found; pass --mlims-pkl or set REFERENCE_DIR. Tried:\n  " +
        "\n  ".join(str(c) for c in candidates)
    )


def num(series: pd.Series) -> pd.Series:
    """Coerce a Series to numeric with NaN for unparseable values."""
    return pd.to_numeric(series, errors="coerce")


def to_int(value: Any) -> int | None:
    """Coerce a single cell to int, returning None when missing/unparseable."""
    parsed = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
    return int(parsed) if pd.notna(parsed) else None


def choose_entries(cle: pd.DataFrame, chemical_ids: set[int]) -> pd.DataFrame:
    """Pick one current CHRO-LIB entry per (chemical_id, lib_id).

    Prefers ``version_current`` truthy, then the most recently modified row, so a
    chemical measured in several libraries yields one entry per library.
    """
    work = cle.copy()
    work["_chem"] = num(work["chemical_id"])
    work = work[work["_chem"].isin(chemical_ids)].copy()
    if work.empty:
        return work
    work["_lib"] = num(work["lib_id"]).astype("Int64")
    vc = work.get("version_current", pd.Series(False, index=work.index))
    work["_current"] = vc.astype(str).str.lower().isin({"1", "1.0", "true", "t", "yes", "y"})
    work["_modified"] = pd.to_datetime(work.get("modified_on"), errors="coerce")
    work = work.sort_values(["_chem", "_lib", "_current", "_modified"], ascending=[True, True, False, False])
    return work.drop_duplicates(["_chem", "_lib"], keep="first").reset_index(drop=True)


def precursor_mz(peaks: pd.DataFrame) -> float:
    """Infer the precursor m/z: the ``m-H``/``m+H`` peak, else max non-isotope mass."""
    names = peaks["name"].fillna("").astype(str)
    base = peaks[names.map(lambda v: bool(_BASE_ION_RE.match(v)))]
    if not base.empty:
        return float(num(base["mass"]).max())
    non_iso = peaks[~names.map(lambda v: bool(_ISOTOPE_RE.search(v)))]
    pool = non_iso if not non_iso.empty else peaks
    value = num(pool["mass"]).max()
    return float(value) if pd.notna(value) else float("nan")


def spectrum_bins(peaks: pd.DataFrame, parent_mz: float, lib_id: int) -> dict[str, Any]:
    """Convert reference peaks into verifier-compatible ion/loss bins + pairs."""
    ion_int: dict[int, float] = {}
    loss_int: dict[int, float] = {}
    pairs: list[list[float]] = []
    for _, row in peaks.iterrows():
        mz = float(num(pd.Series([row.get("mass")])).iloc[0]) if row.get("mass") is not None else float("nan")
        if not np.isfinite(mz) or mz <= 0:
            continue
        ratio = float(num(pd.Series([row.get("mass_ratio")])).fillna(1.0).iloc[0])
        ratio = ratio if np.isfinite(ratio) and ratio > 0 else 1.0
        pairs.append([round(mz, 6), ratio])
        ion_bin, loss_bin = ion_and_loss_bins(mz, parent_mz, int(lib_id))
        if ion_bin:
            ion_int[int(ion_bin)] = max(ratio, ion_int.get(int(ion_bin), 0.0))
        if loss_bin:
            loss_int[int(loss_bin)] = max(ratio, loss_int.get(int(loss_bin), 0.0))
    max_ratio = max([r for _, r in pairs], default=1.0) or 1.0
    ions = sorted(ion_int)
    losses = sorted(loss_int)
    return {
        "ion_bins": ions,
        "loss_bins": losses,
        "ion_intensities": [round(min(1.0, ion_int[b] / max_ratio), 6) for b in ions],
        "loss_intensities": [round(min(1.0, loss_int[b] / max_ratio), 6) for b in losses],
        "ms2_pairs": [[mz, round(min(1.0, r / max_ratio), 6)] for mz, r in sorted(pairs)],
    }


def build(args: argparse.Namespace) -> dict[str, Any]:
    """Resolve PCM entries, read reference spectra, and write per-lib inputs."""
    pkl = resolve_mlims_pkl(args.mlims_pkl)
    tables = pickle.load(pkl.open("rb"))
    cle = tables["chro_lib_entry"]
    clm = tables["chro_lib_masses"]
    clm = clm.assign(_comp=num(clm["comp_id"]).astype("Int64"), _lib=num(clm["lib_id"]).astype("Int64"))
    if "chro_lib_entry_id" in clm.columns:
        clm = clm.assign(_entry=num(clm["chro_lib_entry_id"]).astype("Int64"))

    tasks = pd.read_csv(args.tasks)
    diag = pd.read_csv(args.diagnosis)
    tasks["_chem"] = num(tasks["chemical_id"]).astype("Int64")
    diag["_chem"] = num(diag["control_id"]).astype("Int64")
    targeted = set(int(v) for v in diag["_chem"].dropna())
    name_by_chem = dict(zip(tasks["_chem"], tasks["canonical_name"].astype(str)))
    hint_by_chem = dict(zip(diag["_chem"], diag["conjugate_hints"].fillna("").astype(str)))

    entries = choose_entries(cle, set(int(v) for v in tasks["_chem"].dropna()))
    diag_by_chem = diag.set_index("_chem")

    mapping_rows: list[dict[str, Any]] = []
    scope_obs: dict[tuple[int, str], list[dict[str, Any]]] = {}
    scope_feat: dict[tuple[int, str], list[dict[str, Any]]] = {}
    scope_hyp: dict[tuple[int, str], list[dict[str, Any]]] = {}
    scope_diag: dict[tuple[int, str], list[dict[str, Any]]] = {}

    for _, ent in entries.iterrows():
        chem = int(ent["_chem"])
        lib_id = to_int(ent["_lib"])
        if lib_id is None:
            continue
        comp_id = to_int(ent.get("comp_id"))
        entry_id = to_int(ent.get("chro_lib_entry_id"))
        ion = ion_mode_for_lib(lib_id)
        scope = (lib_id, ion)

        peaks = clm[clm["_comp"].eq(comp_id) & clm["_lib"].eq(lib_id)]
        if entry_id is not None and "_entry" in clm.columns:
            scoped = peaks[peaks["_entry"].eq(entry_id)]
            peaks = scoped if not scoped.empty else peaks

        parent = precursor_mz(peaks) if not peaks.empty else float("nan")
        neutral = neutral_mass_from_precursor(parent, lib_id)
        bins = spectrum_bins(peaks, parent, lib_id) if (not peaks.empty and np.isfinite(parent)) else {
            "ion_bins": [], "loss_bins": [], "ion_intensities": [], "loss_intensities": [], "ms2_pairs": []
        }
        has_ms2 = bool(bins["ion_bins"] or bins["loss_bins"])
        feature_key = f"pcm_{lib_id}_{chem}"
        control_id = str(chem)

        mapping_rows.append({
            "chemical_id": chem,
            "canonical_name": name_by_chem.get(chem, ""),
            "comp_id": comp_id if comp_id is not None else "",
            "lib_id": lib_id,
            "ion_mode": ion,
            "rt": ent.get("rt", ""),
            "ri": ent.get("ri", ""),
            "status": ent.get("status", ""),
            "parent_mz": round(parent, 6) if np.isfinite(parent) else "",
            "neutral_mass": round(neutral, 6) if np.isfinite(neutral) else "",
            "n_ref_peaks": int(len(peaks)),
            "n_ion_bins": int(len(bins["ion_bins"])),
            "n_loss_bins": int(len(bins["loss_bins"])),
            "has_ms2": has_ms2,
            "targeted": chem in targeted,
            "conjugate_hint": hint_by_chem.get(chem, ""),
        })

        if chem not in targeted:
            continue  # only the 92 unresolved controls flow through the pipeline

        scope_obs.setdefault(scope, []).append({
            "control_id": control_id,
            "lib_id": lib_id,
            "ion_mode": ion,
            "feature_key": feature_key,
            "comp_id": comp_id if comp_id is not None else "",
            "parent_mz": round(parent, 6) if np.isfinite(parent) else "",
            "neutral_mass": round(neutral, 6) if np.isfinite(neutral) else "",
            "observed_ion_bins_json": json.dumps(bins["ion_bins"]),
            "observed_loss_bins_json": json.dumps(bins["loss_bins"]),
            "observed_ion_intensities_json": json.dumps(bins["ion_intensities"]),
            "observed_loss_intensities_json": json.dumps(bins["loss_intensities"]),
            "ms2_pairs_json": json.dumps(bins["ms2_pairs"]),
            "n_observed_ion_bins": int(len(bins["ion_bins"])),
            "n_observed_loss_bins": int(len(bins["loss_bins"])),
            "n_ms2_pairs": int(len(bins["ms2_pairs"])),
            "ms2_present": has_ms2,
            "spectrum_source": "chro_lib_masses",
        })
        scope_feat.setdefault(scope, []).append({
            "feature_key": feature_key,
            "comp_id": comp_id if comp_id is not None else "",
            "control_id": control_id,
            "lib_id": lib_id,
            "ion_mode": ion,
            "parent_mz": round(parent, 6) if np.isfinite(parent) else "",
            "neutral_mass_fallback": round(neutral, 6) if np.isfinite(neutral) else "",
            "n_ms2_peaks": int(len(bins["ms2_pairs"])),
            "observed_ion_bins_json": json.dumps(bins["ion_bins"]),
            "observed_loss_bins_json": json.dumps(bins["loss_bins"]),
        })
        proton = 1.007276466812
        scope_hyp.setdefault(scope, []).append({
            "feature_key": feature_key,
            "comp_id": comp_id if comp_id is not None else "",
            "lib_id": lib_id,
            "ion_mode": ion,
            "parent_mz": round(parent, 6) if np.isfinite(parent) else "",
            "adduct_label": "m-H" if ion == "negative" else "m+H",
            "charge": -1 if ion == "negative" else 1,
            "multiplier": 1,
            "offset_mass": round(proton if ion == "negative" else -proton, 6),
            "neutral_mass": round(neutral, 6) if np.isfinite(neutral) else "",
            "recomputed_mz": round(parent, 6) if np.isfinite(parent) else "",
            "adduct_score": 1.0,
            "prevalence": "",
            "prevalence_tier": "chro_lib_precursor",
            "source_reference": "chro_lib_masses",
            "neutral_loss_possible": False,
            "note": "pcm_reference_precursor",
            "adduct_rank": 1,
            "control_id": control_id,
        })
        drow = diag_by_chem.loc[chem] if chem in diag_by_chem.index else None
        if drow is not None:
            drow = drow.iloc[0] if isinstance(drow, pd.DataFrame) else drow
            scope_diag.setdefault(scope, []).append({
                "control_id": control_id,
                "primary_stage": drow.get("primary_stage", "pcm_unresolved"),
                "conjugate_hints": drow.get("conjugate_hints", ""),
                "core_formula": drow.get("core_formula", ""),
                "core_monoisotopic_mass": drow.get("core_monoisotopic_mass", ""),
                "isomer_index": drow.get("isomer_index", ""),
                "canonical_name": drow.get("canonical_name", ""),
                "comp_id": comp_id if comp_id is not None else "",
                "lib_id": lib_id,
                "ion_mode": ion,
                "parent_mz": round(parent, 6) if np.isfinite(parent) else "",
                "observed_neutral_mass": round(neutral, 6) if np.isfinite(neutral) else "",
            })

    out_obs = Path(args.out_observed)
    out_feat = Path(args.out_features)
    out_docs = Path(args.out_docs)
    for d in (out_obs, out_feat, out_docs):
        d.mkdir(parents=True, exist_ok=True)

    written: list[str] = []
    for (lib_id, ion), rows in sorted(scope_obs.items()):
        p = out_obs / f"pcm_lib{lib_id}_{ion}_observed_ms2_lookup.csv"
        pd.DataFrame(rows).to_csv(p, index=False)
        written.append(str(p))
    for (lib_id, ion), rows in sorted(scope_feat.items()):
        p = out_feat / f"pcm_lib{lib_id}_{ion}_feature_rows.csv"
        pd.DataFrame(rows).to_csv(p, index=False)
        written.append(str(p))
    for (lib_id, ion), rows in sorted(scope_hyp.items()):
        p = out_feat / f"pcm_lib{lib_id}_{ion}_neutral_mass_hypotheses.csv"
        pd.DataFrame(rows).to_csv(p, index=False)
        written.append(str(p))
    for (lib_id, ion), rows in sorted(scope_diag.items()):
        p = out_docs / f"pcm_diagnosis_lib{lib_id}_{ion}.csv"
        pd.DataFrame(rows).to_csv(p, index=False)
        written.append(str(p))

    mapping = pd.DataFrame(mapping_rows).sort_values(
        ["targeted", "lib_id", "chemical_id"], ascending=[False, True, True]
    )
    mapping_path = out_docs / "pcm_entry_lib_ion_mode.csv"
    mapping.to_csv(mapping_path, index=False)
    written.append(str(mapping_path))

    lib_counts = mapping.groupby(["lib_id", "ion_mode"]).size().to_dict()
    summary = {
        "mlims_pkl": str(pkl),
        "n_pcm_chemicals": int(mapping["chemical_id"].nunique()),
        "n_entries": int(len(mapping)),
        "n_entries_with_ms2": int(mapping["has_ms2"].sum()),
        "n_targeted_entries": int(mapping["targeted"].sum()),
        "n_targeted_with_ms2": int((mapping["targeted"] & mapping["has_ms2"]).sum()),
        "entries_by_scope": {f"lib{k[0]}_{k[1]}": int(v) for k, v in lib_counts.items()},
        "mapping_csv": str(mapping_path),
        "n_files_written": len(written),
    }
    print(json.dumps(summary, indent=2, sort_keys=True))
    for w in written:
        print(f"  wrote {w}")
    return summary


def main() -> None:
    """Parse CLI arguments and build PCM per-lib pipeline inputs."""
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mlims-pkl", default=None, help="Path to mlims_tables.pkl (default: REFERENCE_DIR / fallbacks)")
    ap.add_argument("--tasks", type=Path, default=Path("docs/pcm/pcm_tasks.csv"))
    ap.add_argument("--diagnosis", type=Path, default=Path("docs/pcm/pcm_diagnosis_seed.csv"))
    ap.add_argument("--out-observed", type=Path, default=Path("data/observed_ms2_lookups"))
    ap.add_argument("--out-features", type=Path, default=Path("data/feature_rows"))
    ap.add_argument("--out-docs", type=Path, default=Path("docs/pcm"))
    build(ap.parse_args())


if __name__ == "__main__":
    main()
