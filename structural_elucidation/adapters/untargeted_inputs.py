"""Build structural-elucidation input tables from an untargeted-pipeline handoff.

The untargeted metabolomics pipeline (``untargeted_workflow``) emits a neutral,
dependency-light handoff describing features that were *not* assigned a
Metabolon/CURATION identity but do have consensus MS2 evidence.  This module
turns that handoff into the per-library CSVs the elucidation stack consumes:

* ``lib{lib}_{ion}_feature_rows.csv``          (candidates + audit)
* ``lib{lib}_{ion}_neutral_mass_hypotheses.csv`` (candidates)
* ``lib{lib}_{ion}_gap_manifest.csv``          (candidates diagnosis/control set)
* ``lib{lib}_{ion}_observed_ms2_lookup.csv``   (verify + report MS2 panels)
* ``ions_novel_lib{lib}_{ion}.csv``            (optional encoder / de novo input)
* ``all_x_catalog.csv``                        (forces the pool non-S0 split)
* ``integrated_per_control_best.csv``          (pool/verify control allow-list stub)
* ``novel_feature_map.csv``                    (feature_id <-> control_id mapping)

Synthetic ``control_id`` values are ``X-<NOVEL_CONTROL_OFFSET + feature_id>`` so
they match the ``X-\\d+`` pattern the pool/verify/rerank code expects while never
colliding with real Metabolon X-compounds.

Handoff JSON schema (``--handoff``)::

    {
      "schema_version": 1,
      "source": "untargeted_workflow",
      "features": [
        {"feature_id": 123, "lib_id": 400, "ion_mode": "positive",
         "parent_mz": 250.093, "ri": 2486.0,
         "peaks": [[mz, intensity], ...],
         "compound_name": "", "prior_match_level": "MS2",
         "prior_matched_name": "..."}
      ]
    }
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Optional

import pandas as pd

from structural_elucidation.adapters.nce_bins import (
    ion_mode_for_lib,
    nce_bin_name,
    observed_bins_for_spectrum,
)

# Keep synthetic control ids well clear of real Metabolon X-compound numbers.
NOVEL_CONTROL_OFFSET: int = 900_000_000
PROTON_MASS: float = 1.007276466812
SUPPORTED_LIBS: frozenset[int] = frozenset({209, 305, 400, 402})


def control_id_for_feature(feature_id: int) -> str:
    """Return the synthetic ``X-<n>`` control id for an untargeted feature id."""
    return f"X-{NOVEL_CONTROL_OFFSET + int(feature_id)}"


def feature_id_from_control_id(control_id: str) -> int | None:
    """Invert :func:`control_id_for_feature`; ``None`` when not a novel id."""
    text = str(control_id or "").strip()
    if not text.startswith("X-"):
        return None
    try:
        value = int(text[2:])
    except ValueError:
        return None
    if value < NOVEL_CONTROL_OFFSET:
        return None
    return value - NOVEL_CONTROL_OFFSET


def neutral_mass_from_precursor(parent_mz: float, lib_id: int) -> float:
    """Neutral mass under the [M+/-H] convention (matches ``utils.phase0``)."""
    mz = float(parent_mz)
    if int(lib_id) in {400, 402}:
        return mz - PROTON_MASS
    return mz + PROTON_MASS


def _adduct_hypotheses(parent_mz: float, lib_id: int) -> list[tuple[str, float, float]]:
    """Enumerate ``(adduct_label, neutral_mass, recomputed_mz)`` hypotheses.

    The first entry is always the primary [M+/-H] hypothesis; the remainder are a
    small common-adduct ladder used only for mass matching in ``candidates``.
    """
    mz = float(parent_mz)
    if int(lib_id) in {400, 402}:  # positive
        rows = [
            ("[M+H]+", mz - PROTON_MASS),
            ("[M+NH4]+", mz - 18.033823),
            ("[M+Na]+", mz - 22.989218),
            ("[M+K]+", mz - 38.963158),
            ("[2M+H]+", (mz - PROTON_MASS) / 2.0),
        ]
    else:  # negative
        rows = [
            ("[M-H]-", mz + PROTON_MASS),
            ("[M+Cl]-", mz - 34.969402),
            ("[M+HCOO]-", mz - 44.998201),
            ("[M+CH3COO]-", mz - 59.013851),
            ("[2M-H]-", (mz + PROTON_MASS) / 2.0),
        ]
    out: list[tuple[str, float, float]] = []
    for label, neutral in rows:
        if neutral is None or neutral <= 0:
            continue
        out.append((label, float(neutral), mz))
    return out


@dataclass
class NovelFeature:
    """One un-curated untargeted feature to elucidate."""

    feature_id: int
    lib_id: int
    ion_mode: str
    parent_mz: float
    ri: float
    peaks: list[tuple[float, float]]
    compound_name: str = ""
    prior_match_level: str = ""
    prior_matched_name: str = ""
    # Adduct-corrected neutral mass from the upstream pipeline's peak grouping
    # (e.g. OpenMS ``dc_charge_adduct_mass``). When present this is preferred over
    # the naive [M+/-H] assumption so non-protonated base peaks are anchored right.
    neutral_mass: Optional[float] = None
    adduct_type: str = ""

    @property
    def control_id(self) -> str:
        """Synthetic X-compound id for this feature."""
        return control_id_for_feature(self.feature_id)

    @property
    def effective_neutral_mass(self) -> float:
        """Best neutral-mass anchor: upstream adduct-corrected value or [M+/-H]."""
        nm = self.neutral_mass
        if nm is not None:
            try:
                value = float(nm)
            except (TypeError, ValueError):
                value = float("nan")
            if math.isfinite(value) and value > 0:
                return value
        return neutral_mass_from_precursor(self.parent_mz, self.lib_id)

    @property
    def primary_adduct_label(self) -> str:
        """Adduct label to record for the primary neutral-mass hypothesis."""
        if self.adduct_type:
            return str(self.adduct_type)
        return "[M+H]+" if int(self.lib_id) in {400, 402} else "[M-H]-"

    @property
    def comp_id(self) -> int:
        """Synthetic integer comp id (the untargeted feature id)."""
        return int(self.feature_id)

    @property
    def feature_key(self) -> str:
        """Unique feature key embedding lib, m/z, RI, and control id."""
        return (
            f"feature_X_{int(self.lib_id)}_mz{float(self.parent_mz):.4f}"
            f"_ri{int(round(float(self.ri)))}_{self.control_id}"
        )


@dataclass
class ScopeInputs:
    """Resolved paths for one ``lib{lib}_{ion}`` scope."""

    lib_id: int
    ion_mode: str
    n_features: int
    feature_rows: Path
    neutral_mass_hypotheses: Path
    gap_manifest: Path
    observed_ms2_lookup: Path
    encoder_table: Path


@dataclass
class AdapterManifest:
    """Everything the driver needs to run the stack on the synthesized inputs."""

    out_dir: Path
    scopes: list[ScopeInputs] = field(default_factory=list)
    all_x_catalog: Path | None = None
    integrated_best: Path | None = None
    novel_feature_map: Path | None = None

    def to_json(self) -> dict[str, Any]:
        """Serialize to a plain JSON-able dict."""
        return {
            "out_dir": str(self.out_dir),
            "all_x_catalog": str(self.all_x_catalog) if self.all_x_catalog else None,
            "integrated_best": str(self.integrated_best) if self.integrated_best else None,
            "novel_feature_map": str(self.novel_feature_map) if self.novel_feature_map else None,
            "scopes": [
                {
                    "lib_id": s.lib_id,
                    "ion_mode": s.ion_mode,
                    "n_features": s.n_features,
                    "feature_rows": str(s.feature_rows),
                    "neutral_mass_hypotheses": str(s.neutral_mass_hypotheses),
                    "gap_manifest": str(s.gap_manifest),
                    "observed_ms2_lookup": str(s.observed_ms2_lookup),
                    "encoder_table": str(s.encoder_table),
                }
                for s in self.scopes
            ],
        }


def _coerce_peaks(raw: Any) -> list[tuple[float, float]]:
    """Coerce a handoff ``peaks`` field into ``[(mz, intensity), ...]``."""
    peaks: list[tuple[float, float]] = []
    for item in raw or []:
        try:
            mz = float(item[0])
            inten = float(item[1]) if len(item) > 1 else 1.0
        except (TypeError, ValueError, IndexError):
            continue
        peaks.append((mz, inten))
    return peaks


def load_handoff(path: Path) -> list[NovelFeature]:
    """Parse a handoff JSON file into :class:`NovelFeature` objects."""
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    records = payload.get("features", payload) if isinstance(payload, dict) else payload
    features: list[NovelFeature] = []
    for record in records or []:
        lib_id = int(record["lib_id"])
        ion_mode = str(record.get("ion_mode") or ion_mode_for_lib(lib_id)).strip().lower()
        raw_neutral = record.get("neutral_mass", None)
        try:
            neutral_mass = float(raw_neutral) if raw_neutral is not None else None
        except (TypeError, ValueError):
            neutral_mass = None
        features.append(
            NovelFeature(
                feature_id=int(record["feature_id"]),
                lib_id=lib_id,
                ion_mode=ion_mode,
                parent_mz=float(record["parent_mz"]),
                ri=float(record.get("ri", 0.0) or 0.0),
                peaks=_coerce_peaks(record.get("peaks")),
                compound_name=str(record.get("compound_name", "") or ""),
                prior_match_level=str(record.get("prior_match_level", "") or ""),
                prior_matched_name=str(record.get("prior_matched_name", "") or ""),
                neutral_mass=neutral_mass,
                adduct_type=str(record.get("adduct_type", "") or ""),
            )
        )
    return features


def _feature_rows_frame(features: Iterable[NovelFeature]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for feat in features:
        ion_bins, loss_bins = observed_bins_for_spectrum(feat.peaks, feat.parent_mz, feat.lib_id)
        rows.append(
            {
                "feature_key": feat.feature_key,
                "comp_id": feat.comp_id,
                "control_id": feat.control_id,
                "lib_id": int(feat.lib_id),
                "ion_mode": feat.ion_mode,
                "parent_mz": float(feat.parent_mz),
                "neutral_mass_fallback": feat.effective_neutral_mass,
                "n_ms2_peaks": len(feat.peaks),
                "observed_ion_bins_json": json.dumps(ion_bins),
                "observed_loss_bins_json": json.dumps(loss_bins),
            }
        )
    return pd.DataFrame(rows)


def _neutral_mass_hypotheses_frame(features: Iterable[NovelFeature]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for feat in features:
        # Primary hypothesis uses the upstream adduct-corrected neutral mass so
        # candidate mass matching is anchored correctly even for non-[M+/-H] base
        # peaks; the standard ladder is appended as additional fallbacks.
        eff = feat.effective_neutral_mass
        hypotheses: list[tuple[str, float, float]] = [
            (feat.primary_adduct_label, eff, float(feat.parent_mz))
        ]
        for label, neutral, recomputed in _adduct_hypotheses(feat.parent_mz, feat.lib_id):
            if abs(float(neutral) - float(eff)) <= 1e-6:
                continue
            hypotheses.append((label, neutral, recomputed))
        for rank, (label, neutral, recomputed) in enumerate(hypotheses, start=1):
            rows.append(
                {
                    "feature_key": feat.feature_key,
                    "comp_id": feat.comp_id,
                    "lib_id": int(feat.lib_id),
                    "ion_mode": feat.ion_mode,
                    "parent_mz": float(feat.parent_mz),
                    "adduct_label": label,
                    "charge": 1,
                    "multiplier": 1,
                    "offset_mass": float(neutral) - float(feat.parent_mz),
                    "neutral_mass": float(neutral),
                    "recomputed_mz": float(recomputed),
                    "adduct_score": 1.0 if rank == 1 else 0.5,
                    "prevalence": 1.0 if rank == 1 else 0.25,
                    "prevalence_tier": "primary" if rank == 1 else "ladder",
                    "source_reference": "untargeted_adapter_adduct_ladder",
                    "neutral_loss_possible": False,
                    "note": "primary [M+/-H]" if rank == 1 else "common adduct ladder",
                    "adduct_rank": rank,
                    "control_id": feat.control_id,
                    "formula_prior_score": 0.0,
                    "ms1_buddy_supported": False,
                }
            )
    return pd.DataFrame(rows)


def _gap_manifest_frame(features: Iterable[NovelFeature]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for feat in features:
        rows.append(
            {
                "control_id": feat.control_id,
                "library_id": int(feat.lib_id),
                "lib_id": int(feat.lib_id),
                "ion_mode": feat.ion_mode,
                "level": 4,
                "sublevel": "4",
                "comp_id": feat.comp_id,
                "feature_key": feat.feature_key,
                "parent_mz": float(feat.parent_mz),
                "neutral_mass_estimate": feat.effective_neutral_mass,
                "has_feature_rows": True,
                "has_mass_hypotheses": True,
                "already_in_unified_pool": False,
                "primary_stage": "not_processed",
                "recovery_track": "untargeted_novel_v1",
                "clean_external_info": feat.prior_matched_name or feat.compound_name or "",
            }
        )
    return pd.DataFrame(rows)


def _observed_lookup_frame(features: Iterable[NovelFeature]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for feat in features:
        ion_bins, loss_bins = observed_bins_for_spectrum(feat.peaks, feat.parent_mz, feat.lib_id)
        rows.append(
            {
                "control_id": feat.control_id,
                "comp_id": feat.comp_id,
                "lib_id": int(feat.lib_id),
                "ion_mode": feat.ion_mode,
                "parent_mz": float(feat.parent_mz),
                "neutral_mass": feat.effective_neutral_mass,
                "ms2_present": True,
                "n_ms2_pairs": len(feat.peaks),
                "n_observed_ion_bins": len(ion_bins),
                "n_observed_loss_bins": len(loss_bins),
                "observed_ion_bins_json": json.dumps(ion_bins),
                "observed_loss_bins_json": json.dumps(loss_bins),
                "ms2_pairs_json": json.dumps([[float(mz), float(inten)] for mz, inten in feat.peaks]),
                "spectrum_source": "untargeted_consensus_ms2",
            }
        )
    return pd.DataFrame(rows)


def _encoder_table_frame(features: Iterable[NovelFeature]) -> pd.DataFrame:
    """Long-format encoder input: one row per MS2 peak (optional de novo path)."""
    rows: list[dict[str, Any]] = []
    for feat in features:
        name = feat.compound_name or feat.control_id
        for mz, inten in feat.peaks:
            rows.append(
                {
                    "feature_key": feat.feature_key,
                    "comp_id": feat.comp_id,
                    "compound_name": name,
                    "rt": float(feat.ri),
                    "parent_mz": float(feat.parent_mz),
                    "mass": float(mz),
                    "intensity": float(inten),
                    "source_html": "",
                    "lib_id": int(feat.lib_id),
                    "ion_mode": feat.ion_mode,
                }
            )
    return pd.DataFrame(
        rows,
        columns=[
            "feature_key",
            "comp_id",
            "compound_name",
            "rt",
            "parent_mz",
            "mass",
            "intensity",
            "source_html",
            "lib_id",
            "ion_mode",
        ],
    )


def _all_x_catalog_frame(features: list[NovelFeature]) -> pd.DataFrame:
    """Catalog that (a) lists novel ids as non-S0 and (b) seeds one dummy S0 row.

    ``pipeline.pool`` only writes ``unified_candidates_per_control_non_s0.csv``
    when the S0 control set for the library is non-empty; otherwise both split
    files come out empty and the rerank stage never surfaces un-truthed novel
    compounds.  A single dummy S0 row per library (with an id that never appears
    in the pool) triggers the split so every real novel control lands in the
    non-S0 file.
    """
    rows: list[dict[str, Any]] = []
    libs = sorted({int(f.lib_id) for f in features})
    for lib in libs:
        rows.append(
            {
                "x_compound_id": f"X-{lib}000000001",
                "control_id": f"X-{lib}000000001",
                "lib_id": lib,
                "evidence_tier": "S0",
            }
        )
    for feat in features:
        rows.append(
            {
                "x_compound_id": feat.control_id,
                "control_id": feat.control_id,
                "lib_id": int(feat.lib_id),
                "evidence_tier": "S2",
            }
        )
    return pd.DataFrame(rows, columns=["x_compound_id", "control_id", "lib_id", "evidence_tier"])


def _integrated_best_frame(features: list[NovelFeature]) -> pd.DataFrame:
    """Stub integrated-best so pool/verify keep novel controls in the allow-list."""
    rows = [
        {
            "control_id": feat.control_id,
            "lib_id": int(feat.lib_id),
            "ion_mode": feat.ion_mode,
            "top10_smiles": "",
            "top10_best_tanimoto": 0.0,
        }
        for feat in features
    ]
    return pd.DataFrame(
        rows, columns=["control_id", "lib_id", "ion_mode", "top10_smiles", "top10_best_tanimoto"]
    )


def _novel_feature_map_frame(features: list[NovelFeature]) -> pd.DataFrame:
    rows = [
        {
            "feature_id": feat.feature_id,
            "control_id": feat.control_id,
            "comp_id": feat.comp_id,
            "lib_id": int(feat.lib_id),
            "ion_mode": feat.ion_mode,
            "parent_mz": float(feat.parent_mz),
            "ri": float(feat.ri),
            "feature_key": feat.feature_key,
            "n_ms2_peaks": len(feat.peaks),
            "prior_match_level": feat.prior_match_level,
            "prior_matched_name": feat.prior_matched_name,
        }
        for feat in features
    ]
    return pd.DataFrame(rows)


def build_inputs(
    features: list[NovelFeature],
    out_dir: Path,
    *,
    write_encoder_table: bool = True,
) -> AdapterManifest:
    """Write all elucidation input tables for ``features`` under ``out_dir``.

    Features are grouped by ``(lib_id, ion_mode)``; unsupported libraries (not in
    :data:`SUPPORTED_LIBS`) are skipped.  Returns an :class:`AdapterManifest`
    describing every file written.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    supported = [f for f in features if int(f.lib_id) in SUPPORTED_LIBS]
    manifest = AdapterManifest(out_dir=out_dir)

    by_scope: dict[tuple[int, str], list[NovelFeature]] = {}
    for feat in supported:
        by_scope.setdefault((int(feat.lib_id), feat.ion_mode), []).append(feat)

    for (lib_id, ion_mode), scope_feats in sorted(by_scope.items()):
        scope = f"lib{lib_id}_{ion_mode}"
        feature_rows_path = out_dir / f"{scope}_feature_rows.csv"
        mass_hyp_path = out_dir / f"{scope}_neutral_mass_hypotheses.csv"
        gap_path = out_dir / f"{scope}_gap_manifest.csv"
        observed_path = out_dir / f"{scope}_observed_ms2_lookup.csv"
        encoder_path = out_dir / f"ions_novel_{scope}.csv"

        _feature_rows_frame(scope_feats).to_csv(feature_rows_path, index=False)
        _neutral_mass_hypotheses_frame(scope_feats).to_csv(mass_hyp_path, index=False)
        _gap_manifest_frame(scope_feats).to_csv(gap_path, index=False)
        _observed_lookup_frame(scope_feats).to_csv(observed_path, index=False)
        if write_encoder_table:
            _encoder_table_frame(scope_feats).to_csv(encoder_path, index=False)

        manifest.scopes.append(
            ScopeInputs(
                lib_id=lib_id,
                ion_mode=ion_mode,
                n_features=len(scope_feats),
                feature_rows=feature_rows_path,
                neutral_mass_hypotheses=mass_hyp_path,
                gap_manifest=gap_path,
                observed_ms2_lookup=observed_path,
                encoder_table=encoder_path,
            )
        )

    catalog_path = out_dir / "all_x_catalog.csv"
    integrated_path = out_dir / "integrated_per_control_best.csv"
    map_path = out_dir / "novel_feature_map.csv"
    _all_x_catalog_frame(supported).to_csv(catalog_path, index=False)
    _integrated_best_frame(supported).to_csv(integrated_path, index=False)
    _novel_feature_map_frame(supported).to_csv(map_path, index=False)
    manifest.all_x_catalog = catalog_path
    manifest.integrated_best = integrated_path
    manifest.novel_feature_map = map_path

    (out_dir / "adapter_manifest.json").write_text(
        json.dumps(manifest.to_json(), indent=2), encoding="utf-8"
    )
    return manifest


def build_inputs_from_handoff(
    handoff_path: Path, out_dir: Path, *, write_encoder_table: bool = True
) -> AdapterManifest:
    """Load a handoff JSON and build all elucidation input tables."""
    features = load_handoff(Path(handoff_path))
    return build_inputs(features, Path(out_dir), write_encoder_table=write_encoder_table)


def main() -> None:
    """CLI: ``python -m structural_elucidation.adapters.untargeted_inputs``."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--handoff", type=Path, required=True, help="Handoff JSON path.")
    parser.add_argument("--out-dir", type=Path, required=True, help="Directory for synthesized inputs.")
    parser.add_argument(
        "--no-encoder-table",
        action="store_true",
        help="Skip writing the long-format encoder/de-novo input table.",
    )
    args = parser.parse_args()
    manifest = build_inputs_from_handoff(
        args.handoff, args.out_dir, write_encoder_table=not args.no_encoder_table
    )
    print(json.dumps(manifest.to_json(), indent=2))


if __name__ == "__main__":
    main()


# Re-export for convenience so callers can render model keys without importing
# the low-level module directly.
__all__ = [
    "NovelFeature",
    "AdapterManifest",
    "ScopeInputs",
    "NOVEL_CONTROL_OFFSET",
    "control_id_for_feature",
    "feature_id_from_control_id",
    "load_handoff",
    "build_inputs",
    "build_inputs_from_handoff",
    "nce_bin_name",
]
