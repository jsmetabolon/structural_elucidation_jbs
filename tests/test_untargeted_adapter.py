"""Tests for the untargeted -> structural-elucidation input adapter.

Covers:
* bit-exact NCE bin grid vs. the vendored ``elucidation.binning.bin_for_mass``
  (numpy-only, importable without rdkit/torch),
* the signed ion/loss bin convention and its round-trip through ``nce_bin_name``,
* the synthesized per-library input tables, the non-S0 catalog seeding, and the
  ``feature_id`` <-> ``control_id`` mapping used to fuse results back.

Heavier cross-checks against ``ModelingAssistant`` are attempted but skipped when
optional deps (rdkit) are unavailable.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from structural_elucidation.adapters import nce_bins
from structural_elucidation.adapters.untargeted_inputs import (
    NOVEL_CONTROL_OFFSET,
    NovelFeature,
    build_inputs,
    control_id_for_feature,
    feature_id_from_control_id,
    load_handoff,
)


# ── bin grid bit-exactness ────────────────────────────────────────────────────

def _vendored_bin_for_mass():
    """Import the vendored numpy-only ``bin_for_mass`` or skip."""
    vendor = REPO_ROOT / "vendor" / "nce"
    if str(vendor) not in sys.path:
        sys.path.insert(0, str(vendor))
    try:
        from elucidation.binning import bin_for_mass  # type: ignore
    except Exception as exc:  # pragma: no cover - optional dependency guard
        pytest.skip(f"vendored elucidation.binning unavailable: {exc}")
    return bin_for_mass


@pytest.mark.parametrize("mass", [50.0, 67.05, 100.0, 137.0492, 249.0857, 300.5, 999.9])
def test_bin_for_mass_matches_vendored(mass: float) -> None:
    """Our grid must equal the vendored NCE grid exactly."""
    vendored = _vendored_bin_for_mass()
    assert nce_bins.bin_for_mass(mass) == int(vendored(mass))


def test_ion_and_loss_bin_signs() -> None:
    """Ion bins (m >= 50) are positive; loss bins are always negative."""
    ion = nce_bins.ion_bin_int(137.05, 400)
    loss = nce_bins.loss_bin_int(113.04)
    assert ion > 0
    assert loss < 0
    assert nce_bins.nce_bin_name(ion) == f"bin_{ion}"
    assert nce_bins.nce_bin_name(loss) == f"bin_m{-loss}"


def test_observed_bins_dedup_and_sorted() -> None:
    """Observed bins are unique and sorted; losses only for fragments < parent."""
    parent = 250.09295
    peaks = [(67.0, 0.5), (137.05, 1.0), (200.1, 0.3), (250.05, 0.9)]  # last ~ parent -> no loss
    ion_bins, loss_bins = nce_bins.observed_bins_for_spectrum(peaks, parent, 400)
    assert ion_bins == sorted(set(ion_bins))
    assert loss_bins == sorted(set(loss_bins))
    assert all(b < 0 for b in loss_bins)
    # The near-precursor peak contributes an ion bin but no loss bin.
    assert len(loss_bins) == 3


def test_ion_loss_match_modeling_assistant_when_available() -> None:
    """Optional: signed ints -> names must equal ModelingAssistant bin names."""
    try:
        from structural_elucidation.nce.runtime import prepare_nce_runtime

        prepare_nce_runtime()
        from elucidation.modelingAssistant import ModelingAssistant  # type: ignore
    except Exception as exc:  # pragma: no cover - optional dependency guard
        pytest.skip(f"ModelingAssistant unavailable: {exc}")

    ma = ModelingAssistant.get_instance(400)
    parent = 250.09295
    for mz in (67.0, 137.05, 200.1):
        assert nce_bins.nce_bin_name(nce_bins.ion_bin_int(mz, 400)) == ma.ion_bin(mz)
        if mz < parent - 0.5:
            assert nce_bins.nce_bin_name(nce_bins.loss_bin_int(parent - mz)) == ma.loss_bin(parent - mz)


# ── control id mapping ────────────────────────────────────────────────────────

def test_control_id_roundtrip() -> None:
    """``feature_id`` -> ``control_id`` -> ``feature_id`` is stable and offset-safe."""
    cid = control_id_for_feature(12345)
    assert cid == f"X-{NOVEL_CONTROL_OFFSET + 12345}"
    assert feature_id_from_control_id(cid) == 12345
    # A real Metabolon-style X id is not mistaken for a novel one.
    assert feature_id_from_control_id("X-25003") is None


# ── input synthesis ───────────────────────────────────────────────────────────

def _sample_features() -> list[NovelFeature]:
    return [
        NovelFeature(
            feature_id=101,
            lib_id=400,
            ion_mode="positive",
            parent_mz=250.09295,
            ri=2486.0,
            peaks=[(67.0, 0.5), (137.05, 1.0), (200.1, 0.3)],
            prior_match_level="MS2",
            prior_matched_name="lib_hit",
        ),
        NovelFeature(
            feature_id=102,
            lib_id=209,
            ion_mode="negative",
            parent_mz=300.12,
            ri=1500.0,
            peaks=[(80.0, 0.9), (150.02, 1.0)],
        ),
        # Unsupported library must be dropped.
        NovelFeature(feature_id=103, lib_id=208, ion_mode="positive", parent_mz=120.0, ri=10.0, peaks=[(50.0, 1.0)]),
    ]


def test_build_inputs_writes_expected_scopes(tmp_path: Path) -> None:
    """All per-scope tables + global stubs are written; lib 208 is skipped."""
    import pandas as pd

    manifest = build_inputs(_sample_features(), tmp_path)
    scope_libs = sorted(s.lib_id for s in manifest.scopes)
    assert scope_libs == [209, 400]  # 208 dropped

    for scope in manifest.scopes:
        assert scope.feature_rows.exists()
        assert scope.neutral_mass_hypotheses.exists()
        assert scope.gap_manifest.exists()
        assert scope.observed_ms2_lookup.exists()

    # feature_rows carry the required candidate columns.
    fr = pd.read_csv(next(s.feature_rows for s in manifest.scopes if s.lib_id == 400))
    for col in ("feature_key", "control_id", "parent_mz", "neutral_mass_fallback",
                "observed_ion_bins_json", "observed_loss_bins_json"):
        assert col in fr.columns
    assert fr.iloc[0]["control_id"] == control_id_for_feature(101)

    # observed lookup loss bins are stored as negative (signed) ints.
    lookup = pd.read_csv(next(s.observed_ms2_lookup for s in manifest.scopes if s.lib_id == 400))
    loss_bins = json.loads(lookup.iloc[0]["observed_loss_bins_json"])
    assert loss_bins and all(b < 0 for b in loss_bins)


def test_all_x_catalog_seeds_non_s0_split(tmp_path: Path) -> None:
    """Catalog must contain one dummy S0 row per lib plus novel non-S0 rows."""
    import pandas as pd

    manifest = build_inputs(_sample_features(), tmp_path)
    catalog = pd.read_csv(manifest.all_x_catalog)
    # Every supported lib gets exactly one S0 seed row.
    s0 = catalog[catalog["evidence_tier"].astype(str).str.contains("S0")]
    assert set(s0["lib_id"]) == {209, 400}
    # Novel controls are non-S0 so they land in the non-S0 pool.
    novel = catalog[catalog["control_id"] == control_id_for_feature(101)]
    assert len(novel) == 1
    assert "S0" not in str(novel.iloc[0]["evidence_tier"])


def test_effective_neutral_mass_prefers_provided() -> None:
    """A provided adduct-corrected neutral mass overrides the [M+/-H] fallback."""
    feat = NovelFeature(
        feature_id=9, lib_id=400, ion_mode="positive", parent_mz=200.0, ri=1.0,
        peaks=[(80.0, 1.0)], neutral_mass=176.9, adduct_type="[M+Na]+",
    )
    assert feat.effective_neutral_mass == pytest.approx(176.9)
    assert feat.primary_adduct_label == "[M+Na]+"
    # Falls back to [M+H]+ neutral mass when none provided.
    feat2 = NovelFeature(
        feature_id=9, lib_id=400, ion_mode="positive", parent_mz=200.0, ri=1.0, peaks=[(80.0, 1.0)]
    )
    assert feat2.effective_neutral_mass == pytest.approx(200.0 - 1.007276466812)
    assert feat2.primary_adduct_label == "[M+H]+"


def test_provided_neutral_mass_flows_into_tables(tmp_path: Path) -> None:
    """The adduct-corrected neutral mass anchors feature_rows/hypotheses/lookup."""
    import pandas as pd

    feat = NovelFeature(
        feature_id=42, lib_id=400, ion_mode="positive", parent_mz=250.0, ri=10.0,
        peaks=[(80.0, 1.0), (150.0, 0.5)], neutral_mass=226.9, adduct_type="[M+Na]+",
    )
    manifest = build_inputs([feat], tmp_path)
    scope = manifest.scopes[0]

    fr = pd.read_csv(scope.feature_rows)
    assert fr.iloc[0]["neutral_mass_fallback"] == pytest.approx(226.9)

    lookup = pd.read_csv(scope.observed_ms2_lookup)
    assert lookup.iloc[0]["neutral_mass"] == pytest.approx(226.9)

    hyp = pd.read_csv(scope.neutral_mass_hypotheses)
    primary = hyp[hyp["adduct_rank"] == 1].iloc[0]
    assert primary["neutral_mass"] == pytest.approx(226.9)
    assert primary["adduct_label"] == "[M+Na]+"
    # The naive [M+H]+ mass should NOT be the anchor here.
    assert abs(float(primary["neutral_mass"]) - (250.0 - 1.007276466812)) > 1.0


def test_load_handoff_roundtrip(tmp_path: Path) -> None:
    """``load_handoff`` parses the documented JSON schema."""
    handoff = {
        "features": [
            {"feature_id": 7, "lib_id": 402, "ion_mode": "positive",
             "parent_mz": 180.1, "ri": 900.0, "peaks": [[60.0, 1.0], [120.0, 0.5]]},
        ]
    }
    path = tmp_path / "handoff.json"
    path.write_text(json.dumps(handoff), encoding="utf-8")
    features = load_handoff(path)
    assert len(features) == 1
    assert features[0].control_id == control_id_for_feature(7)
    assert features[0].peaks == [(60.0, 1.0), (120.0, 0.5)]
