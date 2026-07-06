"""Smoke tests for the hybrid identification-confidence model + report patchers."""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "structural_elucidation" / "reports"))
sys.path.insert(0, str(ROOT / "tools"))

import confidence as C  # noqa: E402


def _fields(**over):
    base = dict(
        source="predicted", structure_similarity=0.6, mass_match_score=0.95,
        abs_mass_error_da=0.001, ms2_spectrum_present=True, ms2_score=0.4,
        pathway_vote_share=0.9, pathway_tier="high", observed_mz_available=True,
        candidate_margin=0.2, library_neighbors_ge_050=3, correlations_real=True,
    )
    base.update(over)
    return base


def test_overall_in_unit_interval_and_axes():
    conf = C.compute_confidence(_fields())
    assert 0.0 <= conf.overall <= 1.0
    assert conf.tier in {"Confirmed", "High", "Moderate", "Low", "Tentative"}
    assert {a.key for a in conf.axes} == {"structure", "mass", "ms2", "pathway"}
    assert 0.0 <= conf.completeness <= 1.0


def test_mass_mismatch_forces_tentative():
    conf = C.compute_confidence(_fields(mass_match_score=0.0, abs_mass_error_da=176.0))
    assert conf.tier == "Tentative"
    assert any("mass" in f.lower() for f in conf.flags)


def test_level1_truth_match_is_confirmed():
    conf = C.compute_confidence(_fields(metabolon_level="1", truth_match=True))
    assert conf.tier == "Confirmed"
    assert conf.overall >= 0.95


def test_level1_without_truth_match_not_confirmed():
    conf = C.compute_confidence(_fields(metabolon_level="1", truth_match=False))
    assert conf.tier != "Confirmed"


def test_exact_is_not_pinned_high():
    # exact catalog match is only ~38% correct -> must not be High/Confirmed alone
    conf = C.compute_confidence(_fields(source="exact", structure_similarity=1.0))
    assert conf.tier in {"Moderate", "Low", "Tentative"}


def test_ms2_absent_marks_axis_unavailable():
    conf = C.compute_confidence(_fields(ms2_spectrum_present=False, ms2_score=None))
    ms2 = next(a for a in conf.axes if a.key == "ms2")
    assert ms2.available is False
    assert conf.completeness < 1.0


def test_public_scored_params_omit_ms2():
    conf = C.compute_confidence(_fields())
    pub = C.scored_parameters_section_html(conf, public=True)
    internal = C.scored_parameters_section_html(conf, public=False)
    assert "MS2 fragment fit" not in pub
    assert "MS2 fragment fit" in internal


def test_public_ms2_section_has_no_ion_detail():
    fields = _fields()
    fields["ms2_detail"] = {"present": True, "n_pairs": 12, "n_ion_bins": 10,
                            "n_loss_bins": 3, "forward": 0.5, "predicted_v2": 0.2, "match": 0.1}
    pub = C.ms2_evidence_section_html(fields, public=True)
    internal = C.ms2_evidence_section_html(fields, public=False)
    for marker in ("fragment peaks", "fragment-ion bins", "neutral-loss bins"):
        assert marker not in pub
    assert "Observed fragment peaks" in internal


def test_confidence_badge_and_card_render():
    conf = C.compute_confidence(_fields(metabolon_level="1", truth_match=True))
    assert "Confidence" in C.confidence_badge_html(conf)
    assert "Overall identification confidence" in C.confidence_metric_card_html(conf)


def test_correlations_table_html():
    rows = C.correlations_table_html([
        {"id": "100", "name": "foo", "score": 0.8, "source": "wgcna_HELIOS_module_1"},
    ])
    assert "foo" in rows and "wgcna_HELIOS_module_1" in rows


def test_consensus_decrowder_reduces_overlap():
    import patch_consensus_crowding as PC
    # two atoms 6px apart, radius 16 -> overlap_ratio ~5.3
    block = (
        '<svg><g class="consensus-atom" transform="translate(0.0,0.0)">'
        '<circle r="16.0" stroke-width="2.5"><title>a</title></circle>'
        '<text y="4.5" font-size="13.0">C</text></g>'
        '<g class="consensus-atom" transform="translate(6.0,0.0)">'
        '<circle r="16.0" stroke-width="2.5"><title>b</title></circle>'
        '<text y="4.5" font-size="13.0">O</text></g></svg>'
    )
    ratio_before, _r, _d = PC._overlap_ratio(block)
    assert ratio_before > 1.0
    new_block, changed = PC.patch_block(block, target=0.92)
    assert changed
    ratio_after, _r2, _d2 = PC._overlap_ratio(new_block)
    assert ratio_after <= 1.0


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
