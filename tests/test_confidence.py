"""Smoke tests for the hybrid identification-confidence model + report patchers."""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "structural_elucidation" / "reports"))
sys.path.insert(0, str(ROOT / "tools"))

import confidence as C  # noqa: E402
import patch_report_confidence as PC  # noqa: E402


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
    # Authentic-standard-confirmed reads 100% by policy (see CONFIRMED_OVERALL).
    assert conf.overall == 1.0
    assert conf.overall == C.CONFIRMED_OVERALL


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


def test_public_scored_params_show_aggregate_ms2():
    conf = C.compute_confidence(_fields())
    pub = C.scored_parameters_section_html(conf, public=True)
    internal = C.scored_parameters_section_html(conf, public=False)
    assert "MS2 fragment fit" in pub
    assert "MS2 fragment fit" in internal
    assert "Observed fragment peaks" not in pub


def test_public_ms2_section_has_no_ion_detail():
    fields = _fields()
    fields["ms2_detail"] = {"present": True, "n_pairs": 12, "n_ion_bins": 10,
                            "n_loss_bins": 3, "forward": 0.5, "predicted_v2": 0.2, "match": 0.1}
    pub = C.ms2_evidence_section_html(fields, public=True)
    internal = C.ms2_evidence_section_html(fields, public=False)
    for marker in ("fragment peaks", "fragment-ion bins", "neutral-loss bins"):
        assert marker not in pub
    assert "Observed fragment peaks" in internal


def test_public_fragment_summary_shows_counts_not_raw_metrics():
    fields = _fields()
    fields["ms2_detail"] = {
        "present": True,
        "n_pairs": 45,
        "n_ion_bins": 45,
        "n_loss_bins": 45,
        "forward": 0.02,
        "predicted_v2": 0.0,
        "match": 0.0,
    }
    candidates = [
        {
            "rank": 1,
            "predicted_formula": "C12H19N3O5",
            "predicted_smiles": "CC(=O)N1CCCC1C(=O)NC(CCC(N)=O)C(=O)O",
            "forward_spectrum_score": 0.019,
            "predicted_spectrum_score_v2": 0.0,
            "ms2_match_score": 0.0,
        },
        {
            "rank": 2,
            "predicted_formula": "C11H18N2O5",
            "predicted_smiles": "CC(C)C(NC(=O)C1CCC(=O)N1)C(=O)NCC(=O)O",
            "forward_spectrum_score": 0.5,
            "predicted_spectrum_score_v2": 0.0,
            "ms2_match_score": 0.0,
        },
    ]
    pub = C.fragment_summary_section_html(fields, candidates, public=True)
    internal = C.fragment_summary_section_html(fields, candidates, public=False)
    assert "Observed fragment peaks" in pub
    assert "Forward fit" in pub
    assert "0.500" in pub
    assert "Forward spectral score" not in pub
    assert "Forward spectral score" in internal
    assert "<!--FRAGMENT_SUMMARY-->" in pub


def test_diagnostic_substructure_rows_html_lists_present_diagnostics():
    rank1 = {
        "lactoyl_hypothesis_score": 0.82,
        "adduct_label": "[M-H]-",
        "adduct_score": 0.71,
        "transform_name": "glucuronidation",
        "biotransform_v2_score": 0.66,
        "recursive__recursive_substructure_weighted_score": 0.44,
    }
    rows = C.diagnostic_substructure_rows_html(rank1, {})
    assert "Lactoyl" in rows
    assert "Adduct ([M-H]-)" in rows
    assert "Biotransformation (glucuronidation)" in rows
    assert "Recursive substructure" in rows
    assert "<!--FRAGMENT_DIAGNOSTICS-->" in rows


def test_diagnostic_substructure_rows_html_fallback_when_empty():
    rows = C.diagnostic_substructure_rows_html({}, {})
    assert "No diagnostic substructure hypotheses were scored" in rows


def test_confidence_badge_and_card_render():
    conf = C.compute_confidence(_fields(metabolon_level="1", truth_match=True))
    assert "Confidence" in C.confidence_badge_html(conf)
    assert "Overall identification confidence" in C.confidence_metric_card_html(conf)


def test_correlations_table_html():
    rows = C.correlations_table_html([
        {"id": "100", "name": "foo", "score": 0.8, "source": "wgcna_HELIOS_module_1"},
    ])
    assert "foo" in rows and "wgcna_HELIOS_module_1" in rows


def test_msi_level1_for_confirmed():
    conf = C.compute_confidence(_fields(metabolon_level="1", truth_match=True))
    assert conf.msi_level == 1
    assert conf.tier == "Confirmed"


def test_msi_level2_probable_structure():
    conf = C.compute_confidence(_fields(
        structure_similarity=0.8, ms2_spectrum_present=True, ms2_score=0.5,
        mass_match_score=0.98, metabolon_level="", truth_match=False,
    ))
    assert conf.msi_level == 2


def test_msi_level3_tentative_when_no_ms2():
    conf = C.compute_confidence(_fields(
        structure_similarity=0.6, ms2_spectrum_present=False, ms2_score=None,
        mass_match_score=0.95,
    ))
    assert conf.msi_level == 3


def test_msi_level4_for_mass_mismatch():
    # observed mass known but candidate mass does not fit -> formula only (L4)
    conf = C.compute_confidence(_fields(
        structure_similarity=0.6, mass_match_score=0.0, abs_mass_error_da=176.0,
        observed_mz_available=True,
    ))
    assert conf.msi_level == 4
    assert "reject" in conf.msi_rationale.lower()


def test_msi_level5_without_observed_mass():
    conf = C.compute_confidence(_fields(
        structure_similarity=0.2, mass_match_score=None, abs_mass_error_da=None,
        observed_mz_available=False,
    ))
    assert conf.msi_level == 5


def test_msi_level_always_in_range():
    for over in ({}, {"source": "exact"}, {"structure_similarity": 0.0},
                 {"mass_match_score": None, "observed_mz_available": False}):
        conf = C.compute_confidence(_fields(**over))
        assert conf.msi_level in {1, 2, 3, 4, 5}
        assert conf.msi_label == C.MSI_LABELS[conf.msi_level]


def test_msi_badge_and_card_render():
    conf = C.compute_confidence(_fields())
    assert "MSI Level" in C.msi_badge_html(conf)
    assert "MSI identification level" in C.msi_metric_card_html(conf)
    section = C.scored_parameters_section_html(conf, public=True)
    assert "MSI identification level" in section


def _baseline_conf(**over):
    """Baseline conf with the exact + strong-mass profile used by the identity checks.

    Later kwargs win, so callers can override any baseline field.
    """
    baseline = dict(
        source="exact", structure_similarity=1.0, name_resolved=True,
        mass_match_score=0.999, abs_mass_error_da=0.001,
        ms2_spectrum_present=True, ms2_score=0.5,
    )
    baseline.update(over)
    return C.compute_confidence(_fields(**baseline))


def test_identity_override_promotes_overall_confidence():
    baseline = _baseline_conf()
    overridden = _baseline_conf(known_identity_override=True)
    assert overridden.identity_status == C.IDENTITY_STATUS_KNOWN
    assert overridden.identity_confidence is not None and overridden.identity_confidence >= 0.90
    assert overridden.overall >= 0.90
    # Promoted known-identity must stay strictly below the authentic-standard 100% ceiling.
    assert overridden.overall <= C.IDENTITY_PROMOTED_MAX
    assert overridden.overall < C.CONFIRMED_OVERALL
    assert overridden.tier == "High"
    assert overridden.msi_level == baseline.msi_level  # still Level 2, not fabricated L1


def test_identity_promotion_clamped_below_confirmed():
    # Even an over-eager curator override of 1.0 must not reach the Confirmed 100% floor.
    conf = _baseline_conf(known_identity_override={"confidence": 1.0, "reason": "over-eager"})
    assert conf.identity_status == C.IDENTITY_STATUS_KNOWN
    assert conf.overall == C.IDENTITY_PROMOTED_MAX
    assert conf.overall < C.CONFIRMED_OVERALL


def test_level1_confirmed_outranks_promoted_level2():
    # Ordering invariant: an authentic-standard-confirmed Level 1 must always read
    # above the best possible curator-known Level 2 promotion.
    l1 = C.compute_confidence(_fields(metabolon_level="1", truth_match=True))
    l2 = _baseline_conf(known_identity_override={"confidence": 1.0, "reason": "max curator override"})
    assert l1.overall > l2.overall
    assert l1.tier == "Confirmed"
    assert l2.tier == "High"


def test_identity_override_dict_uses_custom_confidence_and_reason():
    conf = _baseline_conf(known_identity_override={
        "confidence": 0.93, "reason": "Named library metabolite; identity known"
    })
    assert conf.identity_status == C.IDENTITY_STATUS_KNOWN
    assert conf.identity_confidence == 0.93
    assert "Named library metabolite" in conf.identity_basis


def test_identity_auto_rule_off_by_default_stays_putative():
    # source == exact + name_resolved + strong mass, but no override and auto-rule OFF.
    conf = _baseline_conf()
    assert conf.identity_status == C.IDENTITY_STATUS_PUTATIVE
    assert conf.identity_confidence is None


def test_identity_auto_rule_when_enabled_promotes_named_exact():
    C.IDENTITY_AUTO_RULE = True
    try:
        conf_ms2 = _baseline_conf()
        conf_no_ms2 = _baseline_conf(ms2_spectrum_present=False, ms2_score=None)
    finally:
        C.IDENTITY_AUTO_RULE = False
    assert conf_ms2.identity_status == C.IDENTITY_STATUS_KNOWN
    assert conf_ms2.identity_confidence == C.IDENTITY_CONF_AUTO_MS2
    assert conf_no_ms2.identity_status == C.IDENTITY_STATUS_KNOWN
    assert conf_no_ms2.identity_confidence == C.IDENTITY_CONF_AUTO_NO_MS2


def test_identity_unknown_for_predicted_source():
    conf = C.compute_confidence(_fields(source="predicted", structure_similarity=0.3))
    assert conf.identity_status == C.IDENTITY_STATUS_UNKNOWN
    assert conf.identity_confidence is None


def test_identity_putative_when_mass_mismatch():
    conf = _baseline_conf(mass_match_score=0.0, abs_mass_error_da=176.0)
    assert conf.identity_status == C.IDENTITY_STATUS_PUTATIVE


def test_identity_badge_and_card_render():
    conf = _baseline_conf(known_identity_override=True)
    badge = C.identity_badge_html(conf)
    assert "Identity: Known" in badge and "standard pending" in badge
    card = C.identity_metric_card_html(conf)
    assert "Library identity status" in card
    assert (
        "High Confidence ID, Authentic Standard Confirmation "
        "Pending to Promote to Level 1 Library"
    ) in card
    assert "metric-note" not in card


def test_known_identity_overall_copy():
    conf = _baseline_conf(known_identity_override=True)
    card = C.confidence_metric_card_html(conf)
    assert (
        "This is a High Confidence probable structure based on corroborating "
        "evidence derived from Metabolon’s validated structural elucidation workflow."
    ) in card


def test_mass_card_reports_observed_and_predicted_masses():
    conf = C.compute_confidence(
        _fields(
            mass_match_score=0.999,
            observed_neutral_mass=213.009476,
            candidate_exact_mass=213.009579,
            abs_mass_error_da=0.000103,
        )
    )
    cards = C.identification_parameter_cards_html(
        conf,
        candidate_scores=[],
        pathway_info={},
        top_similarity=0.0,
        neighbor_count=0,
    )
    assert "Strong \u00b7 1.00" in cards
    assert "observed neutral mass 213.009476 Da" in cards
    assert "predicted 213.009579 Da" in cards


def test_pathway_agreement_and_parameter_cards():
    candidates = [
        {
            "predicted_pathway": "Xenobiotics" if idx < 8 else "Lipid",
            "predicted_subpathway": "Chemical",
        }
        for idx in range(10)
    ]
    assert C.candidate_pathway_agreement_note(
        candidates, "Xenobiotics", "predicted_pathway", "Pathway"
    ) == "High Confidence: 8/10 top candidate structures map to this Pathway"
    conf = C.compute_confidence(_fields())
    cards = C.identification_parameter_cards_html(
        conf,
        candidate_scores=candidates,
        pathway_info={"super": "Xenobiotics", "sub": "Chemical", "similarity": 1.0},
        top_similarity=1.0,
        neighbor_count=10,
    )
    assert "High Confidence: 10/10 top candidate structures map to this SubPathway" in cards
    assert "similar to 10 Metabolon library standards" in cards
    assert "exactly matches a Metabolon library standard" in cards
    assert "Fragment Match Score" in cards


def test_report_copy_migration_and_relocation_are_idempotent():
    conf = _baseline_conf(known_identity_override=True)
    fields = _fields(
        display_name="Example compound",
        ms2_detail={"present": True, "n_pairs": 12, "n_ion_bins": 4},
    )
    cards = C.identification_parameter_cards_html(
        conf,
        candidate_scores=[],
        pathway_info={"super": "Xenobiotics", "sub": "Chemical"},
        top_similarity=1.0,
        neighbor_count=10,
    )
    metric = (
        '<div class="metric-card"><div class="metric-label">{}</div>'
        '<div class="metric-value">old</div><div class="metric-note">old</div></div>'
    )
    source = (
        "<body><p>Dashboard \u2022 Public</p>"
        '<span class="badge ok">top-1 candidate</span>'
        "<h3>Most Likely Structure</h3><p><strong>X-1</strong></p>"
        '<p class="note">Exact stereochemistry is not overclaimed; this is the '
        "top-ranked candidate structure for this feature.</p>"
        "<h3>Final Consensus Atom-Position Confidence</h3>"
        '<div class="metrics">'
        + metric.format("Top candidate structural similarity")
        + metric.format("Candidate Chemical Pathway")
        + metric.format("Candidate Chemical SubPathway")
        + metric.format("Library neighbors")
        + "</div>"
        "<h2>Parameters Contributing to Identification Confidence</h2><ul><li>old</li></ul>"
        "<h2>Fragment Scores</h2><p class='note'>old</p><div class='table-wrap'><table>"
        "<thead><tr><th>Fragment</th><th>Match Score</th></tr></thead>"
        "<tbody><tr><td>Lactoyl</td><td></td></tr></tbody></table></div>"
        "<h2>MS Fit Scores</h2><p>old detail</p></section></body>"
    )
    candidates = [
        {
            "rank": 1,
            "predicted_formula": "C6H6",
            "predicted_smiles": "c1ccccc1",
            "forward_spectrum_score": 0.4,
            "predicted_spectrum_score_v2": 0.2,
            "ms2_match_score": 0.1,
        }
    ]
    once = PC.patch_html(source, conf, fields, "X-1", cards, candidates)
    twice = PC.patch_html(once, conf, fields, "X-1", cards, candidates)
    assert twice == once
    assert once.count("<!--ID_PARAMETER_CARDS-->") == 1
    assert "Parameters Contributing to ID Confidence" in once
    assert "<h3>Atom-position confidence mapped to top candidate</h3>" in once
    assert "exact stereochemistry is not known from available evidence" in once
    assert "DB structures above similarity floor" not in once
    assert "Observed fragment peaks" in once
    assert "Forward spectral score" not in once
    assert "<!--FRAGMENT_SUMMARY-->" in once
    assert "<!--FRAGMENT_DIAGNOSTICS-->" in once


def test_identity_badge_hidden_when_not_known():
    conf = _baseline_conf()
    assert C.identity_badge_html(conf) == ""


def test_identity_note_public_omits_ms2():
    conf = _baseline_conf(known_identity_override=True)
    conf.ms2_n_pairs = 12  # simulate observed MS2 pair count present
    pub = C.identity_note_html(conf, public=True)
    internal = C.identity_note_html(conf, public=False)
    assert "MS2" not in pub
    assert "MS2" in internal


def test_identity_note_included_in_scored_params_section():
    conf = _baseline_conf(known_identity_override=True)
    section = C.scored_parameters_section_html(conf, public=True)
    assert "Library identity" in section


def test_load_identity_overrides_missing_file_returns_empty():
    assert C.load_identity_overrides(None) == {}
    assert C.load_identity_overrides(Path("/nonexistent/path.csv")) == {}


def test_load_identity_overrides_seeded_allowlist():
    path = ROOT / "data" / "known_identity_overrides.csv"
    if not path.exists():
        pytest.skip("seeded allowlist not present")
    overrides = C.load_identity_overrides(path)
    assert "X-22834" in overrides  # glycodeoxycholate
    assert overrides["X-22834"]["confidence"] >= 0.90
    # Row with a comma in the reason must round-trip because it's quoted in the CSV.
    assert "2,4-dichlorophenol" in overrides["X-21315"]["reason"]


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
