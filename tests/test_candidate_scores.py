"""Tests for the rank-consistent candidate confidence + delivered-set patcher."""
# ruff: noqa: E501  (inline HTML fixtures are intentionally long)
from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "structural_elucidation" / "reports"))
sys.path.insert(0, str(ROOT / "tools"))

import confidence as C  # noqa: E402
import enrich_report_evidence as E  # noqa: E402
import index_render as IR  # noqa: E402
import patch_candidate_scores as P  # noqa: E402

# ── candidate_confidences: monotone with rank, anchored, tie-aware ───────────

def test_candidate_confidences_monotone_and_anchored():
    rows = [
        {"rank": 1, "candidate_report_score": 0.50},
        {"rank": 2, "candidate_report_score": 0.20},
        {"rank": 3, "candidate_report_score": 0.05},
    ]
    out = C.candidate_confidences(rows, feature_overall=0.62)
    vals = [o["candidate_confidence"] for o in out]
    assert vals == sorted(vals, reverse=True)      # monotone non-increasing with rank
    assert vals[0] == 0.62                          # rank 1 anchored to report overall
    assert vals[1] < vals[0] and vals[2] < vals[1]  # discriminates by ordering score
    assert not any(o["tied"] for o in out)


def test_candidate_confidences_fixes_inversion():
    """A saturated/inverted evidence column becomes monotone with rank."""
    rows = [
        {"rank": 1, "candidate_report_score": 0.50, "evidence_confidence": 0.30},
        {"rank": 2, "candidate_report_score": 0.20, "evidence_confidence": 0.99},
        {"rank": 3, "candidate_report_score": 0.05, "evidence_confidence": 0.95},
    ]
    vals = [o["candidate_confidence"] for o in C.candidate_confidences(rows, 0.6)]
    assert vals == sorted(vals, reverse=True)


def test_candidate_confidences_all_equal_scores_are_flagged_tied():
    rows = [{"rank": i + 1, "candidate_report_score": 0.036} for i in range(5)]
    out = C.candidate_confidences(rows, feature_overall=0.44)
    assert all(o["candidate_confidence"] == 0.44 for o in out)   # honest, not saturated to 0.99
    assert all(o["tied"] for o in out)                            # surfaced as one tie group
    assert len({o["tie_group"] for o in out}) == 1


def test_candidate_confidences_enforces_monotonicity_over_noisy_scores():
    """Even if the ordering column is out of order, the result cannot invert."""
    rows = [
        {"rank": 1, "candidate_report_score": 0.20},
        {"rank": 2, "candidate_report_score": 0.50},  # higher than rank 1 (noisy)
        {"rank": 3, "candidate_report_score": 0.10},
    ]
    vals = [o["candidate_confidence"] for o in C.candidate_confidences(rows, 0.6)]
    assert vals == sorted(vals, reverse=True)


def test_compute_confidence_exposes_raw_evidence_anchor():
    fields = {
        "source": "predicted", "structure_similarity": 0.9, "mass_match_score": 1.0,
        "ms2_spectrum_present": True, "ms2_score": 0.5, "observed_mz_available": True,
    }
    conf = C.compute_confidence(fields)
    assert 0.0 <= conf.raw_evidence <= 1.0
    # raw evidence is fuller-range than the compressed Platt overall for strong evidence,
    # so anchoring the featured candidate confidence to it un-deflates strong best matches.
    assert conf.raw_evidence > conf.overall
    assert conf.as_dict()["raw_evidence"] == conf.raw_evidence


def test_msi_sublevel_curated_surface_and_defensible_derivation():
    base = {"source": "predicted", "structure_similarity": 0.6, "mass_match_score": 1.0, "observed_mz_available": True}
    # curated SubLevel is surfaced authoritatively and NOT overwritten by derivation
    c = C.compute_confidence({**base, "metabolon_sublevel": "3B", "n_isomers": 1})
    assert c.metabolon_sublevel == "3B"
    assert c.msi_sublevel == ""  # fine letters are not fabricated
    # positional-isomer (multiple same-formula candidates) is the one derived fine bin
    c2 = C.compute_confidence({**base, "n_isomers": 3})
    assert c2.msi_level == 3 and c2.msi_sublevel == "3A"
    # non-metabolite curated tokens are ignored
    c3 = C.compute_confidence({**base, "metabolon_sublevel": "Exclude", "n_isomers": 1})
    assert c3.metabolon_sublevel == ""
    assert c.as_dict()["metabolon_sublevel"] == "3B"


def test_candidate_confidences_empty_and_degenerate():
    assert C.candidate_confidences([], 0.5) == []
    out = C.candidate_confidences([{"rank": 1, "candidate_report_score": 0.0}], None)
    assert out[0]["candidate_confidence"] == 0.0


# ── candidate_evidence_rows: join + confidence, anchored to overall ──────────

def test_candidate_evidence_rows_join_and_confidence():
    css = [
        {"rank": "1", "predicted_smiles": "CCO", "candidate_report_score": "0.5"},
        {"rank": "2", "predicted_smiles": "CCC", "candidate_report_score": "0.2"},
    ]
    tc = [
        {"predicted_smiles": "CCO", "evidence_confidence": "0.99", "mass_match_score": "0.98"},
        {"predicted_smiles": "CCC", "evidence_confidence": "0.30", "mass_match_score": "0.95"},
    ]
    rows = E.candidate_evidence_rows(css, tc, feature_overall=0.6)
    assert rows[0]["rank"] == 1 and rows[0]["evidence_confidence"] == 0.99
    assert rows[1]["evidence_confidence"] == 0.30
    # candidate_confidence is present, anchored at rank 1, and monotone
    assert rows[0]["candidate_confidence"] == 0.6
    assert rows[0]["candidate_confidence"] >= rows[1]["candidate_confidence"]


def test_candidate_evidence_rows_row_order_fallback():
    css = [{"rank": "1", "predicted_smiles": "X", "candidate_report_score": "0.4"}]
    tc = [{"predicted_smiles": "Y", "evidence_confidence": "0.5", "mass_match_score": "0.6"}]
    rows = E.candidate_evidence_rows(css, tc, feature_overall=0.5)
    assert rows[0]["evidence_confidence"] == 0.5
    assert rows[0]["candidate_confidence"] == 0.5


# ── shared score-line HTML (source/patch parity) ────────────────────────────

def test_candidate_score_line_features_confidence_not_rank_score():
    row = {"candidate_confidence": 0.44, "mass_match_score": 0.998, "evidence_confidence": 0.999, "abs_mass_error_da": 0.001}
    line = C.candidate_score_line_html(row)
    assert "Candidate confidence: <strong>0.440</strong>" in line
    assert "Mass match: <strong>0.998</strong>" in line
    assert "Formula/class evidence: <strong>0.999</strong>" in line
    assert "Rank score (relative)" not in line


def test_candidate_score_line_tie_note():
    row = {"candidate_confidence": 0.44, "candidate_confidence_tied": True}
    assert "not resolved" in C.candidate_score_line_html(row)


# ── patcher v2: HTML fixtures ───────────────────────────────────────────────

def _widget(option_suffix: str, note: str) -> str:
    """A candidate widget fragment; option_suffix mimics original ('score X') or v1 ('conf X')."""
    return (
        "<section>\n"
        '      <p class="note">Select a candidate to view its consensus atom-position heatmap. The selected structure is the reference; atom circles show maximum-common-substructure retention across the top-10 candidate pool.</p>\n'
        f"{note}"
        '      <select id="candidate-consensus-select">\n'
        f'        <option value="1" selected>Rank 1 \u00b7 C6H12O6{option_suffix}</option>\n'
        f'<option value="2">Rank 2 \u00b7 C6H12O6{option_suffix}</option>\n'
        "      </select>\n"
        '      <div id="candidate-consensus-summary"><div class=\'candidate-summary\'><p><strong>Rank 1</strong> \u00b7 C6H12O6</p><p>old card</p><p><code>CCO</code></p></div></div>\n'
        '      <div class="consensus-box" id="candidate-consensus-panel"></div>\n'
        '      <script type="application/json" id="candidate-consensus-data">{"maps": {}, "summaries": {"1": {"rank": 1, "predicted_formula": "C6H12O6", "predicted_smiles": "CCO", "candidate_report_score": 0.03, "abs_mass_error_da": 0.001}, "2": {"rank": 2, "predicted_formula": "C6H12O6", "predicted_smiles": "CCC", "candidate_report_score": 0.01, "abs_mass_error_da": 0.001}}}</script>\n'
        "<script>\n"
        "      function renderCandidateSummary(rank) {\n"
        "        const summary = consensusSummaries[String(rank)] || {};\n"
        "        const bits = [];\n"
        "      }\n"
        "</script>\n"
        "</section>"
    )


V1_NOTE = '      <p class="note"><strong>Reading the scores:</strong> "Evidence confidence" (library + mass + spectral corroboration) and "Mass match" are the identification signals. The "Rank score (relative)" is a small, uncalibrated value the reranker uses only to order otherwise-equivalent candidates.</p>\n'


def _ev():
    return {
        1: {"candidate_confidence": 0.44, "candidate_confidence_tied": False, "evidence_confidence": 0.999, "mass_match_score": 0.998},
        2: {"candidate_confidence": 0.18, "candidate_confidence_tied": False, "evidence_confidence": 0.30, "mass_match_score": 0.95},
    }


def test_patch_v1_report_to_v2():
    out = P.patch_html(_widget(" \u00b7 conf 0.999", V1_NOTE), _ev())
    assert out is not None
    # featured candidate confidence (anchored rank-1 value) + relabeled class evidence
    assert "Candidate confidence: <strong>0.440</strong>" in out
    assert "Formula/class evidence: <strong>0.999</strong>" in out
    # dropdown now shows the candidate_confidence, not the saturated evidence value
    assert "conf 0.440" in out
    assert "conf 0.999" not in out
    # the misleading raw rank score label is gone
    assert "Rank score (relative)" not in out
    # v2 sentinel present and the note is replaced (not duplicated)
    assert P.V2_SENTINEL in out
    assert out.count("Reading the scores") == 1
    assert "Candidate confidence" in out.split("Reading the scores", 1)[1]


def test_patch_original_report_to_v2_injects_note():
    out = P.patch_html(_widget(" \u00b7 score 0.030", ""), _ev())
    assert out is not None
    assert "conf 0.440" in out and "score 0.030" not in out
    assert P.V2_SENTINEL in out
    assert out.count("Reading the scores") == 1  # injected after the anchor


def test_patch_is_idempotent():
    once = P.patch_html(_widget(" \u00b7 conf 0.999", V1_NOTE), _ev())
    assert once is not None
    assert P.patch_html(once, _ev()) is None  # v2 sentinel -> no re-patch


def test_patch_adds_count_note_when_fewer_than_10():
    out = P.patch_html(_widget(" \u00b7 conf 0.999", V1_NOTE), _ev())  # _ev has 2 candidates
    assert "of up to 10 candidate structures" in out
    assert "<strong>2</strong> of up to 10" in out


def test_patch_marks_ties():
    ev = {
        1: {"candidate_confidence": 0.44, "candidate_confidence_tied": True, "evidence_confidence": 0.999, "mass_match_score": 0.998},
        2: {"candidate_confidence": 0.44, "candidate_confidence_tied": True, "evidence_confidence": 0.999, "mass_match_score": 0.998},
    }
    out = P.patch_html(_widget(" \u00b7 conf 0.999", V1_NOTE), ev)
    assert "(tied)" in out                # dropdown marks the tie
    assert "not resolved" in out          # card explains it


def test_patch_card_uses_shared_score_line():
    """The patched card is byte-identical to the shared confidence.candidate_score_line_html."""
    out = P.patch_html(_widget(" \u00b7 conf 0.999", V1_NOTE), _ev())
    expected = C.candidate_score_line_html(
        {"candidate_confidence": 0.44, "candidate_confidence_tied": False,
         "evidence_confidence": 0.999, "mass_match_score": 0.998, "abs_mass_error_da": 0.001}
    )
    assert expected in out


# ── master index MSI column (unchanged) ─────────────────────────────────────

def test_index_render_has_msi_column_and_sort():
    rows = [
        IR.IndexRow(control_id="X-1", lib="209", name="foo", formula="C6", msi_level="3", msi_label="Tentative candidate structure"),
        IR.IndexRow(control_id="X-2", lib="209", name="bar", formula="C7", msi_level="1", msi_label="Confirmed (authentic standard)"),
    ]
    html_out = IR.render_master_index(rows, "2026-01-01")
    assert "<th title='MSI identification level (1-5)'>MSI</th>" in html_out
    assert "data-sort='3'" in html_out and ">L3<" in html_out
    assert "data-dir" in html_out and "cellVal" in html_out  # sortable JS


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
