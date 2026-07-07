"""Smoke tests for the candidate-widget evidence surfacing + master-index MSI."""
# ruff: noqa: E501  (inline HTML fixtures are intentionally long)
from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "structural_elucidation" / "reports"))
sys.path.insert(0, str(ROOT / "tools"))

import enrich_report_evidence as E  # noqa: E402
import index_render as IR  # noqa: E402
import patch_candidate_scores as P  # noqa: E402


def test_candidate_evidence_rows_join_by_smiles():
    css = [
        {"rank": "1", "predicted_smiles": "CCO"},
        {"rank": "2", "predicted_smiles": "CCC"},
    ]
    tc = [
        {"predicted_smiles": "CCO", "evidence_confidence": "0.99", "mass_match_score": "0.98"},
        {"predicted_smiles": "CCC", "evidence_confidence": "0.30", "mass_match_score": "0.95"},
    ]
    rows = E.candidate_evidence_rows(css, tc)
    assert rows[0] == {"rank": 1, "evidence_confidence": 0.99, "mass_match_score": 0.98}
    assert rows[1]["evidence_confidence"] == 0.30


def test_candidate_evidence_rows_row_order_fallback():
    css = [{"rank": "1", "predicted_smiles": "X"}]  # smiles not in tc -> fallback to row 0
    tc = [{"predicted_smiles": "Y", "evidence_confidence": "0.5", "mass_match_score": "0.6"}]
    rows = E.candidate_evidence_rows(css, tc)
    assert rows[0]["evidence_confidence"] == 0.5


WIDGET_HTML = """<section>
      <p class="note">Select a candidate to view its consensus atom-position heatmap. The selected structure is the reference; atom circles show maximum-common-substructure retention across the top-10 candidate pool.</p>
      <select id="candidate-consensus-select">
        <option value="1" selected>Rank 1 \u00b7 C6H12O6 \u00b7 score 0.030</option>
<option value="2">Rank 2 \u00b7 C6H12O6 \u00b7 score 0.010</option>
      </select>
      <div id="candidate-consensus-summary"><div class='candidate-summary'><p><strong>Rank 1</strong> \u00b7 C6H12O6</p><p>Report score: <strong>0.0300</strong> \u00b7 Mass error (Da): <strong>0.001000</strong></p><p><code>CCO</code></p></div></div>
      <div class="consensus-box" id="candidate-consensus-panel"></div>
      <script type="application/json" id="candidate-consensus-data">{"maps": {}, "summaries": {"1": {"rank": 1, "predicted_formula": "C6H12O6", "predicted_smiles": "CCO", "candidate_report_score": 0.03, "abs_mass_error_da": 0.001}, "2": {"rank": 2, "predicted_formula": "C6H12O6", "predicted_smiles": "CCC", "candidate_report_score": 0.01, "abs_mass_error_da": 0.001}}}</script>
<script>
      function renderCandidateSummary(rank) {
        const summary = consensusSummaries[String(rank)] || {};
        const bits = [];
        candidateConsensusSummary.innerHTML = bits.join('');
      }
</script>
</section>"""


def _ev():
    return {
        1: {"evidence_confidence": 0.999, "mass_match_score": 0.998},
        2: {"evidence_confidence": 0.30, "mass_match_score": 0.95},
    }


def test_patch_html_surfaces_evidence_and_relabels():
    out = P.patch_html(WIDGET_HTML, _ev())
    assert out is not None
    assert "Evidence confidence: <strong>0.999</strong>" in out
    assert "Mass match: <strong>0.998</strong>" in out
    assert "Rank score (relative)" in out
    # dropdown relabeled from "score" to "conf"
    assert "conf 0.999" in out
    assert "score 0.030" not in out
    # embedded JSON carries the new fields
    assert "\"evidence_confidence\": 0.999" in out or "\"evidence_confidence\":0.999" in out
    # reading-the-scores note injected
    assert "Reading the scores" in out


def test_patch_html_idempotent():
    once = P.patch_html(WIDGET_HTML, _ev())
    assert once is not None
    twice = P.patch_html(once, _ev())
    assert twice is None  # sentinel present -> no re-patch


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
