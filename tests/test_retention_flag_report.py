"""The per-candidate retention screen on the internal report's top-10 list.

Demote-only and display-only: a flagged candidate keeps its rank and scores and
is listed last in the internal view; the public variant never shows the flag or
the reordering. The z math is the calibrated axis the confidence model already
applies to rank 1 (reports/confidence.py), extended to every displayed row.
"""
from __future__ import annotations

import math
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent

pytest.importorskip("pandas", reason="the report renderer needs pandas")
pytest.importorskip("rdkit", reason="the report renderer needs rdkit")

import pandas as pd  # noqa: E402

from structural_elucidation.reports import compound as R  # noqa: E402
from structural_elucidation.reports import confidence as C  # noqa: E402

# Committed lib-209 calibration: slope 5.1967, intercept -590.7, rmse 599.9.
# With the observed index at 5000 these proxies land where each test needs:
#   FLAG -> z = (5000 - (5.1967*300 - 590.7)) / 599.9 = 6.72   (flagged)
#   NEAR -> z = 2.90                                            (kept: bar is 3.0)
#   OK   -> z = 0.004                                           (kept)
PROXIES = {"FLAG": 300.0, "NEAR": 741.0, "OK": 1075.4}
OBSERVED_KEY = "feature_X_209_mz100.0_ri5000_X-9"


def fake_proxy(smiles: str) -> float:
    return PROXIES.get(smiles, math.nan)


# --- candidate_retention_flags -------------------------------------------------


def _frame(rows: list[dict]) -> pd.DataFrame:
    return pd.DataFrame(rows)


def test_flags_fire_on_the_calibrated_library_and_only_past_the_bar(monkeypatch):
    monkeypatch.setattr(R, "ri_proxy_from_smiles", fake_proxy)
    out = R.candidate_retention_flags(
        _frame(
            [
                {"rank": 1, "feature_key": OBSERVED_KEY, "predicted_smiles": "OK", "lib_id": "209"},
                {"rank": 2, "feature_key": OBSERVED_KEY, "predicted_smiles": "FLAG", "lib_id": "209"},
                {"rank": 3, "feature_key": OBSERVED_KEY, "predicted_smiles": "NEAR", "lib_id": "209"},
            ]
        )
    )
    assert list(out["retention_implausible"]) == [False, True, False]
    assert out["retention_z"].iloc[1] == pytest.approx(6.72, abs=0.01)
    # The worst true structure in the feasibility audit sits at z ~ 2.9, so the
    # near-threshold case must stay unflagged or the bar has moved.
    assert out["retention_z"].iloc[2] == pytest.approx(2.90, abs=0.01)


def test_rows_without_evidence_stay_unflagged_with_no_z(monkeypatch):
    monkeypatch.setattr(R, "ri_proxy_from_smiles", fake_proxy)
    out = R.candidate_retention_flags(
        _frame(
            [
                # Library 305 never cleared the feasibility audit.
                {"rank": 1, "feature_key": OBSERVED_KEY, "predicted_smiles": "FLAG", "lib_id": "305"},
                # RI 0 means "not measured".
                {"rank": 2, "feature_key": "f_mz100.0_ri0_X-9", "predicted_smiles": "FLAG", "lib_id": "209"},
                # rank_frame's fallback key embeds no index.
                {"rank": 3, "feature_key": "unified_rerank_X-9", "predicted_smiles": "FLAG", "lib_id": "209"},
                # No SMILES / no computable proxy.
                {"rank": 4, "feature_key": OBSERVED_KEY, "predicted_smiles": "", "lib_id": "209"},
                {"rank": 5, "feature_key": OBSERVED_KEY, "predicted_smiles": "MYSTERY", "lib_id": "209"},
            ]
        )
    )
    assert not out["retention_implausible"].any()
    assert out["retention_z"].isna().all()


def test_the_report_level_library_id_backfills_rows_without_one(monkeypatch):
    monkeypatch.setattr(R, "ri_proxy_from_smiles", fake_proxy)
    out = R.candidate_retention_flags(
        _frame([{"rank": 1, "feature_key": OBSERVED_KEY, "predicted_smiles": "FLAG"}]),
        default_lib_id=209.0,  # arrives as a float when it rode a numeric column
    )
    assert bool(out["retention_implausible"].iloc[0]) is True


def test_without_any_calibration_the_screen_disables_rather_than_zeroes(monkeypatch):
    monkeypatch.setattr(R, "ri_proxy_from_smiles", fake_proxy)
    monkeypatch.setattr(C, "load_retention_calibration", lambda path=None: {})
    out = R.candidate_retention_flags(
        _frame([{"rank": 1, "feature_key": OBSERVED_KEY, "predicted_smiles": "FLAG", "lib_id": "209"}])
    )
    assert not out["retention_implausible"].any()
    assert out["retention_z"].isna().all()


def test_an_empty_frame_still_gains_the_columns():
    out = R.candidate_retention_flags(pd.DataFrame())
    assert list(out.columns) == ["retention_z", "retention_implausible"]
    assert out.empty


# --- display-order policy ------------------------------------------------------


def test_flag_truthiness_survives_the_csv_round_trip():
    truthy = [True, 1, 1.0, "1", "1.0", "True", "true"]
    falsy = [False, 0, 0.0, "0", "0.0", "False", "false", "", None, "nan", float("nan")]
    assert all(C.retention_flag_truthy(value) for value in truthy)
    assert not any(C.retention_flag_truthy(value) for value in falsy)


def test_demoted_order_is_stable_and_touches_nothing_else():
    rows = [
        {"rank": 1, "retention_implausible": False},
        {"rank": 2, "retention_implausible": "True"},
        {"rank": 3, "retention_implausible": ""},
        {"rank": 4, "retention_implausible": 1},
        {"rank": 5},
    ]
    ordered = C.retention_demoted_order(rows)
    assert [row["rank"] for row in ordered] == [1, 3, 5, 2, 4]
    # Same objects, same ranks: demote-only means display order and nothing more.
    assert all(any(row is original for original in rows) for row in ordered)
    assert C.retention_demoted_order([]) == []


# --- rendered pages ------------------------------------------------------------

SMILES = "OC(=O)CCC(N)C(=O)O"


def _confidence_level_3() -> C.Confidence:
    fields = dict(
        source="predicted", pathway_vote_share=0.8, observed_mz_available=True,
        structure_similarity=0.60, mass_match_score=0.99,
        ms2_spectrum_present=False, ms2_score=None,
    )
    conf = C.compute_confidence(fields)
    assert conf.msi_level == 3, f"fixture drifted: got L{conf.msi_level}"
    return conf


def _candidate_rows(candidates: int, flags: dict[int, dict]) -> list[dict]:
    rows = []
    for rank in range(1, candidates + 1):
        row = {
            "rank": rank,
            "predicted_formula": "C5H9NO4",
            "predicted_smiles": SMILES,
            "candidate_report_score": 0.04 / rank,
            "forward_spectrum_score": 0.5,
            "predicted_spectrum_score_v2": 0.6,
            "retention_z": "",
            "retention_implausible": False,
        }
        row.update(flags.get(rank, {}))
        rows.append(row)
    return rows


def _render(*, variant: str, candidates: int = 4, flags: dict[int, dict] | None = None) -> str:
    rows = _candidate_rows(candidates, flags or {})
    return R.render_html(
        out_dir=Path("/tmp/does-not-need-to-exist"),
        target={
            "control_id": "X-99999",
            "name": "Glutamate-like conjugate",
            "smiles": SMILES,
            "formula": "C5H9NO4",
            "exact_mass": 147.053158,
            "inchikey": "WHUUTDBJXJRKMK-UHFFFAOYSA-N",
            "tentative_descriptor": {"label": "Glutamate-like conjugate", "confidence": "medium"},
        },
        evidence={"parent_mz": 146.0459, "mass_error_da": 0.0004},
        db_detail={},
        confidence=_confidence_level_3(),
        conf_fields={"ms2_detail": {"present": True, "n_pairs": 12}},
        pathway_info={"source": "predicted", "similarity": 0.6, "super": "Amino Acid"},
        synonyms=[],
        neighbors=[],
        nodes=[],
        edges=[],
        candidate_scores=rows,
        fragment_candidate_rows=rows,
        structure_panel_html="<svg id='the-structure-drawing'></svg>",
        consensus_svg="<svg id='the-consensus-map'></svg>",
        candidate_consensus_maps={rank: "<svg/>" for rank in range(1, candidates + 1)},
        candidate_consensus_summaries={},
        default_consensus_rank=1,
        structural_possibility_tree_svg="",
        structural_possibility_rows=[],
        feature_total_count=12,
        extended_network_svg_html="",
        extended_unknown_rows=[],
        extended_known_rows=[],
        extended_stats={},
        correlation_rows=[],
        correlation_blocker=None,
        variant=variant,
    )


FLAGGED_RANK_2 = {2: {"retention_z": 3.6, "retention_implausible": True}}


def _option_ranks(page: str) -> list[str]:
    return re.findall(r'<option value="(\d+)"', page)


def _fragment_table_ranks(page: str) -> list[str]:
    body = page.split("Per-candidate spectral fit scores", 1)[1].split("</table>", 1)[0]
    return re.findall(r"<tr><td>(\d+)</td>", body)


def test_internal_lists_a_flagged_candidate_last_and_says_why():
    page = _render(variant="internal", flags=FLAGGED_RANK_2)
    assert "Retention screen (demote-only)" in page
    assert "rank 2 (z=3.6)" in page
    assert "retention-implausible" in page
    # Demoted in display order; the rank label itself is untouched.
    assert _option_ranks(page) == ["1", "3", "4", "2"]
    assert "Rank 2" in page
    assert _fragment_table_ranks(page) == ["1", "3", "4"]


def test_public_never_shows_the_flag_or_the_reordering():
    page = _render(variant="public", flags=FLAGGED_RANK_2)
    assert "retention-implausible" not in page
    assert "Retention screen" not in page
    assert _option_ranks(page) == ["1", "2", "3"]
    assert _fragment_table_ranks(page) == ["1", "2", "3"]


def test_a_near_threshold_candidate_is_not_flagged_anywhere():
    page = _render(
        variant="internal",
        flags={2: {"retention_z": 2.9, "retention_implausible": False}},
    )
    assert "retention-implausible" not in page
    assert "Retention screen" not in page
    assert _option_ranks(page) == ["1", "2", "3", "4"]


def test_the_flag_still_demotes_after_a_csv_string_round_trip():
    page = _render(
        variant="internal",
        flags={2: {"retention_z": "3.6", "retention_implausible": "True"}},
    )
    assert _option_ranks(page) == ["1", "3", "4", "2"]
    assert "rank 2 (z=3.6)" in page


def test_an_unflagged_report_renders_without_the_retention_note():
    page = _render(variant="internal")
    assert "Retention screen" not in page
    assert _option_ranks(page) == ["1", "2", "3", "4"]
