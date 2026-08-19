"""Tests for the canonical evidence_confidence definition and its pool backfill."""

from __future__ import annotations

import pandas as pd
import pytest

from structural_elucidation.chemistry.evidence import (
    BASIS_BACKFILLED,
    BASIS_COLUMN,
    BASIS_EMITTED,
    BASIS_UNAVAILABLE,
    MS1_SUPPORT_BONUS,
    backfill_evidence_confidence,
    evidence_confidence,
)
from structural_elucidation.pipeline.pool_v1 import apply_evidence_backfill


def test_evidence_confidence_takes_strongest_term():
    assert evidence_confidence(formula_prior=0.4, ms2_match=0.9, nce_pipeline=0.2) == pytest.approx(0.9)


def test_evidence_confidence_adds_ms1_bonus_and_clips():
    assert evidence_confidence(formula_prior=0.5, ms1_supported=True) == pytest.approx(0.5 + MS1_SUPPORT_BONUS)
    assert evidence_confidence(formula_prior=1.0, ms1_supported=True) == pytest.approx(1.0)


def test_evidence_confidence_tolerates_blank_and_string_inputs():
    assert evidence_confidence(formula_prior="", ms2_match="0.75", ms1_supported="false") == pytest.approx(0.75)


def test_emitted_values_are_never_overwritten():
    df = pd.DataFrame({"evidence_confidence": [0.5], "formula_prior_score": [1.0]})
    score, basis = backfill_evidence_confidence(df)
    assert score.tolist() == [0.5]
    assert basis.tolist() == [BASIS_EMITTED]


def test_backfills_from_existing_formula_prior():
    """The truth-blind operator enumerators emit a prior but no evidence column."""
    df = pd.DataFrame({"evidence_confidence": [0.0], "formula_prior_score": [0.968185]})
    score, basis = backfill_evidence_confidence(df)
    assert score.tolist() == pytest.approx([0.968185])
    assert basis.tolist() == [BASIS_BACKFILLED]


def test_backfills_from_missing_column():
    df = pd.DataFrame({"formula_prior_score": [0.8]})
    score, basis = backfill_evidence_confidence(df)
    assert score.tolist() == pytest.approx([0.8])
    assert basis.tolist() == [BASIS_BACKFILLED]


def test_recomputes_prior_from_formula_and_target_mass():
    df = pd.DataFrame(
        {
            "evidence_confidence": [0.0],
            "formula_prior_score": [0.0],
            "candidate_formula": ["C9H11N3O4"],
            "target_neutral_mass": [225.074955832],
        }
    )
    score, basis = backfill_evidence_confidence(df)
    assert score.iloc[0] > 0.99
    assert basis.tolist() == [BASIS_BACKFILLED]


def test_no_prior_without_a_mass_anchor():
    """A formula alone scores 0.65 on element ratios; that is not evidence."""
    df = pd.DataFrame(
        {
            "evidence_confidence": [0.0],
            "formula_prior_score": [0.0],
            "candidate_formula": ["C9H11N3O4"],
        }
    )
    score, basis = backfill_evidence_confidence(df)
    assert score.tolist() == [0.0]
    assert basis.tolist() == [BASIS_UNAVAILABLE]


def test_spectral_terms_carry_rows_with_no_formula_evidence():
    df = pd.DataFrame({"evidence_confidence": [0.0], "ms2_match_score": [0.42], "ms1_buddy_supported": ["True"]})
    score, basis = backfill_evidence_confidence(df)
    assert score.tolist() == pytest.approx([0.42 + MS1_SUPPORT_BONUS])
    assert basis.tolist() == [BASIS_BACKFILLED]


def test_backfill_is_idempotent():
    df = pd.DataFrame({"evidence_confidence": [0.0, 0.0], "formula_prior_score": [0.9, 0.0]})
    once, _ = backfill_evidence_confidence(df)
    df["evidence_confidence"] = once
    twice, basis = backfill_evidence_confidence(df)
    assert twice.tolist() == pytest.approx(once.tolist())
    assert basis.tolist() == [BASIS_EMITTED, BASIS_UNAVAILABLE]


def test_pool_stage_reports_backfill_counts():
    pool = pd.DataFrame(
        {
            "evidence_confidence": [0.997, 0.0, 0.0],
            "formula_prior_score": [0.997, 0.968, 0.0],
        }
    )
    out, stats = apply_evidence_backfill(pool)
    assert stats == {"n_evidence_backfilled": 1, "n_evidence_unavailable": 1}
    assert out[BASIS_COLUMN].tolist() == [BASIS_EMITTED, BASIS_BACKFILLED, BASIS_UNAVAILABLE]
    assert out["evidence_confidence"].tolist() == pytest.approx([0.997, 0.968, 0.0])


def test_pool_stage_handles_empty_pool():
    out, stats = apply_evidence_backfill(pd.DataFrame())
    assert out.empty
    assert stats == {"n_evidence_backfilled": 0, "n_evidence_unavailable": 0}
