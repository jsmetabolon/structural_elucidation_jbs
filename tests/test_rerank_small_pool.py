"""A tiny labeled pool must still produce scores."""

from __future__ import annotations

import pandas as pd
import pytest

pytest.importorskip("joblib")

from structural_elucidation.pipeline.rerank import select_model


def test_fewer_than_five_training_rows_uses_the_heuristic() -> None:
    pool = pd.DataFrame(
        {
            "control_id": ["X-1", "X-1", "X-2"],
            "label_ge_0p7": [1, 0, 1],
            "predicted_spectrum_score_v2": [0.2, 0.0, 0.4],
            "mass_match_score": [0.9, 0.9, 0.8],
            "formula_prior_score": [0.1, 0.1, 0.1],
            "evidence_confidence": [0.5, 0.4, 0.6],
            "is_library_mlims": [0, 0, 0],
            "is_peptide": [0, 0, 0],
            "min_source_rank_score": [1, 1, 1],
        }
    )
    selected = select_model(pool)
    assert selected["name"] == "heuristic"
    assert selected["model"] is None
    assert len(selected["scores"]) == 3
