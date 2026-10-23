"""Unit tests for the top-K pathway ensemble (T3)."""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parent.parent


def load_tool(name: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "tools" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


E = load_tool("pathway_topk_ensemble")


def _rows(control: str, entries: list[tuple[int, str, str, float]]) -> list[dict]:
    return [
        {
            "control_id": control,
            "rank": rank,
            "predicted_pathway": super_,
            "predicted_subpathway": sub,
            "rerank_score": score,
        }
        for rank, super_, sub, score in entries
    ]


def test_uniform_majority_overrides_rank1_minority():
    frame = pd.DataFrame(
        _rows(
            "X-1",
            [
                (1, "Lipid", "Sterol", 0.9),
                (2, "Amino Acid", "Glycine Metabolism", 0.8),
                (3, "Amino Acid", "Alanine Metabolism", 0.7),
                (4, "Amino Acid", "Glycine Metabolism", 0.6),
            ],
        )
    )
    out = E.build_ensembles(frame, ("rank1", "uniform"), score_col="rerank_score")
    assert out["rank1"].iloc[0]["predicted_super"] == "Lipid"
    uniform = out["uniform"].iloc[0]
    assert uniform["predicted_super"] == "Amino Acid"
    assert uniform["predicted_sub"] == "Glycine Metabolism"  # voted within winning super
    assert uniform["vote_share"] == pytest.approx(0.75)


def test_score_weighted_respects_scores_and_handles_negatives():
    frame = pd.DataFrame(
        _rows(
            "X-1",
            [
                (1, "Lipid", "", 5.0),
                (2, "Amino Acid", "", -1.0),
                (3, "Amino Acid", "", -2.0),
            ],
        )
    )
    out = E.build_ensembles(frame, ("score_weighted",), score_col="rerank_score")
    # Shift-normalized weights: Lipid 7, AA 1 + 0 (+eps) -> Lipid wins.
    assert out["score_weighted"].iloc[0]["predicted_super"] == "Lipid"


def test_agreement_gate_abstains_on_disagreement():
    frame = pd.DataFrame(
        _rows(
            "X-1",
            [
                (1, "Lipid", "", 1.0),
                (2, "Amino Acid", "", 1.0),
                (3, "Xenobiotics", "", 1.0),
            ],
        )
    )
    out = E.build_ensembles(frame, ("agreement_gate",), score_col=None, agreement_threshold=0.6)
    gated = out["agreement_gate"].iloc[0]
    assert gated["predicted_super"] == ""  # abstained: rank1 share 1/3 < 0.6
    concordant = pd.DataFrame(
        _rows("X-2", [(1, "Lipid", "", 1.0), (2, "Lipid", "", 1.0), (3, "Amino Acid", "", 1.0)])
    )
    out2 = E.build_ensembles(concordant, ("agreement_gate",), score_col=None, agreement_threshold=0.6)
    assert out2["agreement_gate"].iloc[0]["predicted_super"] == "Lipid"


def test_unlabeled_candidates_are_ignored_and_empty_pool_abstains():
    frame = pd.DataFrame(
        _rows("X-1", [(1, "", "", 1.0), (2, "Lipid", "", 0.5)])
        + _rows("X-2", [(1, "", "", 1.0), (2, "", "", 0.5)])
    )
    out = E.build_ensembles(frame, ("uniform",), score_col=None)
    by_control = out["uniform"].set_index("control_id")
    assert by_control.loc["X-1", "predicted_super"] == "Lipid"
    assert by_control.loc["X-2", "predicted_super"] == ""
    assert by_control.loc["X-2", "n_labeled"] == 0


def test_rank1_uses_best_labeled_rank_and_k_caps_the_pool():
    frame = pd.DataFrame(
        _rows(
            "X-1",
            [(1, "", "", 1.0), (2, "Lipid", "", 0.9), (3, "Amino Acid", "", 0.8), (4, "Amino Acid", "", 0.7)],
        )
    )
    out = E.build_ensembles(frame, ("rank1", "uniform"), k=3, score_col=None)
    # Rank-1 falls to the best labeled candidate (rank 2).
    assert out["rank1"].iloc[0]["predicted_super"] == "Lipid"
    # k=3 keeps ranks 1-3 only: Lipid 1 vs Amino Acid 1 -> alphabetical tie-break.
    assert out["uniform"].iloc[0]["predicted_super"] == "Amino Acid"


def test_pick_score_column_prefers_request_then_candidates():
    frame = pd.DataFrame({"rerank_score": [1.0], "score": [2.0]})
    assert E.pick_score_column(frame, None) == "rerank_score"
    assert E.pick_score_column(frame, "score") == "score"
    with pytest.raises(ValueError):
        E.pick_score_column(frame, "missing_col")
    assert E.pick_score_column(pd.DataFrame({"rank": [1]}), None) is None
