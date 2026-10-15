"""Reasons a compound is only in the all-adduct-derived summary."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))

from compare_fallback_summary import reason  # noqa: E402


def row(**overrides) -> pd.Series:
    base = {
        "eligible": True,
        "no_inputs": False,
        "fallback_n_le10": 0,
        "earlier_base_pool_only": False,
        "fallback_best_ppm": float("nan"),
    }
    base.update(overrides)
    return pd.Series(base)


def test_reason_order() -> None:
    assert reason(row(eligible=False, no_inputs=True)) == "ineligible"
    assert reason(row(no_inputs=True, fallback_n_le10=3)) == "no_inputs"
    assert reason(row(fallback_n_le10=2, earlier_base_pool_only=True)) == "fallback_mass_only"
    assert reason(row(earlier_base_pool_only=True, fallback_best_ppm=40.0)) == "fallback_base_pool_not_rerun"
    assert reason(row(fallback_best_ppm=40.0)) == "fallback_over_10ppm"
    assert reason(row()) == "fallback_no_candidate"
