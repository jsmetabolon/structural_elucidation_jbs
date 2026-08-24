"""Parity tests for candidates-stage mass matching and memoized RDKit helpers.

The sorted-index ppm lookup and the memo caches must reproduce the original
full-column scan and per-row RDKit calls exactly (same rows, same order,
same values), including tolerance-boundary and NaN edge cases.
"""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest

from structural_elucidation.pipeline import candidates as cand
from structural_elucidation.utils.denovo import canonicalize_smiles, exact_mass_from_smiles, formula_from_smiles
from structural_elucidation.chemistry.formula import formula_prior_score


def _reference_positions(df: pd.DataFrame, target: float, ppm: float, limit: int) -> list[int]:
    """Row positions selected by the original mask + head implementation."""
    mask = cand.ppm_hit_mask(df["library_exact_mass"], target, ppm)
    selected = df[mask].head(limit)
    positional = {label: pos for pos, label in enumerate(df.index)}
    return [positional[label] for label in selected.index]


def test_ppm_match_positions_matches_mask_on_random_library():
    rng = np.random.RandomState(3)
    masses = np.concatenate(
        [
            rng.uniform(50, 1500, size=4000),
            rng.choice(rng.uniform(50, 1500, size=50), size=200),  # duplicates
            [np.nan] * 25,
        ]
    )
    rng.shuffle(masses)
    df = pd.DataFrame({"library_exact_mass": masses})
    index = cand.build_mass_index(df)
    targets = list(rng.uniform(50, 1500, size=200))
    # Exact-boundary targets: value +/- exactly its ppm tolerance.
    for value in masses[:20]:
        if np.isfinite(value):
            tol = abs(float(value)) * 5.0 / 1_000_000.0
            targets.extend([float(value), float(value) + tol, float(value) - tol])
    for target in targets:
        for limit in (5, 120):
            fast = cand.ppm_match_positions(index, target, 5.0, limit).tolist()
            reference = _reference_positions(df, target, 5.0, limit)
            assert fast == reference, f"mismatch at target={target!r} limit={limit}"


def test_ppm_match_positions_with_nonrange_index():
    df = pd.DataFrame(
        {"library_exact_mass": [100.0, 100.0000004, 200.0, np.nan, 100.0000002]},
        index=[10, 20, 30, 40, 50],
    )
    index = cand.build_mass_index(df)
    fast = cand.ppm_match_positions(index, 100.0, 5.0, 120).tolist()
    assert fast == _reference_positions(df, 100.0, 5.0, 120)
    # iloc-selected rows equal the mask-selected rows.
    mask_rows = df[cand.ppm_hit_mask(df["library_exact_mass"], 100.0, 5.0)].head(120)
    pd.testing.assert_frame_equal(df.iloc[fast], mask_rows)


def test_ppm_match_positions_all_nan_masses():
    df = pd.DataFrame({"library_exact_mass": [np.nan, np.nan]})
    index = cand.build_mass_index(df)
    assert cand.ppm_match_positions(index, 100.0, 5.0, 10).tolist() == []


@pytest.mark.parametrize("value", ["CCO", "OCC", "", "not_a_smiles", None, float("nan"), 0])
def test_cached_canonical_matches_direct(value):
    assert cand.cached_canonical_smiles(value) == canonicalize_smiles(value)
    # Second call returns the memoized result unchanged.
    assert cand.cached_canonical_smiles(value) == canonicalize_smiles(value)


def test_cached_mass_and_formula_match_direct():
    for smiles in ["CCO", "c1ccccc1", "C(=O)O"]:
        canonical = canonicalize_smiles(smiles)
        direct_mass = exact_mass_from_smiles(canonical)
        cached_mass = cand.cached_exact_mass(canonical)
        assert (cached_mass == direct_mass) or (math.isnan(cached_mass) and math.isnan(direct_mass))
        assert cand.cached_formula(canonical) == formula_from_smiles(canonical)


def test_cached_formula_prior_matches_direct():
    for formula, mass in [("C2H6O", 46.0418648), ("C2H6O", float("nan")), ("C6H6", 78.04695)]:
        assert cand.cached_formula_prior(formula, mass, 0.05) == pytest.approx(
            formula_prior_score(formula, mass, tolerance_da=0.05), nan_ok=True
        )


def test_canonicalize_many_serial_and_parallel_match():
    values = ["CCO", "OCC", "not_a_smiles", "", "c1ccccc1O"] * 600  # crosses pool threshold
    serial = cand._canonicalize_many(values, workers=1)
    parallel = cand._canonicalize_many(values, workers=2)
    assert serial == parallel
    assert serial[0] == (canonicalize_smiles("CCO") or "")
