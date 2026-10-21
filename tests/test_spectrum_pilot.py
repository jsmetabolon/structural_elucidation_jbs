"""Cohort labels and same-formula isomer enumeration."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))

from spectrum_pilot import (  # noqa: E402
    NO_CANDIDATE,
    NO_OBSERVED,
    OFF_MASS,
    ON_MASS,
    assign_cohort,
    enumerate_constitutional_isomers,
    seed_pool,
)


def test_no_inputs_stay_unpublished_even_with_an_on_mass_seed() -> None:
    assert assign_cohort(no_inputs=True, has_ion_bins=True, on_mass=True, has_candidate=True) == NO_CANDIDATE


def test_ions_and_on_mass_seed_are_the_pilot_cohort() -> None:
    assert assign_cohort(no_inputs=False, has_ion_bins=True, on_mass=True, has_candidate=True) == ON_MASS
    assert assign_cohort(no_inputs=False, has_ion_bins=False, on_mass=True, has_candidate=True) == NO_OBSERVED
    assert assign_cohort(no_inputs=False, has_ion_bins=True, on_mass=False, has_candidate=True) == OFF_MASS
    assert assign_cohort(no_inputs=False, has_ion_bins=True, on_mass=False, has_candidate=False) == NO_CANDIDATE


def test_stored_mass_seed_wins_inside_10_ppm() -> None:
    assert seed_pool(4.0, 1.0) == "stored"
    assert seed_pool(12.0, 3.0) == "adduct_fallback"
    assert seed_pool(None, None) == ""


def test_enumerated_isomers_keep_the_seed_formula_and_drop_the_seed() -> None:
    pytest.importorskip("rdkit")
    isomers = enumerate_constitutional_isomers("CC(=O)O", limit=10)
    assert isomers
    assert "CC(=O)O" not in isomers
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors

    seed = rdMolDescriptors.CalcMolFormula(Chem.MolFromSmiles("CC(=O)O"))
    for smiles in isomers:
        assert rdMolDescriptors.CalcMolFormula(Chem.MolFromSmiles(smiles)) == seed
