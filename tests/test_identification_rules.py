"""Tests for the identification-rule audit harness.

The audit's whole claim is that a stored neutral mass does not follow from the ion
its own label names. That claim is only as good as the inversion used to check it,
so what is tested here is the arithmetic and the label grammar: a parser that
silently mis-read ``(m-2H)/2`` or ``+62.9959`` would manufacture the finding.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent


def load_tool(name: str):
    """Load a tool module by path.

    Registered in ``sys.modules`` before execution because ``@dataclass`` resolves
    annotations through ``sys.modules[cls.__module__]``, which is absent for a
    module loaded from a spec and never registered.
    """

    spec = importlib.util.spec_from_file_location(name, ROOT / "tools" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


A = load_tool("audit_identification_rules")

PROTON = A.ELEMENT_MASS["H"] - A.ELECTRON_MASS


# --- label grammar -----------------------------------------------------------


def test_plain_deprotonation_parses_to_a_single_hydrogen_loss():
    adduct = A.parse_adduct_label("CHRO_LIB:m-H")
    assert adduct is not None
    assert adduct.multiplier == 1
    assert adduct.charge == 1
    assert adduct.isotopes == ()
    assert adduct.is_plain_protonation
    assert adduct.loss_mass == pytest.approx(A.ELEMENT_MASS["H"])


def test_protonation_and_deprotonation_are_the_same_label_shape():
    """Polarity lives in the ion mode column, not in the label."""

    assert A.parse_adduct_label("m+H").is_plain_protonation
    assert A.parse_adduct_label("m-H").is_plain_protonation


def test_multimer_prefix_is_read_as_the_molecule_count():
    assert A.parse_adduct_label("CHRO_LIB:2m-H").multiplier == 2


def test_charge_suffix_is_read_from_the_parenthesised_form():
    adduct = A.parse_adduct_label("CHRO_LIB:(m-2H)/2")
    assert adduct.charge == 2
    assert adduct.loss_mass == pytest.approx(2 * A.ELEMENT_MASS["H"])


def test_isotope_tokens_survive_the_parenthesised_form():
    adduct = A.parse_adduct_label("CHRO_LIB:(m-2H[13C1])/2")
    assert adduct.charge == 2
    assert adduct.isotopes == (("13C", 1),)
    assert adduct.element_hints == {"C": 1}


def test_repeated_isotope_tokens_accumulate_per_element():
    adduct = A.parse_adduct_label("CHRO_LIB:m-H[13C1][18O1]")
    assert dict(adduct.isotopes) == {"13C": 1, "18O": 1}
    assert adduct.element_hints == {"C": 1, "O": 1}


def test_halogen_isotope_token_names_the_element_the_formula_must_carry():
    assert A.parse_adduct_label("CHRO_LIB:m-H[37Cl2]").element_hints == {"Cl": 2}
    assert A.parse_adduct_label("CHRO_LIB:m-H[81Br1]").element_hints == {"Br": 1}


def test_named_neutral_gains_are_summed_as_formulas():
    adduct = A.parse_adduct_label("CHRO_LIB:m-H+CO2")
    expected = A.ELEMENT_MASS["C"] + 2 * A.ELEMENT_MASS["O"]
    assert adduct.gain_mass == pytest.approx(expected)
    assert not adduct.is_plain_protonation


def test_a_bare_numeric_offset_is_not_eaten_by_the_multiplicity_digits():
    """``+62.9959`` is one 62.9959 Da gain, not sixty-two of something."""

    adduct = A.parse_adduct_label("CHRO_LIB:m-H+62.9959")
    assert adduct.numeric_terms == (62.9959,)
    assert adduct.gain_mass == pytest.approx(62.9959)


def test_sodium_adduct_uses_the_neutral_atom_mass():
    adduct = A.parse_adduct_label("CHRO_LIB:m-2H+Na")
    assert adduct.gain_mass == pytest.approx(A.ELEMENT_MASS["Na"])
    assert adduct.loss_mass == pytest.approx(2 * A.ELEMENT_MASS["H"])


def test_rows_that_state_no_ion_are_excluded_rather_than_failed():
    assert A.is_derivationless_label("CHRO_LIB:observed_consensus")
    assert A.is_derivationless_label("feature_neutral_mass_fallback")
    assert A.is_derivationless_label("")
    assert not A.is_derivationless_label("CHRO_LIB:m-H")


def test_the_uppercase_m_convention_is_reported_unparsed_not_guessed():
    """``M+21`` is a different labelling convention; guessing at it would invent data."""

    assert A.parse_adduct_label("CHRO_LIB:M+21") is None
    assert A.label_family(None, "CHRO_LIB:M+21") == "unparsed"


def test_unknown_elements_do_not_parse():
    assert A.token_mass("Xx") is None
    assert A.parse_adduct_label("CHRO_LIB:m+Xx") is None


# --- neutral-mass arithmetic -------------------------------------------------


def test_deprotonation_inverts_to_the_proton_convention():
    adduct = A.parse_adduct_label("m-H")
    assert A.neutral_mass_from_adduct(100.0, adduct, "negative") == pytest.approx(100.0 + PROTON)


def test_protonation_inverts_the_other_way():
    adduct = A.parse_adduct_label("m+H")
    assert A.neutral_mass_from_adduct(100.0, adduct, "positive") == pytest.approx(100.0 - PROTON)


def test_a_dimer_is_half_the_mass_a_monomer_would_be():
    adduct = A.parse_adduct_label("2m-H")
    monomer = A.naive_neutral_mass(153.1034, "negative")
    assert A.neutral_mass_from_adduct(153.1034, adduct, "negative") == pytest.approx(monomer / 2.0)


def test_a_doubly_charged_anion_is_twice_the_mass_a_singly_charged_one_would_be():
    """The X-21470 case: a disulfate read as z=1 lands at exactly half its mass."""

    adduct = A.parse_adduct_label("(m-2H)/2")
    singly = A.naive_neutral_mass(201.11318, "negative")
    assert A.neutral_mass_from_adduct(201.11318, adduct, "negative") == pytest.approx(2.0 * singly)


def test_an_isotopologue_is_lighter_than_the_peak_it_was_measured_at():
    adduct = A.parse_adduct_label("m-H[13C1]")
    plain = A.naive_neutral_mass(100.04035, "negative")
    shifted = A.neutral_mass_from_adduct(100.04035, adduct, "negative")
    assert plain - shifted == pytest.approx(A.ISOTOPE_DELTA["13C"])
    assert A.ISOTOPE_DELTA["13C"] == pytest.approx(1.00336, abs=1e-5)


def test_the_heavy_halogen_shifts_are_two_daltons_not_one():
    """37Cl and 81Br sit at M+2, which is what makes them visible at all."""

    assert A.ISOTOPE_DELTA["37Cl"] == pytest.approx(1.99705, abs=1e-5)
    assert A.ISOTOPE_DELTA["81Br"] == pytest.approx(1.99795, abs=1e-5)


def test_sodiated_loss_adduct_lands_where_hand_calculation_puts_it():
    adduct = A.parse_adduct_label("m-2H+Na")
    got = A.neutral_mass_from_adduct(158.08227, adduct, "negative")
    assert got == pytest.approx(137.10760, abs=1e-5)


@pytest.mark.parametrize(
    "label,mode",
    [
        ("m-H", "negative"),
        ("m+H", "positive"),
        ("2m-H", "negative"),
        ("(m-2H)/2", "negative"),
        ("m-H[13C1]", "negative"),
        ("(m-2H[13C1][18O1])/2", "negative"),
        ("m-2H+Na", "negative"),
        ("2m+H+NH3", "positive"),
    ],
)
def test_the_inversion_round_trips_through_the_forward_equation(label, mode):
    """Inverting m/z to a mass and pushing it back must return the same m/z."""

    adduct = A.parse_adduct_label(label)
    assert adduct is not None
    mz = 321.0987
    neutral = A.neutral_mass_from_adduct(mz, adduct, mode)
    sign = 1.0 if mode == "positive" else -1.0
    forward = (
        adduct.multiplier * neutral
        + adduct.gain_mass
        - adduct.loss_mass
        + adduct.isotope_shift
        - sign * adduct.charge * A.ELECTRON_MASS
    ) / adduct.charge
    assert forward == pytest.approx(mz, abs=1e-9)


# --- families and hints ------------------------------------------------------


@pytest.mark.parametrize(
    "label,expected",
    [
        ("CHRO_LIB:m-H", "protonated"),
        ("CHRO_LIB:m+H", "protonated"),
        ("CHRO_LIB:m-H[13C1]", "isotopologue"),
        ("CHRO_LIB:2m-H", "multimer"),
        ("CHRO_LIB:(m-2H)/2", "multiply_charged"),
        ("CHRO_LIB:m-2H+Na", "adduct_or_cluster"),
        ("CHRO_LIB:m-H+62.9959", "numeric_offset"),
        ("CHRO_LIB:observed_consensus", "consensus_or_fallback"),
    ],
)
def test_families_sort_by_the_term_whose_error_is_largest(label, expected):
    assert A.label_family(A.parse_adduct_label(label), label) == expected


def test_a_doubly_charged_isotopologue_is_filed_under_charge():
    label = "CHRO_LIB:(m-2H[13C1])/2"
    assert A.label_family(A.parse_adduct_label(label), label) == "multiply_charged"


def test_hints_are_merged_from_the_consensus_note_as_well_as_the_label():
    note = "consensus_observed_adducts=m-H,m-H[13C2],m-H[37Cl1]"
    assert A.merge_element_hints("CHRO_LIB:m-H[13C1]", note) == {"C": 2, "Cl": 1}


def test_hint_merging_keeps_the_largest_count_per_element():
    assert A.merge_element_hints("m-H[13C1]", "m-H[13C3]") == {"C": 3}


def test_unrecognised_isotope_tokens_are_dropped_rather_than_mapped():
    assert A.merge_element_hints("m-H[99Zz1]") == {}


# --- end to end over the committed data --------------------------------------


def test_the_audit_runs_over_the_committed_feature_rows():
    paths = sorted((ROOT / "data" / "feature_rows").glob("*_neutral_mass_hypotheses.csv"))
    assert paths, "committed neutral-mass hypotheses are missing"
    result = A.audit(paths)
    totals = result["totals"]
    assert totals["hypothesis_rows"] > 0

    arithmetic = result["adduct_arithmetic"]
    # The finding this audit exists to state: the ladder collapses onto [M+/-H].
    assert arithmetic["share_following_plain_protonation"] > 0.99

    # And the one family that is right is the one the pipeline actually applies.
    protonated = arithmetic["by_family"]["protonated"]
    assert protonated["vs_label"].get("agrees_within_1mda") == protonated["rows"]

    # Halogen evidence exists and is small, which is why it has to be used rather
    # than sampled: losing it loses the only chlorine signal in the pipeline.
    isotopes = result["isotope_evidence"]
    assert isotopes["features_with_halogen_evidence"] > 0
    assert isotopes["features_with_isotope_evidence"] > isotopes["features_with_halogen_evidence"]


def test_the_committed_artifact_matches_a_fresh_run():
    """The JSON in docs/ has to follow from the CSVs in data/, or it is decoration."""

    import json

    paths = sorted((ROOT / "data" / "feature_rows").glob("*_neutral_mass_hypotheses.csv"))
    fresh = A.audit(paths)
    committed = json.loads(
        (ROOT / "docs" / "paper" / "data" / "identification" / "audit_identification_rules.json").read_text()
    )
    assert fresh == committed
