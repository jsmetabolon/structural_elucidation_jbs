"""Tests for the canonical adduct ladder shared by phase0 and analog."""

from __future__ import annotations

import pandas as pd

from structural_elucidation.chemistry import adduct_ladder
from structural_elucidation.evaluation import analog
from structural_elucidation.utils import phase0

PROTON_MASS = 1.007276
MZ = 200.0


def test_default_rows_are_first() -> None:
    pos = adduct_ladder.hypotheses(MZ, "positive")
    neg = adduct_ladder.hypotheses(MZ, "negative")
    assert pos[0][0] == "[M+H]+"
    assert neg[0][0] == "[M-H]-"


def test_offset_rows_invert_as_mz_minus_mz_minus_M() -> None:
    for row in adduct_ladder.load_ladder():
        if row.formula_kind != "offset":
            continue
        hyps = dict(adduct_ladder.hypotheses(MZ, row.ion_mode))
        implied = hyps[row.adduct_label]
        assert abs((MZ - row.mz_minus_M) - implied) < 1e-9


def test_legacy_dimer_and_z2_match_phase0_closed_form() -> None:
    pos = dict(adduct_ladder.hypotheses(MZ, "positive"))
    neg = dict(adduct_ladder.hypotheses(MZ, "negative"))
    assert abs(pos["[M+2H]2+"] - (2.0 * MZ - 2.0 * PROTON_MASS)) < 1e-12
    assert abs(pos["[2M+H]+"] - ((MZ - PROTON_MASS) / 2.0)) < 1e-12
    assert abs(neg["[M-2H]2-"] - (2.0 * MZ + 2.0 * PROTON_MASS)) < 1e-12
    assert abs(neg["[2M-H]-"] - ((MZ + PROTON_MASS) / 2.0)) < 1e-12


def test_phase0_hypotheses_match_ladder() -> None:
    for mode in ("positive", "negative"):
        assert phase0._adduct_neutral_mass_hypotheses(
            MZ, mode
        ) == adduct_ladder.hypotheses(MZ, mode)


def test_analog_mass_match_best_hypothesis_is_on_the_ladder() -> None:
    for mode in ("positive", "negative"):
        labels = {label for label, _ in adduct_ladder.hypotheses(MZ, mode)}
        match = analog.mass_match(180.0, MZ, mode)
        assert match["hypothesis"] in labels
        assert match["adduct"] in labels


def test_acetate_is_not_on_the_ladder() -> None:
    labels = {row.adduct_label for row in adduct_ladder.load_ladder()}
    assert "[M+CH3COO]-" not in labels


def test_phase0_import_smoke() -> None:
    assert callable(phase0.best_adduct_neutral_mass_match)
    match = phase0.best_adduct_neutral_mass_match(
        MZ - PROTON_MASS, MZ, ion_mode="positive"
    )
    assert match["adduct"] == "[M+H]+"


def test_drop_polar_omits_cl_and_formate() -> None:
    all_neg = {label for label, _ in adduct_ladder.hypotheses(MZ, "negative")}
    filtered = {
        label
        for label, _ in adduct_ladder.hypotheses(MZ, "negative", include_polar=False)
    }
    assert "[M+Cl]-" in all_neg
    assert "[M+HCOO]-" in all_neg
    assert "[M+Cl]-" not in filtered
    assert "[M+HCOO]-" not in filtered
    assert "[M-H]-" in filtered
    assert adduct_ladder.include_polar_for_lib(209) is False
    assert adduct_ladder.include_polar_for_lib(400) is True


def test_derive_keeps_protonated_stored_and_inverts_sodium() -> None:
    stored = MZ + PROTON_MASS
    mass, source, shift = adduct_ladder.derive_target_neutral_mass(
        MZ, "negative", "m-H", stored
    )
    assert source == "stored"
    assert mass == stored
    assert shift == 0.0
    na_stored = stored
    mass_na, source_na, shift_na = adduct_ladder.derive_target_neutral_mass(
        MZ, "positive", "m+Na", na_stored
    )
    assert source_na == "label"
    row = adduct_ladder.row_for_observed_label("m+Na", "positive")
    assert row is not None
    expected = MZ - row.mz_minus_M
    assert abs(mass_na - expected) < 1e-9
    assert shift_na is not None and abs(shift_na) > 0.5


def test_annotate_adduct_does_not_change_mass_match_score() -> None:
    from structural_elucidation.reports.compound import annotate_adduct_hypotheses

    df = pd.DataFrame(
        [{"predicted_smiles": "CCO", "mass_match_score": 0.42}]
    )
    out = annotate_adduct_hypotheses(df, 47.049141, "positive")
    assert float(out.iloc[0]["mass_match_score"]) == 0.42
    assert "hypothesized_adduct" in out.columns
