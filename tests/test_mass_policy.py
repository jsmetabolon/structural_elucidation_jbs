"""Tests for the mass-policy audit harness.

The tool's job is to decide whether a scoring change could move a result, so the
pieces that must be right are the ones a wrong answer would hide behind: the ppm
conversion, where the band edges fall, and the paired test that says whether a
difference is distinguishable from noise on 92 controls.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

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


M = load_tool("audit_mass_policy")


# --- ppm conversion ----------------------------------------------------------


def test_ppm_of_a_milli_dalton_at_a_thousand_daltons_is_one():
    assert M.ppm_of(0.001, 1000.0) == 1.0


def test_ppm_scales_inversely_with_mass():
    """The same absolute error is a larger relative error on a lighter compound."""

    light = M.ppm_of(0.001, 250.0)
    heavy = M.ppm_of(0.001, 500.0)
    assert light == 4.0
    assert heavy == 2.0


def test_ppm_is_undefined_without_a_mass_reference():
    assert M.ppm_of(0.001, None) is None
    assert M.ppm_of(0.001, 0.0) is None


# --- policy scoring ----------------------------------------------------------


def test_absolute_policy_reproduces_the_shipped_windows():
    policy = M.CURRENT
    assert policy.score(0.0005, 300.0) == 1.0        # inside the plateau
    assert policy.score(0.001, 300.0) == 1.0         # on the plateau edge
    assert policy.score(0.005, 300.0) == 0.0         # at the tolerance
    assert policy.score(0.010, 300.0) == 0.0         # past it
    assert policy.score(0.0025, 300.0) == 0.5        # halfway down the linear arm


def test_the_absolute_policy_cannot_return_a_score_between_0_8_and_1_0():
    """The plateau and the linear arm do not meet, so that interval is unreachable.

    It matters because MASS_STRONG sits at 0.90, inside the gap, which makes the
    strong-mass condition an exact test for 1.00 rather than the 0.90 bar it reads
    as.
    """

    just_past_plateau = M.CURRENT.score(0.001 + 1e-9, 300.0)
    assert abs(just_past_plateau - 0.8) < 1e-6
    reachable = {M.CURRENT.score(err / 1e6, 300.0) for err in range(0, 6000, 7)}
    assert not any(0.8 < score < 1.0 for score in reachable)


def test_ppm_policy_edges_track_the_compound_mass():
    policy = M.DRAFTED  # plateau 1 ppm, tolerance 5 ppm
    assert policy.edges_da(1000.0) == (0.001, 0.005)
    assert policy.edges_da(500.0) == (0.0005, 0.0025)
    # A 2 ppm error is inside the window at any mass, so it scores the same either way.
    assert policy.score(0.002, 1000.0) == policy.score(0.001, 500.0)


def test_ppm_policy_has_no_score_without_a_mass_reference():
    assert M.DRAFTED.score(0.001, None) is None


# --- bands -------------------------------------------------------------------


def test_bands_use_the_confidence_models_own_cutoffs():
    assert M.band_of(1.0) == M.STRONG
    assert M.band_of(0.90) == M.STRONG
    assert M.band_of(0.89) == M.PARTIAL
    assert M.band_of(0.50) == M.PARTIAL
    assert M.band_of(0.49) == M.FAIL
    assert M.band_of(None) is None


# --- exact McNemar -----------------------------------------------------------


def test_mcnemar_is_one_when_nothing_is_discordant():
    assert M.mcnemar_exact(0, 0) == 1.0


def test_mcnemar_is_one_when_discordance_is_balanced():
    assert M.mcnemar_exact(5, 5) == 1.0


def test_mcnemar_matches_hand_computed_binomial_tails():
    # b=0, c=1: two-sided exact = 2 * P(X <= 0) with X ~ Binom(1, 1/2) = 2 * 1/2 = 1.
    assert M.mcnemar_exact(0, 1) == 1.0
    # b=0, c=5: 2 * (1/32) = 1/16.
    assert abs(M.mcnemar_exact(0, 5) - 0.0625) < 1e-12
    # b=1, c=9: 2 * (1 + 10)/1024 = 22/1024.
    assert abs(M.mcnemar_exact(1, 9) - 22 / 1024) < 1e-12


def test_mcnemar_is_symmetric_in_its_arguments():
    """Gains and losses swap without changing how surprising the split is."""

    assert M.mcnemar_exact(2, 11) == M.mcnemar_exact(11, 2)


def test_mcnemar_never_exceeds_one():
    assert M.mcnemar_exact(3, 4) <= 1.0
    assert M.mcnemar_exact(50, 50) <= 1.0


def test_mcnemar_resolves_a_shift_too_small_for_independent_intervals():
    """Eight flips one way out of ten is significant; the CI comparison is not.

    This is the reason the plan insists the comparison be paired: rank-1 recovery
    on 92 controls has a twenty-point interval, so two independent intervals
    cannot see a change this size at all.
    """

    assert M.mcnemar_exact(1, 9) < 0.05


# --- the gate ----------------------------------------------------------------


def candidate(rank: int, error_da: float | None, truth: float = 0.0) -> "M.Candidate":
    return M.Candidate(rank=rank, smiles=f"C{rank}", error_da=error_da, stored_score=None,
                       truth_tanimoto=truth)


def report_of(*candidates, neutral_mass: float = 300.0, control_id: str = "X-1",
              msi_level: int | None = None, candidate_exact_mass: float | None = None,
              observed_neutral_mass: float | None = None) -> "M.Report":
    return M.Report(control_id=control_id, neutral_mass=neutral_mass, msi_level=msi_level,
                    source="predicted", candidates=list(candidates),
                    candidate_exact_mass=candidate_exact_mass,
                    observed_neutral_mass=observed_neutral_mass)


def test_the_gate_demotes_a_mass_failing_candidate_below_a_passing_one():
    report = report_of(candidate(1, 0.010), candidate(2, 0.0005))
    assert [c.rank for c in M.gated_order(report, M.CURRENT)] == [2, 1]


def test_the_gate_is_stable_where_nothing_fails():
    report = report_of(candidate(1, 0.0005), candidate(2, 0.0008), candidate(3, 0.0002))
    assert [c.rank for c in M.gated_order(report, M.CURRENT)] == [1, 2, 3]


def test_the_gate_preserves_order_within_the_failing_group():
    report = report_of(candidate(1, 0.02), candidate(2, 0.03), candidate(3, 0.0001))
    assert [c.rank for c in M.gated_order(report, M.CURRENT)] == [3, 1, 2]


def test_a_candidate_with_no_mass_error_is_never_demoted():
    """Missing evidence is not failing evidence, which is how the axes treat it."""

    report = report_of(candidate(1, None), candidate(2, 0.0001))
    assert [c.rank for c in M.gated_order(report, M.CURRENT)] == [1, 2]


def test_recovery_is_scored_at_the_truth_threshold():
    report = report_of(candidate(1, 0.0001, truth=0.69), candidate(2, 0.0001, truth=0.71))
    assert report.correct_at(1) is False
    assert report.correct_at(3) is True


# --- mass accuracy split -----------------------------------------------------


def test_accuracy_split_separates_correct_from_incorrect():
    """A tight window should pass the accurate candidate and not the wild one."""

    report = report_of(
        candidate(1, 0.0003, truth=0.95),   # 1 ppm at 300 Da, correct
        candidate(2, 0.0300, truth=0.10),   # 100 ppm, incorrect
    )
    result = M.mass_accuracy_by_correctness([report])
    assert result["correct"]["n"] == 1
    assert result["incorrect"]["n"] == 1
    within_2 = next(w for w in result["windows"] if w["window_ppm"] == 2.0)
    assert within_2["share_of_correct_within"] == 1.0
    assert within_2["share_of_incorrect_within"] == 0.0


def test_accuracy_split_ignores_reports_without_truth():
    """Without a known answer there is nothing to split on."""

    report = report_of(candidate(1, 0.0003, truth=0.0))
    assert M.mass_accuracy_by_correctness([report])["correct"]["n"] == 0


def test_a_window_wide_enough_stops_separating():
    report = report_of(
        candidate(1, 0.0003, truth=0.95),
        candidate(2, 0.0300, truth=0.10),
    )
    widest = M.mass_accuracy_by_correctness([report])["windows"][-1]
    assert widest["share_of_correct_within"] == widest["share_of_incorrect_within"] == 1.0


# --- the conjugate tail ------------------------------------------------------


def test_conjugate_tail_names_a_glucuronide_shifted_candidate():
    """Correct by Tanimoto, wrong by exactly one glucuronide."""

    report = report_of(candidate(1, 176.032, truth=0.80))
    tail = M.conjugate_tail([report])
    assert tail["n_with_implausible_mass"] == 1
    assert tail["by_conjugate"] == {"glucuronide": 1}
    assert tail["n_distinct_controls"] == 1


def test_conjugate_tail_leaves_an_accurate_correct_candidate_alone():
    report = report_of(candidate(1, 0.0003, truth=0.80))
    tail = M.conjugate_tail([report])
    assert tail["n_correct_candidates"] == 1
    assert tail["n_with_implausible_mass"] == 0


def test_conjugate_tail_flags_an_identical_structure_that_cannot_be_a_conjugate():
    """A perfect structural match hundreds of Daltons out is not a conjugate.

    Identical fingerprints cannot differ by a conjugate, so folding it into the
    conjugate count would explain away the wrong thing. In the delivered data the
    real cause was charge state, which the next section tests for.
    """

    report = report_of(candidate(1, 233.066, truth=1.0))
    tail = M.conjugate_tail([report])
    assert tail["identical_structure_but_mass_mismatched"] == ["X-1"]
    assert tail["by_conjugate"] == {"unmatched": 1}


def test_conjugate_tail_reports_modal_masses():
    report = report_of(
        candidate(1, 176.032, truth=0.80),
        candidate(2, 176.032, truth=0.75),
        candidate(3, 162.053, truth=0.72),
    )
    modal = M.conjugate_tail([report])["modal_error_masses_da"]
    assert list(modal.items())[0] == ("176.0320", 2)


# --- charge state ------------------------------------------------------------


def test_a_doubly_charged_precursor_read_as_singly_charged_is_flagged():
    """The X-21470 case: a disulfate whose candidate mass is exactly twice observed."""

    report = report_of(candidate(1, 233.066, truth=1.0), control_id="X-21470", msi_level=1,
                       observed_neutral_mass=233.0668, candidate_exact_mass=466.1331)
    found = M.charge_state_suspects([report])
    assert found["n_suspects"] == 1
    assert found["suspects"][0]["ratio"] == "2x"
    assert found["suspects"][0]["control_id"] == "X-21470"


def test_a_singly_charged_feature_is_not_flagged():
    """The guard that keeps this from being noise: an ordinary report must pass."""

    report = report_of(candidate(1, 0.0003, truth=0.9),
                       observed_neutral_mass=300.0000, candidate_exact_mass=300.0003)
    assert M.charge_state_suspects([report])["n_suspects"] == 0


def test_a_merely_heavy_candidate_is_not_flagged():
    """Being far off is not the signal. Being off by a clean integer factor is.

    A candidate at 2.05x is further from the observed mass than one at exactly 2x,
    and is the ordinary bad match the mass window already handles.
    """

    report = report_of(candidate(1, 315.0, truth=0.0),
                       observed_neutral_mass=300.0, candidate_exact_mass=615.0)
    assert M.charge_state_suspects([report])["n_suspects"] == 0


def test_the_reciprocal_direction_is_flagged_too():
    report = report_of(candidate(1, 306.0, truth=0.0), control_id="X-25107", msi_level=4,
                       observed_neutral_mass=612.0505, candidate_exact_mass=306.0253)
    found = M.charge_state_suspects([report])
    assert found["suspects"][0]["ratio"] == "1/2x"


def test_a_triply_charged_precursor_is_flagged():
    report = report_of(candidate(1, 600.0, truth=0.0),
                       observed_neutral_mass=300.0, candidate_exact_mass=900.0)
    found = M.charge_state_suspects([report])
    assert found["suspects"][0]["ratio"] == "3x"


def test_the_msi_level_of_a_suspect_is_reported():
    """A clean factor of two in the Confirmed tier is the finding, so it must surface."""

    reports = [
        report_of(candidate(1, 233.0, truth=1.0), control_id="X-21470", msi_level=1,
                  observed_neutral_mass=233.0668, candidate_exact_mass=466.1331),
        report_of(candidate(1, 232.0, truth=0.0), control_id="X-24546", msi_level=3,
                  observed_neutral_mass=232.0581, candidate_exact_mass=464.1175),
    ]
    found = M.charge_state_suspects(reports)
    assert found["by_msi_level"] == {"1": 1, "3": 1}
    assert found["reaches_msi_level_1"] is True


def test_a_report_without_both_masses_is_not_considered():
    """The comparison needs an observed mass that is not the candidate's own."""

    report = report_of(candidate(1, 0.0003, truth=0.9), candidate_exact_mass=466.1331)
    found = M.charge_state_suspects([report])
    assert found["n_reports_with_both_masses"] == 0
    assert found["n_suspects"] == 0


# --- a strong claim the evidence does not support ----------------------------


def test_a_strong_score_against_a_huge_error_is_contradictory():
    report = report_of(
        M.Candidate(rank=1, smiles="C", error_da=233.0664, stored_score=0.95, truth_tanimoto=1.0),
        control_id="X-21470", msi_level=1,
    )
    found = M.strong_score_without_support([report])
    assert found["n_contradictory"] == 1
    assert found["contradictory"][0]["control_id"] == "X-21470"


def test_a_strong_score_against_a_tiny_error_is_not_contradictory():
    """Unreproducible is not the same as indefensible.

    A sub-milli-Dalton error with a high score is a fair reading of the evidence
    by some rule this audit has not catalogued, and lumping it in with a 233 Da
    contradiction would bury the one case that matters.
    """

    report = report_of(
        M.Candidate(rank=1, smiles="C", error_da=0.0012, stored_score=0.999, truth_tanimoto=0.9),
    )
    found = M.strong_score_without_support([report])
    assert found["n_contradictory"] == 0


def test_a_reproducible_score_is_not_flagged_at_all():
    """A score a known rule explains is not a disagreement, however it was made."""

    error = 0.002
    report = report_of(
        M.Candidate(rank=1, smiles="C", error_da=error, stored_score=1.0 - error / 1.0,
                    truth_tanimoto=0.9),
    )
    assert M.strong_score_without_support([report])["n"] == 0


# --- truth-definition sensitivity --------------------------------------------


def test_a_correct_hit_that_misses_the_mass_is_lost_under_the_strict_definition():
    """Tanimoto can clear the bar while the molecular formula is wrong.

    A candidate missing a glucuronide is structurally very similar to the truth
    and is 176 Da away from being the compound.
    """

    report = report_of(candidate(1, 176.032, truth=0.80))
    result = M.mass_plausible_recovery([report])
    top1 = result["depths"]["top_1"]
    assert top1["published_n"] == 1
    assert top1["mass_plausible_n"] == 0
    assert top1["controls_lost"] == ["X-1"]


def test_a_correct_hit_that_fits_the_mass_survives():
    report = report_of(candidate(1, 0.0003, truth=0.80))
    result = M.mass_plausible_recovery([report])
    assert result["depths"]["top_1"]["mass_plausible_n"] == 1
    assert result["depths"]["top_1"]["controls_lost"] == []


def test_a_candidate_with_no_recorded_error_cannot_be_shown_to_fit():
    """The strict variant is a bound, so absent evidence counts against it.

    This is deliberately the opposite of how the gate treats a missing error,
    where missing evidence is not failing evidence. Here the question is what
    recovery would be if a fit were required, and an unmeasurable candidate
    cannot demonstrate one.
    """

    report = report_of(candidate(1, None, truth=0.80))
    result = M.mass_plausible_recovery([report])
    assert result["depths"]["top_1"]["published_n"] == 1
    assert result["depths"]["top_1"]["mass_plausible_n"] == 0


def test_the_strict_definition_can_promote_a_deeper_candidate():
    """Depth still counts: if rank 1 fails the mass, a correct rank-2 can carry top-3."""

    report = report_of(candidate(1, 176.032, truth=0.80), candidate(2, 0.0003, truth=0.85))
    result = M.mass_plausible_recovery([report])
    assert result["depths"]["top_1"]["mass_plausible_n"] == 0
    assert result["depths"]["top_3"]["mass_plausible_n"] == 1


def test_reports_without_truth_are_not_counted_as_controls():
    reports = [
        report_of(candidate(1, 0.0003, truth=0.80)),
        report_of(candidate(1, 0.0003, truth=0.0), control_id="X-2"),
    ]
    assert M.mass_plausible_recovery(reports)["n_controls"] == 1

