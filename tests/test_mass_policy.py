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


def report_of(*candidates, neutral_mass: float = 300.0) -> "M.Report":
    return M.Report(control_id="X-1", neutral_mass=neutral_mass, msi_level=None,
                    source="predicted", candidates=list(candidates))


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
