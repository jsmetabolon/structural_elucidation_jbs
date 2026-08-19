"""Tests for the Platt calibration harness.

The load-bearing one is
:func:`test_committed_controls_reproduce_the_shipped_constants`. It pins the fact
that ``docs/paper/data/structure/audit_structure_controls.csv`` carries the exact
inputs the shipped ``PLATT_A`` / ``PLATT_B`` were fitted from, which is what lets
the calibration be re-run and re-measured without a report bundle.
"""
from __future__ import annotations

import csv
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


M = load_tool("calibrate_confidence")
C = M.C

CONTROLS = ROOT / "docs" / "paper" / "data" / "structure" / "audit_structure_controls.csv"


def controls():
    rows = list(csv.DictReader(CONTROLS.open(newline="")))
    raw = [float(r["raw_full"]) for r in rows]
    y = [1.0 if str(r["top1_correct_0p70"]).strip().lower() in {"true", "1"} else 0.0 for r in rows]
    return rows, raw, y


# --- the fit reproduces what shipped ------------------------------------------


def test_committed_controls_reproduce_the_shipped_constants():
    """The committed controls are the calibration input, so no bundle is needed."""
    _, raw, y = controls()
    fit = M.fit_platt(raw, y)
    assert abs(fit.a - C.PLATT_A) < 1e-4, f"{fit.a} vs shipped {C.PLATT_A}"
    assert abs(fit.b - C.PLATT_B) < 1e-4, f"{fit.b} vs shipped {C.PLATT_B}"
    assert fit.n == C.CALIBRATION_N
    assert abs(fit.base_rate - C.CALIBRATION_BASE_RATE) < 1e-3


def test_the_pooled_fit_reproduces_the_unreachable_ceiling():
    """The reachable range is what makes TIER_HIGH unreachable; report it, do not infer it."""
    _, raw, y = controls()
    fit = M.fit_platt(raw, y)
    assert fit.reachable[1] < C.TIER_HIGH
    assert abs(fit.reachable[1] - 0.5907) < 1e-3


def test_the_fit_converges_well_inside_the_iteration_cap():
    _, raw, y = controls()
    fit = M.fit_platt(raw, y)
    assert fit.converged is True
    assert fit.iterations < M.FIT_ITERATIONS
    assert fit.max_gradient < M.FIT_TOLERANCE


# --- degenerate strata --------------------------------------------------------


def test_a_single_class_stratum_is_reported_not_fitted():
    """Logistic regression has no finite solution here; do not return the iteration cap."""
    fit = M.fit_platt([0.1, 0.5, 0.9], [1.0, 1.0, 1.0])
    assert fit.degenerate is True
    assert fit.converged is False
    assert "same class" in fit.note


def test_a_constant_score_is_reported_not_fitted():
    fit = M.fit_platt([0.5, 0.5, 0.5], [1.0, 0.0, 1.0])
    assert fit.degenerate is True
    assert "constant" in fit.note


def test_an_empty_stratum_is_reported_not_fitted():
    fit = M.fit_platt([], [])
    assert fit.degenerate is True
    assert fit.n == 0


# --- per-source ---------------------------------------------------------------


def records(*pairs):
    return [{"raw": raw, "correct": correct, "source": source} for raw, correct, source in pairs]


def test_sources_are_discovered_rather_than_hardcoded():
    fits = M.fit_per_source(records(
        (0.2, 0, "exact"), (0.8, 1, "exact"),
        (0.3, 0, "predicted"), (0.7, 1, "predicted"),
        (0.4, 0, "novel_source"), (0.6, 1, "novel_source"),
    ))
    assert set(fits) == {"exact", "predicted", "novel_source"}


def test_a_blank_source_does_not_become_a_stratum():
    fits = M.fit_per_source(records((0.2, 0, "exact"), (0.8, 1, "exact"), (0.5, 1, "")))
    assert set(fits) == {"exact"}


def test_per_source_fits_the_committed_controls_into_separated_ranges():
    """The result that motivates splitting: the two sources calibrate differently.

    Catalog hits are 38.1% correct and de novo predictions 56.0%, so a shared fit
    has to compromise between them. Fitted apart, the catalog stratum caps well
    below the de novo one instead.
    """
    rows, _, _ = controls()
    labeled = [
        {
            "raw": float(r["raw_full"]),
            "correct": 1 if str(r["top1_correct_0p70"]).strip().lower() in {"true", "1"} else 0,
            "source": r["source"].strip(),
        }
        for r in rows
    ]
    fits = M.fit_per_source(labeled)
    assert set(fits) == {"exact", "predicted"}
    assert all(not f.degenerate for f in fits.values())
    assert fits["exact"].base_rate < fits["predicted"].base_rate
    # The catalog stratum cannot reach where the de novo one can.
    assert fits["exact"].reachable[1] < fits["predicted"].reachable[1]


def test_fit_serializes_with_the_context_needed_to_judge_it():
    fit = M.fit_platt([0.2, 0.4, 0.6, 0.8], [0.0, 0.0, 1.0, 1.0])
    blob = fit.as_dict()
    for key in ("platt_A", "platt_B", "n", "base_rate", "converged", "iterations",
                "max_gradient", "degenerate", "reachable_score_range"):
        assert key in blob, key
    assert len(blob["reachable_score_range"]) == 2
