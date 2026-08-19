#!/usr/bin/env python3
"""Derive + validate calibration constants for the hybrid confidence model.

Uses the shared logic in ``structural_elucidation/reports/confidence.py`` so the
fit matches exactly what the renderer/patcher compute. For each report it builds
the confidence ``fields`` from sidecars, reads the rank-1 ``truth_tanimoto`` label
from ``<slug>_top_candidates.csv``, fits a 1-D Platt scaling of the
completeness-weighted raw evidence score to P(correct), and validates the rubric
tiers against observed correctness.

Ground truth: ``truth_tanimoto`` == 0.0 / empty encodes "no reliable truth" (a
zero fingerprint Tanimoto is chemically implausible for a real comparison), so
only ``truth_tanimoto`` > 0 is a genuine validation point. Correct == truth >= 0.70.

``--per-source`` fits one Platt scaling per candidate source instead of one over
the pool. That is the intervention ``tools/audit_thresholds.py`` recommends when
its ``source_stratification`` block reports ``confounded``, which it currently
does: exact catalog hits and de novo predictions have opposite base rates (0.381
against 0.560), and pooling them produces a rank AUC *below* either stratum's, so
the pooled score orders provenance rather than correctness. Fitting the pool is
what produced the constants in ``confidence.py``.

Usage:
    python3 tools/calibrate_confidence.py <reports_root>
    python3 tools/calibrate_confidence.py <reports_root> --per-source
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "structural_elucidation" / "reports"))
import confidence as C  # noqa: E402

TRUTH_CORRECT = 0.70

# Gradient descent settings. The original ran a fixed 8000 iterations at lr 0.5
# and reported nothing, which is survivable when the fit has 92 points and is
# checked by hand. A per-source fit has 42 or 50, so whether it actually
# converged has to be part of the output rather than an assumption.
FIT_ITERATIONS = 8000
FIT_LEARNING_RATE = 0.5
FIT_TOLERANCE = 1e-8


def _f(v):
    try:
        x = float(v)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(x) else x


def _rank1_truth(path: Path):
    if not path.exists():
        return None
    best, best_rank = None, None
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            r = _f(row.get("rank"))
            r = 1e9 if r is None else r
            if best_rank is None or r < best_rank:
                best_rank, best = r, row
    return _f(best.get("truth_tanimoto")) if best else None


@dataclass(frozen=True)
class PlattFit:
    """One Platt scaling, with enough context to judge whether to trust it."""

    a: float
    b: float
    n: int
    base_rate: float
    converged: bool
    iterations: int
    max_gradient: float
    degenerate: bool = False
    note: str = ""
    # The interval the fit can actually produce, for raw in [0, 1]. A tier cutoff
    # above the upper bound names a tier no report can enter, which is how
    # TIER_HIGH ended up unreachable; reporting it beside the fit means the next
    # person choosing a cutoff sees the ceiling rather than discovering it later.
    reachable: tuple[float, float] = field(default=(0.0, 1.0))

    def as_dict(self) -> dict:
        return {
            "platt_A": round(self.a, 6),
            "platt_B": round(self.b, 6),
            "n": self.n,
            "base_rate": round(self.base_rate, 4),
            "converged": self.converged,
            "iterations": self.iterations,
            "max_gradient": float(f"{self.max_gradient:.3e}"),
            "degenerate": self.degenerate,
            "note": self.note,
            "reachable_score_range": [round(self.reachable[0], 4), round(self.reachable[1], 4)],
        }


def _sigmoid(z: float) -> float:
    if z < -60:
        return 0.0
    if z > 60:
        return 1.0
    return 1.0 / (1.0 + math.exp(-z))


def fit_platt(raw: np.ndarray, y: np.ndarray) -> PlattFit:
    """Fit ``P(correct) = sigmoid(A * raw + B)`` by gradient descent.

    Returns a :class:`PlattFit` rather than a bare pair so the caller can tell a
    converged fit from one that ran out of iterations, and a real fit from a
    stratum that cannot be fitted at all.

    A stratum whose labels are all the same class has no finite maximum-likelihood
    solution - the logistic keeps pushing the slope outward - so it is reported as
    degenerate instead of returning whatever the iteration cap happened to reach.
    """
    raw = np.asarray(raw, dtype=float)
    y = np.asarray(y, dtype=float)
    n = int(raw.size)
    base = float(y.mean()) if n else 0.0
    if n == 0:
        return PlattFit(0.0, 0.0, 0, 0.0, False, 0, float("nan"), True, "no labeled records")
    if y.min() == y.max():
        return PlattFit(
            0.0, 0.0, n, base, False, 0, float("nan"), True,
            f"all {n} labels are the same class; no finite fit exists",
        )
    if raw.std() == 0.0:
        return PlattFit(
            0.0, 0.0, n, base, False, 0, float("nan"), True,
            f"raw score is constant across all {n} records; nothing to order",
        )

    x = (raw - raw.mean()) / (raw.std() + 1e-9)
    a, b = 1.0, 0.0
    grad = float("nan")
    used = FIT_ITERATIONS
    converged = False
    for step in range(FIT_ITERATIONS):
        p = 1.0 / (1.0 + np.exp(-(a * x + b)))
        ga = float(np.mean((p - y) * x))
        gb = float(np.mean(p - y))
        grad = max(abs(ga), abs(gb))
        if grad < FIT_TOLERANCE:
            used, converged = step, True
            break
        a -= FIT_LEARNING_RATE * ga
        b -= FIT_LEARNING_RATE * gb
    mu, sd = float(raw.mean()), float(raw.std() + 1e-9)
    fit_a, fit_b = a / sd, b - a * mu / sd
    return PlattFit(
        a=float(fit_a),
        b=float(fit_b),
        n=n,
        base_rate=base,
        converged=converged,
        iterations=used,
        max_gradient=grad,
        reachable=(_sigmoid(fit_b), _sigmoid(fit_a + fit_b)),
    )


def collect_records(root: Path) -> list[dict]:
    recs = []
    for d in sorted(p for p in root.iterdir() if p.is_dir()):
        data = d / "data"
        if not data.is_dir():
            continue
        slug = d.name.lower()
        fields = C.read_fields_from_sidecars(data, slug)
        if not fields:
            continue
        axes = C.build_axes(fields)
        raw, completeness = C.raw_score(axes)
        conf = C.compute_confidence(fields)
        truth = _rank1_truth(data / f"{slug}_top_candidates.csv")
        recs.append({
            "raw": raw,
            "completeness": completeness,
            "tier": conf.tier,
            "source": fields["source"],
            "truth": truth,
            "labeled": truth is not None and truth > 0.0,
            "correct": 1 if (truth is not None and truth >= TRUTH_CORRECT) else 0,
        })
    return recs


def fit_records(records: list[dict]) -> PlattFit:
    """Fit one Platt scaling over the labeled records handed to it."""
    return fit_platt(
        np.array([r["raw"] for r in records], dtype=float),
        np.array([float(r["correct"]) for r in records], dtype=float),
    )


def fit_per_source(labeled: list[dict]) -> dict[str, PlattFit]:
    """One fit per candidate source, keyed by source name.

    Sources are discovered rather than hard-coded, so a bundle that grows a third
    provenance is fitted rather than silently folded into the pool.
    """
    sources = sorted({str(r.get("source") or "").strip() for r in labeled} - {""})
    return {
        source: fit_records([r for r in labeled if str(r.get("source") or "").strip() == source])
        for source in sources
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("reports_root", nargs="?", type=Path, default=Path("_audit_pull/reports"),
                        help="root of a delivered report bundle (default: _audit_pull/reports)")
    parser.add_argument("--per-source", action="store_true",
                        help="also fit one Platt scaling per candidate source, not just the pool")
    parser.add_argument("--out", type=Path, default=None,
                        help="write the JSON result here as well as to stdout")
    args = parser.parse_args(argv)

    root = args.reports_root
    if not root.is_dir():
        parser.error(f"reports root not found: {root}")
    recs = collect_records(root)
    print(f"[calib] reports: {len(recs)}", file=sys.stderr)
    lab = [r for r in recs if r["labeled"]]
    if not lab:
        parser.error(f"no truth-verified reports under {root}")
    n_exact = sum(1 for r in lab if r["source"] == "exact")
    print(f"[calib] truth-verified (truth>0): {len(lab)} "
          f"(exact {n_exact}, predicted {len(lab) - n_exact})", file=sys.stderr)

    pooled = fit_records(lab)
    A, B = pooled.a, pooled.b
    print(f"[calib] base correctness (truth>0): {pooled.base_rate:.3f}", file=sys.stderr)
    if not pooled.converged:
        print(f"[calib] WARNING pooled fit did not converge "
              f"(max gradient {pooled.max_gradient:.2e} after {pooled.iterations} iterations)",
              file=sys.stderr)
    print(f"[calib] pooled reachable score range: "
          f"[{pooled.reachable[0]:.4f}, {pooled.reachable[1]:.4f}]", file=sys.stderr)

    def prob(rw: float) -> float:
        return 1.0 / (1.0 + math.exp(-(A * rw + B)))

    order = sorted(lab, key=lambda r: prob(r["raw"]))
    k = max(1, len(order) // 4)
    print("[calib] quartile bins by calibrated prob:", file=sys.stderr)
    for i in range(0, len(order), k):
        grp = order[i:i + k]
        if grp:
            ps = [prob(r["raw"]) for r in grp]
            acc = sum(r["correct"] for r in grp) / len(grp)
            print(f"   prob[{min(ps):.2f}-{max(ps):.2f}] n={len(grp)} correct={acc:.2f}", file=sys.stderr)

    print("[calib] rubric tier validation (truth-verified):", file=sys.stderr)
    for tier in ("High", "Moderate", "Low", "Tentative"):
        grp = [r for r in lab if r["tier"] == tier]
        if grp:
            acc = sum(r["correct"] for r in grp) / len(grp)
            print(f"   {tier}: n={len(grp)} observed_correct={acc:.2f}", file=sys.stderr)
    # Full-set tier distribution
    dist = {}
    for r in recs:
        dist[r["tier"]] = dist.get(r["tier"], 0) + 1
    print(f"[calib] full-set tier distribution: {dist}", file=sys.stderr)

    # The pooled keys keep their original names and shape so anything reading
    # this output by hand still finds what it expects.
    result: dict = {
        "weights": C.WEIGHTS,
        "platt_A": round(A, 6),
        "platt_B": round(B, 6),
        "n_calibration": len(lab),
        "base_rate": round(pooled.base_rate, 4),
        "pooled": pooled.as_dict(),
    }

    if args.per_source:
        per_source = fit_per_source(lab)
        result["per_source"] = {name: fit.as_dict() for name, fit in per_source.items()}
        print("[calib] per-source fits:", file=sys.stderr)
        for name, fit in per_source.items():
            if fit.degenerate:
                print(f"   {name}: NOT FITTED - {fit.note}", file=sys.stderr)
                continue
            flag = "" if fit.converged else "  (DID NOT CONVERGE)"
            print(f"   {name}: n={fit.n} base={fit.base_rate:.3f} "
                  f"A={fit.a:.4f} B={fit.b:.4f} "
                  f"reachable=[{fit.reachable[0]:.4f}, {fit.reachable[1]:.4f}]{flag}",
                  file=sys.stderr)
        # The whole point of splitting is that the strata disagree. If they do
        # not, say so, because then pooling was never the problem.
        spreads = [f.base_rate for f in per_source.values() if not f.degenerate]
        if len(spreads) > 1:
            result["base_rate_spread"] = round(max(spreads) - min(spreads), 4)
            print(f"[calib] base-rate spread across sources: {result['base_rate_spread']}",
                  file=sys.stderr)

    payload = json.dumps(result, indent=2, sort_keys=True)
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(payload + "\n")
    print(payload)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
