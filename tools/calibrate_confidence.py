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

Usage:
    python3 tools/calibrate_confidence.py <reports_root>
"""
from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "structural_elucidation" / "reports"))
import confidence as C  # noqa: E402

TRUTH_CORRECT = 0.70


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


def fit_platt(raw: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    x = (raw - raw.mean()) / (raw.std() + 1e-9)
    a, b, lr = 1.0, 0.0, 0.5
    for _ in range(8000):
        p = 1.0 / (1.0 + np.exp(-(a * x + b)))
        a -= lr * np.mean((p - y) * x)
        b -= lr * np.mean(p - y)
    mu, sd = raw.mean(), raw.std() + 1e-9
    return float(a / sd), float(b - a * mu / sd)


def main() -> None:
    root = Path(sys.argv[1] if len(sys.argv) > 1 else "_audit_pull/reports")
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
    n = len(recs)
    print(f"[calib] reports: {n}", file=sys.stderr)
    for axis in ("mass", "ms2"):
        pass
    lab = [r for r in recs if r["labeled"]]
    n_exact = sum(1 for r in lab if r["source"] == "exact")
    print(f"[calib] truth-verified (truth>0): {len(lab)} "
          f"(exact {n_exact}, predicted {len(lab) - n_exact})", file=sys.stderr)

    raw = np.array([r["raw"] for r in lab])
    y = np.array([float(r["correct"]) for r in lab])
    A, B = fit_platt(raw, y)
    print(f"[calib] base correctness (truth>0): {y.mean():.3f}", file=sys.stderr)

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

    print(json.dumps({
        "weights": C.WEIGHTS,
        "platt_A": round(A, 6),
        "platt_B": round(B, 6),
        "n_calibration": len(lab),
        "base_rate": round(float(y.mean()), 4),
    }, indent=2))


if __name__ == "__main__":
    main()
