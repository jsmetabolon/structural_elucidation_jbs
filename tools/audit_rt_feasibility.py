#!/usr/bin/env python3
"""Test whether retention behavior can rule out chromatographically infeasible candidates.

The proposal was to use logP to eliminate candidates whose predicted retention
disagrees with what was observed. Before wiring that into the reports it has to
earn its place, because a filter built on an uncalibrated retention model
removes correct answers as readily as wrong ones.

Two questions, in order:

1. **Is the model calibrated at all?** Correlate the logP-derived retention-index
   proxy in ``structural_elucidation/nce/evidence.py`` against the observed index
   for the controls whose true structure is known. Reported per library, because
   the four libraries are different chromatographic methods and a single fit
   across them would be meaningless.
2. **Does the residual separate right from wrong?** A candidate pool is
   mass-constrained, so the honest test is against mass-matched decoys: hold the
   observed index fixed and ask whether the true structure's residual is smaller
   than that of other real metabolites of similar mass on the same column. That
   is the discrimination an eliminating filter would actually rely on.

Recalibration is fitted leave-one-out so a control never helps predict itself.

The gate is stated up front: wire the axis in only at decoy AUC >= 0.60. Below
that the deliverable is this audit and a written negative result.

Note on inputs: the observed retention index is not a column anywhere. It is
embedded in the feature key, which does survive into the reranked output, so
``ri_from_feature_key`` recovers it without needing the full ion table.

Usage:
    python3 tools/audit_rt_feasibility.py
"""

from __future__ import annotations

import argparse
import csv
import glob
import json
import math
import re
import sys
from pathlib import Path
from typing import Any, Sequence

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from structural_elucidation.utils.domain import ri_from_feature_key  # noqa: E402

try:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Descriptors

    RDLogger.DisableLog("rdApp.*")
    HAS_RDKIT = True
except ImportError:  # pragma: no cover - optional dependency
    Chem = None
    Descriptors = None
    HAS_RDKIT = False

# Wire the axis in only if the residual separates true structures from
# mass-matched decoys at least this well. 0.5 is no information.
DECOY_AUC_GATE = 0.60
# Decoys must be plausible members of the same candidate pool, so they are drawn
# from the same library and within this much of the true neutral mass.
DECOY_MASS_TOLERANCE_DA = 50.0
# A library needs enough labeled controls for its fit to mean anything.
MIN_LIBRARY_CONTROLS = 10

# Per-library bar. Pooled discrimination can be carried by one method, and an
# axis is only defensible where it is working.
MIN_LIBRARY_AUC = 0.60

LIB_FROM_PATH = re.compile(r"lib(\d+)_(\w+?)_feature_rows\.csv$")


def _mean(values: Sequence[float]) -> float:
    return sum(values) / len(values) if values else math.nan


def pearson(xs: Sequence[float], ys: Sequence[float]) -> float | None:
    if len(xs) < 3:
        return None
    mx, my = _mean(xs), _mean(ys)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    if sxx <= 0 or syy <= 0:
        return None
    return sxy / math.sqrt(sxx * syy)


def fit_linear(xs: Sequence[float], ys: Sequence[float]) -> tuple[float, float] | None:
    """Ordinary least squares ``y = slope * x + intercept``."""
    if len(xs) < 3:
        return None
    mx, my = _mean(xs), _mean(ys)
    sxx = sum((x - mx) ** 2 for x in xs)
    if sxx <= 0:
        return None
    slope = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sxx
    return slope, my - slope * mx


def rmse(residuals: Sequence[float]) -> float | None:
    if not residuals:
        return None
    return math.sqrt(sum(value * value for value in residuals) / len(residuals))


def auc(scores: Sequence[float], labels: Sequence[bool]) -> float | None:
    """Rank AUC with ties at 0.5; 0.5 means the score carries no information."""
    n_pos = sum(1 for label in labels if label)
    n_neg = len(labels) - n_pos
    if not n_pos or not n_neg:
        return None
    order = sorted(range(len(scores)), key=lambda index: scores[index])
    ranks = [0.0] * len(scores)
    start = 0
    while start < len(order):
        stop = start
        while stop + 1 < len(order) and scores[order[stop + 1]] == scores[order[start]]:
            stop += 1
        average = (start + stop) / 2.0 + 1.0
        for index in range(start, stop + 1):
            ranks[order[index]] = average
        start = stop + 1
    positive_ranks = sum(rank for rank, label in zip(ranks, labels) if label)
    return (positive_ranks - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)


def load_truth_structures(path: Path) -> dict[str, str]:
    with path.open(newline="") as handle:
        return {
            str(row.get("x_compound_id") or "").strip(): str(row.get("smiles") or "").strip()
            for row in csv.DictReader(handle)
            if str(row.get("smiles") or "").strip()
        }


def load_observed_ri(feature_rows_dir: Path) -> list[dict[str, Any]]:
    """Recover observed retention indices from the feature keys, per library."""
    out: list[dict[str, Any]] = []
    for path in sorted(glob.glob(str(feature_rows_dir / "*_feature_rows.csv"))):
        match = LIB_FROM_PATH.search(path)
        if not match:
            continue
        lib_id, ion_mode = int(match.group(1)), match.group(2)
        with open(path, newline="") as handle:
            for row in csv.DictReader(handle):
                observed = ri_from_feature_key(row.get("feature_key"))
                if not math.isfinite(observed):
                    continue
                control_id = str(row.get("control_id") or "").strip()
                if not control_id:
                    continue
                out.append(
                    {
                        "control_id": control_id,
                        "lib_id": lib_id,
                        "ion_mode": ion_mode,
                        "observed_ri": observed,
                    }
                )
    return out


def build_records(truth: dict[str, str], observed: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Join observed retention indices to the controls whose structure is known."""
    from structural_elucidation.nce.evidence import ri_proxy_from_smiles

    records = []
    for row in observed:
        smiles = truth.get(row["control_id"])
        if not smiles:
            continue
        predicted = ri_proxy_from_smiles(smiles)
        if not math.isfinite(predicted):
            continue
        mol = Chem.MolFromSmiles(smiles) if HAS_RDKIT else None
        if mol is None:
            continue
        records.append(
            {
                **row,
                "smiles": smiles,
                "exact_mass": float(Descriptors.ExactMolWt(mol)),
                "proxy_ri": predicted,
            }
        )
    # One control can appear in several libraries; each is its own measurement.
    return sorted(records, key=lambda record: (record["lib_id"], record["control_id"]))


def calibration_by_library(records: list[dict[str, Any]]) -> dict[str, Any]:
    """How well the proxy tracks the observed index, per chromatographic method."""
    out: dict[str, Any] = {}
    for lib_id in sorted({record["lib_id"] for record in records}):
        subset = [record for record in records if record["lib_id"] == lib_id]
        observed = [record["observed_ri"] for record in subset]
        proxy = [record["proxy_ri"] for record in subset]
        fit = fit_linear(proxy, observed)
        recalibrated_residuals = []
        for index, record in enumerate(subset):
            others = subset[:index] + subset[index + 1 :]
            loo = fit_linear([r["proxy_ri"] for r in others], [r["observed_ri"] for r in others])
            if loo is None:
                continue
            slope, intercept = loo
            record["recalibrated_ri"] = slope * record["proxy_ri"] + intercept
            recalibrated_residuals.append(record["observed_ri"] - record["recalibrated_ri"])
        out[str(lib_id)] = {
            "n": len(subset),
            "ion_mode": subset[0]["ion_mode"],
            "observed_ri_range": [round(min(observed), 1), round(max(observed), 1)],
            "pearson_r": round(pearson(proxy, observed), 4)
            if pearson(proxy, observed) is not None
            else None,
            "raw_rmse": round(rmse([o - p for o, p in zip(observed, proxy)]) or 0.0, 1) or None,
            "recalibrated_loo_rmse": round(rmse(recalibrated_residuals) or 0.0, 1) or None,
            "fit": {"slope": round(fit[0], 4), "intercept": round(fit[1], 1)} if fit else None,
            "sufficient_controls": len(subset) >= MIN_LIBRARY_CONTROLS,
        }
    return out


def decoy_discrimination(
    records: list[dict[str, Any]],
    *,
    mass_tolerance: float = DECOY_MASS_TOLERANCE_DA,
    use_recalibrated: bool = True,
) -> dict[str, Any]:
    """Can the residual tell a control's true structure from a plausible impostor?

    The decoys are other controls' true structures from the same library within
    ``mass_tolerance`` of the true neutral mass, which is the population a
    mass-gated candidate pool would actually contain. Scoring is ``-|residual|``
    so that a better retention fit scores higher.
    """
    key = "recalibrated_ri" if use_recalibrated else "proxy_ri"
    scores: list[float] = []
    labels: list[bool] = []
    per_control: list[dict[str, Any]] = []
    for record in records:
        if key not in record:
            continue
        pool = [
            other
            for other in records
            if other["lib_id"] == record["lib_id"]
            and other["control_id"] != record["control_id"]
            and key in other
            and abs(other["exact_mass"] - record["exact_mass"]) <= mass_tolerance
        ]
        if not pool:
            continue
        true_residual = abs(record["observed_ri"] - record[key])
        scores.append(-true_residual)
        labels.append(True)
        decoy_residuals = []
        for decoy in pool:
            # Same observed index, a different structure's predicted index.
            residual = abs(record["observed_ri"] - decoy[key])
            scores.append(-residual)
            labels.append(False)
            decoy_residuals.append(residual)
        per_control.append(
            {
                "control_id": record["control_id"],
                "lib_id": record["lib_id"],
                "n_decoys": len(pool),
                "true_residual": round(true_residual, 1),
                "median_decoy_residual": round(sorted(decoy_residuals)[len(decoy_residuals) // 2], 1),
                "beats_median_decoy": true_residual
                < sorted(decoy_residuals)[len(decoy_residuals) // 2],
            }
        )
    area = auc(scores, labels)
    wins = sum(1 for row in per_control if row["beats_median_decoy"])
    return {
        "scored_by": f"-|observed RI - {'recalibrated' if use_recalibrated else 'raw'} proxy RI|",
        "mass_tolerance_da": mass_tolerance,
        "n_controls": len(per_control),
        "n_comparisons": len(scores),
        "auc": round(area, 4) if area is not None else None,
        "controls_beating_median_decoy": wins,
        "controls_beating_median_decoy_share": round(wins / len(per_control), 4)
        if per_control
        else None,
        "by_library": _decoy_by_library(scores, labels, per_control),
        "per_control": per_control,
    }


def _decoy_by_library(
    scores: list[float], labels: list[bool], per_control: list[dict[str, Any]]
) -> dict[str, Any]:
    """Split the decoy comparisons back out by library.

    Pooled discrimination can be carried entirely by one method, and a retention
    axis is only defensible on the libraries where it is actually working.
    """
    grouped: dict[int, tuple[list[float], list[bool]]] = {}
    cursor = 0
    for row in per_control:
        span = 1 + row["n_decoys"]
        bucket = grouped.setdefault(row["lib_id"], ([], []))
        bucket[0].extend(scores[cursor : cursor + span])
        bucket[1].extend(labels[cursor : cursor + span])
        cursor += span
    out = {}
    for lib_id, (lib_scores, lib_labels) in sorted(grouped.items()):
        area = auc(lib_scores, lib_labels)
        out[str(lib_id)] = {
            "n_controls": sum(1 for row in per_control if row["lib_id"] == lib_id),
            "n_comparisons": len(lib_scores),
            "auc": round(area, 4) if area is not None else None,
        }
    return out


def _library_phrase(library_ids: Sequence[str]) -> str:
    if not library_ids:
        return "no library"
    if len(library_ids) == 1:
        return f"library {library_ids[0]}"
    return f"libraries {', '.join(library_ids[:-1])} and {library_ids[-1]}"


def qualifying_libraries(calibration: dict[str, Any], decoy: dict[str, Any]) -> dict[str, Any]:
    """Decide which libraries may carry a retention axis, and say why for the rest.

    Three conditions, all of which have to hold. Enough labeled controls for the
    fit to mean anything; a positive slope, because a model saying that greasier
    molecules elute earlier on a reversed-phase column is fitting noise rather
    than chemistry; and discrimination above chance on that library's own decoys.
    """
    out: dict[str, Any] = {}
    for lib_id, block in calibration.items():
        area = (decoy.get("by_library", {}).get(lib_id) or {}).get("auc")
        fit = block.get("fit") or {}
        slope = fit.get("slope")
        reasons = []
        if not block["sufficient_controls"]:
            reasons.append(f"only {block['n']} labeled controls")
        if slope is None:
            reasons.append("no fit could be made")
        elif slope <= 0:
            reasons.append(
                f"fitted slope {slope} is negative, which reverses the expected "
                "relationship between lipophilicity and retention"
            )
        if area is None:
            reasons.append("no decoy comparisons available")
        elif area < MIN_LIBRARY_AUC:
            reasons.append(f"decoy AUC {area} is below the {MIN_LIBRARY_AUC} bar")
        out[lib_id] = {
            "ion_mode": block["ion_mode"],
            "n": block["n"],
            "decoy_auc": area,
            "slope": slope,
            "intercept": fit.get("intercept"),
            "residual_rmse": block["recalibrated_loo_rmse"],
            "qualifies": not reasons,
            "excluded_because": reasons,
        }
    return out


def render_markdown(summary: dict[str, Any]) -> str:
    lines = [
        "# Retention-index feasibility audit",
        "",
        "Generated by `tools/audit_rt_feasibility.py`. Whether a logP-derived",
        "retention model can eliminate infeasible candidates, measured before",
        "anything is wired into the reports.",
        "",
    ]
    if summary.get("status") == "not_evaluated":
        return "\n".join(lines + [f"Not evaluated: {summary['reason']}", ""])
    lines += [
        f"Cohort: {summary['cohort']['n_measurements']} measurements over "
        f"{summary['cohort']['n_controls']} controls with both a known structure and an "
        "observed retention index.",
        "",
        "## Is the proxy calibrated?",
        "",
        "| Library | Ion mode | n | Observed RI range | Pearson r | Raw RMSE | Recalibrated RMSE (LOO) |",
        "|---|---|---|---|---|---|---|",
    ]
    for lib_id, block in summary["calibration"].items():
        marker = "" if block["sufficient_controls"] else " *"
        correlation = "n/a" if block["pearson_r"] is None else f"{block['pearson_r']:.3f}"
        lines.append(
            f"| {lib_id}{marker} | {block['ion_mode']} | {block['n']} | "
            f"{block['observed_ri_range'][0]:.0f}-{block['observed_ri_range'][1]:.0f} | "
            f"{correlation} | {block['raw_rmse']} | {block['recalibrated_loo_rmse']} |"
        )
    lines += [
        "",
        f"`*` fewer than {MIN_LIBRARY_CONTROLS} labeled controls; the fit is not interpretable.",
        "",
        "## Does the residual separate right from wrong?",
        "",
    ]
    for label, block in (
        ("Recalibrated", summary["decoy_recalibrated"]),
        ("Raw proxy", summary["decoy_raw"]),
    ):
        area = "n/a" if block["auc"] is None else f"{block['auc']:.3f}"
        lines.append(
            f"- **{label}**: AUC {area} over {block['n_comparisons']} comparisons on "
            f"{block['n_controls']} controls; the true structure beats the median "
            f"mass-matched decoy for {block['controls_beating_median_decoy']} of them."
        )
    lines += [
        "",
        "| Library | Decoy AUC (recalibrated) | Decoy AUC (raw) |",
        "|---|---|---|",
    ]
    for lib_id in summary["calibration"]:
        recal = (summary["decoy_recalibrated"]["by_library"].get(lib_id) or {}).get("auc")
        raw_auc = (summary["decoy_raw"]["by_library"].get(lib_id) or {}).get("auc")
        lines.append(
            f"| {lib_id} | {'n/a' if recal is None else f'{recal:.3f}'} "
            f"| {'n/a' if raw_auc is None else f'{raw_auc:.3f}'} |"
        )
    lines += ["", "## Which libraries may carry a retention axis", ""]
    for lib_id, block in summary["libraries"].items():
        if block["qualifies"]:
            lines.append(
                f"- **Library {lib_id}** ({block['ion_mode']}): qualifies. Slope "
                f"{block['slope']}, residual RMSE {block['residual_rmse']} index units, "
                f"decoy AUC {block['decoy_auc']}."
            )
        else:
            lines.append(
                f"- **Library {lib_id}** ({block['ion_mode']}): excluded - "
                f"{'; '.join(block['excluded_because'])}."
            )
    verdict = summary["gate"]
    lines += [
        "",
        "## Gate",
        "",
        f"Threshold: decoy AUC >= {verdict['threshold']}. Best observed: "
        f"{verdict['best_auc']}. **{verdict['verdict']}**",
        "",
        verdict["reasoning"],
        "",
    ]
    return "\n".join(lines)


def write_calibration(summary: dict[str, Any], path: Path) -> None:
    """Emit the per-library fit the confidence model loads.

    Written as data rather than baked into the module so the numbers stay
    traceable to the audit that produced them and cannot drift from it.
    """
    libraries = {
        lib_id: {
            "ion_mode": block["ion_mode"],
            "slope": block["slope"],
            "intercept": block["intercept"],
            "residual_rmse": block["residual_rmse"],
            "n_controls": block["n"],
            "decoy_auc": block["decoy_auc"],
        }
        for lib_id, block in summary.get("libraries", {}).items()
        if block["qualifies"]
    }
    path.write_text(
        json.dumps(
            {
                "note": (
                    "Per-library recalibration of the logP-derived retention-index proxy, "
                    "fitted by tools/audit_rt_feasibility.py. Only libraries that cleared the "
                    "audit gate appear here; a library that is absent gets no retention axis, "
                    "which the confidence model reports as unavailable evidence rather than as "
                    "a zero score."
                ),
                "proxy": "structural_elucidation.nce.evidence.ri_proxy_from_smiles",
                "predicted_ri": "slope * proxy_ri + intercept",
                "libraries": libraries,
            },
            indent=2,
        )
        + "\n"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--feature-rows-dir", type=Path, default=Path("data/feature_rows"))
    parser.add_argument(
        "--truth-structures",
        type=Path,
        default=Path("data/known_X_all_available_structures_simple.csv"),
    )
    parser.add_argument("--out-dir", type=Path, default=Path("docs/paper/data/retention"))
    parser.add_argument(
        "--mass-tolerance-da", type=float, default=DECOY_MASS_TOLERANCE_DA
    )
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    if not HAS_RDKIT:
        summary: dict[str, Any] = {
            "status": "not_evaluated",
            "reason": "RDKit is unavailable, so no retention proxy can be computed.",
        }
    else:
        records = build_records(
            load_truth_structures(args.truth_structures),
            load_observed_ri(args.feature_rows_dir),
        )
        if not records:
            summary = {
                "status": "not_evaluated",
                "reason": "No control has both a known structure and an observed retention index.",
            }
        else:
            calibration = calibration_by_library(records)
            recalibrated = decoy_discrimination(
                records, mass_tolerance=args.mass_tolerance_da, use_recalibrated=True
            )
            raw = decoy_discrimination(
                records, mass_tolerance=args.mass_tolerance_da, use_recalibrated=False
            )
            best = max(
                (block["auc"] for block in (recalibrated, raw) if block["auc"] is not None),
                default=None,
            )
            libraries = qualifying_libraries(calibration, recalibrated)
            qualified = [lib for lib, block in libraries.items() if block["qualifies"]]
            passed = best is not None and best >= DECOY_AUC_GATE and bool(qualified)
            summary = {
                "question": (
                    "Can a logP-derived retention model eliminate chromatographically "
                    "infeasible candidates?"
                ),
                "proxy": (
                    "structural_elucidation.nce.evidence.ri_proxy_from_smiles: a linear "
                    "combination of exact mass, MolLogP, ring count and TPSA"
                ),
                "cohort": {
                    "n_measurements": len(records),
                    "n_controls": len({record["control_id"] for record in records}),
                    "source": (
                        "observed retention index parsed from the feature key; true "
                        "structures from the known-X structure set"
                    ),
                },
                "calibration": calibration,
                "decoy_recalibrated": recalibrated,
                "decoy_raw": raw,
                "libraries": libraries,
                "gate": {
                    "threshold": DECOY_AUC_GATE,
                    "best_auc": best,
                    "verdict": "PASS" if passed else "FAIL",
                    "qualifying_libraries": qualified,
                    "reasoning": (
                        "The residual separates true structures from mass-matched decoys well "
                        f"enough on {_library_phrase(qualified)} to justify a retention axis "
                        "there. It carries no information at all before recalibration, so the "
                        "published proxy constants cannot be used as they stand, and the axis is "
                        "soft and demote-only: retention can withdraw support from a candidate "
                        "but never supply it."
                        if passed
                        else "The residual does not separate true structures from mass-matched "
                        "decoys well enough to filter on. A retention axis built on this proxy "
                        "would remove correct answers about as often as wrong ones, so it stays "
                        "out of the reports until a chromatography-specific model trained on "
                        "library compounds with measured indices is available."
                    ),
                },
                "scope_note": (
                    "Decoys are other controls' true structures, which are real metabolites of "
                    "similar mass on the same column but not the actual candidate pool. This "
                    "bounds the discrimination a retention filter could have; it does not "
                    "measure the filter in place."
                ),
            }

    json_path = args.out_dir / "audit_rt_feasibility.json"
    json_path.write_text(json.dumps(summary, indent=2) + "\n")
    markdown_path = args.out_dir / "audit_rt_feasibility.md"
    markdown_path.write_text(render_markdown(summary) + "\n")
    write_calibration(summary, args.out_dir / "ri_calibration.json")

    if summary.get("status") == "not_evaluated":
        print(f"[retention] not evaluated: {summary['reason']}")
    else:
        for lib_id, block in summary["calibration"].items():
            print(
                f"[retention] lib {lib_id} ({block['ion_mode']}): n={block['n']} "
                f"r={block['pearson_r']} raw RMSE={block['raw_rmse']} "
                f"recalibrated RMSE={block['recalibrated_loo_rmse']}"
            )
        gate = summary["gate"]
        print(
            f"[retention] decoy AUC recalibrated={summary['decoy_recalibrated']['auc']} "
            f"raw={summary['decoy_raw']['auc']}; gate {DECOY_AUC_GATE} -> {gate['verdict']}"
        )
    print(f"[retention] wrote {json_path} and {markdown_path}")


if __name__ == "__main__":
    main()
