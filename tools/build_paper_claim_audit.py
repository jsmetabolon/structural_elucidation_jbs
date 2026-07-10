#!/usr/bin/env python3
"""Build the manuscript claim-to-evidence ledger from delivered report artifacts.

The ledger freezes every headline number with its cohort, denominator, metric,
interpretation, and source artifact.  It also reconciles two commonly conflated
mass-mismatch cohorts:

* all reports with rank-1 mass_match_score < 0.5; and
* the subset assigned MSI L4 (the historical 590-case "deep audit").

Usage:
    python3 tools/build_paper_claim_audit.py \
        --reports-root _audit_pull/reports \
        --presentation-summary output/_presentation_audit/audit_presentation_summary.json \
        --out-dir docs/paper/supplement
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any


def _num(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return None if number != number else number


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _ranked_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    return sorted(rows, key=lambda row: _num(row.get("rank")) or 1e9)


def _claim(
    claim_id: str,
    value: Any,
    *,
    cohort: str,
    n: int | None,
    metric: str,
    source: str,
    interpretation: str,
) -> dict[str, Any]:
    return {
        "claim_id": claim_id,
        "value": value,
        "cohort": cohort,
        "n": n,
        "metric": metric,
        "source": source,
        "interpretation": interpretation,
    }


def build(args: argparse.Namespace) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    reports_root = args.reports_root
    pathway_path = args.pathway_summary or (reports_root / "audit_pathway_summary.json")
    pathway = json.loads(pathway_path.read_text())
    loo = json.loads(args.loo_summary.read_text()) if args.loo_summary and args.loo_summary.exists() else None
    presentation = json.loads(args.presentation_summary.read_text())

    total_mass_mismatch = 0
    l4_mass_mismatch = 0
    total_recoverable = 0
    l4_recoverable = 0
    truth_controls = 0
    top1_correct = 0
    top3_correct = 0
    top10_correct = 0
    exact_n = 0
    exact_correct = 0
    predicted_n = 0
    predicted_correct = 0

    for report_dir in sorted(reports_root.glob("X-*")):
        if not report_dir.is_dir():
            continue
        slug = report_dir.name.lower()
        data_dir = report_dir / "data"
        confidence_path = data_dir / f"{slug}_confidence_inputs.json"
        if not confidence_path.exists():
            continue
        fields = json.loads(confidence_path.read_text())
        mass_score = _num(fields.get("mass_match_score"))
        is_mismatch = mass_score is not None and mass_score < 0.5
        is_l4 = int(_num(fields.get("msi_level")) or 0) == 4
        recoverable = isinstance(fields.get("mass_alt"), dict)
        total_mass_mismatch += int(is_mismatch)
        l4_mass_mismatch += int(is_mismatch and is_l4)
        total_recoverable += int(is_mismatch and recoverable)
        l4_recoverable += int(is_mismatch and is_l4 and recoverable)

        candidates = _ranked_rows(data_dir / f"{slug}_top_candidates.csv")
        if not candidates:
            continue
        rank1_truth = _num(candidates[0].get("truth_tanimoto"))
        if rank1_truth is None or rank1_truth <= 0:
            continue
        truth_controls += 1
        similarities = [_num(row.get("truth_tanimoto")) for row in candidates]
        similarities = [value if value is not None else 0.0 for value in similarities]
        top1_correct += int(similarities[0] >= 0.70)
        top3_correct += int(max(similarities[:3], default=0.0) >= 0.70)
        top10_correct += int(max(similarities[:10], default=0.0) >= 0.70)
        source = str(fields.get("source") or "")
        if source == "exact":
            exact_n += 1
            exact_correct += int(similarities[0] >= 0.70)
        else:
            predicted_n += 1
            predicted_correct += int(similarities[0] >= 0.70)

    ref = pathway["reference_db_evaluation"]
    independent = pathway["accuracy"]
    direct = ref.get("by_assignment_source", {}).get("exact_reference_lookup")
    neighbor = independent.get("by_assignment_source", {}).get("neighbor_consensus")
    msi = presentation["msi_level_distribution"]
    tiers = presentation["tier_distribution"]
    claims = [
        _claim(
            "reports_total",
            presentation["n_reports"],
            cohort="delivered unique X-compounds; first-library-wins",
            n=presentation["n_reports"],
            metric="report coverage",
            source=str(args.presentation_summary),
            interpretation="Total report set; not an accuracy denominator.",
        ),
        _claim(
            "structure_rank1_tanimoto_ge_0p70",
            round(top1_correct / truth_controls, 4),
            cohort="truth-verified controls",
            n=truth_controls,
            metric=f"{top1_correct}/{truth_controls} with rank-1 Tanimoto >= 0.70",
            source="per-report *_top_candidates.csv",
            interpretation="In-sample structural recovery on the available truth set.",
        ),
        _claim(
            "structure_top3_tanimoto_ge_0p70",
            round(top3_correct / truth_controls, 4),
            cohort="truth-verified controls",
            n=truth_controls,
            metric=f"{top3_correct}/{truth_controls} with a hit in ranks 1-3",
            source="per-report *_top_candidates.csv",
            interpretation="CASMI-style cumulative top-3 recovery.",
        ),
        _claim(
            "structure_top10_tanimoto_ge_0p70",
            round(top10_correct / truth_controls, 4),
            cohort="truth-verified controls",
            n=truth_controls,
            metric=f"{top10_correct}/{truth_controls} with a hit in ranks 1-10",
            source="per-report *_top_candidates.csv",
            interpretation="CASMI-style cumulative top-10 recovery.",
        ),
        _claim(
            "structure_rank1_exact_source",
            round(exact_correct / exact_n, 4),
            cohort="truth-verified controls whose top call source is exact",
            n=exact_n,
            metric=f"{exact_correct}/{exact_n} with rank-1 Tanimoto >= 0.70",
            source="confidence_inputs.json + *_top_candidates.csv",
            interpretation="Catalog provenance is not independent confirmation.",
        ),
        _claim(
            "structure_rank1_predicted_source",
            round(predicted_correct / predicted_n, 4),
            cohort="truth-verified controls whose top call source is predicted",
            n=predicted_n,
            metric=f"{predicted_correct}/{predicted_n} with rank-1 Tanimoto >= 0.70",
            source="confidence_inputs.json + *_top_candidates.csv",
            interpretation="De novo/predicted-source rank-1 recovery.",
        ),
        _claim(
            "pathway_reference_taxonomy_consistency_super",
            direct["super"]["accuracy"] if direct else ref["super_accuracy_vs_db"],
            cohort="direct exact-structure reference-taxonomy lookups",
            n=direct["super"]["n"] if direct else ref["n_exact_structure_in_db"],
            metric="predicted super-pathway equals reference-database label",
            source=str(pathway_path),
            interpretation="Reference-taxonomy consistency by construction; not prediction accuracy.",
        ),
        _claim(
            "pathway_reference_taxonomy_consistency_sub",
            direct["sub"]["accuracy"] if direct else ref["sub_accuracy_vs_db"],
            cohort="direct exact-structure reference-taxonomy lookups",
            n=direct["sub"]["n"] if direct else ref["n_exact_structure_in_db"],
            metric="predicted subpathway equals reference-database label",
            source=str(pathway_path),
            interpretation="Reference-taxonomy consistency; includes direct lookups.",
        ),
        _claim(
            "pathway_neighbor_known_x_agreement_super",
            neighbor["super"]["accuracy"] if neighbor else None,
            cohort="neighbor-consensus controls with independent known_X labels",
            n=neighbor["super"]["n"] if neighbor else None,
            metric="neighbor-consensus super-pathway equals known_X label",
            source=str(pathway_path),
            interpretation="Independent-label agreement; affected by taxonomy disagreement.",
        ),
        _claim(
            "pathway_neighbor_known_x_agreement_sub",
            neighbor["sub"]["accuracy"] if neighbor else None,
            cohort="neighbor-consensus controls with independent known_X labels",
            n=neighbor["sub"]["n"] if neighbor else None,
            metric="neighbor-consensus subpathway equals known_X label",
            source=str(pathway_path),
            interpretation="Independent-label agreement; subpathway recovery is low.",
        ),
        _claim(
            "pathway_known_x_agreement_super",
            independent["super_accuracy"],
            cohort="controls with independent known_X super-pathway labels",
            n=independent["n_evaluable_super"],
            metric="predicted super-pathway equals known_X label",
            source=str(pathway_path),
            interpretation="Independent-label agreement; affected by taxonomy disagreement.",
        ),
        _claim(
            "pathway_known_x_agreement_sub",
            independent["sub_accuracy"],
            cohort="controls with independent known_X subpathway labels",
            n=independent["n_evaluable_sub"],
            metric="predicted subpathway equals known_X label",
            source=str(pathway_path),
            interpretation="Independent-label agreement; subpathway is substantially harder.",
        ),
        _claim(
            "reference_db_vs_known_x_agreement",
            ref["db_vs_known_x_super_agreement"]["agreement"],
            cohort="identical structures with both taxonomy labels",
            n=ref["db_vs_known_x_super_agreement"]["n"],
            metric="reference-DB super-pathway equals known_X label",
            source=str(pathway_path),
            interpretation="Annotation harmonization ceiling for known_X comparisons.",
        ),
        *(
            [
                _claim(
                    "pathway_neighbor_vote_loo_super",
                    loo["super"]["accuracy"],
                    cohort="reference metabolites under leave-one-out neighbor voting",
                    n=loo["super"]["n"],
                    metric="LOO super-pathway accuracy",
                    source=str(args.loo_summary),
                    interpretation="Intrinsic neighbor-vote mechanism; not the full production predictor.",
                ),
                _claim(
                    "pathway_neighbor_vote_loo_sub",
                    loo["sub"]["accuracy"],
                    cohort="reference metabolites under leave-one-out neighbor voting",
                    n=loo["sub"]["n"],
                    metric="LOO subpathway accuracy",
                    source=str(args.loo_summary),
                    interpretation="Intrinsic neighbor-vote mechanism; not the full production predictor.",
                ),
            ]
            if loo
            else []
        ),
        _claim(
            "rank1_mass_mismatch_all",
            total_mass_mismatch,
            cohort="all delivered reports",
            n=presentation["n_reports"],
            metric="rank-1 mass_match_score < 0.50",
            source="per-report confidence_inputs.json",
            interpretation="598-case full mismatch cohort used for coverage.",
        ),
        _claim(
            "rank1_mass_mismatch_msi_l4_deep_audit",
            l4_mass_mismatch,
            cohort="mass-mismatched reports assigned MSI L4",
            n=total_mass_mismatch,
            metric="mass mismatch AND msi_level == 4",
            source="per-report confidence_inputs.json",
            interpretation="590-case historical deep-audit subset; excludes eight mismatches outside L4.",
        ),
        _claim(
            "rank1_mass_mismatch_recoverable_all",
            total_recoverable,
            cohort="all rank-1 mass mismatches",
            n=total_mass_mismatch,
            metric="lower-ranked candidate with strong mass match",
            source="per-report confidence_inputs.json mass_alt",
            interpretation="Projected rerank opportunity, not achieved structural correctness.",
        ),
        _claim(
            "rank1_mass_mismatch_recoverable_l4",
            l4_recoverable,
            cohort="MSI-L4 mass-mismatch deep-audit subset",
            n=l4_mass_mismatch,
            metric="lower-ranked candidate with strong mass match",
            source="per-report confidence_inputs.json mass_alt",
            interpretation="Historical 519/590 projection; not a demonstrated post-rerank result.",
        ),
        _claim(
            "msi_distribution",
            msi,
            cohort="all delivered reports",
            n=presentation["n_reports"],
            metric="derived MSI/Schymanski rubric",
            source=str(args.presentation_summary),
            interpretation="Custom evidence-derived levels; thresholds require external validation.",
        ),
        _claim(
            "identification_confidence_tier_distribution",
            tiers,
            cohort="all delivered reports",
            n=presentation["n_reports"],
            metric="calibrated structure-identification tier",
            source=str(args.presentation_summary),
            interpretation="Not the pathway-confidence tier.",
        ),
    ]
    snapshot = {
        "sources": {
            str(pathway_path): _sha256(pathway_path),
            str(args.presentation_summary): _sha256(args.presentation_summary),
            **({str(args.loo_summary): _sha256(args.loo_summary)} if args.loo_summary else {}),
        },
        "claims": claims,
    }
    return claims, snapshot


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reports-root", type=Path, default=Path("_audit_pull/reports"))
    parser.add_argument(
        "--presentation-summary",
        type=Path,
        default=Path("output/_presentation_audit/audit_presentation_summary.json"),
    )
    parser.add_argument("--pathway-summary", type=Path, default=None)
    parser.add_argument("--loo-summary", type=Path, default=None)
    parser.add_argument("--out-dir", type=Path, default=Path("docs/paper/supplement"))
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    claims, snapshot = build(args)

    csv_path = args.out_dir / "claim_evidence.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(claims[0]))
        writer.writeheader()
        for row in claims:
            serialized = dict(row)
            if isinstance(serialized["value"], dict):
                serialized["value"] = json.dumps(serialized["value"], sort_keys=True)
            writer.writerow(serialized)
    json_path = args.out_dir / "baseline_claims.json"
    json_path.write_text(json.dumps(snapshot, indent=2) + "\n")
    print(f"[paper-claim-audit] wrote {csv_path} and {json_path}")


if __name__ == "__main__":
    main()
