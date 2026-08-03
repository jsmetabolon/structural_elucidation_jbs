#!/usr/bin/env python3
"""Generate manuscript Supplementary Information from committed audit artifacts."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


def pct(metric: dict[str, Any]) -> str:
    if metric.get("accuracy") is None:
        return "NA"
    low, high = metric["ci95"]
    return f"{100 * metric['accuracy']:.1f}% ({100 * low:.1f}–{100 * high:.1f}%)"


def claim_map(path: Path) -> dict[str, dict[str, Any]]:
    return {row["claim_id"]: row for row in json.loads(path.read_text())["claims"]}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--structure", type=Path, required=True)
    parser.add_argument("--pathway", type=Path, required=True)
    parser.add_argument("--loo", type=Path, required=True)
    parser.add_argument("--presentation", type=Path, required=True)
    parser.add_argument("--claims", type=Path, required=True)
    parser.add_argument("--out", type=Path, default=Path("docs/paper/supplement/supplementary_information.md"))
    args = parser.parse_args()
    args.out.parent.mkdir(parents=True, exist_ok=True)

    structure = json.loads(args.structure.read_text())
    pathway = json.loads(args.pathway.read_text())
    loo = json.loads(args.loo.read_text())
    presentation = json.loads(args.presentation.read_text())
    claims = claim_map(args.claims)
    known = pathway["accuracy"]
    n_reports = int(presentation["n_reports"])
    coverage = presentation["evidence_coverage"]
    gaps = presentation["presentation_gaps"]
    msi_by_tier = presentation["msi_by_tier"]

    lines = [
        "# Supplementary Information",
        "",
        "## S1. Internal toolkit-benchmarking design (prospective)",
        "",
        "The earlier internal study design compares MS2Query, SIRIUS, Spec2Mol, and "
        "MSNovelist with the legacy workflow using CASMI categories: adduct, elemental "
        "formula, structure class, and 2D structure. Library 1 contains 2,091 known "
        "serum/plasma metabolites spanning 9 super-pathways and 135 subpathways. "
        "Library 2 contains 120 internally confirmed NMD/X-compounds plus 500 CASMI "
        "2022 compounds. No comparative toolkit results were available for this "
        "manuscript; the design is prospective and its counts are not pooled with the "
        "3,897-feature deployed evaluation.",
        "",
        "## S2. Cohorts and confidence vocabularies",
        "",
        "| Term | Definition | Source / status |",
        "|---|---|---|",
        "| Delivered set | 3,897 unique X-compounds; first library wins duplicates | report bundle |",
        "| Structure truth set | 92 controls with nonzero `truth_tanimoto` | candidate sidecars |",
        "| Reference-structure cohort | 1,466 rank-1 structures present in reference DB | 1,459 direct lookups |",
        "| `known_X` cohort | 909 controls with independent super-/subpathway labels | taxonomy-discordant |",
        "| P0 | true unknowns with `has_p0_reference=true` | authoritative P0 flags file unavailable |",
        "| P1 | Metabolon Level 1 plus Level 2/2C | authoritative P1 truth directory unavailable |",
        "| Pathway-confidence tier | high / medium / low neighbor-vote confidence | not structure confidence |",
        "| Identification-confidence tier | Confirmed / High / Moderate / Low / Tentative | exploratory Platt score |",
        "| MSI reporting level | evidence-derived L1–L5 rubric | custom; not externally validated |",
        "",
        "**Table S1. Evidence coverage across all delivered reports.**",
        "",
        "| Evidence | n / 3,897 (%) |",
        "|---|---:|",
        f"| Biochemical taxonomy assignment | {n_reports:,} (100.0%) |",
        f"| Observed parent m/z | {round(n_reports * coverage['observed_mz_available']):,} "
        f"({100 * coverage['observed_mz_available']:.1f}%) |",
        f"| Metabolon Level (any) | {round(n_reports * gaps['has_metabolon_level']):,} "
        f"({100 * gaps['has_metabolon_level']:.1f}%) |",
        f"| Observed MS2 spectrum | {round(n_reports * coverage['ms2_spectrum_present']):,} "
        f"({100 * coverage['ms2_spectrum_present']:.1f}%) |",
        f"| Reference-database name resolved | {round(n_reports * gaps['name_resolved_from_db']):,} "
        f"({100 * gaps['name_resolved_from_db']:.1f}%) |",
        f"| Top-candidate mass mismatch | {round(n_reports * coverage['mass_mismatch']):,} "
        f"({100 * coverage['mass_mismatch']:.1f}%) |",
        f"| Public MS2-ion leaks | {presentation['compliance']['public_ms2_leaks']} (0.0%) |",
        "",
        "**Table S2. MSI-style reporting level × identification-confidence tier (counts).**",
        "",
        "| MSI level | Confirmed | High | Moderate | Low | Tentative |",
        "|---|---:|---:|---:|---:|---:|",
        *[
            "| L{level} | {confirmed} | {high} | {moderate} | {low} | {tentative} |".format(
                level=level,
                confirmed=msi_by_tier.get(f"L{level}|Confirmed", 0),
                high=msi_by_tier.get(f"L{level}|High", 0),
                moderate=msi_by_tier.get(f"L{level}|Moderate", 0),
                low=msi_by_tier.get(f"L{level}|Low", 0),
                tentative=msi_by_tier.get(f"L{level}|Tentative", 0),
            )
            for level in range(1, 6)
        ],
        "",
        "## S3. Structural recovery and threshold sensitivity",
        "",
        "**Table S3. Cumulative structural recovery on n=92 truth controls "
        "(Wilson 95% confidence intervals).**",
        "",
        "| Truth threshold | Rank 1 | Top 3 | Top 10 |",
        "|---|---:|---:|---:|",
    ]
    for threshold, metrics in structure["rank_k"].items():
        label = "Exact 2D" if threshold == "exact" else f"Tanimoto ≥ {threshold}"
        lines.append(
            f"| {label} | {pct(metrics['top1'])} | {pct(metrics['top3'])} | {pct(metrics['top10'])} |"
        )
    lines += [
        "",
        "**Table S4. Rank-1 recovery by provenance.**",
        "",
        "| Source | n | Accuracy (95% CI) | Median truth Tanimoto (bootstrap 95% CI) |",
        "|---|---:|---:|---:|",
    ]
    for source, result in structure["source_subgroups"].items():
        metric = result["rank1_tanimoto_ge_0p70"]
        median_ci = result["median_truth_tanimoto_ci95"]
        lines.append(
            f"| {source} | {metric['n']} | {pct(metric)} | "
            f"{result['median_truth_tanimoto']:.3f} ({median_ci[0]:.3f}–{median_ci[1]:.3f}) |"
        )
    lines += [
        "",
        "### S3.1 Cross-validated calibration",
        "",
        f"Full model leave-one-out Brier score: **{structure['calibration_full']['brier_score']:.4f}**; "
        f"calibration slope {structure['calibration_full']['calibration']['slope']:.3f}, "
        f"intercept {structure['calibration_full']['calibration']['intercept']:.3f}. "
        f"Without the pathway axis, Brier score is "
        f"**{structure['calibration_without_pathway']['brier_score']:.4f}**. "
        "These results indicate weak discrimination/calibration and motivate treating "
        "the score as exploratory.",
        "",
        "![Figure S1. Independent pathway agreement by pathway-confidence tier.](../figures/figS1_pathway_tier.svg)",
        "",
        "## S4. Pathway evaluation details",
        "",
        "**Table S5. Independent-label agreement by assignment status.**",
        "",
        "| Status | Super-pathway (95% CI) | Subpathway (95% CI) |",
        "|---|---:|---:|",
    ]
    for status, metrics in known["by_assignment_status"].items():
        lines.append(f"| {status} | {pct(metrics['super'])} | {pct(metrics['sub'])} |")
    lines += [
        "",
        f"Leave-one-out uses k={loo['k']}, similarity floor={loo['similarity_floor']}, "
        f"and similarity exponent={loo['similarity_power']}. It evaluates the neighbor-vote "
        "mechanism only, not the complete production predictor.",
        "",
        "![Figure S2. Neighbor-vote leave-one-out accuracy by similarity.](../figures/figS2_loo_similarity.svg)",
        "",
        "![Figure S3. Independent-label agreement by assignment status.](../figures/figS3_pathway_status.svg)",
        "",
        "**Table S6. Predicted vs independent `known_X` super-pathway distribution (n=909).**",
        "",
        "| Super-pathway | Predicted | Independent label |",
        "|---|---:|---:|",
    ]
    predicted = pathway["bias"]["predicted_super_distribution"]
    truth = pathway["bias"]["truth_super_distribution"]
    for label in sorted(set(predicted) | set(truth), key=lambda key: -(truth.get(key, 0))):
        lines.append(f"| {label} | {predicted.get(label, 0)} | {truth.get(label, 0)} |")
    lines += [
        "",
        "**Table S7. Independent-label super-pathway confusion for the four most common classes.**",
        "",
        "| Independent label | Lipid | Peptide | Xenobiotics | Amino Acid | Other |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    confusion = pathway["confusion_top_true_super"]
    for truth_label in ("Lipid", "Peptide", "Xenobiotics", "Amino Acid"):
        row = confusion.get(truth_label, {})
        main = [int(row.get(label, 0)) for label in ("Lipid", "Peptide", "Xenobiotics", "Amino Acid")]
        other = sum(int(value) for key, value in row.items() if key not in {
            "Lipid", "Peptide", "Xenobiotics", "Amino Acid"
        })
        lines.append(f"| {truth_label} | {main[0]} | {main[1]} | {main[2]} | {main[3]} | {other} |")
    lines += [
        "",
        "## S5. Mass-mismatch audit",
        "",
        f"There are **{claims['rank1_mass_mismatch_all']['value']}** rank-1 mass mismatches "
        f"among {claims['rank1_mass_mismatch_all']['n']} reports. "
        f"**{claims['rank1_mass_mismatch_msi_l4_deep_audit']['value']}** are in the "
        "MSI-L4 deep-audit subset. A lower-ranked mass-consistent candidate exists for "
        f"**{claims['rank1_mass_mismatch_recoverable_all']['value']}/"
        f"{claims['rank1_mass_mismatch_all']['value']}** overall and "
        f"**{claims['rank1_mass_mismatch_recoverable_l4']['value']}/"
        f"{claims['rank1_mass_mismatch_msi_l4_deep_audit']['value']}** in the L4 subset. "
        "These are projected reranking opportunities, not demonstrated correct structures.",
        "",
        "![Figure S4. Rank-1 mass-mismatch cohorts.](../figures/figS4_mass_mismatch.svg)",
        "",
        "## S6. Per-library delivered cohort",
        "",
        "| Library | Ion mode | n | Percentage |",
        "|---|---|---:|---:|",
        "| 209 | negative | 2,389 | 61.3% |",
        "| 305 | negative | 202 | 5.2% |",
        "| 400 | positive | 875 | 22.5% |",
        "| 402 | positive | 431 | 11.1% |",
        "",
        "## S7. Required provenance before submission",
        "",
        "- Complete author list, corresponding-author email, CRediT roles, funding, and conflict disclosures.",
        "- Source-study matrix/species metadata and ethics/consent/IRB determination.",
        "- Spec2Mol encoder/decoder checkpoint hashes, training snapshot, random seeds, "
        "software/CUDA versions, and EC2 hardware.",
        "- Production neighbor-builder parameters and override rules.",
        "- Authoritative P0 flags and P1 truth files for cohort-specific rank-1/top-10 results.",
        "- Immutable repository release/DOI and controlled-access data instructions.",
        "",
        f"Full-set identification-confidence median: {presentation['overall_confidence']['median']:.4f}; "
        f"MSI coverage: {presentation['msi_level_coverage']:.1%}.",
    ]
    args.out.write_text("\n".join(lines) + "\n")
    print(f"[paper-supplement] wrote {args.out}")


if __name__ == "__main__":
    main()
