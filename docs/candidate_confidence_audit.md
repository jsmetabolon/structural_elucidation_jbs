# Candidate-confidence audit (delivered 3,897-report bundle)

Real-bundle before/after for the rank-consistent per-candidate confidence
(`s3://metabolon-rd-raw-data/Structural_Elucidation/se_report_outputs.tgz`,
processed 2026-07-15). Reproduce with the commands at the end; machine-readable
snapshots are in [audit/candidate_confidence_before.json](audit/candidate_confidence_before.json)
and [audit/candidate_confidence_after.json](audit/candidate_confidence_after.json).

Cohort: 3,897 delivered reports; 3,786 have >=2 top candidates and are evaluable
(111 single-candidate reports skipped).

## Problem (what the delivered reports featured: `evidence_confidence`)

| Metric (top-10 per report) | Value |
|---|---:|
| Reports where all top-N are identical | **28.5%** |
| Reports with >=1 rank-vs-score inversion (higher rank, lower score) | **58.8%** |
| Total inversions across the cohort | **19,088** |
| Median Kendall tau vs rank (1.0 = monotone) | **0.619** |

The ordering score itself (`candidate_report_score` = `unified_final_score`) is
monotone with rank (0% inversions, tau 1.0). The problem was entirely the
*featured* number: `evidence_confidence` saturates (shared across mass-consistent
isomers) and is a reranker input feature, not the ranking objective.

## After the fix (reports now feature `candidate_confidence`)

| Metric (top-10 per report) | Value |
|---|---:|
| Rank-vs-score inversions | **0** (tau 1.0) |
| Reports whose rank-1 is genuinely tied (ambiguous headline) | **6.7%** |
| Median rank-1 candidate confidence | **0.431** |
| Shipped (patched sidecar) == freshly computed | **True (all 3,786)** |

`evidence_confidence` is retained and relabeled "Formula/class evidence" (with
"Mass match") as supporting context. Ties concentrate in the low-confidence tail
where the reranker genuinely cannot separate candidates; the top candidate is
distinct in 93.3% of reports.

## Spot checks

**X-02269** (the shared-isomer example): before, ranks 1-5 all showed
`evidence 0.998612` (identical) and rank 8 (0.986) outscored rank 6 (0.317).
After: `candidate_confidence` 0.585 -> 0.112 -> 0.083 -> 0.072 -> 0.065 -> 0.048
-> 0.025 (tied) -> 0.024 (tied) -> 0.018 -> 0.010 - monotone with rank.

**X-22167** (worst inversion): before, rank 1 = `evidence 0.464` while ranks 7-10
= `evidence 0.996`. After: `candidate_confidence` 0.493 -> 0.482 -> 0.482 ->
0.394 -> ... -> 0.311 - monotone.

## Reproduce

```bash
python3 tools/audit_candidate_scores.py <reports_root> --json before.json   # baseline
python3 tools/refresh_candidate_confidence.py <reports_root>                 # add candidate_confidence to confidence_inputs.json
python3 tools/patch_candidate_scores.py <reports_root>                       # upgrade every *.html widget (idempotent)
python3 tools/audit_candidate_scores.py <reports_root> --json after.json     # candidate_confidence inversions = 0; parity True
```
