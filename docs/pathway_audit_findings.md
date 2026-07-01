# Pathway / subpathway prediction audit - findings

Audit of the neighbor-consensus pathway predictions (`predicted_pathway` /
`predicted_subpathway` + confidence tier / status / similarity) that feed the
compound reports. Tool: `tools/audit_pathway_predictions.py` (raw outputs under
`output/_pathway_audit/`).

## Headline: accuracy vs the reference-DB taxonomy the predictor targets
The predictor votes Metabolon `metabolon_super`/`metabolon_sub` labels from the
reference DB, so it is scored against that taxonomy (exact-structure match by
InChIKey):
- **Super-pathway accuracy: 94.5%** (sub-pathway 89.7%) on 1466 exact-structure controls.
- Leave-one-out over 5095 reference metabolites (novel-compound proxy): super **86.7%** / sub 71.1%; by similarity: 92.7% at >=0.7, 77.3% at 0.5-0.7, 68.6% at 0.3-0.5.

## Why an earlier known_X-based number looked low (~56%)
`data/known_X.csv` `SuperPathway` is a *different, independent* Metabolon
annotation from the reference DB. On 267 X-compounds whose structure is identical
to a reference entry, `known_X` and the reference DB agree only **62.5%** of the
time (e.g. dipeptides labeled "Amino Acid" in the DB vs "Peptide" in known_X;
acyl-glycines "Lipid" vs "Amino Acid"). So scoring the predictor against known_X
is capped near 62.5% and understates it - it measures annotation disagreement,
not predictor error. Against known_X the raw super agreement is 56.2% (909
evaluable); still well-calibrated by tier (high 71.8%, medium 62.9%, low 37.7%).

## Coverage + confidence (3897 controls, report-deduped)
- Status: `strong_similarity_vote` 3452, `weak_similarity_backfill` 445.
- Confidence tier: high 1788, medium 666, low 1443.
- `pathway_best_reference_similarity` p10 / p50 / p90 = 0.375 / 0.721 / 1.0.
- 100% of controls carry a rank-1 prediction (vs ~40% exact-DB coverage before).

## Calibration (trustworthy tiers)
Accuracy rises monotonically with confidence tier and structural similarity, both
vs known_X and in the LOO. So the report's tier labels are meaningful: treat
`high` as reliable, `low` / `weak_similarity_backfill` as tentative. Sub-pathway
is intrinsically harder and should be read as indicative (89.7% on exact
structures but drops off quickly with novelty).

## Bias + vocabulary
- Predicted super distribution tracks the truth distribution closely (no severe
  over-prediction of common classes); rare classes are slightly under-predicted.
- Vocabulary drift is minimal.

## Reproduce
```bash
python tools/audit_pathway_predictions.py --equalization-root output/pipeline_equalization --loo
```

## Reuse / cross-check
Canonical accuracy/holdout audit lives in the local tree
(`assign_predicted_pathways.py --evaluate --regression-check`,
`build_v15_pathway_accuracy_audit_v1.py`) and should be evaluated against the
reference-DB Metabolon taxonomy (not known_X) for a consistent number; reconcile
against those outputs when available.
