# Pathway / subpathway prediction audit - findings

> **Superseded. Retained as a record of the original investigation, not as current
> results.** Several headline claims below were withdrawn by the scientific revision
> recorded in `docs/paper/REVISION_LOG.md`, and several counts predate a later audit
> run. The authoritative figures are in `docs/paper/manuscript.md` and the committed
> artifacts under `docs/paper/data/`. Read the corrections here before the body.
>
> | Claim below | Status |
> |---|---|
> | "Super-pathway accuracy: 94.5% (sub-pathway 89.7%)" | **Withdrawn.** These were reference-taxonomy *consistency* on exact-structure matches, where the label is copied from the same database it is scored against. That is circular as an accuracy test. The manuscript now reports it as 100% consistency by construction for the 1,459 direct lookups and does not present it as prediction accuracy. |
> | Coverage statuses 3452 / 445 | **Stale.** Current audit: exact reference match 1,459, strong-similarity vote 1,995, weak-similarity backfill 443. |
> | Confidence tiers 1788 / 666 / 1443 | **Stale.** Current audit: high 1,863, medium 616, low 1,418. |
> | Similarity p50 = 0.721 | **Stale.** Current audit: 0.733. |
> | known_X super agreement 56.2% | **Stale.** Current audit: 55.3% (503/909). |
> | Tier accuracies 71.8 / 62.9 / 37.7 | **Stale.** Current audit: high 70.3%, medium 61.4%, low 36.9%. |
> | "treat `high` as reliable" | **Overstated.** Agreement is monotonic by tier, which supports prioritization, but the top tier still agrees with the independent labels only about seven times in ten. |
>
> One conclusion below held up and was later confirmed: the 62.5% agreement between
> the two taxonomies on identical structures is a genuine ceiling on any
> `known_X`-scored result. The manuscript now carries it, along with the subpathway
> equivalent of 34.8%, which this document did not measure.

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
