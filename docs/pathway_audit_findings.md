# Pathway / subpathway prediction audit - findings

Audit of the neighbor-consensus pathway predictions (`predicted_pathway` /
`predicted_subpathway` + confidence tier / status / similarity) that feed the
compound reports. Tool: `tools/audit_pathway_predictions.py` (raw outputs under
`output/_pathway_audit/`). Ground truth: `data/known_X.csv`
`SuperPathway` / `Draft Subpathway` (excluding blank / `Unknown` / data-entry
sentinels like "Ask DTS"), joined on `control_id`.

## Coverage + confidence (3897 controls, report-deduped)
- Status: `strong_similarity_vote` 3452, `weak_similarity_backfill` 445.
- Confidence tier: high 1788, medium 666, low 1443.
- `pathway_best_reference_similarity` p10 / p50 / p90 = 0.375 / 0.721 / 1.0.
- 100% of controls carry a rank-1 prediction (vs ~40% exact-DB coverage before).

## Accuracy vs known_X ground truth (909 evaluable)
- **Super-pathway accuracy: 56.2%.**
- Sub-pathway accuracy: 24.2% (43.0% given a correct super-pathway).
- By confidence tier: high **71.8%** (n=365), medium 62.9% (n=175), low 37.7% (n=369).
- By best-similarity bin: [0.0,0.3) 25.6%, [0.3,0.5) 31.9%, [0.5,0.7) 50.3%, [0.7,1.0) **69.6%**.
- By status: `strong_similarity_vote` 58.2%, `weak_similarity_backfill` 36.1%.

## Intrinsic mechanism accuracy - leave-one-out over the reference DB (5095 labeled metabolites)
- Super-pathway **86.7%**, sub-pathway 71.1%.
- By best-similarity bin: [0.3,0.5) 68.6%, [0.5,0.7) 77.3%, [0.7,1.0) **92.7%**.
- (This is the core Tanimoto k-NN consensus, not the full production predictor
  with MW priors / xenobiotic classifier / SMARTS overrides - that is audited by
  the local `build_v15_pathway_accuracy_audit_v1.py`.)

## Confusion (top true super -> predicted)
- Lipid -> mostly Lipid, else Peptide / Xenobiotics (chemically adjacent, e.g. acyl-glycines, lipopeptides).
- Peptide -> Peptide, else Lipid / Amino Acid.
- Xenobiotics -> Xenobiotics, else Lipid / Amino Acid.
- Amino Acid -> Amino Acid, else Xenobiotics / Peptide.

## Bias + vocabulary
- Predicted super distribution tracks the truth distribution closely (no severe
  over-prediction of common classes); rare classes (Cofactors and Vitamins,
  Energy) are slightly under-predicted.
- Vocabulary drift is minimal (only "Partially Characterized Molecules" predicted
  but absent from truth vocab); truth-only labels are data-entry sentinels.

## Interpretation
- The mechanism is strong (86.7% LOO); accuracy on the X-compounds (56.2% super)
  is lower because the unknowns are genuinely more novel - most sit in the
  0.3-0.7 similarity range, where accuracy is 32-50%.
- Accuracy is **strongly gated by structural similarity / confidence tier**, so
  the report's tier labels are trustworthy: treat `high` (~72% here, ~93% at high
  similarity) as reliable, `low`/`weak_similarity_backfill` (~37%) as tentative.
- Sub-pathway prediction is much harder (24%); it should be read as indicative.

## Reproduce
```bash
python tools/audit_pathway_predictions.py --equalization-root output/pipeline_equalization --loo
```

## Reuse / cross-check
Canonical accuracy/holdout audit lives in the local tree
(`assign_predicted_pathways.py --evaluate --regression-check`,
`build_v15_pathway_accuracy_audit_v1.py`). This cloud audit reports accuracy on
the exact delivered predictions; reconcile against those local outputs when
available.
