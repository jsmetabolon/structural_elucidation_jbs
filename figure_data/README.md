# Figure source data

Copied from `docs/paper/` (committed ~2026-07-10 in `1235c80`).
These are the audited aggregate artifacts consumed by `tools/paper_figures.py`.
Upstream proprietary inputs (`_audit_pull/reports`, `_bundle/metabolite_synonyms.sqlite`) are not in this checkout.

| Figure | Source file(s) |
|--------|----------------|
| fig1_workflow.svg | hand-authored SVG (no tabular input) |
| fig2_structure_accuracy.svg | structure/audit_structure_summary.json |
| fig3_confidence_reliability.svg | structure/audit_structure_summary.json |
| fig4_identification_distributions.svg | presentation/audit_presentation_summary.json |
| fig5_pathway_evaluations.svg | pathway/audit_pathway_summary.json + pathway/audit_pathway_loo.json |
| figS1_pathway_tier.svg | pathway/audit_pathway_summary.json |
| figS2_loo_similarity.svg | pathway/audit_pathway_loo.json |
| figS3_pathway_status.svg | pathway/audit_pathway_summary.json |
| figS4_mass_mismatch.svg | claims/baseline_claims.json |
| (supporting) | structure/audit_structure_controls.csv, pathway/audit_pathway_accuracy.csv, claims/claim_evidence.csv |

Regenerate SVGs:
```bash
python3 tools/paper_figures.py \
  --structure docs_figure_data/structure/audit_structure_summary.json \
  --pathway docs_figure_data/pathway/audit_pathway_summary.json \
  --loo docs_figure_data/pathway/audit_pathway_loo.json \
  --presentation docs_figure_data/presentation/audit_presentation_summary.json \
  --claims docs_figure_data/claims/baseline_claims.json \
  --out-dir docs/paper/figures
```
