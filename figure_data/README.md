# Figure source data

**The source of truth is `docs/paper/`.** This directory holds no data of its own; it
is the output location for `tools/build_figures_presentation.py`.

It previously held a copy of the audited artifacts, taken in `1235c80`, plus flattened
CSV extracts and PNG renders. Nothing read any of it, no tool regenerated it, and the
claim files had drifted from `docs/paper/supplement/` by the time it was removed. A
duplicate that cannot be regenerated can only get staler, so it was deleted rather
than resynchronised.

## Where the data lives now

| Figure | Source file(s) |
|--------|----------------|
| fig1_workflow.svg | hand-authored SVG (no tabular input) |
| fig2_structure_accuracy.svg | `docs/paper/data/structure/audit_structure_summary.json` |
| fig3_confidence_reliability.svg | `docs/paper/data/structure/audit_structure_summary.json` |
| fig4_identification_distributions.svg | `docs/paper/data/presentation/audit_presentation_summary.json` |
| fig5_pathway_evaluations.svg | `docs/paper/data/pathway/audit_pathway_summary.json` + `audit_pathway_loo.json` |
| figS1_pathway_tier.svg | `docs/paper/data/pathway/audit_pathway_summary.json` |
| figS2_loo_similarity.svg | `docs/paper/data/pathway/audit_pathway_loo.json` |
| figS3_pathway_status.svg | `docs/paper/data/pathway/audit_pathway_summary.json` |
| figS4_mass_mismatch.svg | `docs/paper/supplement/baseline_claims.json` |
| figS5_selective_prediction.svg | `docs/paper/data/selective_prediction.json` |

Per-control data behind the aggregates is committed at
`docs/paper/data/structure/audit_structure_controls.csv` and
`docs/paper/data/pathway/audit_pathway_accuracy.csv`.

## Regenerate the figures

```bash
python3 tools/paper_figures.py \
  --structure docs/paper/data/structure/audit_structure_summary.json \
  --pathway docs/paper/data/pathway/audit_pathway_summary.json \
  --loo docs/paper/data/pathway/audit_pathway_loo.json \
  --presentation docs/paper/data/presentation/audit_presentation_summary.json \
  --claims docs/paper/supplement/baseline_claims.json \
  --selective docs/paper/data/selective_prediction.json \
  --out-dir docs/paper/figures
```

## If you need PNGs or flattened tables

Render any committed SVG:

```bash
python3 -c "import cairosvg, sys; cairosvg.svg2png(url=sys.argv[1], write_to=sys.argv[2], scale=2)" \
  docs/paper/figures/fig2_structure_accuracy.svg fig2.png
```

The tabular extracts were flattened views of the JSON artifacts above; read those
directly rather than reintroducing a copy.
