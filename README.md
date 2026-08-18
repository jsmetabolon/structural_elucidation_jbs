# Structural Elucidation Inference

Self-contained inference repo for the Metabolon X-compound pipeline: spectrum features in, ranked SMILES predictions out (full rerank path).

## Quick start (EC2 CPU)

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

# Large assets (encoder, decoder weights, embedding banks, NCE modelfiles)
./fetch_assets.sh --source ec2 --host user@your-ec2-host
# or from local parent workspace:
./fetch_assets.sh --source local --parent ../

# Run equalization rerank pipeline (verifier -> pool -> rerank)
chmod +x fetch_assets.sh run_inference.sh
START_AT=verifier RUN_LIBS=305 ./run_inference.sh
```

## Layout

```
structural_elucidation_jbs/
├── structural_elucidation/   # installable package (pipeline, inference, reports, nce, ...)
├── vendor/nce/               # Vendored NewCompoundElucidation package
├── vendor/modelfiles/        # NCE pkl assets (fetched)
├── models/
│   ├── decoder/              # vocab/config JSONs (shipped)
│   ├── encoder/              # best_encoder.pt (fetched)
│   ├── decoder_weights/      # .pt decoder weights (fetched)
│   └── embedding_banks/      # retrieval latents (fetched)
├── data/                     # ions, known_X, base pools, feature rows, SMILES subset
├── assets/fonts/montserrat/  # brand typeface, vendored (SIL OFL 1.1)
├── fetch_assets.sh
├── run_inference.sh
├── pyproject.toml
└── tests/test_smoke.py
```

Install the package in editable mode (optional but recommended):

```bash
pip install -e .
```

## Key module locations

| Role | Module |
|------|--------|
| SMILES inference | `structural_elucidation.inference.predict` |
| Merged inference | `structural_elucidation.inference.merged` |
| Verifier scoring | `structural_elucidation.pipeline.verify` |
| Unified pool | `structural_elucidation.pipeline.pool` |
| Rerank | `structural_elucidation.pipeline.rerank` |
| Post-rerank enrichment | `structural_elucidation.pipeline.enrichment` |
| HTML reports | `structural_elucidation.reports.compound` |

## Pipeline

`run_inference.sh` runs the equalization path per library:

1. **Library-anchored candidates** (optional if `START_AT=verifier`)
2. **Verifier scoring** via NCE `PredictedSpectrumMatcher`
3. **Unified candidate pool** merge with base pool
4. **GBM/logistic rerank** → `output/pipeline_equalization/rerank_lib*/unified_rerank_all_top10.csv`

Environment variables:

| Var | Default | Description |
|-----|---------|-------------|
| `OUT_ROOT` | `output/pipeline_equalization` | Output directory |
| `RUN_LIBS` | `209,305,400,402` | Comma-separated libs |
| `START_AT` | `verifier` | `full`, `candidates`, or `verifier` |

## Dual-lane inference (optional)

For retrieval + SELFIES de novo on new spectra:

```bash
python -m structural_elucidation.inference.merged \
  --feature_table data/ions_unknown_merged_clean.csv \
  --encoder_weights models/encoder/best_encoder.pt \
  --lib_id 400 --ion_mode positive \
  --artifact_prefix output/predict_lib400
```

Requires fetched encoder weights and embedding banks.

## Compound HTML reports

Generate a self-contained report bundle (`index.html`, `internal.html`, `guides.html`, `data/`, `assets/`) for any X compound from rerank predictions:

```bash
# Fetch metabolite synonyms DB (+ optional WGCNA cache) first:
./fetch_assets.sh --source local --parent ../

python -m structural_elucidation.reports.compound \
  --control-id X-25003 \
  --predictions output/pipeline_equalization/rerank_lib209/unified_rerank_all_top10.csv \
  --out-dir output/reports/X-25003
```

Optional flags: `--target-smiles`, `--target-name`, `--candidate-pool`, `--hypotheses`, `--peak-csv`, `--db`, `--known-x`, `--curation-ledger`, `--public-candidate-n`.

Reports degrade gracefully when `vendor/metabolite_synonyms.sqlite` or WGCNA cohort files are missing (synonym/network/correlation sections use placeholders or are omitted).

### What a report is allowed to claim

The MSI level governs how much of a structure the report shows, so the same
pipeline output produces a different page depending on what the evidence
supports:

| MSI level | Structure | Name | Public candidate list |
|---|---|---|---|
| 1-2 | shown | asserted | 1 |
| 3 | shown as a hypothesis | hedged | 3 (`--public-candidate-n`) |
| 4 | withheld | molecular formula | none |
| 5 | withheld | observed mass | none |

Internal reports keep ten candidates for triage at Levels 1-3. Human and
laboratory verdicts come from `data/curation_ledger.csv`, the only place a
reviewer's finding enters the model; a failed authentic-standard run recorded
there caps the level at 4 and withdraws the name.

The evidence rules behind that level — what counts as an isotope peak, when a
halogen may be claimed, which molecular formulas are reachable, and how a
precursor m/z becomes a neutral mass — are written down with their constants and
their current enforcement status in
[docs/identification_rules.md](docs/identification_rules.md).

## Batch reports for all X compounds

Generate one report bundle per `control_id` found in the rerank CSVs, plus a master `index.html`:

```bash
# Quick scope check (no generation): how many compounds per lib?
DRY_RUN=1 OUT_ROOT=output/pipeline_equalization ./generate_all_reports.sh

# Generate everything (resumes by default; one isolated subprocess per compound)
OUT_ROOT=output/pipeline_equalization JOBS=8 ./generate_all_reports.sh
```

Or call the module directly:

```bash
python -m structural_elucidation.reports.batch \
  --equalization-root output/pipeline_equalization \
  --out-dir output/reports --jobs 8
```

Outputs land under `output/reports/`:

- `output/reports/<control_id>/` — full per-compound bundle (same as the single-compound builder)
- `output/reports/index.html` — searchable/sortable master table linking every compound
- `output/reports/reports_index.csv`, `output/reports/batch_summary.json`

Key flags: `--jobs N` (parallelism), `--overwrite` (rebuild existing), `--limit N` (quick test), `--dry-run`, `--control-ids X-1,X-2` / `@file`, `--index-only` (rebuild master index from existing manifests), `--observed-lookup-dir` (per-lib peak CSVs for MS-fit scoring). One compound failing never aborts the batch; failures are recorded in `batch_summary.json` with a `_build.log` in the compound dir.

## Verification

```bash
pytest tests/test_smoke.py -m smoke
# or
python -m pytest tests/test_smoke.py -m smoke
```

## Expected outputs

After `./run_inference.sh`:

- `output/pipeline_equalization/rerank_lib{209,305,400,402}/unified_rerank_all_top10.csv`
- Per-lib verifier and unified pool summaries under `output/pipeline_equalization/`
