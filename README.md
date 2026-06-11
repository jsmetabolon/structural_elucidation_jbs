# Spec2Mol v21 Inference

Self-contained inference repo for the v21 clean Metabolon pipeline: spectrum features in, ranked SMILES predictions out (full rerank path).

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
spec2mol_v21_inference/
├── src/                    # Python modules (predict_smiles.py, rerank, verifier, ...)
├── vendor/nce/             # Vendored NewCompoundElucidation package
├── vendor/modelfiles/      # NCE pkl assets (fetched)
├── models/
│   ├── decoder/            # vocab/config JSONs (shipped)
│   ├── encoder/            # best_encoder.pt (fetched)
│   ├── decoder_weights/    # .pt decoder weights (fetched)
│   └── embedding_banks/    # retrieval latents (fetched)
├── data/                   # ions, known_X, base pools, feature rows, SMILES subset
├── fetch_assets.sh
├── run_inference.sh
├── src/build_compound_report.py  # per-compound HTML report builder
└── tests/smoke_test.py
```

## Key renames from parent repo

| Original | This repo |
|----------|-----------|
| `predict_x209_smiles.py` | `src/predict_smiles.py` |
| `x209_retrieval_scoring.py` | `src/retrieval_scoring.py` |
| `x209_merge_policy.py` | `src/merge_policy.py` |
| `x209_phase0_utils.py` | `src/phase0_utils.py` |
| `x400_ri_audit_utils.py` | `src/ri_audit_utils.py` |
| `merged_encoder/predict_merged.py` | `src/predict_merged.py` |

## Pipeline

`run_inference.sh` runs the v21 equalization path per library:

1. **Library-anchored candidates** (optional if `START_AT=verifier`)
2. **Verifier scoring** via NCE `PredictedSpectrumMatcher`
3. **Unified candidate pool** merge with v20 base pool
4. **GBM/logistic rerank** → `output/pipeline_equalization_v21/rerank_lib*/unified_rerank_all_top10.csv`

Environment variables:

| Var | Default | Description |
|-----|---------|-------------|
| `OUT_ROOT` | `output/pipeline_equalization_v21` | Output directory |
| `RUN_LIBS` | `209,305,400,402` | Comma-separated libs |
| `START_AT` | `verifier` | `full`, `candidates`, or `verifier` |

## Dual-lane inference (optional)

For retrieval + SELFIES de novo on new spectra:

```bash
export PYTHONPATH="$PWD/src:$PYTHONPATH"
python src/predict_merged.py \
  --feature_table data/ions_unknown_merged_v21_clean.csv \
  --encoder_weights models/encoder/best_encoder.pt \
  --lib_id 400 --ion_mode positive \
  --artifact_prefix output/predict_lib400
```

Requires fetched encoder weights and embedding banks.

## Compound HTML reports

Generate a self-contained report bundle (`index.html`, `internal.html`, `guides.html`, `data/`, `assets/`) for any X compound from rerank predictions:

```bash
export PYTHONPATH="$PWD/src:$PYTHONPATH"

# Fetch metabolite synonyms DB (+ optional WGCNA cache) first:
./fetch_assets.sh --source local --parent ../

python src/build_compound_report.py \
  --control-id X-25003 \
  --predictions ../models/pipeline_equalization_v21/rerank_lib209/unified_rerank_all_top10.csv \
  --out-dir output/reports/X-25003
```

Optional flags: `--target-smiles`, `--target-name`, `--candidate-pool`, `--hypotheses`, `--peak-csv`, `--db`, `--known-x`.

Reports degrade gracefully when `vendor/metabolite_synonyms.sqlite` or WGCNA cohort files are missing (synonym/network/correlation sections use placeholders or are omitted).

## Regenerating src/ from parent

```bash
python tools/package_from_parent.py
```

Copies and renames modules from `../spec2mol_merged_fullrun_metabolon_only`.

## Verification

```bash
python tests/smoke_test.py
```

## Expected outputs

After `./run_inference.sh`:

- `output/pipeline_equalization_v21/rerank_lib{209,305,400,402}/unified_rerank_all_top10.csv`
- Per-lib verifier and unified pool summaries under `output/pipeline_equalization_v21/`
