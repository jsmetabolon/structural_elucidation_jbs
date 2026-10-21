# Untargeted-pipeline integration

This document describes how the structural elucidation stack is attached to the
[`untargeted-pipeline`](https://github.com/MetabolonInternal/untargeted-pipeline)
so that it runs automatically on features that were **not** assigned a
Metabolon/CURATION identity but do have consensus MS2 evidence.

The consolidated integration lives on `codex/structural-elucidation-latest`.
It combines the October development baseline (`a8ea4d5`) with the untargeted
adapter and the latest publication-level, identity-confidence, and MSI tooling.
Use the exact companion commit documented by the untargeted pipeline checkout
when reproducing a run; the driver never fetches or switches branches.

## Overview

The two repositories intentionally do **not** share a Python environment
(elucidation pins `torch==1.12.0` / `pandas<2`; the untargeted pipeline needs
`pandas>=2` / Python 3.11+). Integration therefore runs elucidation
**out-of-process** in its own virtual environment:

```mermaid
flowchart LR
  A["untargeted_workflow\nStructuralElucidationStep\n(pandas 2.x venv)"] -->|handoff.json| B["run_novel.py\n(elucidation venv)"]
  B --> C["adapters.untargeted_inputs\nsynthesize per-lib tables"]
  C --> D["candidates -> verify -> pool -> rerank -> scoring -> reports"]
  D -->|novel_elucidation_summary.csv| A
```

The untargeted step writes a neutral **handoff JSON**, invokes
`python -m structural_elucidation.pipeline.run_novel` in this repo's venv, then
merges the resulting `novel_elucidation_summary.csv` back into its feature table.

## Entry point

```bash
python -m structural_elucidation.pipeline.run_novel \
    --handoff /path/to/handoff.json \
    --out-root output/novel_run \
    --jobs 8
```

Key flags:

| Flag | Meaning |
| --- | --- |
| `--handoff` | Handoff JSON produced by the untargeted step (required). |
| `--out-root` | Root for all outputs (required). |
| `--inputs-dir` | Where synthesized input tables go (default `<out-root>/inputs`). |
| `--assets-root` | Root holding fetched model/reference assets (default: repo root). |
| `--stages` | Subset of `candidates,verify,pool,rerank,reports` (default all). Scoring always runs for available predictions, including when reports are omitted. |
| `--libs` | Optional comma-separated lib filter (e.g. `400,402`). |
| `--enable-denovo` | Also run encoder de-novo/retrieval inference (needs encoder assets). |
| `--top-k` | Predictions kept per feature (default 10). |
| `--dry-run` | Synthesize inputs and print the planned commands only. |

## Handoff schema

```json
{
  "schema_version": 1,
  "source": "untargeted_workflow",
  "features": [
    {
      "feature_id": 123,
      "lib_id": 400,
      "ion_mode": "positive",
      "parent_mz": 250.09295,
      "ri": 2486.0,
      "peaks": [[67.0, 0.5], [137.05, 1.0]],
      "compound_name": "",
      "prior_match_level": "MS2",
      "prior_matched_name": "some_lib_hit",
      "neutral_mass": 278.1517,
      "adduct_type": "[M+Na]+"
    }
  ]
}
```

`neutral_mass` and `adduct_type` are optional but recommended: when the upstream
pipeline has resolved the adduct (e.g. via OpenMS peak grouping), pass the
adduct-corrected neutral mass so candidate mass matching is anchored correctly
even when the representative (base-peak) ion is not `[M+H]+`/`[M-H]-`. When
`neutral_mass` is omitted the adapter falls back to the naive `[M+/-H]`
assumption from `parent_mz`.

Only libraries `209, 305, 400, 402` are supported (others are skipped). Each
feature is given a synthetic control id `X-<900000000 + feature_id>` that matches
the `X-\d+` pattern the pool/verify/rerank stages require while never colliding
with real Metabolon X-compounds.

## What `run_novel` does that `run_inference.sh` does not

`run_inference.sh` reuses pre-built `data/base_pools/` for de-novo/retrieval
candidates and assumes registered Metabolon X-compounds. Brand-new untargeted
features have neither, so `run_novel`:

1. Synthesizes `feature_rows`, `neutral_mass_hypotheses`, `gap_manifest`, and a
   bit-exact `observed_ms2_lookup` (NCE ion/loss bins) via
   `structural_elucidation.adapters.untargeted_inputs`.
2. Seeds a synthetic `all_x_catalog.csv` with one dummy `S0` row per library so
   `pipeline.pool` actually writes `unified_candidates_per_control_non_s0.csv`
   (otherwise the split comes out empty and un-truthed novel controls never reach
   the reranker).
3. Tolerates the verifier's non-zero acceptance exit code (expected when there
   are no known winners) as long as `scored_candidates.csv` is produced.
4. Optionally runs `inference.merged` for de-novo/retrieval candidates.
5. Joins `rerank_lib*/unified_rerank_all_top10.csv` back to the originating
   `feature_id` in `novel_elucidation_summary.csv`.
6. Scores each feature before rendering and saves canonical JSON. HTML consumes
   those same confidence fields. A report failure cannot discard completed scores.

Only features in the current handoff are summarized. Features without predicted
structures have no invented MSI or confidence values. Prior putative MS1/MS2
names are provenance; they cannot establish authentic-standard confirmation.

## Outputs

Under `--out-root`:

- `inputs/` — synthesized per-library input tables + `novel_feature_map.csv`.
- `library_anchored_lib*/`, `verifier_lib*/`, `unified_pool_lib*/`,
  `rerank_lib*/unified_rerank_all_top10.csv` — standard stage artifacts.
- `reports/` — per-compound HTML bundles + master `index.html` (when reports run).
- `scores/<control_id>.json` — canonical evidence fields, selected and raw
  candidates, naming decisions, confidence values, publication calls, and errors.
- `novel_elucidation_summary.csv` — one row per elucidated feature with
  `feature_id, control_id, lib_id, se_rank1_smiles, se_rank1_score,
  se_topk_smiles_json, se_report_path`.
- `run_novel_summary.json` — summary schema version 2, per-scope stage status,
  companion Git commit and dirty state, assets root, calibration SHA256, scoring
  errors, and report errors. Source archives without Git record a null commit.

The summary adds the following columns without changing the handoff schema:

| Columns | Meaning |
| --- | --- |
| `se_msi_level`, `se_msi_label`, `se_msi_rationale` | Evidence-based MSI level 1–5, label, and reasoning. |
| `se_confidence`, `se_confidence_tier` | Report headline confidence (0–1) and tier; can include explicit laboratory/curator assurance. |
| `se_evidence_score` | Model-only evidence score (0–1), kept separate from assurance. |
| `se_reported_level`, `se_reported_identification` | Publication call and level; the level is text and preserves `2*`. |
| `se_candidate_identity_confidence` | Candidate identity score (0–100), using the shipped publication calibration. |
| `se_confidence_flags_json` | JSON list of unavailable evidence and scoring limitations. |

Evidence MSI and publication levels use their respective existing rules and can
differ. These columns also differ from reranking scores and the candidate
`evidence_confidence` match score. Pure NCE scoring is preferred for the candidate
identity margin: a successful zero stays zero; a model error or absent spectrum
is unavailable and uses the existing evidence-score fallback. Calibration tiers
without sufficient labels retain their recorded prior-based treatment. This
integration does not refit the supplied calibration or promote exploratory scores
to externally validated probabilities.

The standalone `tools/build_exec_summary.py`, `tools/build_exec_summary_workbook.py`,
and `tools/build_msi_level_figure.py` retain their existing CLIs. The optional
`generous` chart interpretation does not replace the pipeline's evidence MSI.

## Environment setup

`run_novel` (and every stage it spawns) runs in this repo's venv:

```bash
cd structural_elucidation_jbs
python3.10 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt && pip install -e .
```

The pinned `torch==1.12.0` stack uses **Python 3.10**; the untargeted pipeline
uses its separate Python 3.11+ environment. Pass `--se-python` explicitly when
using an existing elucidation environment outside this checkout.

## Required assets

Fetch the model/reference assets into their canonical repo locations with
`fetch_assets.sh` (see its `--help`), or point `--assets-root` at a tree with the
same layout.

`--assets-root` routes candidate reference tables, NCE models, and the optional
report reference database. Source code always comes from this checkout. The
small static adduct ladder and identity-confidence calibration are committed
with the source; model checkpoints and cohort outputs are external assets.

Required for the default `candidates,verify,pool,rerank` path:

- `vendor/modelfiles/` — NCE predicted-spectrum models (`smiles2ms2_models_{lib}*.pkl`).
  Without them the verifier falls back to a portable RDKit fragment-overlap score.
- `vendor/nce/elucidation/table_builders/hmdb_minimal.csv` (+ optional
  `HMDB_36_classyfire_21_annotations.csv`) — library-anchored candidate search.
- `data/nce_full_smiles_candidates.csv` — NCE snapshot candidates (optional).

Additionally required only when `--enable-denovo` is set:

- `models/encoder/best_encoder.pt`, `models/decoder_weights/*.pt`,
  `models/embedding_banks/` and the `SMILES_DB/` retrieval tables.

Optional:

- `vendor/metabolite_synonyms.sqlite` — richer pathway/synonym report detail.
- Open Babel (`obabel`) — nicer 2D structure rendering in reports.
