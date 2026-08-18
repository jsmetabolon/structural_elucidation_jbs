# X-Compound Structure Prediction Pipeline

This repository contains the Metabolon X-compound structure prediction and reporting pipeline. It starts from discovery-derived feature tables and observed MS2 evidence for unknowns such as `X-25003`, proposes and reranks candidate molecular structures, adds biochemical context, and produces evidence-tagged reports.

The tool is best understood as an augmentation layer around the normal Metabolon Curate/Sally and structural-confirmation workflows. It does not replace production compound calls, authentic-standard confirmation, or expert chemistry review. Instead, it expands the candidate space beyond the existing library and makes the evidence, uncertainty, and unresolved ambiguity visible.

The pipeline is modular and resumable. Candidate generation, spectral verification, pool assembly, reranking, post-rerank enrichment, report generation, and audits each write explicit artifacts. A full run is designed for the EC2 elucidation environment, while rerank-only and report-only workflows can reuse existing intermediate files.

## What This Adds

Compared with library-only lookup and manual structure review, this project adds:

- Two complementary prediction lanes: mass-gated library/analogue retrieval and constrained de novo SMILES generation.
- Mass-, adduct-, and optional retention-index-aware filtering of candidate structures.
- Library-anchored and conjugate-aware candidates for formula and biotransformation rescue.
- Candidate-to-spectrum verification using NCE predicted-spectrum scoring and observed MS2 evidence.
- A unified cross-source pool and GBM/logistic reranker that returns the top 10 structures per X-compound.
- Post-rerank super-pathway/subpathway, lipid-class, and peptide-composition enrichment.
- Evidence-aware reports that separate raw rank score, mass fit, spectral support, biochemical context, and missing evidence.
- MSI-style reporting levels that present a structure as a hypothesis or withhold it when support is insufficient.
- Public and internal report variants, plus a sortable index across all X-compounds.

The intended result is not an automatic final identity. It is a ranked, traceable set of hypotheses that makes the next curation, standard-acquisition, or synthesis decision more efficient.

## Why This Is Complementary

Curate/Sally remains the production path for released compound identities and statuses. It uses library mass/retention-index windows, curated standards, calling rules, and expert review. The X-compound pipeline asks a different question: given the available MS1/MS2 evidence for an unresolved feature, which structures are plausible, how strongly is each supported, and what biochemical context is consistent with the proposal?

The main benefits are:

- Broader candidate coverage: de novo generation can propose structures outside existing spectral libraries.
- Continuity with Metabolon data: outputs retain X-compound, feature, library, ion-mode, and candidate provenance.
- Evidence-tagged interpretation: retrieval, generation, mass, formula, MS2, pathway, and source evidence remain separate.
- Rank-aware prioritization: reports expose the top 10 candidates rather than only a single opaque call.
- Controlled claims: MSI-style levels distinguish confirmed, probable, tentative, formula-only, and mass-only outputs.
- Graceful degradation: missing MS2, database, or WGCNA evidence is shown explicitly instead of being silently imputed.
- Direct comparison with known references: exact reference lookup is labeled as provenance, not authentic-standard proof.

In short, Curate/Sally provides the production interpretation, while this tool provides a broader structure-hypothesis and evidence-review layer for unresolved X-compounds.

## Relationship To Curate/Sally

The normal Curate/Sally flow is optimized for final production calls. The X-compound pipeline is optimized for discovery and prioritization. It starts from feature-level evidence prepared by the existing discovery process rather than from raw vendor files.

When the rank-1 candidate exactly matches a reference-database structure, the report can reuse its name and biochemical taxonomy. This is a direct lookup and is labeled as such; it is not treated as independent proof that the candidate is correct for the X-compound.

For non-reference structures, the pipeline uses retrieval/de novo candidate generation, observed-mass and MS2 scoring, and neighbor-based biochemical context. The result remains a hypothesis until it is independently confirmed through curation and, where appropriate, analysis of an authentic standard.

The report layer therefore complements production curation by showing what was predicted, why it ranked highly, what evidence is missing, and what alternative structures remain plausible.

## Prediction And Evidence Strategy

Each X-compound moves through several distinct evidence layers:

- `predicted_smiles` and `rank`: the candidate structure and its final ordering.
- `candidate_report_score`: the raw relative rerank score used to order candidates; it is not a probability of correctness.
- `evidence_confidence`: a pipeline evidence summary based on library, mass, and spectral support.
- `mass_match_score`, formula fields, and mass error: consistency with the observed neutral-mass hypothesis.
- Predicted-spectrum, forward-spectrum, NCE, and MS2-match scores: candidate support from fragmentation evidence.
- Candidate source fields: retrieval, library anchor, de novo generation, peptide, conjugate, or biotransformation origin.
- Pathway source and pathway-confidence tier: direct reference lookup or neighbor-consensus biochemical context.
- Identification-confidence tier and MSI-style level: report-layer summaries controlling how the structural claim is presented.

MS2 and authentic-standard evidence are more specific than formula or library proximity. A high raw rank score or exact database source must not be read as confirmation by itself. The report keeps these quantities separate so reviewers can inspect the complete evidence stack.

The delivered reports require a specific caveat on the mass axis. Their
`mass_match_score` values were not produced by the 0.001–0.005 Da function
formerly documented in `inference.predict`. Most were produced by
`chemistry.conjugates.linear_mass_match_score` with a 1.0 Da default, while some
candidate generators used 5.0 or 0.005 Da. These scales differ by as much as
1,000-fold, so delivered mass scores are internally consistent within a
generator but are not directly comparable across all candidate sources. The
raw `abs_mass_error_da` is the safer quantity to inspect. The discrepancy,
measured impact, and recommended correction are recorded in
`docs/mass_axis_decisions_needed.md` and
`docs/paper/data/thresholds/audit_mass_policy.json`.

## Installation

Clone the repository, create an environment, install dependencies, and fetch the large model/data assets:

```bash
git clone git@github.com:jsmetabolon/structural_elucidation_jbs.git
cd structural_elucidation_jbs
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
pip install -e .
./fetch_assets.sh --source local --parent ../
```

The full EC2 workflow additionally expects encoder/decoder weights, embedding banks, NCE model files, observed-MS2 lookups, base candidate pools, integrated per-library candidates, and the metabolite reference database. `run_full_pipeline_ec2.sh` performs a preflight check before running.

## Basic Usage

Check the full EC2 environment without starting an expensive run:

```bash
PREFLIGHT_ONLY=1 \
OUT_ROOT=models/pipeline_equalization_$(date -u +%Y%m%d) \
./run_full_pipeline_ec2.sh
```

Run the full prediction, equalization, and enrichment workflow into a fresh output directory:

```bash
RUN_ID=$(date -u +%Y%m%d_%H%M%S)
RUN_ID="$RUN_ID" \
OUT_ROOT="models/pipeline_equalization_${RUN_ID}" \
RUN_LIBS=209,305,400,402 \
./run_full_pipeline_ec2.sh
```

Reuse existing candidates and run only verifier, pool, and rerank stages:

```bash
START_AT=verifier \
RUN_LIBS=209,305,400,402 \
OUT_ROOT=output/pipeline_equalization \
./run_inference.sh
```

Generate a self-contained report for one X-compound:

```bash
python -m structural_elucidation.reports.compound \
  --control-id X-25003 \
  --predictions output/pipeline_equalization/rerank_lib209/unified_rerank_all_top10_with_pathways.csv \
  --peak-csv data/observed_ms2_lookups/lib209_negative_observed_ms2_lookup.csv \
  --out-dir output/reports/X-25003
```

Generate all per-compound reports and the master index:

```bash
python -m structural_elucidation.reports.batch \
  --equalization-root output/pipeline_equalization \
  --predictions-name unified_rerank_all_top10_with_pathways.csv \
  --observed-lookup-dir data/observed_ms2_lookups \
  --out-dir output/reports \
  --jobs 8
```

## Pipeline Flow

The complete conceptual flow is:

1. Build or load per-library X-compound feature rows, precursor/adduct hypotheses, and observed MS2 lookups.
2. Encode the feature spectrum and run Track A mass-gated retrieval against the embedding bank.
3. Run Track B constrained de novo decoding, with mass/adduct and optional retention-index conditioning.
4. Add library-anchored, peptide, conjugate, and biotransformation candidates where supported.
5. Score candidates against observed MS2 using the NCE/predicted-spectrum verifier.
6. Merge all candidate sources with the base unified pool.
7. Train/apply the unified reranker and retain the top 10 candidates per control.
8. Add direct-reference or neighbor-consensus biochemical taxonomy, lipid class, and peptide composition.
9. Build per-compound public/internal reports, structure graphics, candidate consensus maps, and evidence sidecars.
10. Audit structural recovery, pathway provenance, evidence coverage, MSI presentation, and public-data compliance.

Each stage writes explicit CSV/JSON outputs, allowing expensive steps to be reused and downstream stages to be rerun independently.

## Inputs

The main inputs are:

- Feature tables with `control_id`, library, ion mode, precursor m/z, retention index, and binned MS1/MS2 evidence.
- Encoder and decoder checkpoints plus SELFIES vocabularies.
- Retrieval embedding banks and structure/SMILES lookup tables.
- Per-library base candidate pools and integrated candidate tables.
- Observed-MS2 lookup tables and neutral-mass/adduct hypotheses.
- `known_X.csv` and known-structure truth tables for evaluation and report metadata.
- Metabolite reference/synonym SQLite database for names, structures, and biochemical taxonomy.
- Optional WGCNA/co-abundance evidence for cohort-specific biological context.

Runs should use immutable input snapshots and a fresh output directory. The model/data hashes, software environment, random seed, and hardware should be recorded with production runs.

## Outputs

The primary prediction outputs are:

- `rerank_lib*/unified_rerank_all_top10.csv`: final ranked candidate structures and evidence fields.
- `rerank_lib*/unified_rerank_all_top10_with_pathways.csv`: candidate table with biochemical taxonomy enrichment.
- `verifier_lib*/scored_candidates.csv`: candidate-level spectral verification.
- `unified_pool_lib*/unified_candidates_per_control.csv`: merged cross-source candidate pool.
- Per-stage summary JSON/Markdown diagnostics and model-selection reports.

The report outputs are:

- `output/reports/<control_id>/index.html`: public per-compound report.
- `output/reports/<control_id>/internal.html`: internal report with aggregate MS2 detail.
- `output/reports/<control_id>/data/`: candidate, evidence, pathway, network, and audit sidecars.
- `output/reports/<control_id>/assets/`: structure and similarity graphics.
- `output/reports/index.html`: searchable/sortable master report index.
- `output/reports/reports_index.csv` and `batch_summary.json`: flat inventory and build status.

The delivered-set audit and summary outputs are:

- `docs/paper/data/exec_summary/exec_summary_top10.csv`: all top-ten candidate
  rows for 3,897 X-compounds, including rank, candidate scores, mass error,
  truth basis, Tanimoto, and an indication of which candidate best matched the
  available answer.
- `docs/paper/data/exec_summary/exec_summary_by_compound.csv`: one-row
  compound-level rollup of the same export.
- `docs/paper/data/parameters/evidence_parameters.csv`: the scored evidence
  axes, weights, strength bands, definitions, report language, and
  consequences (Table S10).
- `docs/paper/data/thresholds/audit_mass_policy.json`: candidate-level replay
  of the delivered mass scorers and proposed mass-policy interventions.
- `docs/paper/data/feces_pathway/feces_pathway_update_weighted.xlsx`: expected
  metabolites detected by pathway and subpathway, with known and
  pathway-annotated X-compound contributions on separate workbook tabs.

## How To Interpret Results

All structural outputs should be treated according to their evidence:

- Rank 1 is the current preferred hypothesis, not automatically the true structure.
- Top 3 and top 10 capture useful alternatives, including isomers that may be indistinguishable without stronger MS2 or retention evidence.
- The raw rank score is a relative ordering value, not a match probability.
- The aggregate identification-confidence score is exploratory; leave-one-out calibration is weak.
- Delivered `mass_match_score` values come from mixed tolerances and must not
  be interpreted as one common instrumental-accuracy scale; use
  `abs_mass_error_da` and the mass-policy audit for comparisons.
- MSI-style L1–L3 reports may show a structure, while L4–L5 deliberately withhold a structure and report formula/mass only.
- Direct pathway lookup is taxonomy consistency by construction; neighbor-consensus assignments should be interpreted using their pathway-confidence tier and provenance.
- Public reports omit observed fragment-ion detail; internal reports retain aggregate spectral evidence.

## Current Validation Snapshot

On the available 92 truth-verified controls, structural recovery at Tanimoto ≥0.70 is:

- Rank 1: 44/92 = 47.8% (95% CI 37.9–57.9%).
- Top 3: 59/92 = 64.1% (95% CI 54.0–73.2%).
- Top 10: 68/92 = 73.9% (95% CI 64.1–81.8%).

Across the delivered 3,897-report set, 2,726 reports (70.0%) present a structure at MSI-style L1–L3 and 1,171 (30.0%) withhold one at L4–L5. The exploratory confidence model has leave-one-out Brier score 0.255 and should not be interpreted as an externally calibrated probability.

The complete top-ten export contains 34,296 candidate rows. It reproduces the
same truth-control recovery figures above. Among the 163 Level 2C compounds,
only 26 have a reference formula; 23 of those 26 have a top-ten candidate with
the same formula. The remaining 3,642 compounds have neither structural truth
nor an assessable Level 2C formula and are left blank rather than counted as
false.

The mass-policy audit found that the delivered candidate sources were scored
under incompatible mass tolerances. Correct candidates have median mass error
1.07 ppm, and 81.4% are within 5 ppm; correct-versus-incorrect separation is
about 22 percentage points at 2–5 ppm but reverses at the shipped 1.0 Da
window. Replaying a hard mass gate does not improve rank-1 recovery, however:
the useful measurement and the proposed intervention are different questions.
The recommended next step is to unify candidate sources onto one
instrument-scale scorer, then rerun and recalibrate rather than tuning the
current mixed-scale gate.

For biochemical taxonomy, 1,459 exact-reference cases are direct lookups and are not independent accuracy tests. Among 644 neighbor-consensus assignments with independent labels, agreement is 52.2% for super-pathway and 18.8% for subpathway. Read those against the annotation ceiling rather than against a perfect score: on 267 controls whose structure is identical to a reference entry, where the label is copied and the assignment is right by construction, the reference database and the independent labels agree only 62.5% at super-pathway and 34.8% at subpathway. Intrinsic leave-one-out validation of the neighbor-vote mechanism is 86.7% super-pathway and 71.1% subpathway, but does not include all production priors and overrides.

The feces pathway update combines 5,189 named Global-panel metabolites with
705 X-compounds in the feces detection history (5,894 records in the full
comparison universe). Of those, 2,265 named metabolites are detected in feces
and 704 X-compounds have a report and usable pathway annotation, so 2,969
records drive the weighted averages. Each pathway value is the sum of
per-metabolite detection rates, rounded to the nearest integer. For example,
Peptide rises from 42 expected known-metabolite detections to 134 after adding
the pathway-annotated X-compounds; at subpathway level, Dipeptide rises 27 to 57
and Polypeptide 3 to 49. These are expected detection counts, not compositional
percentages.

The authoritative P0 and P1 cohort source files are separate from this 92-control set and must be evaluated independently when available.

## Documentation

- `README.md`: installation, module layout, quick start, outputs, and report commands.
- `docs/compound_report_workflow.md`: report-builder internals, evidence sections, networks, and batch generation.
- `docs/pipeline_refinement_handoff.md`: current upstream limitations and prioritized refinements.
- `docs/mass_axis_decisions_needed.md`: delivered mass-scorer discrepancy,
  measured ppm behavior, and decisions for correction and validation.
- `docs/exec_summary_and_feces_pathway.md`: complete top-ten export and
  known-plus-X feces pathway/subpathway detection analysis.
- `docs/paper/data/parameters/evidence_parameters.csv`: evidence ranges,
  definitions, report language, and consequences.
- `docs/paper/manuscript.md`: scientific validation, limitations, and interpretation.
- `docs/paper/supplement/supplementary_information.md`: cohort definitions, detailed tables, and provenance gaps.
- `docs/pcm_prediction_handoff.md`: constrained prediction workflow for Partially Characterized Molecules.
- `run_full_pipeline_ec2.sh`: complete EC2 preflight and orchestration.
- `run_inference.sh`: verifier/pool/rerank orchestration by library.

