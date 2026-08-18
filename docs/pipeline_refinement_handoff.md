# Pipeline-Refinement Handoff (Structural-Elucidation Reports)

Self-contained handoff for the next agent/engineer. The report layer
(`structural_elucidation/reports/`) and its delivered bundle now surface the best
available evidence and an honest, calibrated identification confidence. The
items below cannot be fixed in the report layer alone - they require changes in
the upstream inference/reranking pipeline (not runnable in the report VM) or new
data. Each item lists the goal, current evidence, where it plugs in, ordered
steps, expected output, verification, and open questions.

Context you need:
- Reports are one dir per `control_id` under a `reports/` tree; each carries
  `data/<slug>_*.csv/json` sidecars. The confidence model is
  `structural_elucidation/reports/confidence.py` (axes: structure 0.40, mass 0.30,
  MS2 0.20, pathway 0.10; Platt-scaled on 92 truth-verified controls; tiers
  Confirmed/High/Moderate/Low/Tentative). Evidence is folded in by
  `tools/enrich_report_evidence.py` and injected by `tools/patch_report_confidence.py`.
- Compliance constraint: observed MS2 ions are internal-only; public reports must
  never contain per-fragment/ion detail. Keep this in any pipeline output that
  reaches the public bundle.
- Current audited coverage: observed m/z 97%, MS2 present 44%, mass-mismatch 15%,
  names DB-resolved 37%, WGCNA correlations 13% (placeholder), Level 62%; tiers
  Confirmed 44 / Moderate 1416 / Low 1344 / Tentative 1093.

## 1. Real client-specific co-abundance (WGCNA) - highest impact
- Goal: replace the reference-cohort WGCNA placeholder in "Top Correlating
  Metabolites" with per-dataset measured co-abundance, and (only once real) add it
  as an independent confidence axis.
- Current: `correlations` are placeholder cohort neighbors for 13% of controls
  (labeled illustrative); not fed into confidence.
- Plug-in: produce a per-control `parent_neighbors`-style table (columns:
  `control_id, neighbor_name, neighbor_id, correlation_r, source`) for the client
  dataset; `tools/enrich_report_evidence.py::load_wgcna` already consumes that
  shape. To make it a confidence axis, extend `confidence.py::build_axes` with a
  fifth axis (re-fit the Platt constants via `tools/calibrate_confidence.py`).
- Steps: (a) compute the client co-abundance matrix -> per-control top-k real
  neighbors; (b) drop into the bundle; (c) re-enrich + re-patch; (d) if adding as
  an axis, re-calibrate + re-validate tiers.
- Expected output: real correlations for ~all controls; "placeholder" labeling
  removed. Verify: `tools/audit_presentation.py` `has_wgcna_correlations` -> ~1.0.
- Open question: correlation-strength -> confidence weight (needs truth-labeled
  validation before it can move the score).

## 2. Hard mass/formula gate in reranking - highest-payoff, quantified
- Goal: stop a candidate whose formula/mass disagrees with the observed feature
  from taking rank 1. This is the single largest recoverable-structure lever.
- Evidence (audit of the 590 rank-1 mass-mismatch reports = the entire MSI Level-4
  "mass rejected" group; `mass_match_score < 0.5`):
  - **519 / 590 (88%) already have a mass-consistent candidate at a lower rank**
    (`mass_match_score >= 0.90`). The reranker (`unified_rerank_v1`) is placing a
    mass-inconsistent candidate at #1 while the right-mass one sits at #2..N.
  - It is **not** a tolerance/strictness issue: rank-1 errors are large, not
    borderline - 308 are >= 5 Da and the **median nonzero error is ~176 Da**
    (glucuronide-scale). Deltas cluster on chemically real conjugates/adducts:
    hexose 162.053, glucuronide 176.032, sulfate 79.957, H2O 18.011, NH3, CH2, ...
    A looser mass tolerance cannot rescue a 162/176 Da gap.
  - **248 / 590 have a *blank* formula on the rank-1 SMILES**, so no mass could be
    computed at all (`mass_match = 0` by default) - a formula-assignment gap, not a
    mismatch.
  - Concrete: X-25984 rank-1 aglycone `C10H12N4O5` (err 162.05 = hexose) vs
    **rank-2 `C16H22N4O10` fitting to 3.4e-5 Da**; X-21434 rank-1 mass_match 0.0 vs
    rank-2 0.996; X-12811 a glucuronide 176 Da heavier than the feature.
- Current (report layer): these are demoted (Tentative confidence, MSI Level 4) and
  each such report now **names the specific recoverable alternative** - rank,
  formula, mass match, delta-mass, and SMILES - via the upgraded "Mass-consistent
  alternative (recoverable structure)" callout (`confidence.py` + the enriched
  `mass_alt` in `tools/enrich_report_evidence.py`). The report layer cannot re-rank.
- Plug-in (three ordered upstream fixes, by payoff):
  1. **Mass-aware rerank / hard gate**: before final ordering, if rank-1's
     `abs_mass_error_da` to the observed neutral mass exceeds tolerance (e.g. > 0.02
     Da / few ppm) AND a lower-ranked candidate passes, the failing candidate cannot
     hold rank 1. Recovers a structure for ~519 reports (~13% of all 3,897):
     MSI L4 -> at least L3 (mass-consistent hypothesis), higher where structure/MS2
     agree.
  2. **Assign a formula to every predicted SMILES** (RDKit) so the 248 blank-formula
     rank-1 cases become gate-able instead of defaulting to `mass_match = 0`.
  3. **Conjugate-aware candidate generation** for the residual ~71 with no
     right-mass candidate anywhere in the list: enumerate common conjugates
     (glucuronide / sulfate / hexose / etc.) of the predicted aglycone, since the
     observed feature is frequently a conjugate not in the current candidate set.
- Steps: (a) ensure observed neutral mass/adduct is available at rerank time
  (see item 4 / feature rows); (b) add fix 1, then 2, then 3; (c) regenerate rerank
  CSVs; (d) re-enrich + re-patch (or re-render).
- Expected output: mass-mismatch rank-1 rate -> ~0; MSI L4 count drops by up to
  ~519. Verify: `tools/audit_presentation.py` `mass_mismatch` -> near 0; the
  mass-alt callout becomes rare; `msi_level_distribution` L4 falls, L2/L3 rise.
- Open question: exact tolerance + how to treat adduct ambiguity (the same feature
  may be explainable as multiple adducts/conjugates).
- Reference prototype: `tools/prototype_mass_gate_rerank.py` (read-only) applies this
  gate to the delivered `top_candidates.csv` and confirms the projected impact:
  598 rank-1 mass-rejected today, **523 (87%) recoverable** by promoting an in-list
  mass-consistent candidate (~523 MSI L4 -> L3+), 75 with no in-list alt (need
  conjugate enumeration). Use it as the spec for the reranker change.

## 3. Broader MS2 capture + a better candidate-vs-spectrum scorer
- Goal: raise MS2 coverage above 44% and make the spectral score more
  discriminative (MS2 is the strongest independent corroboration and the main
  lever for higher confidence).
- Current: MS2 present for 44%; the report uses the pipeline
  `predicted_spectrum_score_v2` (the report's own raw-m/z matcher could not consume
  the binned observed lookups). Observed MS2 is stored as bins, not decodable m/z.
- Plug-in: (a) propagate an observed MS2 spectrum for every feature that has one;
  (b) expose raw peak m/z (internal-only) so a proper candidate-fragment match can
  run; (c) improve the scorer (fragment + neutral-loss coverage, intensity
  weighting).
- Steps: (a) fix MS2 capture/propagation upstream; (b) emit raw peaks alongside the
  binned form; (c) re-score; (d) re-enrich reports.
- Expected output: MS2 coverage up; confidence better calibrated (see item 4).
  Verify: `ms2_spectrum_present` fraction; calibration decile spread widens.
- Open question: keep raw MS2 ions strictly internal (compliance).

## 4. Expand the authentic-standard / truth set (calibration)
- Goal: sharpen the confidence calibration (currently only n=92 truth-verified
  controls; structure-correctness base rate ~48%) and let more IDs reach
  Confirmed/High.
- Current: 44 Confirmed = Level-1 controls whose top candidate matches the known
  structure; calibration Platt fit is intentionally conservative on n=92.
- Plug-in: add more confirmed structures (authentic-standard IDs) with SMILES to
  `known_X_all_available_structures_simple.csv`; `tools/calibrate_confidence.py`
  and the Level-1 Confirmed path pick them up automatically.
- Steps: (a) collect more Level-1/standard-confirmed structures; (b) re-run
  `tools/calibrate_confidence.py` to refit Platt + re-validate tiers; (c) bake the
  new constants into `confidence.py`; (d) re-enrich + re-patch.
- Expected output: tighter probability calibration; more Confirmed. Verify:
  calibration decile-bin monotonicity + per-tier observed correctness improve.
- Open question: whether to also fit per-library or per-pathway calibration once
  the labeled set is large enough.
- **Second route, no new chemistry: hold out the library against itself.** The
  steps above buy truth with acquisition lead time. The reference library is
  already thousands of compounds with both a known structure and an acquired
  spectrum, so running those spectra through the pipeline with each compound
  suppressed from retrieval makes every one of them a held-out control. Precedent
  is item-adjacent and already in the repo: `tools/audit_pathway_predictions.py
  --loo` does this for pathway labels over 5,095 reference metabolites and reports
  it as mechanism validation rather than accuracy.
  - What it needs: the library standard spectra and structures (both exist), a way
    to suppress one compound from retrieval and library matching for its own
    query without disturbing the rest, and an inference run. Unlike the pathway
    version it cannot be replayed offline from the database, because candidate
    generation needs the models and the spectra.
  - What it gives: thousands of paired comparisons instead of 92, which is what
    any question of the form "does this reranking change help" actually requires.
    It is the missing prerequisite for validating item 2 before shipping it.
  - What it does not give: a rate transferable to real unknowns. Catalogued
    compounds are a different population, and this work already measures the gap
    at 38.1% exact-source rank-1 recovery against 56.0% predicted-source. Report
    it as measured on catalogued chemistry, the way the pathway LOO is reported.
  - Sequencing: worth doing before, not after, the authentic-standard expansion.
    It is cheaper, it arrives sooner, and it tells you whether a proposed rerank
    change is worth the standards budget in the first place.
  - Costed in `docs/library_holdout_validation_scope.md`, including the three
    places a held-out compound could leak into its own answer and the one number
    that has to be established before any of it is worth scoping further: how many
    library standards have a usable acquired MS2 spectrum.

## 5. Evidence-aware headline candidate score (fold evidence_confidence into rerank)
- Goal: make the number the reranker emits as the candidate's headline score reflect
  the actual identification evidence, not just generative likelihood.
- Evidence: `candidate_report_score` (= `unified_final_score` = `rerank_score`) is a
  raw, uncalibrated generative-rerank value. Across all 3,897 rank-1 candidates its
  median is **0.035** while the pipeline's own `evidence_confidence` median is
  **0.996**; 57% of reports have `evidence_confidence >= 0.90` but a headline score
  < 0.10, and 91% of library-hit + mass-perfect rank-1s show a headline score < 0.10
  (e.g. X-02269 phenyl glucopyranoside: `evidence_confidence 0.999`, `mass_match
  0.9998`, headline score 0.036). The score is small partly because it spreads over
  many mass-consistent isomers.
- Current (report layer, stopgap): reports now surface `evidence_confidence` +
  `mass_match` in the top-10 widget and relabel the raw value "Rank score (relative)"
  (see `tools/patch_candidate_scores.py` + `compound.py`), so readers stop reading
  the 0.036 as a match quality. This does not change ranking.
- Plug-in: in the reranker, combine `evidence_confidence` (library + mass + spectral)
  with the generative score into the final ordering/score, or expose a calibrated
  0-1 `candidate_confidence` column and use it as the headline. Keep the raw
  generative score available as a secondary isomer tie-breaker.
- Caveat (do not overclaim): `evidence_confidence` is shared across mass-consistent
  isomers (all five C12H16O6 isomers of X-02269 score 0.9986), so it certifies the
  formula/library class, not the specific isomer. Retain an isomer-margin term (the
  small generative-score gap) so a confident class does not become a false-confident
  single structure - this mirrors why the report's overall confidence stays Moderate
  / MSI Level 3 for such cases.
- Verify: headline candidate score median rises and tracks correctness on the truth
  set; the report note about "relative" rank score can then be removed.
- Reference prototype: `tools/prototype_mass_gate_rerank.py` shows the headline shift
  on the delivered data (median raw `candidate_report_score` 0.035 ->
  `evidence_confidence` 0.997).

## 6. Multiply charged precursors are converted as though singly charged
- Goal: stop assigning a neutral mass that is out by a clean integer factor, and
  stop such a feature reaching the Confirmed tier.
- Evidence: `X-21470` has `observed_neutral_mass` 233.0668 and a rank-1
  `candidate_exact_mass` of 466.1331, a ratio of exactly 2.0000. The compound is
  5-androstenetriol 3,17-disulfate (`C19H30O9S2`), a disulfate, which readily
  forms `[M-2H]2-` in negative mode. The structure, its Tanimoto of 1.0 and the
  annotation are all correct; the m/z-to-neutral-mass conversion assumed one
  charge. Across the 2,471 delivered reports carrying both masses, five show an
  exact integer ratio: `X-21470`, `X-24546`, `X-24637`, `X-25814` at 2x and
  `X-25107` at 0.5x.
- Why it is worth fixing at 0.2%: **`X-21470` is MSI Level 1, the Confirmed
  tier.** A feature whose neutral mass is wrong by a factor of two carries the
  pipeline's highest-confidence label, and nothing flagged it. The same report
  also states a `mass_match_score` of 0.95 against a recorded error of 233 Da, so
  two of its own fields disagree and the one a reader sees is the reassuring one.
  The exposure class -- disulfates, diphosphates, dicarboxylates -- is common
  enough that the rate may be higher in other matrices than in this cohort.
- Cheap detector, already running: `charge_state_suspects` in
  `tools/audit_mass_policy.py` tests whether a candidate's mass is a near-exact
  integer multiple or reciprocal of the observed neutral mass (relative tolerance
  1e-4). No mass window catches these, because a factor of two is orders of
  magnitude too large to read as a tolerance problem and far too exact to be
  noise. Fixing item 2's mass gate will not address this.
- Plug-in: infer charge at feature extraction from isotope-envelope spacing
  (`1/z` Da between peaks) rather than assuming `z=1`, and carry `z` through to
  the neutral-mass conversion. Failing that, when a candidate matches at an exact
  integer multiple, re-evaluate under that charge assumption instead of scoring
  it as a huge mass error.
- Verify: `charge_state_suspects` reports zero, and no report reaches MSI Level 1
  while carrying an unexplained integer mass ratio.

## 7. `evidence_confidence` is 0-filled for generators that never emit it
- Goal: stop reporting an absent measurement as a measured zero, and stop the
  reranker penalising candidates for which the score was simply never computed.
  Code for this landed in this repo (see "What is already done" below); what
  remains is a re-run, part of which is possible on the report VM and part of
  which is not (see "Where each step can run").
- Evidence, measured on the delivered 3,897-report bundle
  (`reports/*/data/<slug>_top_candidates.csv`):
  - **3,213 of 34,306 candidate rows (9.4%) carry `evidence_confidence = 0.0`,
    504 of them at rank 1.** Whether a row is zero is decided entirely by which
    generator produced it, not by the evidence available for it.
  - Generators that always populate it: `library_anchored_candidates_v1`
    (18,433 rows, 0 zeros), `peptide_first_denovo_v1` (9,533, 0), `balanced`
    (2,584, 0), `wgcna_seeded_v1`, `mlims_anchored_candidates_v1`,
    `mass_aware_combinatorial_v1`, `mlims_strict_overlay_v1`, and the other
    targeted-chemistry tracks.
  - Generators that never populate it: the retrieval/decoder tracks
    (`ec2_bundle_retrieval_selfies_v1` 1,497 zeros, `retrieval_fallback` 924,
    `selfies_fallback` 227) and the single-operator chemistry enumerators
    (`dehydro_v1` 218, `hydroxyl_v1` 124, `methylation_v1` 113, `demethyl_v1` 49,
    `n_acetyl_v1`, `nucleoside_v1`, `glucuronide_v1`, `carnitine_acyl_v1`,
    `glycine_v1`, `peptide_first_v1`).
  - **All 527 enumerator zeros carry a nonzero `formula_prior_score`**, so the
    canonical definition would score them ~0.999. The zero is contradicted by the
    row's own other columns.
  - Mechanism: `pool_v1.dedup` coerced the score columns with
    `numeric(work, col, default=0.0)`, so a generator that omits the column is
    indistinguishable from one that measured nothing. Nothing downstream could
    tell the two apart.
  - Not an MS2-availability artefact: `X-18345` has an MS2 spectrum scoring 0.98
    and still reports 0.0, while 269 of 600 sampled compounds with no MS2 at all
    have a nonzero rank-1 score.
  - Consequence beyond display: `evidence_confidence` is a trained reranker
    feature (`rerank.FEATURE_COLUMNS`), carries 0.12 in `rerank.heuristic_score`,
    and is an F2 promotion-gate input (`chemistry/diagnostics.py`). Affected rows
    were docked up to ~0.12 of rerank score against an otherwise identical
    library-anchored candidate.
  - User-visible symptom that surfaced this: 16 truth controls at MSI Level 1
    ("Confirmed (authentic standard)", `candidate_confidence` 0.97) display an
    evidence confidence of 0.0 — e.g. `X-12127`, whose rank-1 candidate is the
    correct structure at Tanimoto 1.0 with `formula_prior_score` 1.0,
    `mass_match_score` 1.0 and zero mass error.
- Second defect found alongside it, **not** fixed: the column has two different
  definitions in-tree. `pipeline/candidates.py` uses
  `max(formula_prior, forward, ms2_match, nce_pipeline) + 0.05 * ms1_supported`,
  while `chemistry/conjugates.py` uses `0.35 + 0.35 * formula_prior + 0.15 * ms1`,
  and `chemistry/diagnostics.py` documents a third meaning (chemistry-class
  diagnostic `evidence_score` joined onto truth-blind rows). The conjugate
  formula is live in the delivered data: every nonzero `n_acetyl_v1` row reads
  0.699 = `0.35 + 0.35 * 0.997`. So the column is not on one scale even where it
  is populated. Unifying it changes existing nonzero values and therefore needs
  the same retrain, which is why it was left as a decision rather than folded in.
- What is already done in this repo (no pipeline run required to review):
  - `structural_elucidation/chemistry/evidence.py` — the canonical definition
    (`evidence_confidence`) plus `backfill_evidence_confidence`, which returns a
    completed score and a per-row basis of `emitted` / `backfilled` /
    `unavailable`. `pipeline/candidates.py` now calls the shared helper instead of
    its own inline copy.
  - `pipeline/pool_v1.py` — `apply_evidence_backfill` runs on the final pool
    (after `dedup`, so the per-term maxima it aggregates are the inputs the
    definition needs) on both the standard and `--base-unified-pool-csv` paths.
    Writes `evidence_confidence_basis` and reports `n_evidence_backfilled` /
    `n_evidence_unavailable` in `unified_candidate_pool_summary.json`. Toggle with
    `--no-evidence-backfill`.
  - `pipeline/rerank.py` — `add_features` repeats the backfill (idempotent) so a
    pool CSV written before this change cannot feed 0.0 into the model. Toggle
    with `--no-evidence-backfill`.
  - `tests/test_evidence_backfill.py` — 12 tests, pass locally without
    rdkit/sklearn.
  - Deliberate guard, do not "simplify" it away: a formula prior is only
    recomputed when a **finite target neutral mass** is available.
    `formula_prior_score` falls back to element-ratio plausibility alone when the
    mass is missing and returns 0.65 for any well-formed formula, which would
    manufacture evidence for the ~1,151 retrieval/fallback rows that carry no
    mass at all. Those rows are labelled `unavailable` and left at 0.0.
- Measured effect on the real pools, which are on the report VM under
  `../spec2mol_merged_fullrun_metabolon_only/models/` (v14 and v20, all four
  scopes). **Every zero row recovers a score; none is `unavailable`**, because the
  pool carries `target_neutral_mass` on 78.6% of rows plus three alternates that
  cover most of the rest:
  - `unified_candidate_pool_v14_lib209_negative`: 16,032 of 40,179 rows were zero
    (39.9%), all 16,032 backfilled, 1,607 controls. The pool stage ran in 9.5 s
    via `--base-unified-pool-csv` and reported `acceptance_pass: true`.
  - `unified_pool_v20_lib209_negative`: 1,836 of 10,392 zero (17.7%), all
    recovered. Largest zero source is
    `peptide_first_v1_1|peptide_first_v1_1_verifier` (1,629 rows, median new
    score 0.931), then `nce_direct_v1` (128, median 0.467) and `hydroxyl_v1`
    (66, median 0.976).
  - Other v20 scopes carry zeros in the same proportion and are ready to run:
    lib305 953/2,190, lib400 4,105/9,043, lib402 534/999.
  - A dry run over the delivered report sidecars instead recovers only 977 of
    3,213 (30.4%), because those sidecars drop `target_neutral_mass` on the
    retrieval rows. Quote the pool figure, not the sidecar one.
- Where each step can run. The pool stage and a rerank A/B are runnable on the
  report VM against the pools above, using the pipeline's own environment
  (`../spec2mol_merged_fullrun_metabolon_only/.venv/bin/python`: Python 3.10,
  pandas 1.5.3, sklearn 1.6.1, rdkit 2026.03.3, torch 1.12). The 12 tests in
  `tests/test_evidence_backfill.py` pass under that interpreter as well as under
  the repo's own Python 3.14, so the helper is safe on pandas 1.5. What is *not*
  reproducible locally is the delivered 3,897-compound bundle: the local pools
  cover 2,169 controls (v14, four scopes) and 624 (v20), and
  `run_full_pipeline_ec2.sh` preflight fails here on `merged_encoder/`,
  `decoder/`, `references/*.pkl` and the 19 per-source `candidate_metrics_all.csv`
  inputs. Published recovery figures can only be re-measured on EC2.
- Ordered steps:
  1. Local, done: read-only dry run per scope, recording recovered vs
     `unavailable` per source. If a source recovers 0, check its column names
     against `FORMULA_COLUMNS` / `NEUTRAL_MASS_COLUMNS` before adding an alias.
  2. Local: re-run `pool_v1` per scope with `--base-unified-pool-csv` and check
     `n_evidence_backfilled` > 0 and `n_evidence_unavailable` in the summary JSON.
     Done for v14 lib209; the other three scopes remain.
  3. Local A/B: `tools/ab_evidence_backfill_rerank.py` trains `rerank` twice on
     the same pool, with and without `--evidence-backfill`, and reports the
     rank-1/top-10 counts, which controls crossed Tanimoto 0.70 in either
     direction, and how much the model leans on the feature in each arm. Feed it
     the **unpopulated** pool so the arms differ only in feature values. It
     symlinks `known_X.csv` and `known_X_all_available_structures_simple.csv`
     from `data/` into the repo root for the run (`rerank.run` hardcodes them
     there) and removes the links afterwards. `data/known_X.csv` is
     byte-identical to the pipeline tree's copy; the structures file is **not**,
     so both arms deliberately use the `data/` version. Scope caveat: only 55 of
     the 58 lib209 truth controls appear in the v14 pool, so this shows
     direction, not a publishable rate.
  4. EC2: retrain `rerank` per scope on the full inputs. `evidence_confidence` is
     a model feature, so the model must be refit, not merely re-applied: LOOCV,
     the G3 gate and the per-source isotonic calibrators all regenerate.
  5. EC2: compare recovery against the published baseline — 44/92 rank 1 (47.8%),
     64.1% top 3, 73.9% top 10. Decide before running whether a drop is
     acceptable; the fix is correctness-motivated, so ranks may move both ways.
  6. EC2: re-run the F2 promotion gate; candidates that previously failed it now
     carry evidence and may promote.
  7. Regenerate reports and re-export `tools/build_exec_summary.py` and the paper
     artifacts.
  8. Display half, still open: carry `evidence_confidence_basis` into the report
     sidecars and have `tools/patch_candidate_scores.py` / `compound.py` and
     `build_exec_summary.py` render blank rather than `0.0` for `unavailable`
     rows, so an uncomputed score never again reads as a measured zero.
- Verify: no row in the regenerated pool has `evidence_confidence = 0.0` with a
  nonzero `formula_prior_score`; every `0.0` carries basis `unavailable`; the 16
  MSI Level-1 controls listed above no longer show 0.0; recovery figures are
  re-measured and the manuscript numbers updated if they moved.
- Open questions: (a) is the `candidates.py` max-based definition the one to
  standardise on, and should `conjugates.py` be migrated to it (changes existing
  nonzero values); (b) should the retrieval/fallback tracks compute a formula
  prior at generation time so they stop being structurally unscoreable; (c) does
  the F2 gate's numeric-to-tier mapping in `chemistry/diagnostics.py` still hold
  once ~1,000 more rows carry a score.

## Cross-cutting: from-source regeneration + render_html parity task
Once items 1-4 land upstream, prefer regenerating reports from source
(`python -m structural_elucidation.reports.batch` with `--observed-lookup-dir`,
`--db`, and the pathway-enriched rerank CSVs) rather than patching, so the
delivered bundle is a single clean render.

Parity status (updated). The following now render identically from source
(`compound.py::render_html`/`build_report`) and from the enrich->patch flow, verified
on lib209 by `tools/diff_source_vs_patched.py` (public MS2-ion leaks = 0):
the overall-confidence badge/card, scored-parameters table, MSI badge/card/line,
the top-10 candidate-widget evidence (`evidence_confidence` + `mass_match` +
relabeled "Rank score (relative)"), the Metabolon Level line, the
mass-consistent-alternative callout, the internal MS2-evidence panel
(`ms2_evidence_section_html`), and the master-index MSI column + column sorting
(`index_render.py`). `mass_consistent_alt` is shared in `confidence.py`.

Two DATA-SOURCE gaps remain before a from-source regen can replace the patch chain
(the render code is ready; the inputs differ):
1. Observed-MS2 source mismatch: the two paths key MS2 presence off *different*
   observed files. `build_report` (via `--peak-csv`, auto-derived from
   `--observed-lookup-dir`) uses the per-lib `*_observed_ms2_lookup.csv`, while the
   enrich flow uses `per_control_all_x.csv`. They disagree for ~47% of lib209
   controls (e.g. X-11640 is in the per-lib lookup but not in `per_control_all_x`,
   so from-source shows the MS2 panel while the patched report does not). Fix: pick
   one canonical observed-MS2 source + presence criterion and use it for both the
   MS2 panel and the MS2 confidence axis in both paths.
2. `mass_match_score`: the current `_bundle/rerank/*.csv` lacks this column, so
   `mass_consistent_alt` cannot fire from source (the enrich flow had it from the
   originally-generated `candidate_score_summary.csv`). Either emit
   `mass_match_score` in the rerank CSV or derive it from `abs_mass_error_da` + the
   observed neutral mass in `compact_candidate_scores`.

A single-compound and a 93-compound (lib209) from-source build both succeed in-VM
(~1 s/compound; ~15-20 min for all 3,897 at `--jobs 8`). DECISION: do not switch
delivery yet - close the two gaps above and re-run `tools/diff_source_vs_patched.py`
to confirm full parity first. Until then, the enrich -> patch flow
(`tools/enrich_report_evidence.py`, `tools/patch_report_confidence.py`,
`tools/patch_candidate_scores.py`, `tools/patch_consensus_crowding.py`, and
`tools/rebuild_master_index.py`) remains the source of truth for delivery.
