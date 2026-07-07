# Pipeline-Refinement Handoff (Structural-Elucidation Reports)

Self-contained handoff for the next agent/engineer. The report layer
(`structural_elucidation/reports/`) and its delivered bundle now surface the best
available evidence and an honest, calibrated identification confidence. The four
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
1. Observed MS2: the from-source MS2 panel + the MS2 confidence axis must be fed the
   observed spectrum (run `batch` with `--peak-csv` / the observed lookup). Without
   it, `build_report` cannot know observed-MS2 presence, so the panel is
   (correctly) omitted whereas the enrich flow shows it from
   `observed_ms2_lookups`.
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
