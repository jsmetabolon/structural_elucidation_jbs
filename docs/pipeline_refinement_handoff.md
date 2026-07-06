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

## 2. Hard mass/formula gate in reranking
- Goal: stop a candidate whose formula/mass disagrees with the observed feature
  from taking rank 1 (currently 15% of reports; e.g. X-12811 - a glucuronide
  176 Da heavier than the feature).
- Current: the report demotes these to Tentative + shows a "mass-consistent
  alternative candidate" callout, but cannot re-rank.
- Plug-in: in the reranker, before final ordering, apply a hard gate: if
  `abs_mass_error_da` to the observed neutral mass exceeds tolerance (e.g. > 0.02
  Da / a few ppm) AND a lower-ranked candidate passes, the failing candidate
  cannot hold rank 1.
- Steps: (a) ensure observed neutral mass/adduct is available at rerank time
  (see item 4 / feature rows); (b) add the gate; (c) regenerate rerank CSVs.
- Expected output: mass-mismatch rank-1 rate -> ~0. Verify:
  `tools/audit_presentation.py` `mass_mismatch` -> near 0; the mass-alt callout
  becomes rare.
- Open question: exact tolerance + how to treat adduct ambiguity.

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

## Cross-cutting: from-source regeneration + render_html parity task
Once items 1-4 land upstream, prefer regenerating reports from source
(`python -m structural_elucidation.reports.batch` with `--observed-lookup-dir`,
`--db`, and the pathway-enriched rerank CSVs) rather than patching, so the
delivered bundle is a single clean render.

Parity status (important): the overall-confidence badge/card + scored-parameters
table are already wired into `render_html` (via the shared `confidence.py`
builders), so those match a from-source render. The remaining gap-fix DISPLAY
elements are currently applied only by `tools/patch_report_confidence.py` and are
NOT yet in `render_html`: the WGCNA correlations table, the MS Fit MS2-evidence
panel, the resolved display name, the Metabolon Level line, and the
mass-consistent-alternative callout. Before a from-source regen will match the
delivered (patched) bundle, wire these into `compound.py::render_html` (build the
same normalized fields in `build_report` from the bundle inputs, then call the
existing `confidence.py` builders: `correlations_table_html`,
`ms2_evidence_section_html`, etc.). Until then, the enrich -> patch flow
(`tools/enrich_report_evidence.py` + `tools/patch_report_confidence.py` +
`tools/patch_consensus_crowding.py`) is the source of truth for delivery.
