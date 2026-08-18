---
name: MSI sub-levels and their confidence bands
overview: 'Ship the part of the Metabolon-style sub-level table that the delivered evidence supports - the curated SubLevel surfaced rather than derived, a generated level-by-sub-level count table, and a definition CSV that records what each letter requires - and hold the score-derived 2a/2b/2c banding until the score it would be cut from is monotone with correctness. Neither candidate score is today: evidence_confidence is flat above 0.98 and candidate_confidence is inverted (AUC 0.378) on its own output.'
todos:
  - id: definition-csv
    content: 'Write docs/paper/data/identification/msi_sublevel_definitions.csv recording each level and sub-level, the evidence it requires, whether that evidence is computable in this repo today, and the source it comes from - definitions only, no assignment'
    status: pending
  - id: counts-generator
    content: 'Write tools/build_msi_level_counts.py to generate docs/paper/data/msi_level_counts.csv from data/known_X.csv plus the delivered exec summary, replacing the hand-committed copy on cursor/evidence-confidence-calculation-012b, and add it to the CI regenerate-and-diff step'
    status: pending
  - id: surface-curated
    content: 'Surface the curated SubLevel on the report card and master index as a labeled cross-reference next to metabolon_level, withheld for the 1,403 compounds that have none rather than defaulted'
    status: pending
  - id: availability-flag
    content: 'Replace the "Confidence Score" reading of evidence_confidence in the matrix workbook with the one distinction it actually makes - evidence present versus absent at the zero boundary - and relabel the column accordingly'
    status: pending
  - id: recalibrate-platt
    content: 'Re-fit the Platt (A, B) against the retention axis on a bundle that carries axis-level columns, so the evidence score spans a usable range and TIER_HIGH stops being unreachable; measure the tier and MSI distribution shift before adopting'
    status: pending
  - id: monotone-score
    content: 'Establish a per-candidate probability that is monotone with correctness out-of-sample - candidate the isotonic rerank_score_calibrated from cursor/candidate-confidence-rework-f8a2 - and gate it on the sub-level band audit reporting ordered bands'
    status: pending
  - id: derive-bands
    content: 'Only once monotone-score passes: choose cutoffs from the calibrated score by equal-precision bands rather than round numbers, and assign MSI Sub-level from them'
    status: pending
  - id: evidence-handoff
    content: 'Write the upstream handoff for the three inputs the Schymanski sub-letters actually name: a reference MS2 spectral library, observed retention time per X-compound, and candidate-versus-observed diagnostic fragment matching'
    status: pending
isProject: false
---
# MSI sub-levels and their confidence bands

Status: proposal, not started. The measurement it rests on *is* done and
committed — [tools/audit_sublevel_bands.py](../tools/audit_sublevel_bands.py),
regenerating [docs/paper/data/identification/audit_sublevel_bands.json](../docs/paper/data/identification/audit_sublevel_bands.json),
with the conclusions written up in
[docs/msi_sublevel_band_findings.md](../docs/msi_sublevel_band_findings.md).
Read that findings document first; this plan assumes its four results and does
not re-derive them.

## Goal

A Metabolon-style evidence table wants MSI levels 1–5 split into sub-levels
(`2a`/`2b`/`2c`, `3a`/`3b`/`3c`, `5a`) with an associated confidence range on
each. Today this repo has neither the sub-levels nor the ranges: sub-levels exist
only as Metabolon's curated `SubLevel` metadata read from `data/known_X.csv`, and
the pipeline's MSI level and the report's score tiers are separate mechanisms
that never meet.

The goal is to ship the sub-level table with each cell sourced from something
that supports it, and to leave visibly empty the cells that nothing does. Not to
populate every cell. A sub-level column filled by banding a score that does not
track correctness is worse than an absent one, because it reads as a graded
claim and the grading runs backwards.

## The decision this unblocks, and why it could not be answered as asked

The request was for cutoffs and for which score they apply to. Neither can be
chosen yet, for reasons that are properties of the scores rather than of any
particular threshold:

- **`evidence_confidence`** (the matrix *Confidence Score*) is flat where bands
  would go. 62.4% of compounds sit at or above 0.98; across its top seven
  deciles P(rank-1 correct) moves only between 13.6% and 19.1%.
- **`candidate_confidence`** (the report overall) is inverted once the 44
  policy-pinned Level 1 reports are held out — AUC 0.378 against truth, lowest
  decile most accurate. Those 44 are correct by construction and carry a display
  value rather than model output, which is why the same score reads AUC 0.7865
  when they are left in.
- **`TIER_HIGH = 0.65` is unreachable.** Platt scaling bounds the evidence score
  to `[0.1794, 0.5907]`, so a top band cut at 0.65 holds only reports that are
  already Level 1, and `2a`/`3a` come out empty at any cutoff above 0.5907.
- **`2a`/`2b` are not score bands in the cited scheme.** In Schymanski 2014 they
  are evidence-type distinctions — library spectrum match versus diagnostic
  evidence — and this repo can compute neither.

So the plan splits in two. Stages 1–3 ship now and do not depend on the decision.
Stages 4–6 are what makes the decision answerable, in order. Stage 7 is the
original request, gated on stage 5 passing.

## What this does not plan

**Changing `assign_msi_level`.** The coarse level rule stays exactly as it is.
Every stage here adds a sub-level *alongside* the level; none of them may move a
compound between levels. The audit in
[tools/audit_thresholds.py](../tools/audit_thresholds.py) owns the coarse
boundaries and a change there is a different piece of work.

**Overwriting `evidence_score`.** `confidence.py` deliberately keeps two numbers
— the frozen model output and the headline that folds in human verdicts (see the
`Confidence` docstring). A sub-level must never write back into `evidence_score`,
or every downstream audit becomes a statement about policy instead of evidence.

**Re-opening the mass policy.** Owned by
[plans/mass_policy_impact_measurement.plan.md](mass_policy_impact_measurement.plan.md).

## Constraints the receiving agent needs

1. **Tools are standard library only, read-only, and reproducible from committed
   inputs.** CI runs the audit tools and then `git diff --exit-code`
   (`.github/workflows/publication-checks.yml`), so any new generator must be
   deterministic — sort keys, round floats, never embed a timestamp.
2. **Two work streams live on unmerged branches** and stages 4 and 5 depend on
   them. `candidate_confidence` is *read* by `tools/build_exec_summary.py:242`
   on `main` but produced by nothing on `main`; its producer
   (`candidate_confidence_column`, which prefers an isotonic
   `rerank_score_calibrated`) is on `cursor/candidate-confidence-rework-f8a2`.
   The matrix workbook `tools/build_all_metabolites_matrix_workbook.py` and the
   hand-committed `docs/paper/data/msi_level_counts.csv` are on
   `cursor/evidence-confidence-calculation-012b`. Confirm the merge state before
   starting either stage.
3. **The truth cohort is small and selected.** 255 truth-verified compounds, 211
   once the pinned ones are held out, 23 of those correct. Any cutoff fitted on
   it must be validated out-of-sample; in-sample fitting on 92 controls is how
   the current Platt constants ended up with an unreachable tier.
4. **Do not fabricate a letter for an uncurated unknown.** 1,403 of the 3,897
   delivered compounds have no curated `SubLevel`. Blank is the correct value.

## Stage 1 — Definition CSV (`definition-csv`)

The one output that is correct regardless of how the banding question resolves,
because a definition is not an assignment. It also forces the `2a`/`2b`
semantics to be written down as the cited scheme defines them rather than as
score bins.

Write `docs/paper/data/identification/msi_sublevel_definitions.csv` with one row
per level and sub-level and these columns:

| Column | Meaning |
|---|---|
| `level` | 1–5 |
| `sublevel` | `2a`, `2b`, `2c`, `3a`… or blank for the bare level |
| `label` | short display name |
| `evidence_required` | what the scheme demands, in evidence terms not score terms |
| `computable_today` | `yes` / `no` / `curated_only` |
| `blocking_input` | for `no`, the missing input named in stage 8 |
| `source` | `schymanski_2014`, `metabolon_curated`, or `repo_rubric` |

Populate `metabolon_curated` rows from the values actually present in
`data/known_X.csv` (`2C`, `3A`, `3B`, `3C`, `3D`, plus the bare levels and the
curation states `KNOWN`, `Released`, `Exclude`, `Artifact`, `4?`). Note in the
file header that `2A`, `2B` and `5A` appear nowhere in the curated data — the
`2a`/`2b`/`5a` vocabulary in the request comes from the slide, not from a
Metabolon export in this repo.

**Expected output:** one CSV, no code changes elsewhere.
**Verification:** every `sublevel` value present in `data/known_X.csv` has a row;
no row claims `computable_today: yes` for a sub-level whose `evidence_required`
names MS² library matching or retention time.

## Stage 2 — Generated level counts (`counts-generator`)

`docs/paper/data/msi_level_counts.csv` on
`cursor/evidence-confidence-calculation-012b` is a hand-committed table with no
generator behind it, which means it cannot be refreshed and cannot be checked.
Replace it with `tools/build_msi_level_counts.py`.

Inputs: `data/known_X.csv` for the curated level and sub-level, and
`docs/paper/data/exec_summary/exec_summary_by_compound.csv` for the delivered
MSI level and the Known/X split. Emit the same columns the hand-made file used —
`MSI Level`, `MSI Sub-level`, `Known`, `X`, `Total` — so anything already
consuming it keeps working.

Two corrections to make while regenerating, both visible in the hand-made copy:

- It mixes case (`2c` and `2C`, `3a` and `3A`) as if they were different
  buckets. Normalize, and keep the curated casing.
- It carries `Exclude`, `Artifact` and `Internal Standard` as if they were
  sub-levels. They are curation states; put them in a separate section or a
  `is_curation_state` column rather than in the level rollup.

**Expected output:** the tool, the regenerated CSV, a test, and a line in the CI
regenerate-and-diff block.
**Verification:** `Total` equals `Known` + `X` on every row; the level totals
reconcile against `msi_level_distribution` in
`docs/paper/data/presentation/audit_presentation_summary.json`
(L1 44, L2 193, L3 2,489, L4 1,165, L5 6); re-running the tool produces no diff.

## Stage 3 — Surface the curated sub-level (`surface-curated`)

2,494 of 3,897 delivered compounds carry a curated `SubLevel`, 949 of them a
fine letter. It is expert-labeled and authoritative, and it needs no model. It
is already read by `tools/build_exec_summary.py::load_metabolon_levels` and
already lands in both exec-summary CSVs as `metabolon_sublevel`; the missing
piece is the report surface.

Render it next to the existing Metabolon-level cross-reference in
`scored_parameters_section_html` and as a column in the master index
(`structural_elucidation/reports/index_render.py`). Label it as Metabolon's
curated assignment, not as this pipeline's output — the same way
`metabolon_level` is already labeled. Withhold it where absent.

**Verification:** a report whose compound has no curated sub-level shows no
sub-level row at all; `tools/diff_source_vs_patched.py` parity markers still
match; the count of reports showing a sub-level equals 2,494.

## Stage 4 — Relabel the availability flag (`availability-flag`)

The matrix workbook presents `rank1_evidence_confidence` as *Confidence Score*.
On the delivered set that column is 62.4% saturated above 0.98 and flat against
truth above 0.9955, so as a confidence it is not one. The single real
distinction it makes is at zero: deciles 1–2 sit at exactly 0.0 and score 0%
correct, everything above scores a flat 13.6–19.1%.

Keep the number if a consumer depends on it, but stop calling it confidence.
Either rename to something that describes it (`evidence_available`,
`pipeline_evidence_score`) or add a boolean beside it and demote the float to a
diagnostic. Update the workbook's `Data_Dictionary` sheet, which currently reads
"Canonical rank-1 evidence_confidence".

**Verification:** no artifact presents a score labeled *confidence* whose band
precision in `audit_sublevel_bands.json` is unordered.

## Stage 5 — Recalibrate Platt (`recalibrate-platt`)

`PLATT_A = 1.887084`, `PLATT_B = -1.520154` were fitted in-sample on 92 controls
and predate the retention axis, and they are the direct cause of the unreachable
tier: `sigmoid(A + B) = 0.5907` is the ceiling. Re-fit on a bundle that carries
axis-level columns — the same bundle `tools/audit_thresholds.py` needs to
populate its currently `not_evaluated` `msi_rule` block — using
`tools/calibrate_confidence.py`.

This fixes reachability. **It will not fix ordering.** Platt scaling is monotone,
so rescaling an inverted score leaves it inverted; stage 6 is not optional.
Adopting a new `(A, B)` moves every tier boundary and therefore the reported tier
distribution, so measure before adopting.

**Verification:** re-run `tools/audit_sublevel_bands.py` and confirm
`unreachable_cutoffs` is empty; diff the tier and MSI distributions against
`audit_presentation_summary.json` and record the shift in the findings doc.

## Stage 6 — Establish a monotone score (`monotone-score`)

The gate for the whole original request. Produce a per-candidate probability
that is monotone with correctness *out-of-sample*, and prove it.

Leading candidate: the isotonic `rerank_score_calibrated` that
`candidate_confidence_column` on `cursor/candidate-confidence-rework-f8a2`
already prefers when the column exists. Isotonic regression is the right shape
for this — it is monotone by construction — but it is fitted per source and the
current cohort is small enough to overfit, so validation must be out-of-sample
or leave-one-out, matching what `tools/audit_structure_predictions.py` does for
the structure axis.

**Exit criterion, and it is a hard one:** `tools/audit_sublevel_bands.py` reports
`monotone: true` for at least one band scheme on the model-scored cohort, with
each band separated by more than `MONOTONE_MARGIN` (0.02). Until that passes,
stage 7 does not start.

## Stage 7 — Derive the bands (`derive-bands`)

Only after stage 6. Then, and only then, the original question has an answer, and
the answer should not be round numbers:

- Choose cutoffs by **equal precision**, not equal width or equal population.
  Pick the score values at which P(correct) crosses chosen levels, so the band a
  reader sees corresponds to a stated accuracy.
- Anchor the top band inside the reachable range. After stage 5 that range is
  whatever the new `(A, B)` gives; do not reuse 0.65 without re-checking it.
- Assign into a **new** field. Do not overwrite `metabolon_sublevel`, which is
  curated ground truth and the only thing available to validate a derived letter
  against. Name the derived one distinctly (`derived_sublevel`) and carry a
  source column, the way the matrix workbook already carries
  `MSI Sub-level Source`.
- Report agreement with the curated letter as a headline number. The prior
  attempt on `cursor/candidate-confidence-rework-f8a2` reached 9.5% exact.

## Stage 8 — Upstream evidence handoff (`evidence-handoff`)

Stages 5–7 can produce a defensible *graded* sub-level. They cannot produce
`2a`/`2b` as Schymanski defines them, because those name evidence types this
repo does not have. Write the handoff for the three missing inputs:

1. A **reference MS² spectral library** — the thing that distinguishes `2a` from
   `2b`, and the 2-series from the 3-series.
2. **Observed retention time per X-compound.** Feature tables carry only m/z.
   `docs/paper/data/retention/ri_calibration.json` covers only the libraries that
   cleared `tools/audit_rt_feasibility.py`.
3. **Candidate-versus-observed diagnostic fragment matching**, computed from the
   observed spectra rather than predicted ones.

The same three are already named in section 5 of `docs/msi_sublevel_findings.md`
on `cursor/candidate-confidence-rework-f8a2` (unmerged, so there is nothing to
link to yet); reconcile with it rather than duplicating it.

## Verification for the whole plan

- `python3 tools/audit_sublevel_bands.py` runs clean and produces no diff.
- `python3 -m pytest tests/test_sublevel_bands.py` passes. Note that
  `test_committed_audit_records_the_two_blocking_findings` asserts the bands are
  *not* ordered — when stage 6 lands, that test is supposed to fail, and the
  correct response is to re-read the findings and update both, not to delete the
  assertion.
- Level totals reconcile against `audit_presentation_summary.json` everywhere
  they appear.
- No new artifact labels an unordered score as confidence.

## Unresolved questions

1. **Is the sub-level scheme in the screenshot Metabolon's or a slide's?** The
   requested vocabulary (`2a`, `2b`, `5a`, and ranges like `2a (0.9–1.0)`) does
   not appear in `data/known_X.csv`, which has `2C` and `3A`–`3D` only. If a real
   Metabolon sub-level specification exists, it should replace stage 1's
   reconstruction — that would change what the definition CSV says and possibly
   make some cells computable that this plan marks as blocked.
2. **What is the sub-level for, commercially?** An include/exclude filter for a
   client data set and a display grade have different requirements: a filter
   needs a defensible cut point, a display grade needs an ordered scale. The
   screenshot's `INCLUDE IN DATA SET` / `INCLUDE IN PATHWAYS` columns suggest the
   former, which stage 3 plus a curated-letter filter would already serve
   without stages 5–7.
3. **Does `5a` mean anything here?** No `5A` exists in the curated data and MSI
   Level 5 is 6 delivered reports. If it is a real bucket, it needs a definition
   before it needs a range.
4. **Which `(A, B)` bundle?** Stage 5 needs a report bundle carrying axis-level
   columns, and this checkout does not have one — the same gap that leaves
   `msi_rule` unevaluated in `docs/paper/data/thresholds/audit_thresholds_summary.json`.
   Locating it is a prerequisite, not part of the stage.
