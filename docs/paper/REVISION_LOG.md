# Scientific Revision Log

## Why the manuscript changed

A read-only scientific/editorial review identified cohort conflation, circular
pathway evaluation, missing uncertainty, in-sample confidence calibration, and
incomplete reproducibility. This revision prioritizes validity over preserving
earlier headline claims.

## Editorial revision (white-paper writing guidance)

A later editorial pass applied the internal white-paper writing guidance to both
publication surfaces. In this manuscript the changes are presentational rather
than evidentiary: the summary was reshaped into problem, approach, evidence and
consequence; figure captions now lead with the message each figure supports;
Applications states what an analyst can do that was previously impractical, with
the supporting result and the bounding limitation named alongside it; and a
terminology paragraph in Methods fixes recovery, accuracy, agreement and
consistency as distinct terms. Comparative statements are limited to defensible
wording because no head-to-head benchmark has been run. Section order, figure
count, table count, citations and audited values are otherwise unchanged, apart
from the coverage reconciliation recorded below.

## Editorial revision (library-matching framing, report figure, glossary)

A later pass addressed how both surfaces described spectral-library matching,
added a figure showing the delivered report, and defined the vocabulary each
document depends on.

- **Library matching is described by its coverage rather than as a shortfall.**
  The Introduction previously opened on "Library-retrieval-only workflows cannot
  propose genuinely novel structures" and called the manual route that follows a
  "bottleneck"; Section 1.1 said MS2Query was "bounded by its reference library"
  and described the generative tools as "not restricted to database membership"
  and "independent of spectral libraries". Read together, that framed the method
  which produces the named fraction of every deliverable, and which this
  workflow uses as an evidence axis, as a deficient one. The boundary is now
  located in catalogue coverage rather than capability, the two routes are
  described as answering different questions, and the Discussion states that a
  library match remains the most direct route wherever the compound is
  catalogued. The brief's background section, headed "Why the established routes
  fall short", was renamed and reordered the same way.
- **No audited value moved, and no equivalence claim was added.** Exact-source
  rank-1 recovery of 38.1% (16/42) and predicted-source 56.0% (28/50) both
  remain, with their overlapping intervals and the conclusion that catalog
  provenance is not confirmation. Two sentences changed sense rather than
  content: Results previously asserted that provenance "confers no recovery
  advantage" and now says the cohort does not resolve a difference in either
  direction, and the Discussion's "provenance is not confidence" is now
  "provenance is not by itself confirmation", scoped to how this workflow's own
  source field should be read. The retracted phrase "essentially
  indistinguishable" was not reintroduced.
- **Results gained Section 3.6 and Figure 6.** Claim control was described in
  3.1 and 3.3 but never shown. The figure pairs the public and internal
  renderings of one feature, which differ in how much of the candidate list each
  may expose. It is a screen capture of a demonstration control built from data
  committed here, labelled as such in the caption, and not a delivered customer
  result. The brief carries the public panel alone.
- **Methods gained Section 2.7 as a glossary**, and the former Section 2.7
  became 2.8. The four evaluation senses moved there out of the Methods
  paragraph that carried them, joined by the provenance vocabulary — reference
  library, retrieval, analogue search, de novo generation, catalog provenance,
  authentic standard. Supplementary S2 carries the same provenance terms. The
  brief carries a plain-language key-terms table. Both were subsequently moved
  to the back of their documents and expanded; see below.

## The glossaries moved to the back and grew

Placing the definitions inside Methods 2.7 put them in the one section a reader
looking up a word would not search, and the same was true of the brief's
key-terms table inside How to read a result. Both are now sections in their own
right at the back of their documents, and both were expanded from roughly a
clause per term to a definition.

- **Manuscript.** A `## Glossary` section closes the document, with 53 entries
  in seven groups: the measurement, where a candidate comes from, comparing
  structures, evaluation vocabulary, what a result claims, biochemical taxonomy,
  and failure modes. Section 2.7 reverted to *Software and reproducibility*, and
  the pointer in 2.6 now names the four evaluation senses rather than only citing
  a section number. Entries carry their audited numbers where a term is
  meaningless without them, so catalog provenance is defined next to its 38.1%
  rank-1 recovery, the annotation ceiling next to the two rates that bound it,
  and sparse MS2 next to its 43.6% coverage.
- **Brief.** A `## Key terms` section sits between the conclusion and the
  availability block, with 24 entries against the previous 10. The reporting
  outcome table took the number the glossary vacated, so it is now Table 2 and
  the glossary is Table 3.
- **Two mechanical traps worth recording.** The citation-parity check reads any
  line starting with a digit and a period as a reference entry, so the manuscript
  glossary has to be bullets rather than a numbered list. In the brief, table
  rows cost nothing against the word budget, because the word counter skips lines
  beginning with a pipe, but the scanner that checks statistics against the claim
  ledger does read them, so a glossary definition may not carry a bare numeral.
- **The brief now runs to ten pages of ten.** The key-terms table needs a page.
  Shrinking the report capture to recover it was measured at three heights and
  does not work, so the figure keeps its size and roughly half of the last page
  is what remains for the availability and legal text.
- `Glossary` was added to the end of `SECTION_ORDER` in
  `tools/build_manuscript_docx.py`. Unrecognised headings are appended after that
  tuple, so the section would have landed in the same place regardless; naming it
  makes the placement a decision. In the Word build it therefore follows the
  declarations, which is worth the corresponding author confirming.

## Review pass on inflated annotation counts and rigor

A comment pass on both surfaces asked for the competitive argument to address a
practice rather than a product: matching every mass-spec signal against a broad
reference database and reporting many thousands of metabolites, including weak
matches and including duplicates, adducts and isotopologues of one compound. No
company is named in either document, and comparison to other workflows is kept
to the manuscript.

- **The Introduction gained the opposite failure mode.** It described only the
  coverage boundary of library matching. It now also names the three mechanisms
  that inflate an annotation list without adding identifications: adducts and
  isotopologues annotated separately, one compound counted again per
  chromatographic method, and matches admitted with no similarity floor. The
  consequence stated is that annotation totals are not comparable between
  workflows unless the redundancy and admission criteria behind them are given.
- **The Discussion gained a false-positive control paragraph.** Every mechanism
  in it was already in the paper, spread across the mass-mismatch section, the
  level rubric and the coverage result, where they could not be seen as a
  system. Nothing new is claimed: it closes on the same limit as before, that a
  mass-consistent alternative sits lower in 523 of 598 demoted lists, so claim
  limitation is containment rather than repair.
- **Confidence wording.** The comment asked us to state that we provide a
  validated score. That would be false, given the leave-one-out Brier of 0.255
  already reported here. The manuscript instead states that both signals were
  tested, publishes what each supports and what it does not, and notes that a
  confidence value published without that statement cannot be acted on
  differently from one that was never tested. The brief carries the same point
  in plainer words and does not use "validated" anywhere.
- **Limitations reframed, not softened.** Every entry is unchanged. The section
  now opens by naming them as the bounds within which the results hold, because
  the same facts read the other way round are the rigor the competitive argument
  rests on.
- **Blocked on missing inputs.** Two comments asked for matrix-specific content
  and one for a table of improved percentage coverage. Specimen-level matrix and
  species metadata were never retained in the reporting bundle, and no audited
  baseline exists for the named fraction of a study, so there is nothing to
  improve on. Both documents carry `[REQUIRED: ...]` placeholders rather than
  estimates.
- **The brief moved up a level.** The pathway section lost four bolded
  sub-paragraphs and six percentages, and the confidence-ordering paragraph lost
  half its length, in both cases keeping the caveat and dropping the arithmetic.
  That freed a page, which is what paid for the new competitive material; the
  brief is back to nine pages of ten.

## Material claim changes

- **Subpathway annotation ceiling reported.** The manuscript gave the super-pathway
  harmonization ceiling, 62.5% on 267 identical structures, but not the subpathway
  one, although the same audit block already carried it: the reference database and
  `known_X` agree only 34.8% (93/267) at that level. Without the bound, the 18.8%
  neighbor-consensus subpathway agreement read as model error alone. Results,
  Discussion and Limitations now carry both bounds and note that agreement
  corresponds to roughly 83% and 54% of the achievable rates, with the assumption
  behind that normalization stated rather than assumed.
- **Pathway coverage counts reconciled with the audit artifact.** Results
  previously reported assignment statuses as 3,452 strong-similarity vote and 445
  weak-similarity backfill, pathway-confidence tiers of high 1,788 / medium 666 /
  low 1,443, and a best-reference similarity median of 0.721. The committed audit
  summary gives exact reference match 1,459, strong-similarity vote 1,995 and
  weak-similarity backfill 443; tiers high 1,863 / medium 616 / low 1,418; and a
  median of 0.733. The manuscript now states the audited values and reports the
  three assignment statuses separately rather than folding exact matches into the
  strong-similarity count.
- **Pathway 94.5%/89.7% headline removed.** In the delivered reports, 1,459
  exact-structure cases are direct reference-taxonomy lookups. They are 100%
  internally consistent by construction and are no longer described as prediction
  accuracy.
- **Independent neighbor-consensus result added.** Among 644 controls with
  independent `known_X` labels, agreement is 52.2% super-pathway (95% CI
  48.3–56.0%) and 18.8% subpathway (15.9–22.0%).
- **LOO claim narrowed.** The 86.7%/71.1% result is explicitly an intrinsic test of
  the neighbor-vote mechanism (k=10, similarity floor 0.30, squared weights), not
  the full production predictor.
- **Structural top-k results added.** On n=92: rank 1 47.8%, top 3 64.1%, top 10
  73.9%, all with Wilson intervals.
- **Exact-vs-predicted language qualified.** Exact-source 38.1% and
  predicted-source 56.0% have overlapping intervals; “indistinguishable” was
  removed and no equivalence claim is made.
- **Confidence downgraded to exploratory.** Leave-one-out Brier score is 0.255 and
  calibration slope 0.056; removing the pathway axis slightly improves Brier to
  0.253. Absolute probabilities are no longer presented as externally calibrated;
  report-facing explanatory text was updated to match.
- **Mass-mismatch cohorts reconciled.** There are 598 total mismatches; 590 are in
  the MSI-L4 deep-audit subset. Recoverable lower-rank alternatives are 523/598
  overall and 519/590 in the subset. These are reranking opportunities, not proven
  correct structures.

## Packaging changes

- Added machine-readable pathway, structure, presentation, and LOO audit artifacts.
- Added a claim-to-evidence ledger with source hashes.
- Rebuilt all figures from committed JSON data with cohort/n/CI labeling.
- Added Supplementary Information with cohort glossary and extended tables.
- Added a one-command regeneration/verification entry point.
- Added a Nature initial-submission DOCX: 67-character title, 136-word summary,
  12-point body, double spacing, line/page numbering, five embedded
  figures plus one main table, and Nature-ordered References/Methods/declarations.
  The body was Times New Roman until the brand pass recorded below.
- Added explicit required human inputs for authorship, funding, disclosures,
  ethics, run-environment provenance, and data-access details.
- Extended the package verifier to assert the Results pathway coverage counts and
  the best-reference similarity median against
  `data/pathway/audit_pathway_summary.json`. This is the check that would have
  caught the drift recorded above.
- Added a `reference_db_vs_known_x_agreement_sub` row to the claim ledger and to
  `tools/build_paper_claim_audit.py`, drawn from the same audited cohort as the
  existing super-pathway ceiling, so the bound is citable in both documents.
- Deleted four superseded figures that were no longer regenerated by
  `tools/paper_figures.py` and were referenced by nothing:
  `fig2_pathway_by_similarity.svg`, which still rendered the retracted 94.5%
  reference-DB headline; `fig5_confidence_tiers.svg`, which still described the
  aggregate score as a calibrated probability; and `fig3_pathway_by_tier.svg` and
  `fig4_msi_distribution.svg`, both superseded by current figures. The package
  verifier now scans every file under `figures/` for the stale phrases, not just
  the manuscript text, since an unreferenced figure is invisible to the link checks.
- Added `tools/analyze_selective_prediction.py` and Supplementary Figure S5,
  measuring accuracy against the fraction of each cohort retained. It runs from the
  committed per-control audit outputs and needs no controlled bundle.
- Added `tools/analyze_taxonomy_harmonization.py` and Supplementary Table S8,
  separating vocabulary disagreement from prediction error by applying the same
  correction to the predictor and to the annotation ceiling. Both new analyses are
  wired into the full rebuild.
- Added `tools/analyze_library_stratification.py` and Supplementary Table S9,
  measuring the library-imbalance limitation the manuscript states but never
  quantified. Inferred agreement is consistent across the three libraries with
  enough controls to score, with every interval overlapping at both levels. The
  structural cohort cannot be stratified: two libraries contribute fewer than ten
  truth-verified controls each.
- Supplementary Information now states both annotation ceilings in section S4,
  carries the selective-prediction and harmonization analyses as section S8, and the
  library stratification as section S9.
- Fixed a nondeterministic sort in the Supplementary Table S6 generator that
  reordered equal-count rows between runs and produced spurious diffs.
- Added `tools/audit_thresholds.py` and `tools/audit_rt_feasibility.py` with their
  committed artifacts under `data/thresholds/` and `data/retention/`. Both run
  from committed inputs and are wired into the rebuild and CI, so the numbers the
  reporting rules rest on are reproducible from this repository.

## Packaging changes for the report figure

- Both verifiers parsed every linked figure as XML, which was right while every
  figure was a generated SVG and made a screenshot impossible to reference. The
  check now branches on format: SVGs still have to parse, and PNGs are validated
  against their signature, IHDR dimensions and a minimum width, so renaming a
  file's extension does not get it past the gate. `tools/figure_formats.py` holds
  the shared check and `tests/test_figure_formats.py` asserts each way it fails.
- Display-item budgets moved from five figures to six in the manuscript and from
  four figures and two tables to five and three in the brief.
- `tools/capture_report_screens.py` produces both figures. It can screenshot a
  delivered bundle, and it can build a demonstration bundle from committed data,
  so the figures are regenerable without the controlled report bundle or the
  reference database. It is wired into neither `rebuild_paper_artifacts.py` nor
  CI, because capture needs a headless browser that the documentation
  dependencies leave out; `fig1_workflow.svg` and the two hand-authored brief
  diagrams are the existing precedent for a figure the generators do not own.
- The public capture is redacted with the brief's own disclosure term list and
  then re-scanned before the image is written. `tools/verify_public_brief.py`
  exposes that matching as `disclosed_terms()` so the guard and the gate cannot
  drift apart.
- The manuscript DOCX builder scales a figure down when its aspect ratio would
  run it off the page, and the brief builder applies a much tighter height cap:
  a report screenshot printed any taller cannot share a page, and the brief has
  ten pages to hold everything including the legal text still to be written.
- The report template gave inline `code` chips a near-white background while
  they inherited the header's white text, so the molecular formula in every
  report subtitle was invisible. Fixed in `structural_elucidation/reports/compound.py`.

## Reporting-rule changes awaiting a report bundle

The reporting rules changed after the delivered set was produced. The manuscript
and the Supplementary Information still describe the rules that generated the
committed distributions, which keeps each document internally consistent, but
both the rule descriptions and the distributions have to be updated together the
next time the reports are regenerated. Whoever does that should re-read Methods
(the MSI level definitions), Results (the tier and level distributions) and
Supplementary Table S2 against the regenerated artifacts.

What changed, and what it will move:

- **Levels 4-5 now withhold the structure in the rendered HTML.** The published
  description of the workflow already said they did, and the 2,726 shown against
  1,171 withheld were counts of MSI assignment rather than of anything the
  renderer did: it drew the structure panel unconditionally. The claim is now
  true of the artifact as well as of the rubric. The counts themselves do not
  change for this reason, but they will move for the ones below.
- **A catalog hit needs observed-MS2 fit >= 0.50 to reach Level 2**, against 0.30
  for a de novo prediction. Methods currently gives a single 0.30 bar. The
  evidence is in the manuscript already: exact-source rank-1 recovery is 38.1%
  against 56.0% for predicted, while a catalog hit pins the structure axis to
  1.0. Expect some Level 2 reports to move to Level 3.
- **A fifth evidence axis, retention-index fit, at weight 0.05.** It is available
  only on libraries 209 and 400, whose retention models cleared
  `tools/audit_rt_feasibility.py`, and it can lower the evidence score without
  raising it. A three-sigma residual demotes the tier and blocks Levels 2 and 3.
  This will move both distributions and makes `PLATT_A`/`PLATT_B` stale: re-run
  `tools/calibrate_confidence.py` on the bundle before quoting confidence values.
- **A new Refuted tier**, entered when the curation ledger records a failed
  authentic-standard run or a curator rejection. No delivered report carries one
  yet, so the tier distribution gains a category with a count of zero until such
  a verdict is recorded.
- **The headline confidence and the model's evidence score are now separate
  fields.** Audits that read `overall` were measuring a blend of evidence and
  display policy; `evidence_score` is the quantity to calibrate and to report.
  The committed structure audit already uses the raw axis score, so its numbers
  are unaffected, but anything new should read `evidence_score`.
- **The pathway axis now scores a real vote share.** `pathway_vote_share` was in
  the confidence model's input contract from the start, but nothing ever wrote
  it, so the axis had been falling back to a three-valued tier map of 0.90 /
  0.60 / 0.30 for every report. The reports now reconstruct the neighbor vote
  and supply the actual share. The effect is bounded: the axis carries weight
  0.10 of 1.05, so moving it across its full range shifts `evidence_score` by
  0.047, and the largest single change is roughly 0.04 for a report whose tier
  said high while its vote was in fact split. That is a correction rather than
  drift, since the tier was only ever a proxy for the share. It is one more
  reason `PLATT_A`/`PLATT_B` need re-fitting on the next bundle.
- **The pathway audit's `vote_share_pctiles` were all NaN** for the delivered
  set, because the report sidecar never carried the field the audit reads. They
  will populate on the next regeneration, so that block of
  `audit_pathway_summary.json` becomes reportable for the first time.

## Presentational revision (brand palette and typeface)

Both surfaces were restyled to the Metabolon brand. No number, claim, cohort,
caption message, or section order changed; this pass is presentational in full.

- **Colour and chart furniture now come from `tools/brand.py`.** The palette, the
  categorical series order, and the axis, gridline, label and whisker treatment
  were read off the brand guide and the chart formatting in
  `figure_examples.xlsx` rather than chosen. Both documents' figures read from
  that one module, so neither can be restyled without the other following.
- **Both documents are set in Montserrat**, including this manuscript, which was
  12-point Times New Roman. The double spacing, 12-point body and line numbering
  a Nature initial submission asks for are unchanged. A journal that requires a
  serif will re-set the manuscript in its own template regardless.
- **Montserrat is vendored at `assets/fonts/montserrat/`** under the SIL OFL, and
  `tools/brand.py` registers it with fontconfig before cairosvg loads. This is
  not convenience packaging: figures are rasterised at build time, so a checkout
  without the family produces documents whose figures are set in a fallback with
  every other check still passing. Both release gates now fail on that condition.
- **Two workflow diagrams were re-laid out, not just recoloured.** Montserrat is
  wider than the Helvetica they were drawn against and the labels overflowed
  their boxes. The new box widths and line breaks were measured against the real
  font.
- **The public brief's body dropped from 10.5 to 10 point** to pay for the wider
  face inside an unchanged ten-page budget. Its headings and captions were
  measured at their original sizes and kept them.

## Documentation revision (evidence-score parameter table)

Supplementary Table S10 now states what each evidence parameter measures, what
it weighs, and where its bands fall. It is generated by
`tools/build_parameter_table.py` from the constants that apply those bands, so
it cannot describe a model other than the one that ran; `tests/test_parameter_table.py`
fails if a weight or cutoff moves without the export being regenerated. No
scoring behaviour changed in this pass.

Writing it down required reconciling a drafted parameter sheet against the
implementation, and four things did not match. They are recorded here rather
than quietly resolved, because the draft has circulated.

- **The drafted mass bands of 0-1, 1-3 and 3-5 ppm are not implemented and never
  were.** `mass_match_score` scores an absolute mass error: 1.00 at or below
  0.001 Da, 0.00 at or above 0.005 Da, linear between. Because the window is
  absolute, its ppm equivalent moves with the compound — near 3.3 ppm at 300 Da,
  near 1.7 ppm at 600 Da — so no single ppm band can describe it. Table S10
  quotes Da with the ppm equivalent at a stated reference mass.
- **The mass sub-score cannot take a value between 0.80 and 1.00.** The plateau
  returns 1.00 at 0.001 Da while the linear arm returns 0.80 at the same point,
  so the interval between them is unreachable. `MASS_STRONG = 0.90`, which
  governs the strong-mass condition for MSI Level 2, therefore behaves as an
  exact test for 1.00 rather than as the 0.90 bar it reads as. This is a
  property of the scorer, not of the table; changing it would move every
  affected report and is left for a pass that can re-calibrate and re-audit.
- **Mass Match and Mass Formula Match are one axis**, not two. Both are the
  `mass` axis, whose report label is already "Mass / formula match".
- **Library neighbors carries no weight.** It is a count shown on the report
  card; `build_axes` never reads it, and the drafted five-neighbor threshold
  does not exist. It appears in S10 marked unscored so it is documented without
  being mistaken for an axis. Conversely the draft omitted retention-index fit,
  which is a real axis at weight 0.05 and the only one that can demote a report
  to Tentative and block MSI Levels 2 and 3 by itself; S10 includes it.

## Presentational revision (the manuscript joins the brand)

The manuscript was branded in part during the earlier palette pass and finished
here. No claim, number, cohort or section changed.

- **The headings and captions were still black.** The face, the figures and the
  table fills had already moved; the title, both heading levels and the captions
  had not. They are now Dark Mode Blue and brand muted, matching the brief.
- **The typeface question is closed.** It had been open since the palette pass and
  was listed as needing a person. The decision is that this is the internal
  document and it is branded. `BRANDED` in `tools/build_manuscript_docx.py` now
  governs the body face and the heading colour together, and both the release gate
  and `tests/test_manuscript_docx.py` read `BODY_FONT` from the builder rather than
  naming a face, so producing a journal-plain submission is one flag rather than
  four style blocks and two assertions. Part of why the question stayed open was
  that reversing it looked expensive; it no longer is.
- **The Nature mechanics are kept deliberately.** Twelve-point double-spaced body,
  line numbering, seven display items. Those are structural requirements rather
  than styling, and keeping them is what makes the flag cheap.
- **IMRaD was retained.** The writing guidance describes a white-paper arc, and
  this document is Abstract, Introduction, Methods, Results, Discussion,
  Limitations, Conclusion. Those are different documents rather than one document
  formatted two ways, and the guidance's "separate explanation from proof" rule is
  already satisfied by the section structure, with Results paragraphs leading on
  their claim rather than their evidence.
- **One review comment was closed.** The direction to highlight that other
  approaches' confidence scores are of unclear meaning and untested accuracy had
  been answered only for this workflow's own score. The related-work section now
  separates score interpretability from accuracy, frames it as a reporting gap
  rather than a defect, and points at Section 3.3, which reports that the score
  here orders candidates within a feature and does not separate features. No
  accuracy comparison was added and no tool is criticised.
- **Both guidance worksheets were refreshed** to carry review 6, which they had
  missed, and the review log's note that they were stale is retired.

## Decision record: the mass axis

Measuring the drafted ppm bands turned up something larger than the question
that prompted it. The audit is `tools/audit_mass_policy.py`, run over the 3,897
delivered reports; its output is `docs/paper/data/thresholds/audit_mass_policy.json`.
The reader was validated by reproducing the committed 598 rank-1 mass
mismatches, the 15.35% rate, and the 523 with an in-list alternative exactly.

No scoring behaviour changed in this pass. What changed is what is written down.

### Settled

- **The delivered reports were not scored by the function the documentation
  described.** Two functions were named `mass_match_score`, with different
  signatures and tolerances 200-fold apart. `pipeline/candidates.py` imported the
  one in `chemistry/conjugates.py`, which runs at a 1.0 Da default; the
  documentation described the one in `inference/predict.py`, which has a 0.001 Da
  plateau and a 0.005 Da tolerance. Nothing failed, because the names matched and
  the signatures were never checked against each other. The conjugates function
  is now `linear_mass_match_score` and each definition points at the other.
- **The scale varies by generator, and systematically.** Inverting the score
  recovers the tolerance each candidate was scored under: library-anchored and
  peptide de novo candidates at 1.0 Da, the balanced assembly path at 5.0 Da, and
  the biotransformation families at 0.005 Da. Each generator is consistent with
  itself. So the mass sub-score is not one quantity: candidates on scales that
  differ by up to 200-fold are ranked against each other, averaged with the other
  axes, and the Platt coefficients were fitted on the mixture. This is the most
  consequential finding here and it was not previously recorded anywhere.
- **The drafted 0-1 / 1-3 / 3-5 ppm bands are rejected.** They shrink rank-1
  strong-mass calls from 2,436 to 1,424 and made recovery worse on every
  comparison. The question is closed.
- **The mass tolerance is not a useful dial.** Swept from 0.005 Da to 5.0 Da, a
  thousand-fold range, rank-1 recovery moves by two controls out of 92 and no
  setting gains a single one. Losses fall as the tolerance loosens, so the gate
  does least harm where it fires least. Meanwhile rank-1 strong-mass claims range
  from 2,057 to 3,105 across the same sweep: the axis relabels reports without
  discriminating between them.
- **Handoff item 2 should not ship as specified.** Its projection of roughly 523
  recoverable reports counts cases where a mass-consistent alternative exists,
  which is a claim about mass consistency rather than about correctness. Applied
  to the truth controls the gate gained nothing at rank 1 and lost three to five.

### Open, and why they cannot be closed here

- **Whether a 1.0 Da tolerance was intended.** At 300 Da that is roughly 3,300
  ppm, which admits essentially every candidate the generator proposes. The
  profile carrying it is named `strict_1da`, which suggests deliberate, but one
  Dalton is not strict for high-resolution data. Someone with the pipeline's
  history should say which.
- **Whether to unify the scorers, rescore and recalibrate.** The technical work
  is clear and the cost is not: a pipeline re-run, a Platt refit, and re-auditing
  every number in both documents. It should be decided together with the question
  above rather than slipped in as cleanup.
- **Whether the truth set can support any of this.** Only 9 of the 598
  mass-mismatch reports have a known correct answer. That single number is why
  the two decisions above rest on judgement rather than measurement, and it
  gates any future rerank change as well. Expanding it is handoff item 4, and it
  now looks like a prerequisite rather than a parallel improvement.

