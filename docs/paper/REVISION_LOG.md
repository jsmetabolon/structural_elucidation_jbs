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
  12-point Times New Roman, double spacing, line/page numbering, five embedded
  figures plus one main table, and Nature-ordered References/Methods/declarations.
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

