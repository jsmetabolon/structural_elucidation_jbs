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
- Supplementary Information now states both annotation ceilings in section S4 and
  carries the two new analyses as section S8.
- Fixed a nondeterministic sort in the Supplementary Table S6 generator that
  reordered equal-count rows between runs and produced spurious diffs.

