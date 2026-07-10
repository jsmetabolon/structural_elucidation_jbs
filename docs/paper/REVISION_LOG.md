# Scientific Revision Log

## Why the manuscript changed

A read-only scientific/editorial review identified cohort conflation, circular
pathway evaluation, missing uncertainty, in-sample confidence calibration, and
incomplete reproducibility. This revision prioritizes validity over preserving
earlier headline claims.

## Material claim changes

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
- Added explicit required human inputs for authorship, funding, disclosures,
  ethics, run-environment provenance, and data-access details.

