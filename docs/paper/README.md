# Manuscript: Evidence-Aware Reporting of Metabolite Structure Hypotheses

This folder holds the manuscript, committed aggregate audit artifacts, vector
figures, and Supplementary Information.

- `manuscript.docx` — primary Nature initial-submission manuscript (figures embedded).
- `manuscript.md` — editable source.
- `supplement/supplementary_information.md` — cohort glossary, extended tables,
  calibration, pathway, and mass-mismatch details.
- `supplement/claim_evidence.csv` / `baseline_claims.json` — claim-to-source ledger.
- `data/` — compact pathway, structure, presentation, and LOO audit artifacts.
- `figures/` — main and supplementary SVG figures generated from `data/`.

A separate customer-facing summary of the same work lives in `docs/public_brief/`.
It shares this folder's audited data but deliberately withholds all methodology;
see that folder's README before editing either document.

Both documents follow the internal white-paper writing guidance. The outline that
governs their structure and figure messages is `docs/publication_outline.md`, and
the record of the four review passes is `docs/publication_review.md`. Read the
outline before restructuring either document, since the two are meant to stay
consistent with each other.

## Verify from a clean checkout

```bash
python3 tools/rebuild_paper_artifacts.py --verify-only
```

This validates the abstract word count, required files, audited counts, stale
claims, citations, figure links/order, and SVG syntax, then runs the public
brief's release gate. Human-supplied submission fields (authors, email, funding,
disclosures, ethics, repository DOI) are reported as warnings unless
`tools/verify_paper_package.py --strict-human-inputs` is used.

`.github/workflows/publication-checks.yml` runs this on every pull request that
touches `docs/`, `tools/` or `tests/`, together with the document tests and a check
that the derived artifacts still reproduce from their committed inputs. The page
count is not in CI: `--check-pages` needs LibreOffice, which is too large to install
on every run, so it stays a manual check before release.

The stale-claim scan covers every markdown file under `docs/` and every figure.
Three records are exempt because their purpose is to document what was withdrawn:
`REVISION_LOG.md`, `pathway_audit_findings.md`, and `publication_review.md`.

## Fully regenerate with controlled inputs

```bash
python3 tools/rebuild_paper_artifacts.py \
  --reports-root _audit_pull/reports \
  --db _bundle/metabolite_synonyms.sqlite \
  --known-x data/known_X.csv
```

This one command regenerates the presentation summary, stratified pathway audit,
neighbor-vote LOO artifact, structural rank-k/calibration audit, claim ledger,
all SVGs, Supplementary Information, and then verifies the package. The report
bundle and reference database are proprietary/controlled inputs and are not
stored in Git.

Authoritative P0/P1 results are intentionally absent until
`v20_all_libraries_top10_predictions_with_p0_flags.csv` and
`models/p1_reference_v1/` are supplied; the 92-control calibration set is not
substituted for those cohorts.

## Build the Nature Word manuscript

The Word generator enforces the Nature-oriented title/summary/display-item limits,
uses 12-point Times New Roman with double spacing, adds page and continuous line
numbers, converts SVGs to high-resolution embedded PNGs, and preserves Word-native
tables and references.

```bash
python3 tools/build_manuscript_docx.py \
  --source docs/paper/manuscript.md \
  --output docs/paper/manuscript.docx
```

The initial-submission DOCX contains the five main figures inline and one main
table (six display items total). Detailed evidence/MSI tables remain in
Supplementary Information.

## Current audited results

- 3,897 delivered reports; 2,726 show a structure (L1–L3), 1,171 withhold it.
- Pathway assignment statuses: exact reference match 1,459, strong-similarity vote
  1,995, weak-similarity backfill 443; tiers high 1,863 / medium 616 / low 1,418.
- Structure truth set (n=92): rank 1 47.8%, top 3 64.1%, top 10 73.9%.
- Direct pathway lookups (n=1,459): 100% reference consistency by construction,
  not independent prediction accuracy.
- Neighbor consensus vs independent labels (n=644): super-pathway 52.2%,
  subpathway 18.8%. Read against the annotation ceiling below, that is roughly 83%
  and 54% of the achievable rates.
- Annotation harmonization ceiling (n=267 identical structures): the reference
  database and the independent labels agree only 62.5% at super-pathway and 34.8%
  at subpathway, which bounds any independently scored pathway result.
- Neighbor-vote mechanism LOO (n=5,095): super-pathway 86.7%, subpathway 71.1%.
- Leave-one-out identification-confidence Brier score 0.255; the aggregate score
  is exploratory rather than externally calibrated.
