# Manuscript: Evidence-Aware Reporting of Metabolite Structure Hypotheses

This folder holds the manuscript, committed aggregate audit artifacts, vector
figures, and Supplementary Information.

- `manuscript.docx` — the internal manuscript, branded, with figures embedded. It
  keeps Nature initial-submission mechanics (12-point double-spaced body, line
  numbering, seven display items) so a submission stays cheap to produce; set
  `BRANDED = False` in `tools/build_manuscript_docx.py` to build that version.
- `manuscript.md` — editable source.
- `supplement/supplementary_information.md` — cohort glossary, extended tables,
  calibration, pathway, and mass-mismatch details.
- `supplement/claim_evidence.csv` / `baseline_claims.json` — claim-to-source ledger.
- `data/` — compact pathway, structure, presentation, and LOO audit artifacts.
- `figures/` — main and supplementary SVG figures generated from `data/`, plus
  `fig6_report_example.png`, a screen capture of the rendered report.

A separate customer-facing summary of the same work lives in `docs/public_brief/`.
It shares this folder's audited data but deliberately withholds all methodology;
see that folder's README before editing either document.

Both documents follow the internal white-paper writing guidance. The outline that
governs their structure and figure messages is `docs/publication_outline.md`, and
the record of the review passes is `docs/publication_review.md`. Read the outline
before restructuring either document, since the two are meant to stay consistent
with each other.

`docs/writing_guidance_manuscript.md` answers every prompt in that guidance for
this manuscript, including how its IMRaD sections map onto the guidance's section
names and why limitations are collected here rather than only distributed. It is
built to Word by `tools/build_guidance_answers_docx.py`.

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

## Build the Word manuscript

The Word generator enforces the Nature-oriented title/summary/display-item limits,
uses a 12-point double-spaced body, adds page and continuous line numbers,
converts SVGs to high-resolution embedded PNGs, and preserves Word-native tables
and references.

**Branded or journal-plain, on one flag.** `BRANDED` in
`tools/build_manuscript_docx.py` governs the body face and the heading colour
together: true gives the internal document this repository ships, Montserrat with
Dark Mode Blue headings; false gives a journal submission, Times New Roman with
black headings. Nothing else differs, because the mechanics a submission actually
requires are kept under both settings. The release gate and
`tests/test_manuscript_docx.py` read `BODY_FONT` from the builder rather than
naming a face, so flipping the flag does not mean editing assertions too.

Montserrat is vendored at `assets/fonts/montserrat/` and registered by
`tools/brand.py` regardless of the flag, because the figures are branded either
way and are rasterised through cairosvg at build time; an unresolved family would
put them into a fallback face without failing anything. Both release gates check
that it resolves.

```bash
python3 tools/build_manuscript_docx.py \
  --source docs/paper/manuscript.md \
  --output docs/paper/manuscript.docx
```

The initial-submission DOCX contains the six main figures inline and one main
table (seven display items total). Detailed evidence/MSI tables remain in
Supplementary Information. A figure taller than the page is scaled down rather
than clipped, which matters only for the report capture.

## Regenerate the report capture

Figure 6 is a screenshot rather than a chart, so `tools/paper_figures.py` does
not produce it and neither does CI: capture needs a headless browser that the
documentation dependencies leave out. Regenerate it deliberately.

```bash
python3 -m pip install playwright pillow && python3 -m playwright install chromium
python3 tools/capture_report_screens.py --demo
```

`--demo` builds its own report bundle from data committed here, so it works in
any checkout. Point it at a delivered bundle with `--report-dir` when one is
available. The same command rewrites the brief's capture, which is redacted and
re-scanned with the brief's disclosure term list before it is written.

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
