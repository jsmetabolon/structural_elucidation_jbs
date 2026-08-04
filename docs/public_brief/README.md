# Public capability brief

A customer-facing, high-level companion to the scientific manuscript in
`docs/paper/`. It describes what the structural-elucidation capability delivers,
and how it works only at the level of concept, never of parameters.

- `public_brief.md` — editable source.
- `public_brief.docx` — the distributable Word file (currently nine pages,
  against a ten-page budget).
- `figures/` — four SVGs: two hand-authored diagrams plus two charts
  regenerated from the same audited JSON the manuscript uses.

## Structure

The document follows the internal white-paper writing guidance: an executive
summary carrying the whole argument, the problem, why the established routes fall
short, what the capability does, how biochemical context is assigned, the
evidence, what a customer receives, how to read a result, where it fits, and a
conclusion. The outline governing that
structure — including the message each figure must support and the
capability-to-value pairs behind Table 1 — is `docs/publication_outline.md`, and
the four review passes are recorded in `docs/publication_review.md`. Read the
outline before restructuring this document or the manuscript, since the two are
meant to stay consistent with each other.

Two habits from the guidance are load-bearing here and easy to undo by accident.
Substantive paragraphs end on the consequence for the reader rather than on the
mechanism, and comparative statements use defensible wording rather than implied
superiority, because no head-to-head benchmark exists.

## Two rules this document exists under

**It must not be reproducible.** No model families, feature lists, fingerprint or
scoring parameters, level thresholds, internal library numbers, pipeline
version strings, or controlled input filenames. This is enforced by a scanner,
not by reviewer attention.

What the brief may state is the concept: what a stage does, not how it does it.
Pathway assignment is described at that level, as a classification either looked
up from a catalogued compound or inferred from chemically similar ones, because
pathway context appears on every feature including those where no structure is
shown, and a customer cannot judge such a call without knowing which of the two
produced it. Candidate ranking is deliberately not described to the same depth.
The asymmetry is intentional rather than an oversight.

**It must not overstate.** The brief omits detail symmetrically. It drops the
full validation analysis in both directions rather than keeping the flattering
half, and it never presents the confidence output as a probability of
correctness, because leave-one-out calibration does not support that reading.
Every statistic it asserts is checked against the audited claim ledger at
`docs/paper/supplement/claim_evidence.csv`.

## Rebuild and verify

```bash
python3 tools/rebuild_paper_artifacts.py --brief-only
```

This regenerates the figures from committed audit artifacts, rebuilds the Word
file, and runs the release gate. It needs no controlled inputs, so it works in
any checkout.

To verify without rebuilding:

```bash
python3 tools/verify_public_brief.py --check-pages
```

`--check-pages` renders through LibreOffice to confirm the page count and is
skipped with a warning when LibreOffice is absent. The gate checks disclosure,
claim accuracy, overclaim language, and the word, figure, table, and page
budgets.

## Before external release

The brief carries `[REQUIRED: ...]` placeholders for contact details, commercial
availability, platform positioning, trademark usage, legal disclaimer, and the
document version and approval reference. The verifier reports these as warnings;
run it with `--strict-human-inputs` to make them blocking once they have been
filled in.

The positioning placeholder is deliberate. The writing guidance frames this work
as a platform enhancement, but no platform positioning is documented in this
repository, and whether existing deliverables may be named externally is a
marketing decision rather than a technical one.

Watch the page budget while filling these in. The document renders to nine pages
of ten, and the last page currently holds only the closing lines, so replacing
the placeholders with real text should absorb that page rather than add one. The
margin is about a page in total, so long legal boilerplate may require trimming
elsewhere. Re-run the gate with `--check-pages` afterwards.

Automated checks cannot clear brand, trademark, or claim-substantiation review.
Marketing and legal sign-off is still required.
