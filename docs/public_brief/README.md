# Public capability brief

A customer-facing, high-level companion to the scientific manuscript in
`docs/paper/`. It describes what the structural-elucidation capability delivers,
and how it works only at the level of concept, never of parameters.

- `public_brief.md` — editable source.
- `public_brief.docx` — the distributable Word file, against a ten-page budget.
  Measure it rather than trusting a number written down here; see the note on
  pages below.
- `figures/` — four SVGs (two hand-authored diagrams plus two charts regenerated
  from the same audited JSON the manuscript uses) and one PNG screen capture of
  a rendered report.

## Structure

The document follows the internal white-paper writing guidance: an executive
summary carrying the whole argument, the problem, why the established routes fall
short, what the capability does, how biochemical context is assigned, the
evidence, what a customer receives, how to read a result, where it fits, and a
conclusion. The outline governing that
structure — including the message each figure must support and the
capability-to-value pairs behind Table 1 — is `docs/publication_outline.md`, and
the review passes are recorded in `docs/publication_review.md`. Read the outline
before restructuring this document or the manuscript, since the two are meant to
stay consistent with each other.

`docs/writing_guidance_public_brief.md` answers every prompt in that guidance for
this brief. It is an internal working record rather than a release artifact: it
names what the brief withholds and why, which is exactly the material this
folder's release gate exists to keep out. It is built to Word by
`tools/build_guidance_answers_docx.py`.

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

Both rules apply to the report screenshot as well as to the prose, and a scanner
cannot read pixels. `tools/capture_report_screens.py` therefore redacts the page
before capturing it, by removing every block whose text matches the same term
list the release gate applies, and then re-reads the rendered text and refuses
to write the image if anything survives. Redaction removes elements and never
rewrites them: editing the words inside a screenshot would make it a drawing of
a report rather than a picture of one. Regenerate with:

```bash
python3 -m pip install playwright pillow && python3 -m playwright install chromium
python3 tools/capture_report_screens.py --demo
```

The capture is deliberately short. A report screenshot printed taller than the
cap in `tools/build_public_brief_docx.py` cannot share a page, so Word gives it
one of its own and a tenth of the page budget goes to white space.

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

Watch the budgets while filling these in. The word count is close to its 2,600
ceiling, leaving room for tens of words rather than hundreds, and that is the
looser of the two constraints.

Pages are the one to watch, and the count is not portable. Different LibreOffice
versions and different sets of installed fonts disagree by about a page on the
same file, so re-run the gate with `--check-pages` in your own environment after
every edit, without exception, and treat any figure quoted here as a rough guide
rather than a fact about the document.

The key-terms table is what pushes the length: it needs a page, and no amount of
shrinking the report capture buys that page back, which was measured rather than
assumed. Roughly half of the last page is the room the availability and legal
text has to fit into. Two things have each cost a whole page before: a caption
gaining a single line, and the report capture being printed tall enough that
Word stops letting it share a page.

The body size in `tools/build_public_brief_docx.py` is the knob to reach for
first if the count goes over. It is set a full point below where this document
turns a page, which is deliberate headroom, not slack to be spent.

If the placeholder text will not fit, shorten the key-terms definitions before
touching prose. Table rows cost no words at all, since the word counter skips
lines beginning with a pipe, but they do cost lines.

Automated checks cannot clear brand, trademark, or claim-substantiation review.
Marketing and legal sign-off is still required.
