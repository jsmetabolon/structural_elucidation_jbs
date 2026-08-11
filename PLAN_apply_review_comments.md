# Plan: apply the sixteen review comments to both documents

Scratch file, deliberately at the repository root and deliberately without YAML
frontmatter, so nothing can swallow it on the way to the screen. Delete it once
the work is done; it is not part of the repository record.

**Goal.** Apply the sixteen comments to both publication surfaces: raise the
brief to a higher level, add an unnamed contrast against indiscriminate database
matching that inflates annotation counts, make rigor and validation the spine of
both documents, and keep named workflow comparisons in the manuscript only.

---

## 1. What the comments resolve to

**The competitive argument, unnamed.** The competitor is the practice of matching
every mass-spec signal against a broad reference database and reporting many
thousands of metabolites, including bad matches, and including duplicates,
adducts and isotopes of the same molecule. No company is named anywhere. The
point to land is that our approach is more rigorous and yields more biological
insight, without overstating. Named workflow comparisons stay in the manuscript.

**Rigor as the spine.** Validation and accuracy are central; limitations are
framed as restraint rather than shortfall; false-positive control is stated
explicitly; the score we publish is transparent and its meaning is stated.

**Level 1 as the foundation, not the ceiling.** Give Level 1 identifications
wherever possible and do not stop there; a large Level 1 library is what makes
predicting new structures possible; the capability supplements what Level 1
teaches and extends signatures into new biology.

**Structure and impact.** Problem in the first two sections, solution after.
Allude to impact in the executive summary.

Three comments cannot be answered from committed data. They are handled in
section 5.

---

## 2. The brief: raise the level

The brief reads at manuscript register in places. Two passages carry most of
that, and they are where the "higher level" instruction bites hardest.

- **How biochemical context is assigned** runs four bolded sub-paragraphs and six
  percentages, including a comparison of two label vocabularies disagreeing with
  each other. Compress to the reader-facing point: context reaches a feature
  either by lookup or by inference from similar compounds, the two do not carry
  equal weight, each report says which, the broad class is dependable and the
  finer one is a lead. Keep the honesty, drop the arithmetic.
- **The ordering paragraph in How to read a result** explains what the confidence
  output can and cannot rank. Keep the substance, since it is exactly the "do not
  overstate" instruction, but say it in half the words.

Raising the level also relieves the page budget. The brief sits at ten pages of
ten, so the words this frees are what pays for the new competitive material.

---

## 3. The brief: the new argument

### 3.1 Add a third route to "What the established routes cover, and what they leave"

The section currently covers library matching and expert interpretation. Add,
without naming anyone, the practice of matching every detected signal against a
broad database and reporting whatever comes back. State what is wrong with it in
the reader's own terms: a list of several thousand metabolite names is not
several thousand findings when the same molecule appears more than once, when
adducts and isotopes of one compound are counted separately, and when weak
matches are reported with the same weight as strong ones. A long list transfers
the work of deciding what to trust onto the reader.

Putting it here keeps the section skeleton intact, which matters because the
outline and the review passes govern it, and it satisfies the comment about
establishing the problem across the first two sections before proposing the
solution.

### 3.2 Answer it with what we actually do

All of this is in the audited record rather than asserted:

- Features, not raw signals: reproducible precursor features with their
  fragmentation scans aggregated into a consensus spectrum.
- Adduct-aware mass gating, so a candidate is scored against the neutral mass an
  adduct hypothesis implies rather than against the raw measurement.
- Isotope and carbon agreement among the evidence the ranking reads.
- Deduplication: a feature appearing in more than one library is counted once,
  which is why the delivered set is 3,897 distinct features rather than a larger
  number counted repeatedly.
- Mass-mismatch demotion: 598 top candidates whose mass did not fit were demoted
  rather than shown.
- Withholding: 1,171 of 3,897 features are delivered with no structure at all.

That last pair is the strongest version of "we don't inflate counts", because it
is a number we deliberately gave up.

### 3.3 Rigor, Level 1, impact

- **Rigor and restraint.** The brief already says withholding is designed
  behaviour. Strengthen it into the document's spine: the restraint is what makes
  the rest usable, and because the same criteria run over every feature, results
  are comparable across studies rather than only within one.
- **Level 1 framing**, in Where this fits and the executive summary: Level 1
  identifications are given wherever the evidence supports them, the capability
  begins where that stops, and the size of the Level 1 library is precisely what
  makes proposing new structures possible, since every proposal is judged against
  catalogued chemistry. This reinforces the library-matching rebalance already in
  place rather than cutting across it.
- **Impact in the executive summary.** It currently ends on completeness and
  comparability. Add the consequence: the unnamed fraction becomes usable
  biology, and confirmation effort goes to the right molecules.
- **Keep the no-benchmark disclaimer.** Counter-intuitive but important: now that
  the brief contrasts our rigor against another practice, the line saying no
  head-to-head comparison has been run becomes more necessary, not less. Shorten
  it, do not remove it.

---

## 4. The manuscript: where the comparisons live

- **Introduction.** Add the indiscriminate-matching practice as a distinct failure
  mode alongside the coverage boundary already described, unnamed and framed as a
  class of practice. The paper can be more specific about the mechanism than the
  brief: annotation counts inflated by redundant adduct and isotope features, and
  by matches accepted below any useful similarity bar.
- **Section 1.1 keeps naming published tools.** MS2Query, SIRIUS, CANOPUS,
  Spec2Mol and MSNovelist are cited academic software, not companies, and this is
  where workflow comparison belongs per the instruction. Sharpen the descriptions;
  do not add a company anywhere.
- **A false-positive control paragraph in the Discussion.** The mechanisms are
  already in the paper but scattered across the mass-mismatch section, the level
  rubric and the coverage result. Pull them into one statement of how a wrong call
  is prevented from reaching a report, and state the limit honestly: claim
  limitation contains the error, it does not correct it, and 523 of 598
  mass-mismatched features have a better alternative sitting lower in the list.
- **Rigor framing in Limitations.** Keep every limitation, but lead the section by
  saying these are the bounds within which the results hold, rather than
  presenting them as apology.

---

## 5. One comment I plan to answer differently from how it is written

"Set up saying that we transparently provide a validated meaningful score."

Taken literally this would be false. The aggregate confidence score has a
leave-one-out Brier of 0.255, and both documents currently say it is exploratory
and not externally calibrated. The brief's release gate also blocks the phrase
"validated accuracy".

What is true, and is a stronger claim because it can be defended:

- We publish the calibration test and what it showed, which is more than is
  usually shown.
- The biochemical confidence signal does order assignments, gaining 28 percentage
  points of agreement between the whole cohort and its most confident fifth, so it
  can be used for prioritization.
- The structural grade limits what a report may claim rather than ranking one
  feature against another, and we say so.

Both documents will therefore say what the score means and how it was tested, and
will not call it validated. Flagged because it is a deliberate departure from the
comment as written.

---

## 6. Blocked on data

Three comments need inputs the repository does not contain. Unless the data
arrives first, each gets a `[REQUIRED: ...]` placeholder in the style both
documents already use, so they are visibly parked rather than silently dropped.

- **Matrix-specific comments, and matrix-specific impact.** The manuscript states
  twice that specimen-level matrix and species metadata were not retained in the
  reporting bundle, and the supplement already lists that metadata as required
  before submission.
- **A table of improved percentage coverage.** All 22 claims in the ledger concern
  the 3,897 unnamed features. There is no audited figure for the fraction of a
  study that gets named today, so there is no "before" to improve on, and the
  brief's gate rejects any statistic that does not resolve to the ledger.

---

## 7. Guardrails

- No company names in either document. The brief's scanner already blocks the tool
  names and stays as is; no relaxation is needed now.
- Every number in the brief must resolve to
  `docs/paper/supplement/claim_evidence.csv`. Cutting statistics is safe, adding
  is not. Table rows are free against the word budget but are still scanned for
  statistics.
- Avoid the banned overclaim phrases, "validated accuracy" in particular.
- The manuscript's audited counts must survive verbatim, and the retracted
  equivalence phrasing must not return.
- Watch pages, not words. The brief is at ten of ten, and level-raising should give
  a page back rather than consume one. Re-check with `--check-pages` after every
  pass.

---

## 8. Ordered steps

1. **brief-level** — compress the two over-detailed passages; measure the pages freed.
2. **brief-competitor** — add the unnamed third route and answer it with the audited controls.
3. **brief-rigor-l1** — restraint as spine, Level 1 framing, impact in the summary, shortened disclaimer.
4. **manuscript-competitor** — indiscriminate matching in the Introduction; sharpen 1.1.
5. **manuscript-rigor** — false-positive control paragraph; reframe the Limitations lead.
6. **score-wording** — say what the score means and how it was tested, in both documents.
7. **blocked-placeholders** — `[REQUIRED: ...]` for the matrix comments and the coverage table.
8. **records-verify** — revision log, review record, both guidance worksheets, all Word rebuilds, all gates.

---

## 9. Verification

```bash
python3 tools/build_manuscript_docx.py --source docs/paper/manuscript.md --output docs/paper/manuscript.docx
python3 tools/rebuild_paper_artifacts.py --brief-only
python3 tools/rebuild_paper_artifacts.py --verify-only
python3 tools/verify_public_brief.py --check-pages
python3 -m pytest -q tests/test_public_brief.py tests/test_manuscript_docx.py tests/test_figure_formats.py tests/test_x_compound_overview.py
ruff check tools/
```

Then reproduce the CI derived-artifact regeneration and confirm `git diff --exit-code`
is clean apart from the rebuilt Word files.
