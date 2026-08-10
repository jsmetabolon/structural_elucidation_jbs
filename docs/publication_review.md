# Publication review record

The white-paper writing guidance specifies four review passes. This file records what
each pass checked across both publication surfaces, what it found, and what changed as
a result. Items under "Carried" are open and need a human decision; automated checks
cannot close them.

Documents reviewed: `docs/public_brief/public_brief.md` and `docs/paper/manuscript.md`.

## Review 1 — Argument

*Is the problem clear? Is the solution differentiated? Is the evidence sufficient? Is
the value explicit?*

**White paper.** The problem is stated twice at different depths: compressed in the
executive summary and priced out in the section that follows, in both the biological
direction (discarded features often carry the differentiating signal) and the
operational one (manual naming reaches a handful of features per study). Differentiation
now rests on a stated baseline rather than assertion: the Background section names the
two established routes and the specific place each stops short, so the capability is
positioned against something concrete. Value is explicit in three places — the
capability-to-value table, the closing sentence of each substantive paragraph, and the
conclusion.

*Found and changed.* The document previously had no Background section and no
Conclusion, so the argument opened without a baseline and ended on an availability
notice. Both were added.

**Manuscript.** The argument is a scientific one and was already sound. The gap was that
the introduction ended on the contributions themselves rather than on their consequence,
and Applications listed uses without connecting them to results. Both now close on what
an analyst can do that was previously impractical.

## Review 2 — Audience

*Can a broad scientific reader understand it? Are terms defined? Is unnecessary detail
removed?*

**White paper.** Two terms were reaching the reader undefined.

- "Leave-one-out test" was replaced with a plain description of the procedure: hiding a
  reference metabolite's own labels and asking the capability to recover them.
- "Super-pathway" is now defined in place as the broad biochemical class, with
  "subpathway" positioned as the finer level.

Unnecessary detail is constrained structurally rather than editorially: the release gate
blocks model families, descriptors, scoring parameters, thresholds, internal library
numbers, and infrastructure names, so the white paper describes approach classes rather
than named tools.

**Manuscript.** The audience is specialist and the detail level is appropriate to it. No
changes were needed beyond the terminology paragraph noted below, which serves both
audience clarity and rigor.

## Review 3 — Scientific rigor

*Are claims supported? Are limitations stated? Are accuracy, confidence, identification
and annotation used precisely?*

**Claims supported.** Every numeral in the white paper resolves to a row in
`docs/paper/supplement/claim_evidence.csv`; the gate fails the build otherwise. The
manuscript's audited fractions are asserted by `tools/verify_paper_package.py`.

*Found and changed.* The manuscript's coverage counts had drifted from the committed
audit artifact — assignment statuses, pathway-confidence tiers, and the best-reference
similarity median all disagreed with `data/pathway/audit_pathway_summary.json`. The
prose now states the audited values and reports the three assignment statuses
separately instead of folding exact matches into the strong-similarity count. The
verifier was extended to assert these counts, since nothing had been checking them.

**Limitations stated.** In the white paper, limitations sit next to the claims they
qualify rather than in a separate section: control-set size and selection appear with
the recovery result, the ranked-list caveat with the top-ten figure, taxonomy
disagreement with the pathway result, and the refusal to present confidence as a
probability of correctness with the reporting table. The manuscript retains its
dedicated Limitations section.

**Terminology.** A paragraph in Methods now fixes four evaluation terms as distinct:
*recovery* against a known structure, *accuracy* only where a label was withheld from
the mechanism that produced the prediction, *agreement* with an imperfect independent
label source, and *consistency* with the database a value was copied from. It also
separates *identification*, *annotation*, and *structure hypothesis*, and keeps
identification-confidence tiers distinct from pathway-confidence tiers.

*Found and changed.* Auditing the manuscript against its own new definitions caught two
violations, both in figure captions: "accuracy advantage" where structural recovery was
meant, and "a single pathway accuracy" where a performance figure was meant. Both were
corrected.

*Accepted simplification.* The white paper reports a threshold-based structural match as
"correct". That simplification is now made explicit in the text — a proposal counts as
correct when it matches the known answer within a fixed structural-similarity threshold
— with the precise definition carried by the manuscript.

*Later addition — disclosing the pathway mechanism.* The brief was subsequently extended
to describe how biochemical context is assigned, and doing so changed what the evidence
owed the reader. The brief had reported only the held-out validation of the inference
mechanism, 86.7% and 71.1%, which is the friendlier of the two available measurements. It
now also reports agreement measured on the inferred path alone against a separate internal
label source: 52.2% super-pathway and 18.8% subpathway across 644 controls. Presenting the
mechanism while quoting only the flattering figure would have breached the second rule the
document exists under.

*Found while doing so — a missing denominator.* Reporting those figures raised the question
of why the two measurements differ so much, and the answer turned out to be partly
arithmetic. They are scored against different biochemical vocabularies, and the audit
already recorded how far those vocabularies agree with each other on identical structures:
62.5% at the super-pathway level, which both documents carried, and 34.8% at the subpathway
level, which neither did. Those rates bound what any predictor targeting the reference
taxonomy can score against the independent labels. Against the bound, the inferred path
reaches roughly 83% of the achievable super-pathway rate and 54% of the achievable
subpathway rate, so the apparent collapse at the finer level is substantially a vocabulary
mismatch. Both documents now carry both bounds; the manuscript also states the assumption
the normalization rests on, namely that the disagreement rate measured on reference
structures transfers to the non-reference cohort, which is untested. The missing figure was
added to the claim ledger so the brief could cite it at all.

That addition also amended the document's governing disclosure rule from "never how" to
"the conceptual how, never the parameters". The rationale, and the deliberate asymmetry
with candidate ranking, are recorded in `docs/public_brief/README.md`.

## Review 4 — Style

*Sentence length, repetition, transitions, figures and captions, final polish.*

- **Repetition.** "Worth acting on" appeared three times in the white paper across the
  executive summary, the reporting table discussion, and the conclusion. It now appears
  once; the other two make the same point differently.
- **Polish.** "Adding compounds to a reference library, which compounds the value" put a
  noun and a verb sense of the same word in adjacent clauses; rewritten. "Two points of
  context matter when reading those figures" sat directly beneath the recovery figure and
  read as a reference to the figure rather than to the percentages; changed to
  "percentages".
- **Fragile cross-reference.** "The most important line in that table is the fourth"
  pointed at a row by position and would silently become wrong if a row were added. It
  now refers to the last two rows, which is also more accurate, since both withhold a
  structure.
- **Captions.** All nine captions across the two documents now lead with the message the
  figure supports and follow with the methodological detail, rather than describing
  content only.

## Review 5 — Fairness to the incumbent method

*A later pass, prompted by how both documents treated spectral-library matching.*

- **The finding.** Neither document said anything untrue about library matching, but the
  cumulative framing was dismissive of it. The manuscript opened its second paragraph on
  what library-retrieval-only workflows "cannot" do and called the route that follows a
  "bottleneck"; it described MS2Query as "bounded by" its library while the generative
  tools were "not restricted to database membership" and "independent of spectral
  libraries", which reads as escape from a defect. The white paper's background section
  was headed "Why the established routes fall short" and told the reader what a library
  search "cannot see". Library matching produces the named fraction of every deliverable,
  is the first step of this workflow, and supplies one of the evidence axes the reranker
  reads, so the framing undersold a method the work depends on.
- **The correction.** The boundary is now coverage rather than capability, stated once in
  each document and not repeated. The two routes are described as answering different
  questions. The Discussion and Where this fits both say that a library match remains the
  most direct answer wherever the compound is catalogued, and that well-corroborated
  hypotheses feed back into the library, so the two compound rather than compete.
- **What was deliberately not softened.** Exact-source rank-1 recovery of 38.1% against
  56.0% for predicted-source stays, with its overlapping intervals, and so does the
  conclusion that catalog provenance is not confirmation. Two sentences were rescoped
  rather than removed: Results no longer asserts that provenance "confers no recovery
  advantage", which claimed more than 42 controls can support, and says instead that the
  cohort does not resolve a difference in either direction. That is a narrower claim, not
  a friendlier one, and it was checked against the standing ban on equivalence language.
- **Vocabulary as part of the fix.** A neutral definition of catalog provenance,
  retrieval, analogue search and de novo generation is what keeps "exact source" reading
  as a record of where a candidate came from rather than as a verdict on the route. The
  glossary in Methods and the key-terms table in the white paper carry that work.
- **Captions.** Eleven captions now, with the two report screenshots following the same
  message-first rule and each stating that the example is a demonstration rather than a
  delivered customer result.

## Carried

These need a human and are not closable by review or by the gates.

- **Release placeholders.** The white paper carries seven `[REQUIRED: ...]` placeholders
  and the manuscript eight, covering contact details, commercial availability,
  trademark and legal disclaimer, platform positioning, authorship, funding,
  disclosures, ethics, and data access. Run the verifiers with `--strict-human-inputs`
  to make them blocking once filled.
- **Platform positioning.** The writing guidance frames the work as a tool or platform
  enhancement, but no platform positioning is documented in this repository. A
  placeholder was added rather than a guess, and it also asks whether existing
  deliverables may be named externally.
- **Comparative claims.** No head-to-head benchmark exists, so both documents restrict
  comparative language to defensible wording. If quantitative comparison data becomes
  available, those passages can be strengthened.
- **Brand, trademark, and claim-substantiation review.** Marketing and legal sign-off
  remains required regardless of gate status.
