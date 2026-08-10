# White Paper Writing Guidance: answers for the public capability brief

This worksheet answers every prompt in the internal white-paper writing guidance for one document: the customer-facing capability brief at `docs/public_brief/public_brief.md`. A companion worksheet, `docs/writing_guidance_manuscript.md`, answers the same prompts for the scientific manuscript. The two overlap deliberately, so each can be read on its own.

This is an internal working record about a customer-facing document. It is not itself for release: it names the things the brief withholds, and explains why, which is precisely the material the brief's release gate exists to keep out. The governing structure it records is shared with `docs/publication_outline.md`, and the review passes with `docs/publication_review.md`.

---

## White paper objective

Explain that a single untargeted study delivers thousands of reproducible features that no laboratory can name and that most analyses therefore set aside; show that this capability turns each one into a documented structure hypothesis carrying an explicit grade; give the measured evidence for how often the correct structure is found and how completely the delivered set is covered; and explain that the result is a prioritization instrument that tells a research team which confirmation experiment is worth running.

The argument runs problem, solution, evidence, impact, in that order, and the executive summary carries all four in miniature so a reader who stops there still leaves with it.

---

## The four-question test, section by section

The guidance asks four questions of every section. The brief has ten.

**Executive summary**

- *Understand:* the unnamed fraction of an untargeted dataset is real, quantifiable biology that is currently discarded.
- *Believe:* it can be turned into graded structure hypotheses without overclaiming.
- *Evidence:* 3,897 delivered features, every one reported, 2,726 with a structure and 1,171 with one deliberately withheld.
- *Why it matters:* nothing is silently dropped from the deliverable, so the unnamed fraction stays inside the biological interpretation.

**The unnamed fraction of an untargeted dataset**

- *Understand:* unnamed features are reproducible signals, not noise, and they frequently include the compounds that separate the study groups.
- *Believe:* setting them aside is expensive both biologically and operationally.
- *Evidence:* the conventional route to naming one feature reaches a handful per study, so the choice of which few to pursue is made without a consistent basis for comparison.
- *Why it matters:* the finding that would have mattered most can be the one never examined.

**What the established routes cover, and what they leave**

- *Understand:* library matching and expert interpretation both work, and each covers a different part of the problem.
- *Believe:* the gap between them is coverage and throughput, not quality.
- *Evidence:* a molecule that has not been catalogued is outside the space a catalogue search covers; expert confirmation is the scarcest resource in the workflow.
- *Why it matters:* what neither supplies on its own is a comparable statement of evidence for every feature, which is what triage requires.

**What the capability does**

- *Understand:* four conceptual stages, detect, propose, weigh, report.
- *Believe:* the value sits in the last two, because proposing plausible molecules is comparatively easy and deciding which proposals have earned the right to be shown is not.
- *Evidence:* the grade is produced by the workflow rather than attached afterwards, which Figure 1 shows as a stage rather than a label.
- *Why it matters:* the capability-to-value table converts each stage into a consequence for a research team.

**How biochemical context is assigned**

- *Understand:* pathway context reaches a feature by one of two paths, a lookup against a catalogued compound or an inference from chemically similar ones.
- *Believe:* the two do not carry equal weight, and every report says which produced the call.
- *Evidence:* 1,459 of the delivered features were looked up; the inferred path recovers the correct broad class 86.7% of the time when a reference metabolite's own labels are hidden, and agrees with a separate internal label source 52.2% of the time, against a ceiling of 62.5% on identical structures.
- *Why it matters:* the broad class is worth acting on, the finer subpathway is a lead, and a reader cannot judge either without knowing which path produced it.

**Evidence that it works**

- *Understand:* the capability was scored against controls whose true structures are known independently.
- *Believe:* the correct structure is usually near the top of a short ranked list.
- *Evidence:* 47.8% at rank one, 64.1% within the top three, 73.9% within the top ten, on 92 controls.
- *Why it matters:* a shortlist with a measured hit rate changes the arithmetic of following a feature up.

**What you receive**

- *Understand:* the deliverable is a self-contained report per feature plus a searchable index.
- *Believe:* it can be browsed by a biologist and audited by a chemist.
- *Evidence:* Figure 4 shows a report; Figure 5 shows that coverage and restraint arrive together across the delivered set.
- *Why it matters:* the reader can picture what lands, rather than taking a description of it on trust.

**How to read a result**

- *Understand:* the grades follow the identification-level framework the field already uses, and the key terms mean specific things.
- *Believe:* withholding is designed behaviour rather than a gap in coverage.
- *Evidence:* the last two rows of the reporting-outcome table are the ones that carry the argument, and the confidence output limits claims rather than stating a probability of correctness.
- *Why it matters:* the reader gets an operating manual, including where the ordering can be trusted and where it cannot.

**Where this fits**

- *Understand:* this complements library matching and authentic-standard confirmation rather than replacing either.
- *Believe:* the two routes reinforce each other, because well-supported hypotheses are candidates for entering the library.
- *Evidence:* comparative statements are limited to defensible wording, since no head-to-head benchmark exists.
- *Why it matters:* it tells a prospective customer where the capability sits in work they already do.

**Conclusion**

- *Understand:* the unnamed fraction is unexamined biology, and this makes it examinable.
- *Believe:* the three takeaways, restated as impact.
- *Evidence:* everything already established, restated without new numbers.
- *Why it matters:* the capability does not remove the final experiment, it indicates which one to run.

---

## Before writing: capture in a sentence or two each

### The problem or unmet need

A single untargeted LC-MS/MS study yields thousands of reproducible signals, but only a fraction can be matched to an authentic chemical standard, so the rest are delivered as unnamed features and set aside. That discarded remainder is reproducible and quantifiable, and it frequently contains the compounds that differentiate the study groups.

### Why existing approaches are insufficient

Matching against a spectral library is the most precise route short of an authentic standard and produces the named part of every deliverable, but its reach is set by the size of the catalogue, so genuinely novel and study-specific molecules sit outside it. Expert manual interpretation followed by purchase or synthesis is authoritative but scales to a handful of features per study. What neither supplies on its own is a comparable, explicit statement of how much evidence stands behind each call, for every feature.

### What the tool / GDP enhancement does

For every unnamed feature it assembles candidate structures that reach beyond the reference catalogue, orders them by how well several independent lines of evidence agree, assigns biochemical pathway context, and issues a graded evidence report that proposes a specific structure only where the evidence supports one.

### What is technically differentiated

Three things, in the order they matter: candidate generation is not bounded by catalogue membership; ranking draws on several independent evidence axes rather than a single similarity measure; and the reporting layer is an explicit claim-control mechanism that withholds a structure when corroboration is insufficient. The grade is therefore produced by the process rather than asserted after it.

### What evidence supports the claims

On 92 controls whose structures are known independently, the correct structure was the single top-ranked proposal 47.8% of the time, appeared within the top three 64.1% of the time, and within the top ten 73.9% of the time. A leave-one-out test of the pathway assignment mechanism across 5,095 reference metabolites recovered the correct broad class 86.7% of the time. Across 3,897 delivered features, every feature received a report and pathway context, a structure was proposed for 2,726, and a structure was withheld for 1,171.

### What users can now do that they could not do before

They can act on the unnamed fraction of their data instead of setting it aside: rank it, place it in a biochemical pathway, and take a short, chemically specific shortlist into confirmation rather than choosing a few features by intuition. Because the capability is computational it is applied uniformly to every feature, so the grades are comparable across a dataset and across studies.

### The three most important takeaways

1. Every unnamed feature receives a documented, graded result. Nothing is silently dropped from the deliverable.
2. What the workflow declines to claim matters as much as what it proposes. Withholding is designed behaviour, and it is what makes the results that are shown worth acting on.
3. This is a prioritization instrument, not a replacement for authentic-standard confirmation. It makes confirmation cheaper by pointing it at the right molecules.

---

## Audience and objective

### Primary audience

Study leads and data users, in biopharma, biotech, and academic research, who receive an untargeted metabolomics deliverable and have to decide what to do about the unnamed features in it. They are scientifically literate but are not mass-spectrometry method developers, and their decision is where to spend confirmation effort.

### Secondary audience

Analytical chemists, metabolomics core directors, and computational biologists who will scrutinize how the confidence grades were derived and how far they can be trusted before recommending the capability internally.

### The objective sentence

After reading this paper, the reader should understand **that the unnamed fraction of an untargeted dataset can be turned into ranked, pathway-placed, evidence-graded structure hypotheses**, believe **that the grades are earned rather than asserted, because the workflow withholds a structure whenever the evidence does not support one**, and feel confident that **they can use the output to decide which features are worth the cost of confirmation, with the limits of each claim stated plainly**.

---

## The job of each section

The guidance names six sections. The brief has ten, so the mapping is given alongside.

| Guidance section | Section in the brief | What it needs to accomplish |
|---|---|---|
| Executive summary | Executive summary | Deliver the whole problem-solution-evidence-impact arc in miniature and close on the three takeaways, so a reader who stops here still leaves with the argument. |
| Introduction | The unnamed fraction of an untargeted dataset | Establish the unnamed-feature problem and price it in both directions: biologically, discarded features often carry the differentiating signal; operationally, the conventional route to naming one is too slow to apply at scale. |
| Background | What the established routes cover, and what they leave | Explain how identification confidence is communicated in the field today and what each existing route covers, locating the gap in catalogue coverage and expert throughput rather than in the quality of either method. Ends by naming the requirement the rest of the paper answers. |
| Tool description | What the capability does; How biochemical context is assigned | Describe the four stages so a non-specialist can follow the logic, each paragraph ending on the consequence for the reader; then explain the two paths by which a feature gets its pathway, attaching to each the evidence that bounds it. |
| Validation and evidence | Evidence that it works | Keep what was done, what was observed, what it means, and why it matters as four distinct moves, and state the limits of the evidence in the same breath as the result. |
| Conclusion | Conclusion | Restate the three takeaways as impact, define where the capability fits against authentic-standard confirmation, and avoid any comparative claim the evidence cannot carry. |
| (additional) | What you receive | Turn the deliverable into something concrete a reader can picture, and show that coverage and restraint arrive together. |
| (additional) | How to read a result | Give the reader an operating manual: the key terms in plain language, then what each reporting outcome asserts and what the reasonable next step is. |
| (additional) | Where this fits | Position the capability against the routes a customer already uses, without a comparative claim the evidence cannot carry. |
| (additional) | Availability | Carry the required commercial, legal, and contact placeholders through to marketing and legal review. |

### User value, kept in mind throughout

User value is the closing note of each substantive paragraph rather than a section of its own. The habit is easy to undo by accident, so it is worth stating: paragraphs end on the consequence for the reader, not on the mechanism that produced it.

### Limitations, sprinkled in to remain trustworthy

Limitations are placed next to the claims they qualify rather than collected into a list the reader can skip.

- The control set is modest at 92 compounds and was assembled from the material available rather than sampled at random, so the percentages describe that set rather than a rate that can be assumed for new data. Stated in the same section as the recovery figures.
- A ranked list is not a confirmed answer; a top-ten result reduces a boundless chemical space to a short list, which is not the same as a single settled answer.
- The confidence output limits claims rather than stating a probability of correctness, and it orders candidates within a feature but does not reliably order one feature against another.
- Pathway context orients rather than settles a classification: the broad class is worth acting on, the finer subpathway is a lead worth checking against local conventions.
- A looked-up pathway assignment is firmer than an inferred one, and part of any disagreement with another label source is vocabulary rather than chemistry.
- No head-to-head benchmark against other computational approaches exists, and none is claimed.

---

## Translating features into value

| Technical capability | Why it matters to the target audience |
|---|---|
| A report for every unnamed feature | Nothing disappears between the instrument and the deliverable, so the unnamed fraction stays inside the biological interpretation instead of being written off. |
| Candidate generation that reaches beyond the reference library | Novel and study-specific molecules can be described rather than skipped, which is where the differentiating biology often sits. |
| Several independent evidence axes combined into one ordering | Reduces dependence on a single evidence source, so a candidate that looks good on one measure alone does not reach the top of the list unchallenged. |
| Graded reporting aligned to community identification levels | Increases transparency and supports more informed prioritization, because the output is interpretable without learning a bespoke vocabulary. |
| Deliberate withholding when evidence is insufficient | What is shown has earned the right to be shown, so confirmation effort is not spent chasing a structure the data never supported. |
| Pathway context on every feature, including those without a structure | Extends the biological usefulness of features that fall short of confirmation, because an unnamed feature can still be placed in a narrative. |
| Uniform automated application across the whole study | Enables consistent application at scale, so grades are comparable between the first feature and the three-thousandth, and between studies. |
| The supporting evidence shown alongside the call | A reviewer can see what the call rests on instead of trusting a single score, which is what makes the result defensible to a third party. |

---

## Other guidelines

### Comparative wording

No head-to-head benchmark has been run, so every comparative statement uses defensible wording. All six of the guidance's suggested phrases appear in the brief, and none of them is used loosely:

- *provides additional evidence* — "the capability contributes additional evidence", in Where this fits.
- *increases transparency* and *supports more informed prioritization* — paired in the capability-to-value table.
- *reduces dependence on a single evidence source* — in both the table and Where this fits.
- *enables consistent application at scale* — used of the comparison with manual interpretation, where the claim is about throughput rather than accuracy.
- *extends the biological utility of annotations that fall short of confirmation* — the pathway-context row of the table.

The brief also states outright that no head-to-head comparison is claimed, so a reader cannot mistake the absence of one for an implied result.

### Separating explanation from proof

The four moves are kept distinct, and the Evidence section is written as four labelled paragraphs so the separation is visible on the page rather than merely intended:

- **What we did.** The capability was evaluated against internal controls whose true structures are known independently, so each proposal could be scored against a known answer rather than against another prediction.
- **What we observed.** 47.8% at rank one, 64.1% within the top three, 73.9% within the top ten.
- **What it means.** The control set is modest and was not sampled at random, and a ranked list is not a settled answer.
- **Why the reader should care.** A shortlist with a measured hit rate changes the arithmetic of following a feature up.

### Feature to value, worked

*Weak:* Candidate structures are not restricted to compounds present in the reference library.

*Stronger:* Candidate generation reaches beyond the compounds a reference library already contains, which is what allows genuinely novel and study-specific molecules to be described rather than skipped, and that is often where the differentiating biology sits.

The revision does two things at once: it ends on the consequence rather than the mechanism, and it states the capability positively rather than by negation of the incumbent method.

---

## Outline

### Section outline

1. **Executive summary.** The whole argument in miniature, closing on the three takeaways.
2. **The unnamed fraction of an untargeted dataset.** The problem, priced biologically and operationally, closing on the absence of a consistent basis for comparison.
3. **What the established routes cover, and what they leave.** What library matching and expert interpretation each cover, where the gap is, and the requirement the rest of the paper answers.
4. **What the capability does.** The four stages at the level of concept, with the capability-to-value table.
5. **How biochemical context is assigned.** The two paths, the evidence bounding each, and what that means for using them.
6. **Evidence that it works.** Done, observed, means, matters, as four distinct moves.
7. **What you receive.** The deliverable made concrete, with a report shown rather than described.
8. **How to read a result.** Key terms, then what each grade asserts and the reasonable next step.
9. **Where this fits.** Complementary positioning against library matching and standard confirmation.
10. **Conclusion.** The takeaways restated as impact.
11. **Availability.** Commercial, legal, and contact placeholders for marketing and legal review.

### Figures

Five figures, in the order they appear. Two are hand-authored diagrams, two are charts regenerated from the same audited data the manuscript uses, and one is a screen capture of a rendered report.

- `fig_workflow_public.svg` — This figure demonstrates that *each stage adds a distinct kind of evidence and the final stage decides how much may responsibly be claimed*, which matters because *it shows the confidence grade is a product of the workflow rather than a label attached at the end*.
- `fig_pathway_public.svg` — This figure demonstrates that *biochemical context reaches a feature by one of two paths, one a lookup against a catalogued compound and the other an inference from chemically similar ones*, which matters because *the two carry very different weight, and a reader who cannot tell them apart cannot tell which pathway calls to lean on*.
- `fig_recovery_public.svg` — This figure demonstrates that *the correct structure is usually near the top of a short ranked list on controls with independently known structures*, which matters because *it converts an unbounded chemical space into a shortlist a chemist can actually work through*.
- `fig_report_public.png` — This figure demonstrates that *a delivered report states the proposed structure and the grade limiting it in the same place*, which matters because *it turns an abstract description of the deliverable into something the reader can picture, and shows the grade travelling with the call rather than sitting apart from it*.
- `fig_coverage_public.svg` — This figure demonstrates that *every delivered feature is reported and placed in a pathway while a specific structure appears only where the evidence supports one*, which matters because *it shows completeness and restraint are delivered together rather than traded against each other*.

---

## Review passes

### Review 1 — Argument

**Is the problem clear?** Yes. It is stated in the executive summary, developed in its own section, and priced in two directions rather than asserted: biologically, the discarded features often carry the differentiating signal; operationally, the conventional route reaches a handful of features per study.

**Is the solution differentiated?** Yes, on three specific axes rather than on a general claim of being better: candidate generation is not bounded by catalogue membership, ranking uses several independent evidence axes, and the reporting layer withholds a structure when corroboration is insufficient. The third is the one competitors are least likely to have, and it is given the most space.

**Is the evidence sufficient?** For the claims actually made, yes, and the claims were trimmed to what the evidence carries. Recovery is reported on a named cohort with its size given, coverage is reported across the full delivered set, and the pathway result is reported with the annotation ceiling that bounds it. Nothing is claimed about performance on new data, and no benchmark against other tools is asserted.

**Is the value explicit?** Yes. Every substantive paragraph ends on the consequence for the reader, the capability-to-value table states it row by row, and the conclusion restates the three takeaways as impact rather than as a summary of numbers.

### Review 2 — Audience

**Can a broad scientific reader understand it?** Yes. The document assumes familiarity with untargeted metabolomics but not with mass-spectrometry method development, and it never requires the reader to hold a parameter or a threshold in mind.

**Are terms defined?** Yes, and this was a gap until recently. A key-terms table now fixes the ten words that carry the most weight, including the vocabulary of library matching, before the reporting-outcome table uses them.

**Is unnecessary detail removed?** Yes, and symmetrically. Model families, descriptors, scoring parameters, thresholds, internal library numbers, pipeline versions, and infrastructure names are all absent, and their absence is enforced by a scanner rather than by reviewer attention. Detail is dropped in both directions rather than keeping the flattering half.

### Review 3 — Scientific rigor

**Are claims supported?** Yes. Every statistic in the brief must resolve to a row in the audited claim ledger at `docs/paper/supplement/claim_evidence.csv`, and the release gate fails the build on any number that does not. That check is what keeps marketing copy from drifting away from what was measured.

**Are limitations stated?** Yes, and next to the claims they qualify rather than collected where they can be skipped. The cohort caveat sits in the evidence section, the ranked-list caveat beside the recovery figures, the annotation-vocabulary caveat beside the pathway agreement figures.

**Are accuracy, confidence, identification, and annotation used precisely?** Yes. The brief never calls the confidence output a probability of correctness, because leave-one-out calibration does not support that reading, and banned overclaim phrases are enforced by the gate. It distinguishes a structure hypothesis from a confirmed identity, and it distinguishes the biochemical confidence grade from the structural one, since only the former orders one feature against another.

### Review 4 — Style

**Sentence length.** Varied deliberately, with the load-bearing claims in short sentences so they are not buried in qualification.

**Repetition.** "Worth acting on" appeared three times, in the executive summary, the reporting-table discussion, and the conclusion; it now appears once and the other two make the point differently.

**Transitions.** Each section ends on the consequence that motivates the next, so the argument moves rather than accumulating.

**Figures and captions.** All five captions lead with the message the figure supports and follow with the detail, rather than describing content only. The two screenshots additionally state that the example is a demonstration compound rather than a customer result.

**Final polish.** "Adding compounds to a reference library, which compounds the value" put a noun and a verb sense of the same word in adjacent clauses, and was rewritten; the same collision reappeared later in "the two routes compound rather than compete" and was caught again. A bullet that had run into the following paragraph without a blank line was separated. "Two points of context matter when reading those figures" sat directly beneath a figure and read as a reference to it, so it now says "percentages".

### Review 5 — Fairness to the incumbent method

Added after a reader observed that the brief was dismissive of library matching.

**The finding.** Nothing in the brief was untrue, but the cumulative framing undersold the method that produces the named part of every deliverable and that this capability itself depends on. The background section was headed "Why the established routes fall short" and told the reader what a library search "cannot see".

**The correction.** The heading now describes what each route covers, the paragraph leads with library matching's strengths before locating its limit in catalogue size, and Where this fits says matching stays the first step, with well-supported hypotheses feeding back into the library so the two reinforce each other.

**Vocabulary as part of the fix.** Defining reference library, library match, and authentic standard neutrally in the key-terms table is what keeps the rest of the document from implying a verdict on a route it relies on.

---

## Carried, for a human

- **Release placeholders.** Seven `[REQUIRED: ...]` placeholders remain, covering contact details, commercial availability, platform positioning, trademark usage, legal disclaimer, and the document version and approval reference. Run the gate with `--strict-human-inputs` to make them blocking once filled.
- **Page budget.** The brief renders to nine pages of ten, and pages now bind before words do. A caption gaining a single line has been enough to cost a whole page, so re-run the gate with `--check-pages` after any edit near a figure.
- **Brand, trademark, and claim substantiation.** Marketing and legal sign-off is still required regardless of gate status.
