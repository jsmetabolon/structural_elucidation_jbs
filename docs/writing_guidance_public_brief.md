# White Paper Writing Guidance, completed for the public capability brief

This is the writing-guidance worksheet with the answers filled in. Every heading, prompt, question and table below is reproduced from the guidance document as issued, in its original order and wording. Anything indented under a prompt is the answer to it.

The document being answered for is the customer-facing capability brief, `docs/public_brief/public_brief.md`. A companion worksheet, `docs/writing_guidance_manuscript.md`, answers the same prompts for the scientific manuscript. The two overlap deliberately, so each can be read on its own.

This worksheet is an internal working record about a customer-facing document. It is not itself for release: it names the things the brief withholds, and explains why, which is precisely the material the brief's release gate exists to keep out.

---

## White Paper Objective:

Explain the problem with current approaches to metabolite ID and LC-MS feature use, show how the new capability addresses that problem, provide evidence that it works, and explain why the improvement matters to researchers.

  A single untargeted study delivers thousands of reproducible features that no laboratory can name, and most analyses set them aside along with the biology they carry. The brief shows that each one can be turned into a documented structure hypothesis with an explicit grade, gives the measured recovery and coverage, and closes on the practical consequence: the output tells a research team which confirmation experiment is worth running.

**“problem → solution → evidence → impact”**

  The executive summary carries all four in miniature, so a reader who stops there still leaves with the argument, and the body then takes them in the same order across the unnamed-fraction section, the two established routes, what the capability does, and the evidence.

## A useful test for every section is:

- What does the reader need to understand here?
- What do we want the reader to believe after reading it?
- What evidence supports that belief?
- Why does it matter to the reader?

  Answered below for each of the brief's eleven sections.

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

  - *Understand:* three routes exist; two work and are bounded by what they cover, and the third trades rigor for volume.
  - *Believe:* the first two are limited by coverage and throughput rather than quality, and the third produces a number that is not a measure of coverage at all.
  - *Evidence:* an uncatalogued molecule is outside the space a catalogue search covers; expert confirmation is the scarcest resource in the workflow; and a list built by reporting whatever a broad database returns counts one compound repeatedly and prints weak matches like strong ones.
  - *Why it matters:* none of the three supplies a comparable statement of evidence for every feature, which is what triage requires, and the third actively transfers that work to the reader.

  **What the capability does**

  - *Understand:* four conceptual stages, detect, propose, weigh, report.
  - *Believe:* the value sits in the last two, because proposing plausible molecules is comparatively easy and deciding which proposals have earned the right to be shown is not.
  - *Evidence:* the grade is produced by the workflow rather than attached afterwards, which Figure 1 shows as a stage rather than a label.
  - *Why it matters:* the capability-to-value table converts each stage into a consequence for a research team.

  **How biochemical context is assigned**

  - *Understand:* pathway context reaches a feature by one of two paths, a lookup against a catalogued compound or an inference from chemically similar ones.
  - *Believe:* the two do not carry equal weight, and every report says which produced the call.
  - *Evidence:* hiding a reference metabolite's own labels and asking the capability to recover them shows the broad class coming back right most of the time and the subpathway less reliably. The agreement rates and the annotation ceiling that bounds them are reported in the manuscript rather than here, which is the level distinction the review asked for.
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

  - *Understand:* the grades follow the identification-level framework the field already uses.
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

  **Key terms**

  - *Understand:* what each of the twenty-four words the document leans on actually means.
  - *Believe:* the vocabulary is being used precisely rather than loosely.
  - *Evidence:* the definitions themselves, including the distinction between an identification and an annotation.
  - *Why it matters:* it sits at the back because a glossary is consulted rather than read, and a reader who needs one term should not have to find it mid-argument.

## Before writing the white paper, capture this info in a sentence or two each:

- The problem or unmet need

  A single untargeted LC-MS/MS study yields thousands of reproducible signals, but only a fraction can be matched to an authentic chemical standard, so the rest are delivered as unnamed features and set aside. That discarded remainder is reproducible and quantifiable, and it frequently contains the compounds that differentiate the study groups.

- Why existing approaches are insufficient

  Three routes exist and they fail differently. Matching against a spectral library is the most precise route short of an authentic standard and produces the named part of every deliverable, but its reach is set by the size of the catalogue, so genuinely novel and study-specific molecules sit outside it. Expert manual interpretation followed by purchase or synthesis is authoritative but scales to a handful of features per study. The third fails in the opposite direction: matching every detected signal against a broad database and reporting whatever returns inflates a list without adding identifications, because one compound is counted repeatedly across its ionized and isotopic forms and across chromatographic methods, and weak matches are presented like strong ones. What none of them supplies on its own is a comparable, explicit statement of how much evidence stands behind each call, for every feature.

- What the tool / GDP enhancement does

  For every unnamed feature it assembles candidate structures that reach beyond the reference catalogue, orders them by how well several independent lines of evidence agree, assigns biochemical pathway context, and issues a graded evidence report that proposes a specific structure only where the evidence supports one.

- What is technically differentiated

  Four things, in the order they matter: candidate generation is not bounded by catalogue membership; ranking draws on several independent evidence axes rather than a single similarity measure; the reporting layer is an explicit claim-control mechanism that withholds a structure when corroboration is insufficient; and the count is kept honest before any of that, by scoring against the neutral mass so ionized forms are not read as separate compounds, weighing isotope agreement, and counting a feature once across chromatographic methods. The grade is therefore produced by the process rather than asserted after it, and the deliverable reports compounds rather than returns.

- What evidence supports the claims

  On 92 controls whose structures are known independently, the correct structure was the single top-ranked proposal 47.8% of the time, appeared within the top three 64.1% of the time, and within the top ten 73.9% of the time. A leave-one-out test of the pathway assignment mechanism across 5,095 reference metabolites recovered the correct broad class 86.7% of the time. Across 3,897 delivered features, every feature received a report and pathway context, a structure was proposed for 2,726, and a structure was withheld for 1,171.

- What users can now do that they could not do before

  They can act on the unnamed fraction of their data instead of setting it aside: rank it, place it in a biochemical pathway, and take a short, chemically specific shortlist into confirmation rather than choosing a few features by intuition. Because the capability is computational it is applied uniformly to every feature, so the grades are comparable across a dataset and across studies.

- The three most important takeaways

  First, every unnamed feature receives a documented, graded result, and nothing is silently dropped from the deliverable. Second, what the workflow declines to claim matters as much as what it proposes: withholding is designed behaviour, and it is what makes the results that are shown worth acting on. Third, this is a prioritization instrument rather than a replacement for authentic-standard confirmation, and it makes confirmation cheaper by pointing it at the right molecules.

## Define the audience and objective

Describe the primary audience

  Study leads and data users, in biopharma, biotech, and academic research, who receive an untargeted metabolomics deliverable and have to decide what to do about the unnamed features in it. They are scientifically literate but are not mass-spectrometry method developers, and their decision is where to spend confirmation effort.

Describe a secondary audience

  Analytical chemists, metabolomics core directors, and computational biologists who will scrutinize how the confidence grades were derived and how far they can be trusted before recommending the capability internally.

After reading this paper, the reader should understand ______, believe ______, and feel confident that ______.

  After reading this paper, the reader should understand **that the unnamed fraction of an untargeted dataset can be turned into ranked, pathway-placed, evidence-graded structure hypotheses**, believe **that the grades are earned rather than asserted, because the workflow withholds a structure whenever the evidence does not support one**, and feel confident that **they can use the output to decide which features are worth the cost of confirmation, with the limits of each claim stated plainly**.

## Define the job of each section

The guidance names six sections. The brief has eleven, so each row names the section or sections that carry the job.

| Section | What it needs to accomplish |
|---|---|
| Executive summary | *Executive summary.* Deliver the whole problem-solution-evidence-impact arc in miniature and close on the three takeaways, so a reader who stops here still leaves with the argument. |
| Introduction | *The unnamed fraction of an untargeted dataset.* Establish the unnamed-feature problem and price it in both directions: biologically, discarded features often carry the differentiating signal; operationally, the conventional route to naming one is too slow to apply at scale. |
| Background | *What the established routes cover, and what they leave.* Explain how identification confidence is communicated today and what each existing route covers, locating the gap in catalogue coverage and expert throughput rather than in the quality of either method. Ends by naming the requirement the rest of the paper answers. |
| Tool Description | *What the capability does*, then *How biochemical context is assigned.* Describe the four stages so a non-specialist can follow the logic, each paragraph ending on the consequence for the reader; then explain the two paths by which a feature gets its pathway, attaching to each the evidence that bounds it. |
| Validation/evidence | *Evidence that it works.* Keep what was done, what was observed, what it means, and why it matters as four distinct moves, and state the limits of the evidence in the same breath as the result. |
| Conclusion | *Conclusion.* Restate the three takeaways as impact, define where the capability fits against authentic-standard confirmation, and avoid any comparative claim the evidence cannot carry. |

The brief carries four sections the guidance does not name, and each earns its place:

| Section | What it needs to accomplish |
|---|---|
| What you receive | Turn the deliverable into something concrete a reader can picture, and show that coverage and restraint arrive together. |
| How to read a result | Give the reader an operating manual: what each reporting outcome asserts and what the reasonable next step is. |
| Where this fits | Position the capability against the routes a customer already uses, without a comparative claim the evidence cannot carry. |
| Key terms | Define, at the back where a reader can look one up, the vocabulary that decides how a result should be read. |

| User value | (Always keep this in mind throughout) |
|---|---|
| Limitations | (Sprinkle these in to remain trustworthy and acknowledge not all is perfect) |

**User value, kept in mind throughout.** Value is the closing note of each substantive paragraph rather than a section of its own. The habit is easy to undo by accident, so it is worth stating: paragraphs end on the consequence for the reader, not on the mechanism that produced it.

**Limitations, sprinkled in.** Each sits next to the claim it qualifies rather than collected where a reader can skip them.

- The control set is modest at 92 compounds and was assembled from the material available rather than sampled at random, so the percentages describe that set rather than a rate that can be assumed for new data. Stated in the same section as the recovery figures.
- A ranked list is not a confirmed answer; a top-ten result reduces a boundless chemical space to a short list, which is not the same as a single settled answer.
- The confidence output limits claims rather than stating a probability of correctness, and it orders candidates within a feature but does not reliably order one feature against another.
- Pathway context orients rather than settles a classification: the broad class is worth acting on, the finer subpathway is a lead worth checking against local conventions.
- A looked-up pathway assignment is firmer than an inferred one, and part of any disagreement with another label source is vocabulary rather than chemistry.
- No head-to-head benchmark against other computational approaches exists, and none is claimed.

## Translate features into value

Avoid ending paragraphs on the feature. End on the consequence.

For each section, focus not only on explaining what you built, but also on the problem it addresses, the evidence that it works, and why the improvement matters to the user. The goal is not to remove technical depth, but to organize and translate it so that a broader scientific audience can follow the logic and understand the value.

| Technical capability | Why it matters to target audience |
|---|---|
| Ex. Evidence Scores | Several independent evidence axes are combined into one ordering, which reduces dependence on a single evidence source, so a candidate that looks good on one measure alone does not reach the top of the list unchallenged |
| A report for every unnamed feature | Nothing disappears between the instrument and the deliverable, so the unnamed fraction stays inside the biological interpretation instead of being written off |
| Candidate generation that reaches beyond the reference library | Novel and study-specific molecules can be described rather than skipped, which is where the differentiating biology often sits |
| Graded reporting aligned to community identification levels | Increases transparency and supports more informed prioritization, because the output is interpretable without learning a bespoke vocabulary |
| Deliberate withholding when evidence is insufficient | What is shown has earned the right to be shown, so confirmation effort is not spent chasing a structure the data never supported |
| Pathway context on every feature, including those without a structure | Extends the biological usefulness of features that fall short of confirmation, because an unnamed feature can still be placed in a narrative |
| Uniform automated application across the whole study | Enables consistent application at scale, so grades are comparable between the first feature and the three-thousandth, and between studies |

Weak:

The workflow calculates an evidence score using spectral and structural information.

Stronger:

The workflow combines spectral and structural evidence into an interpretable score, helping researchers prioritize the most strongly supported annotations for biological interpretation and follow-up.

  The same move, applied to this document. *Weak:* candidate structures are not restricted to compounds present in the reference library. *Stronger:* candidate generation reaches beyond the compounds a reference library already contains, which is what allows genuinely novel and study-specific molecules to be described rather than skipped, and that is often where the differentiating biology sits.

  The revision does two things at once: it ends on the consequence rather than the mechanism, and it states the capability positively rather than by negation of the incumbent method.

## Other Guidelines

When comparing to other tools or competitors, either have concrete quantitative comparison data or use defensible wording such as:

- provides additional evidence
- increases transparency
- supports more informed prioritization
- reduces dependence on a single evidence source
- enables consistent application at scale
- extends the biological utility of Level 2+ annotations

  No head-to-head benchmark has been run, so all six of these appear in the brief and none is used loosely. *Provides additional evidence* and *reduces dependence on a single evidence source* are paired in Where this fits. *Increases transparency* and *supports more informed prioritization* are paired in the capability-to-value table. *Enables consistent application at scale* is used only of the comparison with manual interpretation, where the claim is about throughput rather than accuracy. *Extends the biological utility* appears as the pathway-context row, phrased for a customer who does not use level numbers. The brief also states outright that no head-to-head comparison is claimed, so a reader cannot mistake the absence of one for an implied result.

## Separate explanation from proof

Scientists often mix method detail, results, interpretation, and value in one paragraph. Keep them distinct:

- What we did
- What we observed
- What it means
- Why the reader should care

  The Evidence section is written as four labelled paragraphs, so the separation is visible on the page rather than merely intended.

  - **What we did.** The capability was evaluated against internal controls whose true structures are known independently, so each proposal could be scored against a known answer rather than against another prediction.
  - **What we observed.** 47.8% at rank one, 64.1% within the top three, 73.9% within the top ten.
  - **What it means.** The control set is modest and was not sampled at random, and a ranked list is not a settled answer.
  - **Why the reader should care.** A shortlist with a measured hit rate changes the arithmetic of following a feature up.

This makes the argument much easier to follow.

## Make an outline before another full draft

Include:

- A section outline with one or two sentences describing the purpose of each section and how to accomplish that purpose

  1. **Executive summary.** The whole argument in miniature, closing on the three takeaways.
  2. **The unnamed fraction of an untargeted dataset.** The problem, priced biologically and operationally, closing on the absence of a consistent basis for comparison.
  3. **What the established routes cover, and what they leave.** What library matching and expert interpretation each cover, where the gap is, and the requirement the rest of the paper answers.
  4. **What the capability does.** The four stages at the level of concept, with the capability-to-value table.
  5. **How biochemical context is assigned.** The two paths, the evidence bounding each, and what that means for using them.
  6. **Evidence that it works.** Done, observed, means, matters, as four distinct moves.
  7. **What you receive.** The deliverable made concrete, with a report shown rather than described.
  8. **How to read a result.** What each grade asserts and the reasonable next step.
  9. **Where this fits.** Complementary positioning against library matching and standard confirmation.
  10. **Conclusion.** The takeaways restated as impact.
  11. **Key terms.** A consultable glossary of the vocabulary the document depends on, at the back rather than in the flow.
  12. **Availability.** Commercial, legal, and contact placeholders for marketing and legal review.

- Proposed figures and the message each figure supports

  Five figures: two hand-authored diagrams, two charts regenerated from the same audited data the manuscript uses, and one screen capture of a rendered report.

For each figure, complete this sentence:

This figure demonstrates that ______, which matters because ______.

  - `fig_workflow_public.svg` — This figure demonstrates that *each stage adds a distinct kind of evidence and the final stage decides how much may responsibly be claimed*, which matters because *it shows the confidence grade is a product of the workflow rather than a label attached at the end*.
  - `fig_pathway_public.svg` — This figure demonstrates that *biochemical context reaches a feature by one of two paths, one a lookup against a catalogued compound and the other an inference from chemically similar ones*, which matters because *the two carry very different weight, and a reader who cannot tell them apart cannot tell which pathway calls to lean on*.
  - `fig_recovery_public.svg` — This figure demonstrates that *the correct structure is usually near the top of a short ranked list on controls with independently known structures*, which matters because *it converts an unbounded chemical space into a shortlist a chemist can actually work through*.
  - `fig_report_public.png` — This figure demonstrates that *a delivered report states the proposed structure and the grade limiting it in the same place*, which matters because *it turns an abstract description of the deliverable into something the reader can picture, and shows the grade travelling with the call rather than sitting apart from it*.
  - `fig_coverage_public.svg` — This figure demonstrates that *every delivered feature is reported and placed in a pathway while a specific structure appears only where the evidence supports one*, which matters because *it shows completeness and restraint are delivered together rather than traded against each other*.

## When reviewing what’s written:

## Review 1: Argument

- Is the problem clear?

  Yes. It is stated in the executive summary, developed in its own section, and priced in two directions rather than asserted: biologically, the discarded features often carry the differentiating signal; operationally, the conventional route reaches a handful of features per study.

- Is the solution differentiated?

  Yes, on specific axes rather than on a general claim of being better: candidate generation is not bounded by catalogue membership, ranking uses several independent evidence axes, the reporting layer withholds a structure when corroboration is insufficient, and redundancy is removed before anything is counted. The withholding is the one competitors are least likely to have and it gets the most space, because restraint is the whole argument against a long unqualified list.

- Is the evidence sufficient?

  For the claims actually made, yes, and the claims were trimmed to what the evidence carries. Recovery is reported on a named cohort with its size given, coverage is reported across the full delivered set, and the pathway result is reported with the annotation ceiling that bounds it. Nothing is claimed about performance on new data, and no benchmark against other tools is asserted.

- Is the value explicit?

  Yes. Every substantive paragraph ends on the consequence for the reader, the capability-to-value table states it row by row, and the conclusion restates the three takeaways as impact rather than as a summary of numbers.

## Review 2: Audience

- Can a broad scientific reader understand it?

  Yes. The document assumes familiarity with untargeted metabolomics but not with mass-spectrometry method development, and it never requires the reader to hold a parameter or a threshold in mind.

- Are terms defined?

  Yes, and this was a gap until recently. A Key terms section at the back defines twenty-four words, including the vocabulary of library matching and the distinction between an identification and an annotation. It sits at the back deliberately: a glossary is consulted rather than read, and the first version was buried mid-document where nobody would look one up.

- Is unnecessary detail removed?

  Yes, and symmetrically. Model families, descriptors, scoring parameters, thresholds, internal library numbers, pipeline versions, and infrastructure names are all absent, and their absence is enforced by a scanner rather than by reviewer attention. Detail is dropped in both directions rather than keeping the flattering half. A later pass raised the level again: the pathway agreement rates and the annotation ceiling went back to the manuscript, leaving the brief to say what a reader should do with a broad class and a subpathway. Statistics now appear only where they carry the argument, which is recovery, coverage, and the counts behind the restraint claim.

## Review 3: Scientific rigor

- Are claims supported?

  Yes. Every statistic in the brief must resolve to a row in the audited claim ledger at `docs/paper/supplement/claim_evidence.csv`, and the release gate fails the build on any number that does not. That check is what keeps marketing copy from drifting away from what was measured.

- Are limitations stated?

  Yes, and next to the claims they qualify rather than collected where they can be skipped. The cohort caveat sits in the evidence section and the ranked-list caveat beside the recovery figures. They also do double duty: the same restraint that limits what the brief claims is what distinguishes it from an unqualified list, so the withheld count is presented as evidence rather than as apology.

- Are accuracy, confidence, identification, and annotation used precisely?

  Yes. The brief never calls the confidence output a probability of correctness, because leave-one-out calibration does not support that reading, and banned overclaim phrases are enforced by the gate. Identification and annotation are now defined against each other in Key terms, because a feature can carry an annotation without ever receiving an identification, and that is what lets a feature with no structure stay in the biological story. The brief also distinguishes a structure hypothesis from a confirmed identity, and biochemical confidence from structural confidence, since only the former orders one feature against another.

## Review 4: Style

- Sentence length

  Varied deliberately, with the load-bearing claims in short sentences so they are not buried in qualification.

- Repetition

  "Worth acting on" appeared three times, in the executive summary, the reporting-table discussion, and the conclusion; it now appears once and the other two make the point differently.

- Transitions

  Each section ends on the consequence that motivates the next, so the argument moves rather than accumulating.

- Figures and captions

  All five captions lead with the message the figure supports and follow with the detail, rather than describing content only. The screen capture additionally states that the example is a demonstration compound rather than a customer result.

- Final polish

  "Adding compounds to a reference library, which compounds the value" put a noun and a verb sense of the same word in adjacent clauses, and was rewritten; the same collision reappeared later in "the two routes compound rather than compete" and was caught again. A bullet that had run into the following paragraph without a blank line was separated. "Two points of context matter when reading those figures" sat directly beneath a figure and read as a reference to it, so it now says "percentages".

## Review 5: Fairness to the incumbent method

Added after a reader observed that the brief was dismissive of library matching. Not part of the guidance as issued.

- The finding

  Nothing in the brief was untrue, but the cumulative framing undersold the method that produces the named part of every deliverable and that this capability itself depends on. The background section was headed "Why the established routes fall short" and told the reader what a library search "cannot see".

- The correction

  The heading now describes what each route covers, the paragraph leads with library matching's strengths before locating its limit in catalogue size, and Where this fits says matching stays the first step, with well-supported hypotheses feeding back into the library so the two reinforce each other.

- Vocabulary as part of the fix

  Defining reference library, library match, catalogued compound, and authentic standard neutrally in Key terms is what keeps the rest of the document from implying a verdict on a route it relies on.

---

## Carried, for a human

Not part of the guidance. These need a decision and are not closable by review or by the gates.

- **Release placeholders.** Seven `[REQUIRED: ...]` placeholders remain, covering contact details, commercial availability, platform positioning, trademark usage, legal disclaimer, and the document version and approval reference. Run the gate with `--strict-human-inputs` to make them blocking once filled.
- **Page budget.** The brief renders to ten pages of ten, with roughly half of the last page free for the availability and legal text. Pages bind well before words do, so re-run the gate with `--check-pages` after every edit.
- **Brand, trademark, and claim substantiation.** Marketing and legal sign-off is still required regardless of gate status.
