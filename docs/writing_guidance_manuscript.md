# White Paper Writing Guidance, completed for the scientific manuscript

This is the writing-guidance worksheet with the answers filled in. Every heading, prompt, question and table below is reproduced from the guidance document as issued, in its original order and wording. Anything indented under a prompt is the answer to it.

The document being answered for is the scientific manuscript, `docs/paper/manuscript.md`, prepared for Nature initial submission. A companion worksheet, `docs/writing_guidance_public_brief.md`, answers the same prompts for the customer-facing capability brief. The two overlap deliberately, so each can be read on its own.

The guidance was written for a white paper, and this document is not one. Three differences shape the answers throughout, and they are worth stating before the prompts rather than discovering mid-worksheet.

- **The audience is adversarial by design.** A reviewer is looking for the claim that outruns its evidence, so the guidance's instruction to make value explicit is subordinate here to making the bound on each claim explicit.
- **Limitations are not sprinkled.** The guidance asks for them to be distributed to maintain trust. The manuscript does that *and* carries a dedicated Limitations section, because a journal expects one and because several of the limits are structural rather than local.
- **The section names differ.** Tool description maps to Methods, validation and evidence to Results, and the manuscript adds Discussion, Limitations, Future work and a Glossary, which the guidance does not name. The mapping is given in the section-jobs table.

---

## White Paper Objective:

Explain the problem with current approaches to metabolite ID and LC-MS feature use, show how the new capability addresses that problem, provide evidence that it works, and explain why the improvement matters to researchers.

  Establish that untargeted LC-MS/MS detects far more reproducible features than can be identified against authentic standards and that the remainder is set aside for want of a comparable confidence statement; describe a workflow that combines retrieval and de novo candidate generation with evidence-aware reranking, biochemical taxonomy assignment and MSI-style claim control; report what it recovers and where it does not, with uncertainty and provenance attached to every figure; and state plainly that the contribution is a reporting framework that makes the unidentified fraction usable rather than a claim of higher accuracy.

**“problem → solution → evidence → impact”**

  The same order, but the impact statement is deliberately modest: the closing claim is about what an analyst can now do, not about outperforming anything. The abstract carries all four inside 200 words.

## A useful test for every section is:

- What does the reader need to understand here?
- What do we want the reader to believe after reading it?
- What evidence supports that belief?
- Why does it matter to the reader?

  Answered below for each of the manuscript's nine sections.

  **Abstract**

  - *Understand:* what was evaluated, on what, and what was found.
  - *Believe:* the contribution is a reporting framework rather than a higher-accuracy predictor.
  - *Evidence:* recovery at three ranks with confidence intervals, a Brier score, and taxonomy results separated by provenance, all within 200 words.
  - *Why it matters:* it makes the unidentified fraction usable instead of discarded.

  **Introduction**

  - *Understand:* how identification confidence is communicated today, and what the established routes cover.
  - *Believe:* the gap is coverage and throughput, not the quality of a library match.
  - *Evidence:* the identification-level frameworks that formalize the problem, and the toolkit landscape that spans retrieval to de novo generation.
  - *Why it matters:* it names the requirement the rest of the paper answers, and states the three contributions in terms of what an analyst can do.

  **Methods**

  - *Understand:* what the workflow does at each stage and how each result was scored.
  - *Believe:* the evaluation design separates prediction from lookup, so nothing circular is reported as accuracy.
  - *Evidence:* four evaluation terms are fixed as distinct measurements, and each is used only where it applies.
  - *Why it matters:* a reader can attribute any result below to a specific stage rather than to the pipeline as a whole.

  **Results**

  - *Understand:* what was measured, cohort by cohort.
  - *Believe:* the figures mean what their provenance allows and no more.
  - *Evidence:* Wilson intervals throughout, cohort sizes named in every claim, and the pathway result split three ways by evaluation context.
  - *Why it matters:* quoting a single pathway figure without its provenance would misstate the result by tens of percentage points.

  **Discussion**

  - *Understand:* the two results that reframe how the outputs should be read.
  - *Believe:* claim limitation is the defensible mechanism, and it contains error rather than correcting it.
  - *Evidence:* poor cross-validated calibration, and a mass-consistent alternative already sitting lower in 523 of 598 candidate lists.
  - *Why it matters:* it tells an analyst which of the two confidence signals can be used for prioritization and which cannot.

  **Limitations**

  - *Understand:* where the evidence stops.
  - *Believe:* the authors know where it stops.
  - *Evidence:* seven named limits, each tied to the specific analysis it bounds.
  - *Why it matters:* it is what allows the rest of the paper to be taken at face value.

  **Conclusion**

  - *Understand:* what the workflow produces across the full delivered set.
  - *Believe:* the central result is the reporting framework.
  - *Evidence:* coverage and recovery restated without new numbers.
  - *Why it matters:* the unidentified fraction becomes a graded, prioritizable set of hypotheses.

  **Future work**

  - *Understand:* what would settle the questions this study leaves open.
  - *Believe:* the gaps are known and addressable.
  - *Evidence:* named cohorts, a prospective hard mass gate, and a formal benchmark design.
  - *Why it matters:* it tells a reader what the next paper has to show.

  **Glossary**

  - *Understand:* the precise sense of every term the paper depends on, across the measurement, candidate provenance, structural comparison, evaluation, claims, taxonomy and failure modes.
  - *Believe:* the vocabulary is being used consistently, and the four evaluation senses in particular are not interchangeable.
  - *Evidence:* the definitions themselves, with audited numbers attached where a term means nothing without them.
  - *Why it matters:* a reader checking whether a reported figure measures prediction or lookup can settle it in one place.

## Before writing the white paper, capture this info in a sentence or two each:

- The problem or unmet need

  Untargeted LC-MS/MS metabolomics detects far more reproducible features than can be identified with authentic chemical standards, and the unidentified remainder carries substantial biological information; communicating how confident an identification is, and what kind of identification it is, remains the central obstacle to using it.

- Why existing approaches are insufficient

  Matching an observed spectrum against a reference library is the established route and the most precise one short of an authentic standard, but what bounds it is coverage: a structure absent from the library is outside the space the search returns, so extending the named fraction has meant manual chemistry at a throughput that limits both library expansion and the discovery of dataset-specific metabolites. Published computational approaches span retrieval to de novo generation, but none has been benchmarked head to head against this workflow, so no accuracy comparison is made.

- What the tool / GDP enhancement does

  It predicts candidate structures for each unknown feature from its MS1/MS2 evidence, reranks them using library, mass and spectral evidence, assigns a pathway and subpathway by consensus over structurally similar reference metabolites, and renders a per-compound report carrying an aggregate evidence score and an evidence-derived MSI-style reporting level.

- What is technically differentiated

  Three contributions, and none of them is accuracy: rank-aware structural validation with uncertainty and source-stratified recovery; provenance-aware taxonomy evaluation that separates direct database lookup from neighbour consensus and from agreement with an independent label source; and evidence-limited reporting for every feature, with a rubric that presents a structure hypothesis only when minimum corroboration criteria are met.

- What evidence supports the claims

  Among 92 truth-verified controls, structural recovery was 47.8% at rank 1 (44/92; 95% CI 37.9–57.9%), 64.1% within the top 3 (59/92), and 73.9% within the top 10 (68/92). Leave-one-out validation showed weak aggregate-score calibration (Brier 0.255), which is why the score is presented as an ordering device rather than a probability. Taxonomy separates by provenance: 1,459 direct lookups are consistent by construction, while 644 neighbour-consensus assignments with independent labels agree 52.2% at super-pathway and 18.8% at subpathway, against a ceiling of 62.5% and 34.8% measured on identical structures. Coverage is complete across 3,897 features, with 2,726 showing a structure and 1,171 withholding one.

- What users can now do that they could not do before

  Rather than selecting a handful of features for manual interpretation, an analyst can rank every feature, place it in a biochemical context, and triage it against the same evidence standard, so the decision of where to spend confirmation effort rests on a comparable basis across the dataset. Three applications follow: driving standard acquisition from a measured shortlist, identifying which hypotheses are corroborated enough to justify extending the reference library, and keeping the unidentified fraction inside the biological interpretation because taxonomy is assigned even where no structure is shown.

- The three most important takeaways

  First, provenance is not by itself confirmation, and the reporting rubric is built around that: an exact catalog source records how a candidate was proposed, and does not stand in for an authentic standard. Second, apparent pathway performance depends entirely on the evaluation context, so any single figure quoted without its provenance is misleading. Third, the contribution is a reporting framework that grades every feature, exposes the provenance of each claim, and withholds a structure when corroboration is insufficient.

## Define the audience and objective

Describe the primary audience

  Metabolomics method developers, computational chemists, and reviewers who assess whether an identification workflow's claims are supported. They will read the Methods before the Results, will check whether a reported figure measures prediction or lookup, and are the reason the glossary exists.

Describe a secondary audience

  Analytical chemists and study leads who will not reimplement the workflow but need to know what an output level asserts before acting on it, and internal readers deciding whether the capability is ready to underpin a commercial deliverable.

After reading this paper, the reader should understand ______, believe ______, and feel confident that ______.

  After reading this paper, the reader should understand **that every unidentified feature can be graded against one evidence standard, with the provenance of each claim exposed**, believe **that the reported figures mean exactly what their cohorts allow, because circular measurements are named as such and the annotation ceiling bounding the taxonomy result is reported alongside it**, and feel confident that **the workflow's outputs can be used to prioritize confirmation effort, with the calibration and cohort limits stated plainly enough to judge where they cannot**.

## Define the job of each section

| Section | What it needs to accomplish |
|---|---|
| Executive summary | *Abstract.* Deliver problem, approach, evidence and consequence in under 200 words, and close on what the workflow is for rather than on what it computed. |
| Introduction | *1. Introduction, with 1.1 Toolkit landscape.* Establish the confidence-communication problem, describe what the established routes cover and where the coverage boundary falls, place the workflow in the published landscape without an accuracy comparison, and end the contributions list on the consequence for an analyst. |
| Background | *1.1 Toolkit landscape (related work).* Describe each published approach factually and symmetrically, so the retrieval and generative ends read as answering different questions rather than as one superseding the other. |
| Tool Description | *2. Methods, 2.1 to 2.7.* Give dataset, workflow, pathway assignment, confidence model, level rubric, evaluation design and reproducibility in enough detail to be assessed, and state explicitly what provenance was not preserved in the report bundle. |
| Validation/evidence | *3. Results, 3.1 to 3.6.* Report the measurement first and interpret it second, cohort by cohort, with Wilson intervals and the evaluation context attached to every figure. |
| Conclusion | *6. Conclusion.* State what the workflow produces across the full set, name the central result as the reporting framework rather than accuracy, and end on the practical consequence. |

The manuscript carries five sections the guidance does not name, four of them because a journal expects them:

| Section | What it needs to accomplish |
|---|---|
| 4. Discussion, with 4.1 Applications | Separate what was observed from what it means; state for each application what an analyst can now do, the result supporting it, and the limitation bounding it. |
| 5. Limitations | Carry the structural limits a journal expects: cohort selection, calibration set size, sparse MS2, circular taxonomy tests, the custom rubric, placeholder co-abundance evidence, and reproducibility constraints. |
| 7. Future work | Name the cohorts, experiments and benchmark that would close each gap. |
| Declarations | Carry authorship, funding, competing interests, ethics and data-availability placeholders through to submission. |
| Glossary | Define, at the back where a reader can look one up, the vocabulary the paper depends on, including the four evaluation senses that separate a validated result from a circular one. |

| User value | (Always keep this in mind throughout) |
|---|---|
| Limitations | (Sprinkle these in to remain trustworthy and acknowledge not all is perfect) |

**User value, kept in mind throughout.** Value is stated as capability rather than as benefit, which is the register a reviewer will accept. The Applications subsection is where the guidance's instruction lands most directly: each of the three entries names what an analyst can do that was previously impractical, the specific result that supports it, and the limitation that bounds it, in that order.

**Limitations, sprinkled in and also collected.** Distributed next to the claims they qualify:

- Overlapping intervals on the exact-versus-predicted subgroups, stated in the same sentence as the two figures.
- The annotation ceiling, reported alongside every independently scored pathway result.
- The leave-one-out pathway result described as validating the neighbour-vote mechanism only, not the production predictor with its external priors and overrides.
- The mass-mismatch alternatives named as a reranking opportunity rather than as recovered structures.

Collected in the Limitations section: availability-selected single-organization cohorts; a confidence model originally fit on the same 92 controls used to evaluate it; observed MS2 present for only 43.6% of features; direct reference lookups being circular as accuracy tests; a custom MSI-style rubric whose boundaries are unvalidated; placeholder co-abundance evidence; and proprietary spectra, checkpoints and reference database that are not publicly deposited.

## Translate features into value

Avoid ending paragraphs on the feature. End on the consequence.

For each section, focus not only on explaining what you built, but also on the problem it addresses, the evidence that it works, and why the improvement matters to the user. The goal is not to remove technical depth, but to organize and translate it so that a broader scientific audience can follow the logic and understand the value.

| Technical capability | Why it matters to target audience |
|---|---|
| Ex. Evidence Scores | Reranking over several independent axes reduces dependence on a single evidence source, so a candidate strong on one axis alone does not reach rank 1 unchallenged |
| Retrieval and de novo generation run as complementary tracks | Candidate coverage is not bounded by catalogue membership, while a library match still contributes its evidence where the compound is catalogued |
| Rank-aware evaluation with Wilson intervals and source stratification | A reader can see how much of the reported performance is provenance rather than prediction, which a single accuracy figure would conceal |
| Provenance-separated taxonomy evaluation, with the annotation ceiling reported beside it | Prevents a lookup that is consistent by construction from being read as prediction accuracy, and tells a reader how much of an apparent shortfall is vocabulary rather than model error |
| MSI-style claim control enforced in the rendered report | The reporting level governs the artifact and not merely the metadata, so a withheld structure is a designed outcome rather than a missing field |
| Committed audit artifacts, a claim ledger and one-command regeneration | Each reported figure traces to a versioned artifact, so the numbers can be reproduced from the repository rather than taken on trust |

Weak:

The workflow calculates an evidence score using spectral and structural information.

Stronger:

The workflow combines spectral and structural evidence into an interpretable score, helping researchers prioritize the most strongly supported annotations for biological interpretation and follow-up.

  The same move, applied to this document. *Weak:* the reporting layer assigns an MSI-style level from the available evidence and suppresses the structure at levels 4 and 5. *Stronger:* the reporting level is enforced in the rendered page rather than only in the rubric, so the 1,171 features that receive no structure are a designed outcome rather than a gap in coverage, and an analyst can trust that a structure shown to them cleared a stated bar.

  For this audience the stronger version still has to hold its bound, so the sentence stops at what the artifact does rather than claiming the withheld features were correctly withheld.

## Other Guidelines

When comparing to other tools or competitors, either have concrete quantitative comparison data or use defensible wording such as:

- provides additional evidence
- increases transparency
- supports more informed prioritization
- reduces dependence on a single evidence source
- enables consistent application at scale
- extends the biological utility of Level 2+ annotations

  No head-to-head benchmark exists, so the manuscript restricts itself to this wording and says so twice, once in the toolkit landscape and once in the Discussion. Combining the retrieval and generative routes *provides additional evidence per candidate* and *reduces dependence on any single evidence source*; the framework *increases transparency* about what each call rests on and *supports more informed prioritization* of which features merit standard acquisition. Against manual interpretation the claim is throughput, so it *enables consistent application at scale*. Because taxonomy is assigned even where no structure is shown, the workflow *extends the biological utility* of features that fall short of confirmation.

  One comparative statement was removed rather than reworded. An earlier draft described the exact-source and predicted-source recovery rates as effectively the same; the intervals overlap on 42 and 50 controls, which supports neither a difference nor an equivalence, so the claim was withdrawn and the package verifier now fails the build if that phrasing returns.

## Separate explanation from proof

Scientists often mix method detail, results, interpretation, and value in one paragraph. Keep them distinct:

- What we did
- What we observed
- What it means
- Why the reader should care

  Results reports the measurement and Discussion interprets it, and the split was enforced during revision rather than assumed. The clearest instance is the confidence result.

  - **What we did.** Recomputed probabilities by leave-one-control-out Platt fitting and evaluated them by Brier score, calibration slope, reliability bins, and an ablation omitting the pathway axis.
  - **What we observed.** Brier 0.255, slope 0.056, intercept −0.082; ranking the controls by leave-one-out confidence and scoring successively larger prefixes leaves rank-1 recovery essentially flat, at 50.0% for the most confident fifth against 47.8% for the whole cohort.
  - **What it means.** The score does not order one feature against another, which is a stronger and more useful statement than poor calibration alone.
  - **Why the reader should care.** Prioritization guidance can be offered for pathway assignments, where similarity does order the calls, and must be withheld for structural ones.

This makes the argument much easier to follow.

## Make an outline before another full draft

Include:

- A section outline with one or two sentences describing the purpose of each section and how to accomplish that purpose

  1. **Abstract.** Problem, approach, evidence, consequence, in under 200 words.
  2. **Introduction.** The confidence-communication problem, what the established routes cover, the toolkit landscape, and three contributions ending on what an analyst can do.
  3. **Methods.** Dataset and its three non-random evaluation subsets; the two-track workflow and the reranker; pathway assignment by lookup or neighbour consensus; the exploratory confidence model; the MSI-style level rubric; the evaluation design; software and reproducibility, including what provenance was not preserved.
  4. **Results.** Coverage; structural recovery with source stratification; confidence tiers and levels with the calibration result; taxonomy separated by evaluation context; the mass-mismatch failure mode; and the delivered artifact.
  5. **Discussion.** The two reframing results, how the two confidence signals differ in what they can order, and the positioning against published tools without an accuracy claim; then Applications, each with its supporting result and bounding limitation.
  6. **Limitations.** Seven structural limits, each tied to the analysis it bounds.
  7. **Conclusion.** What the workflow produces, what the central result is, and the practical consequence.
  8. **Future work.** External validation cohorts, a prospective hard mass gate, better MS2 coverage and scoring, real co-abundance evidence, and a formal CASMI-style benchmark.
  9. **Glossary.** A consultable reference for the vocabulary the paper depends on, closing the document.

- Proposed figures and the message each figure supports

  Six main figures and one main table. Supplementary figures S1 to S5 and tables S1 to S9 carry the extended analyses.

For each figure, complete this sentence:

This figure demonstrates that ______, which matters because ______.

  - `fig1_workflow.svg` — This figure demonstrates that *evidence generation, ranking, taxonomy assignment and claim control are separable stages with distinct failure modes*, which matters because *each result in the paper can be attributed to a specific stage rather than to the pipeline as a whole*.
  - `fig2_structure_accuracy.svg` — This figure demonstrates that *recovery improves substantially from rank 1 to top 10, and that the two source subgroups have overlapping intervals on a cohort of this size*, which matters because *it justifies reporting ranked hypotheses and refuses the inference that an exact-source call is a confirmed one*.
  - `fig3_confidence_reliability.svg` — This figure demonstrates that *the aggregate score is poorly calibrated under leave-one-out validation*, which matters because *it is the evidence for presenting the score as an ordering device rather than a probability of correctness*.
  - `fig4_identification_distributions.svg` — This figure demonstrates that *the reporting rubric concentrates results in hypothesis-level categories and leaves the strongest tier empty*, which matters because *it shows the claim-control rules bind in practice rather than only in principle*.
  - `fig5_pathway_evaluations.svg` — This figure demonstrates that *apparent pathway performance depends entirely on which evaluation context is used*, which matters because *quoting a single pathway accuracy number without its provenance would misstate the result by tens of percentage points*.
  - `fig6_report_example.png` — This figure demonstrates that *the reporting level is enforced in the rendered page, with the public and internal renderings of one feature differing in how much of the candidate list each may expose*, which matters because *claim control was otherwise only ever described, and this shows it binding on the artifact a customer actually receives*.
  - **Table 1** — The one display item that cannot be a figure, because the point is the four-way comparison of pathway evaluation contexts and their intervals side by side.

## When reviewing what’s written:

## Review 1: Argument

- Is the problem clear?

  Yes. The Introduction states it in one sentence and the abstract in its first, and the problem is confidence communication rather than detection, which is the distinction the rest of the paper turns on.

- Is the solution differentiated?

  Yes, but deliberately not on accuracy. The three contributions are evaluation design, provenance separation and claim control. The manuscript says outright that the contribution is not higher structural recovery, which is a stronger position than an unsupported comparison would be.

- Is the evidence sufficient?

  For the claims made, yes. Where it is not sufficient, the claim was withdrawn rather than hedged: the authoritative P0 and P1 cohort analyses are absent because their source files were unavailable, and the 92-control set was explicitly not substituted for them.

- Is the value explicit?

  Yes, in the Applications subsection, where each entry states the new capability, the result supporting it, and the bound. Elsewhere value is stated as capability rather than benefit, which suits the audience.

## Review 2: Audience

- Can a broad scientific reader understand it?

  Within the intended audience, yes. The paper assumes familiarity with LC-MS/MS identification and with the identification-level frameworks, and it does not assume familiarity with the specific model families, which are described in one place and then referred to by role.

- Are terms defined?

  Yes, in a Glossary that closes the document: 53 entries in seven groups, covering the measurement, where a candidate comes from, comparing structures, the evaluation vocabulary, what a result claims, biochemical taxonomy and the failure modes. It fixes the four evaluation senses a reader could otherwise conflate, separates identification from annotation and from structure hypothesis, keeps the two kinds of confidence tier apart, and defines the provenance vocabulary neutrally. Entries carry their numbers where a term means nothing without them, so catalog provenance is defined alongside the recovery rate measured for it. The glossary moved out of Methods because a reader looking up a word does not search the middle of a Methods section for it; a pointer in 2.6 names the four evaluation senses at the point where the distinction between them starts to bite.

- Is unnecessary detail removed?

  Mostly. Extended tables, threshold sensitivity, harmonization scenarios, selective prediction and library stratification all moved to Supplementary Information, leaving the main text to carry one figure per argument. What remains that a reader might call unnecessary is the record of provenance that was *not* preserved in the report bundle, which is kept because a reviewer needs it to judge reproducibility.

## Review 3: Scientific rigor

- Are claims supported?

  Yes, and the package verifier asserts it mechanically: the audited counts must appear verbatim in the prose, the pathway coverage counts and the best-reference similarity median are checked against the audit artifact, and the headline statistics must agree across the manuscript, the paper README and the operations overview. That check exists because those counts had drifted apart once before.

- Are limitations stated?

  Yes, both distributed and collected. The most consequential are that the confidence model was originally fit on the same controls used to evaluate it, and that direct reference lookups cannot serve as accuracy tests.

- Are accuracy, confidence, identification, and annotation used precisely?

  This is the review that changed the paper most. *Accuracy* is now used only where a prediction is scored against a label withheld from the mechanism that produced it, which excludes the direct lookups an earlier version reported as pathway accuracy. *Recovery*, *agreement* and *consistency* carry the other three senses. *Identification* is a graded claim about identity, *annotation* is any label including taxonomy, and a *structure hypothesis* is explicitly not asserted as confirmed. Identification-confidence tiers and pathway-confidence tiers are never used interchangeably. All of these are now defined in the Glossary rather than left to context.

## Review 4: Style

- Sentence length

  Long where a claim needs its qualification attached, short where the claim is the point. The bolded lead-ins in Discussion and Results carry the argument for a reader skimming.

- Repetition

  The provenance argument is made once in Results, once in Discussion and once in the Glossary, which is intentional, since each is a different move: measurement, interpretation, definition.

- Transitions

  Each Results subsection reports before it interprets, and the Discussion opens by naming the two results that reframe the rest, so a reader knows what the section is for.

- Figures and captions

  All six captions lead with the message the figure supports and follow with the methodological detail. Four superseded figures were deleted rather than left unreferenced, after one of them was found still rendering a retracted headline that nothing was checking.

- Final polish

  Methods numbering was corrected twice, once when the glossary was inserted into that section and again when it was moved out to the back. The glossary is written as bullets rather than a numbered list, because the citation-parity check reads any line beginning with a digit and a period as a reference entry.

## Review 5: Fairness to the incumbent method

Added after a reader observed that the manuscript was dismissive of library matching. Not part of the guidance as issued.

- The finding

  The Introduction opened on what library-retrieval-only workflows "cannot" do and called the route that follows a "bottleneck"; the toolkit landscape described one tool as "bounded by" its library while the generative tools were "not restricted to database membership" and "independent of spectral libraries". Read together, that framed the method which produces the named fraction of every deliverable, and which this workflow uses as an evidence axis, as a deficient one.

- The correction

  The boundary is now coverage rather than capability, the two ends of the toolkit spectrum are described as answering different questions, and the Discussion states that a library match remains the most direct route where the compound is catalogued and that well-corroborated hypotheses are candidates for entering the library.

- What was deliberately not softened

  Exact-source rank-1 recovery of 38.1% (16/42) and predicted-source 56.0% (28/50) both stand, with their overlapping intervals, and so does the conclusion that catalog provenance is not confirmation. Two sentences were rescoped rather than removed: Results no longer asserts that provenance confers no recovery advantage, which claims more than 42 controls can support, and now says the cohort does not resolve a difference in either direction. That is a narrower claim, not a friendlier one.

---

## Carried, for a human

Not part of the guidance. These need a decision and are not closable by review or by the gates.

- **Submission placeholders.** Eight `[REQUIRED: ...]` items remain, covering the author list, corresponding-author email, CRediT roles, funding, competing interests, ethics and consent determination, run-environment provenance, and the repository release or DOI. Run the verifier with `--strict-human-inputs` to make them blocking once filled.
- **Display items.** Six figures and one table is a self-imposed budget rather than a journal rule. If the corresponding author would rather not spend a display item on the report capture, it can move to Supplementary Information at no cost to the argument.
- **Where the glossary sits in the Word build.** It closes the document, which puts it after the declarations, since Nature order runs References, Methods, then the declarations. That is unusual back matter and is a placement decision rather than an accident: `Glossary` is the last entry in `SECTION_ORDER`. Moving it earlier, or into Supplementary Information, is a one-line change.
- **Reporting rules ahead of the artifact.** The reporting rules changed after the delivered set was produced. The manuscript and Supplementary Information describe the rules that generated the committed distributions, which keeps each internally consistent, but both must be updated together the next time the reports are regenerated. `docs/paper/REVISION_LOG.md` lists what will move.
- **Authoritative cohorts.** The P0 and P1 analyses remain absent until their source files are supplied, and the 92-control set must not be substituted for them.
