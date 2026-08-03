# Structural Elucidation for Unknown Metabolites

**Naming the unnamed: evidence-graded structure hypotheses for untargeted LC-MS/MS features**

**Metabolon, Inc.** — [REQUIRED: contact email] — [REQUIRED: public web URL]

---

## Executive summary

Untargeted metabolomics detects far more molecules than any laboratory can name. The unnamed remainder is not noise: it is reproducible, quantifiable, and often the part of the data most closely tied to the biology under study. Metabolon's structural elucidation capability turns those unnamed features into documented structure hypotheses with an explicit, graded statement of how much evidence stands behind each one.

Applied to a delivered set of 3,897 unnamed features, the capability produced a complete evidence report for every feature, proposed a specific structure for 2,726 of them, and deliberately withheld a structure for the remaining 1,171 where the evidence did not justify one. That last behaviour is the point of the system as much as the first.

## The unknown-metabolite problem

A typical untargeted run yields thousands of reproducible signals. Only a fraction correspond to compounds that have been matched to an authentic chemical standard. The rest are reported as unnamed features, and most analyses set them aside.

Setting them aside is expensive in two directions. Biologically, discarded features frequently include the compounds that differentiate the study groups. Operationally, the conventional route to naming one of them, in which a chemist manually interprets the spectra and then buys or synthesizes a candidate to confirm it, is slow enough that it can only ever be applied to a handful of features per study.

The alternative is not to guess. It is to produce a defensible hypothesis and to say plainly how strong that hypothesis is, so that a scientist can decide which features are worth the cost of confirmation.

Scale is part of the answer. Because the capability is computational, it is applied uniformly to every unnamed feature in a study rather than to the few that a chemist has time to examine. The same criteria are applied to the first feature and the three-thousandth, which makes the resulting confidence grades comparable across a whole dataset and across studies.

## What the capability delivers

Every unnamed feature receives a report. There are no silent omissions and no features that simply disappear from the deliverable.

- **A structure hypothesis, where the evidence supports one.** Roughly seven in ten features receive a specific proposed structure.
- **Biochemical pathway context for every feature**, including those where no structure is shown, so that an unnamed feature can still be placed in a biological narrative.
- **A graded confidence statement** aligned to the community-standard identification levels used across the metabolomics field.
- **The underlying evidence**, presented so that a reviewer can see what the call was based on rather than having to trust a single score.

![Figure 1. How the capability works.](figures/fig_workflow_public.svg)

*Figure 1. From detected signal to a graded, documented result. Each stage adds evidence, and the final stage decides how much may responsibly be claimed.*

## How it works

The capability runs in four conceptual stages.

1. **Detect.** Reproducible features and their fragmentation evidence are extracted from the untargeted run.
2. **Propose.** For each feature, a set of candidate structures is assembled that is not limited to compounds already present in a reference library. This is what allows genuinely novel and study-specific molecules to be described rather than skipped.
3. **Weigh.** Candidates are ordered by how well several independent lines of evidence agree with one another, rather than by any single measure.
4. **Report.** Each result is graded, given pathway context, and rendered as a standalone report. Where the supporting evidence falls short, the structure is withheld and the report says so.

The value is concentrated in the third and fourth stages. Proposing plausible molecules is comparatively easy; deciding which proposals have earned the right to be shown to a customer is the hard part, and it is where the system is deliberately conservative.

## Evidence at a glance

The capability has been evaluated against a set of internal controls whose true structures are known independently. On those controls, the correct structure was the single top-ranked proposal 47.8% of the time, appeared within the top three proposals 64.1% of the time, and appeared within the top ten 73.9% of the time.

![Figure 2. Structural recovery on truth-verified controls.](figures/fig_recovery_public.svg)

*Figure 2. How often the correct structure appears at the top of the candidate list, and within the top three and top ten, across 92 controls with independently known structures.*

Two points of context matter when reading those figures. The control set is modest in size at 92 compounds and was assembled from the material available rather than sampled at random, so the percentages describe performance on that set rather than a guaranteed rate on new data. And the results are ranked lists: a top-ten result is genuinely useful, because it reduces a boundless chemical space to a short list a chemist can work through, but it is not the same as a single confirmed answer.

Pathway assignment was evaluated separately. A leave-one-out test of the pathway assignment mechanism across 5,095 reference metabolites recovered the correct super-pathway 86.7% of the time. Assignment at the finer subpathway level is substantially harder, and agreement between independent biochemical taxonomies is itself imperfect, so pathway context is best used to orient and prioritize rather than as a settled classification.

## What you receive

Each feature is delivered as a self-contained report, with a searchable index across the full set.

- **The proposed structure and a short list of alternatives**, so that a reviewer sees the runners-up rather than only the winner.
- **The supporting evidence**, including how well the candidate accounts for the observed mass and fragmentation behaviour.
- **Biochemical pathway and subpathway context.**
- **An identification level** and a plain-language statement of what that level does and does not assert.
- **Cross-references** to the independently assigned Metabolon Level, so the new output can be read alongside existing deliverables.

![Figure 3. Reporting coverage across the delivered set.](figures/fig_coverage_public.svg)

*Figure 3. Across 3,897 delivered features, every one receives a report and pathway context, while a specific structure is shown only where the evidence supports it.*

## How to read a result

Results are graded on the confidence-level framework the metabolomics community already uses, so the output is interpretable without learning a bespoke vocabulary.

**Table 1. What each reporting outcome asserts.**

| Reporting outcome | What it means | Reasonable next step |
|---|---|---|
| Confirmed against an authentic standard | The identity is established | Use directly |
| Probable structure | Multiple independent lines of evidence agree | Strong candidate for confirmation |
| Tentative structure | A specific, evidence-consistent hypothesis | Prioritize if biologically interesting |
| Formula-level only | Composition supported; no structure claimed | Treat as an unnamed but characterized feature |
| Mass-level only | The feature is real and quantifiable | Track across studies |

The most important line in that table is the fourth. When the system cannot support a structure, it does not offer one. Withholding is a designed behaviour, not a gap in coverage, and it is what makes the results that are shown worth acting on.

For the same reason, the confidence output is used to rank and to limit claims, not as a statement of the probability that a structure is correct. Reporting a hypothesis as a hypothesis is a more useful and more honest product than reporting it as a number that implies more certainty than the underlying evidence can support.

## Where this fits

This capability complements, rather than replaces, standard-confirmed identification. It is designed for three situations.

- **Prioritizing discovery.** When a study points to unnamed features, the capability turns them into a ranked, chemically specific shortlist so that confirmation effort goes to the most promising targets.
- **Extending reference coverage.** Well-supported hypotheses provide a defensible starting point for adding compounds to a reference library.
- **Interpreting the unnamed fraction.** Pathway context for every feature, including those without a structure, keeps the unnamed portion of the data inside the biological interpretation.

It is deliberately not a replacement for authentic-standard confirmation, and it is not presented as one. A proposed structure is a hypothesis to be tested, and the reporting is built so that the distinction between a tested identity and a proposed one is never blurred. Customers who need a settled identity should continue to treat confirmation as the final step; the capability makes that step cheaper by pointing it at the right molecules.

## Availability

[REQUIRED: describe commercial availability, engagement model, and how a prospective customer initiates a project.]

[REQUIRED: confirm trademark usage, legal disclaimer, and marketing approval before external release.]

Detailed methodology, the full validation analysis, and the accompanying limitations are documented in a separate scientific manuscript available under appropriate agreement. [REQUIRED: state the disclosure process and any confidentiality terms.]

---

*Prepared by Metabolon, Inc. [REQUIRED: document version, date, and approval reference.]*
