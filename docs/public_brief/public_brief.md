# Structural Elucidation for Unknown Metabolites

**Naming the unnamed: evidence-graded structure hypotheses for untargeted LC-MS/MS features**

**Metabolon, Inc.** — [REQUIRED: contact email] — [REQUIRED: public web URL]

---

## Executive summary

Untargeted metabolomics detects far more molecules than any laboratory can name. The unnamed remainder is not noise: it is reproducible, quantifiable, and often the part of the data most closely tied to the biology under study. But because those features arrive with no structure and no statement of confidence, most analyses set them aside, and the biology they carry is set aside with them.

Metabolon's structural elucidation capability closes that gap. Compounds confirmable against an authentic standard are still reported that way, because that remains the strongest answer; this begins where that stops. It turns unnamed features into documented structure hypotheses, each carrying an explicit, graded statement of how much evidence stands behind it. Applied to a delivered set of 3,897 unnamed features, it produced a complete evidence report for every feature, proposed a specific structure for 2,726 of them, and deliberately withheld a structure for the remaining 1,171 where the evidence did not justify one.

Three consequences follow for a research team. Nothing is silently dropped from the deliverable, so the unnamed fraction stays inside the biological interpretation and can extend a biochemical signature into areas the reference catalogue does not yet cover. What the workflow declines to claim is as informative as what it proposes, which is what makes the results that are shown worth acting on. And because every feature is graded against the same criteria, results stay comparable across studies as well as within one, so confirmation effort can be aimed at the molecules most likely to repay it.

## The unnamed fraction of an untargeted dataset

A typical untargeted run yields thousands of reproducible signals. Only a fraction correspond to compounds that have been matched to an authentic chemical standard. The rest are reported as unnamed features: real, measurable, and anonymous.

Setting them aside is expensive in two directions. Biologically, the discarded features frequently include the compounds that differentiate the study groups, so the finding that would have mattered most can be the one never examined. Operationally, naming one by hand and then buying or synthesizing a candidate to confirm it is slow enough to reach a handful of features per study. The deeper cost is that the choice of which few to pursue is made without any consistent basis for comparison.

The alternative is not to guess. It is to produce a defensible hypothesis for every feature and to say plainly how strong each one is, so that a scientist can decide which are worth the cost of confirmation.

## What the established routes cover, and what they leave

Three routes put a name to a signal today. Two of them work and are bounded by what they cover. The third trades rigor for volume, and it is the one worth understanding before reading any list of metabolite names.

The first is matching an observed spectrum against a reference library. It is fast, well understood, and the most precise answer available short of a standard, which is why it produces the named part of every untargeted deliverable and why it is the first step here too. What bounds it is the size of the catalogue rather than the quality of the match: an uncatalogued molecule is outside the space such a search covers, and novel, study-specific chemistry tends to sit outside it.

The second is expert manual interpretation, followed by acquisition of a candidate compound and confirmation against it. This remains the most authoritative answer available and the proper standard for a settled identity. It is also the scarcest resource in the workflow, which is why it reaches a few features rather than thousands.

The third is to match every detected signal against a broad database and report whatever comes back. It produces an impressive number, and that is the problem. A list of several thousand metabolite names is not several thousand findings when the same molecule is counted more than once, when different ionized and isotopic forms of one compound each return as their own entry, and when a weak match is printed in the same typeface as a strong one. Nothing in the list says which is which, so the work of deciding what to trust passes to the reader, who has less to go on than the workflow did. Volume of that kind is not coverage.

Where such a list does carry a confidence number, the question to ask is what the number was tested against, and the answer is usually not given. A score with no stated test, and no statement of which comparisons it does and does not support, cannot be acted on any differently from a score that was never tested at all. That is the standard we hold ourselves to below, and it is a fair one to apply to anybody's output, including ours.

What none of the three supplies on its own is the thing a study needs in order to triage: a consistent, comparable statement of how much evidence supports each candidate identity. The metabolomics community communicates that through identification levels, which describe the kind of evidence behind a call, from an exact mass alone up to a structure confirmed against an authentic standard. Applying that vocabulary to every unnamed feature in a study, rather than to the few a chemist has time to examine, is the requirement this capability was built to meet.

## What the capability does

Every unnamed feature receives a report. There are no silent omissions and no features that quietly disappear from the deliverable.

![Figure 1. How the capability works.](figures/fig_workflow_public.svg)

*Figure 1. Each stage contributes a distinct kind of evidence, and the final stage decides how much may responsibly be claimed. The confidence grade is produced by the workflow rather than attached to the result afterwards.*

The capability runs in four conceptual stages.

1. **Detect.** Reproducible features and their fragmentation evidence are extracted from the untargeted run, establishing that a signal is real and measurable before anything is claimed about its identity.
2. **Propose.** For each feature, a set of candidate structures is assembled that reaches beyond the compounds a reference library already contains, which is what allows genuinely novel and study-specific molecules to be described rather than skipped.
3. **Weigh.** Candidates are ordered by how well several independent lines of evidence agree with one another rather than by any single measure, so a candidate that looks convincing on one axis alone does not reach the top of the list unchallenged.
4. **Report.** Each result is graded, given biochemical pathway context, and rendered as a standalone report. Where the supporting evidence falls short, the structure is withheld and the report says so.

The value is concentrated in the last two stages. Proposing plausible molecules is comparatively easy; deciding which proposals have earned the right to be shown is the hard part, and it is where the capability is deliberately conservative.

**How a wrong or repeated call is kept out.** Two things have to be prevented: the same molecule being counted several times, and a structure that the evidence does not support being shown at all. The controls for both sit inside the stages above rather than being applied at the end.

Against repetition: what enters is a reproducible feature with its fragmentation evidence combined across the scans it appeared in, not a raw signal; candidates are scored against the mass of the uncharged molecule implied by how it ionized, so different ionized forms of one compound are not read as different compounds; isotope agreement is among the evidence weighed; and a feature seen by more than one chromatographic method is counted once.

Against unsupported structures: a top candidate whose mass does not fit is demoted rather than shown, which happened for 598 features, and where the evidence does not support a structure, none is offered — 1,171 of 3,897 features arrive with no structure at all. That last number is the clearest evidence of the approach, because it is coverage we could have claimed and did not.

None of this makes the workflow incapable of being wrong, and it is worth saying so plainly. What it does is make a wrong call visible as a weak one: an unsupported proposal is either withheld or arrives carrying a grade that tells you not to build on it yet. Containing the consequences of an error is a different thing from never making one, and only the first is claimed here.

**Table 1. What each part of the workflow means for a research team.**

| Capability | Why it matters to you |
|---|---|
| A report for every unnamed feature | The unnamed fraction stays inside the biological interpretation instead of being written off |
| Candidates that reach beyond the reference library | Novel and study-specific molecules can be described, which is often where the differentiating biology sits |
| Several independent lines of evidence combined into one ordering | Reduces dependence on a single evidence source, so one flattering measure cannot carry a weak candidate |
| Grading aligned to community identification levels | Increases transparency and supports more informed prioritization, with no bespoke vocabulary to learn |
| A structure withheld when the evidence is insufficient | Confirmation effort is not spent chasing a structure the data never supported |
| Pathway context on every feature | Extends the biological usefulness of features that fall short of confirmation |
| The same criteria applied to every feature | Enables consistent application at scale, so grades stay comparable across a study and between studies |

## How biochemical context is assigned

Every feature is placed in a biochemical pathway and subpathway, including the features where no structure is shown. That context arrives by one of two paths, and they do not carry equal weight.

![Figure 2. The two paths to a biochemical assignment.](figures/fig_pathway_public.svg)

*Figure 2. Biochemical context reaches a feature either by lookup against a catalogued compound or by inference from chemically similar ones. The two are not equally strong, which is why every report says which of them produced the call.*

**When the proposed structure is already catalogued**, the classification established for that compound is carried across unchanged. It is a lookup rather than a prediction, so it is exactly as dependable as the reference catalogue and no more. **When the structure is new to the catalogue**, the classification comes from the consensus of the catalogued compounds most chemically similar to it, because molecules that are chemically alike usually take part in similar biochemistry.

That second path was tested by hiding a reference metabolite's own labels and asking the capability to recover them, across several thousand metabolites. The broad class comes back right most of the time; the finer subpathway is less reliable, because fine distinctions turn on exactly the small differences that chemical similarity smooths over. So the broad class is worth acting on, the subpathway is a lead worth checking against your own conventions, a looked-up assignment is firmer than an inferred one, and every report tells you which of the two you are reading.

## Evidence that it works

**What was done.** The capability was scored against internal controls whose true structures are known independently, so each proposal was judged against a known answer rather than against another prediction. A proposal counts as correct when it matches that answer within a fixed structural-similarity threshold.

**What was observed.** On those controls, the correct structure was the single top-ranked proposal 47.8% of the time, appeared within the top three proposals 64.1% of the time, and appeared within the top ten 73.9% of the time.

![Figure 3. Structural recovery on truth-verified controls.](figures/fig_recovery_public.svg)

*Figure 3. The correct structure is usually near the top of a short ranked list, across 92 controls with independently known structures. That is what converts an unbounded chemical space into a shortlist a chemist can work through.*

**What it means.** Two caveats travel with those percentages. The control set is modest at 92 compounds and was assembled from the material available rather than sampled at random, so it describes performance on that set rather than a rate to assume for new data. And these are ranked lists: a top-ten result reduces a boundless chemical space to a short one, which is useful but is not a settled answer.

**Why it matters to you.** A shortlist that usually contains the right molecule changes the arithmetic of following a feature up: instead of committing chemistry time on intuition, a team works through a small ranked set whose hit rate has been measured.

## What you receive

Each feature is delivered as a self-contained report, with a searchable index across the full set, so the deliverable can be browsed by a biologist and audited by a chemist.

- **The proposed structure and a short list of alternatives**, so that a reviewer sees the runners-up rather than only the winner.
- **The supporting evidence**, including how well the candidate accounts for the observed mass and fragmentation, so a call can be checked rather than taken on trust.
- **Biochemical pathway and subpathway context**, including where no structure is shown.
- **An identification level** with a plain-language statement of what it does and does not assert.
- **Cross-references** to the independently assigned Metabolon Level, so the new output can be read alongside existing deliverables.

![Figure 4. A delivered report.](figures/fig_report_public.png)

*Figure 4. Every feature arrives as a page like this one, stating the proposed structure and the grade limiting it together, so what a result claims is visible rather than buried in a score. The panel is the top of one report, and the example is a demonstration compound.*

![Figure 5. Reporting coverage across the delivered set.](figures/fig_coverage_public.svg)

*Figure 5. Across 3,897 delivered features, every one receives a report and pathway context while a specific structure appears only where the evidence supports it, so completeness and restraint arrive together.*

## How to read a result

Results are graded on the confidence-level framework the metabolomics community already uses, so the output is interpretable without learning a bespoke vocabulary. The words that carry the most weight in a report are defined in Key terms, at the back of this document.

**Table 2. What each reporting outcome asserts.**

| Reporting outcome | What it means | Reasonable next step |
|---|---|---|
| Confirmed against an authentic standard | The identity is established | Use directly |
| Probable structure | Multiple independent lines of evidence agree | Strong candidate for confirmation |
| Tentative structure | A specific, evidence-consistent hypothesis | Prioritize if biologically interesting |
| Formula-level only | Composition supported; no structure claimed | Treat as an unnamed but characterized feature |
| Mass-level only | The feature is real and quantifiable | Track across studies |

The most important lines are the last two. When the capability cannot support a structure it does not offer one, and that withholding is designed behaviour rather than a gap in coverage. It is the reason a structure that is shown carries weight.

The confidence output works the same way: it limits what a report claims rather than giving the probability that a structure is correct, and we say which of those it is. Within a feature the candidate list is genuinely ordered, so working down a short list pays. Between features, use the biochemical confidence described earlier, which does separate stronger assignments from weaker ones, rather than the structural grade. We tested both, and publishing what each supports is what makes a grade usable rather than decorative.

## Where this fits

This capability complements library matching and authentic-standard confirmation rather than replacing either. Matching against the reference library stays the first step and the most direct answer wherever a feature corresponds to a catalogued compound; what this adds is a graded result for the features that fall outside that coverage. The relationship runs both ways, which is the part most easily missed: every proposal is weighed against catalogued chemistry and placed in biochemistry by the catalogued compounds it most resembles, so a larger confirmed library does not only name more compounds outright, it makes the predictions about everything else better. The capability is designed for two situations in particular.

- **Prioritizing discovery.** When a study points to unnamed features, the capability turns each of them into a short, chemically specific candidate list, so the work goes into a handful of named possibilities rather than starting from nothing.
- **Extending reference coverage.** Well-supported hypotheses give a defensible starting point for adding compounds to a reference library, and every addition raises the share of future studies that can be named outright, so the two routes reinforce each other rather than compete.

Used alongside a reference-library search, the capability contributes additional evidence and reduces dependence on a single evidence source. Used alongside manual interpretation, it enables consistent application at scale. It has not been benchmarked head to head against other computational approaches, and no such comparison is claimed here.

[REQUIRED: add matrix-specific impact and a coverage-improvement table. Both need source-study matrix metadata and a per-study baseline for the named fraction, neither retained in the reporting bundle.]

## Conclusion

The unnamed fraction of an untargeted dataset is not a reporting inconvenience. It is unexamined biology, and this capability makes it examinable.

Three points carry the argument. Every unnamed feature receives a documented, graded result. What the workflow declines to claim matters as much as what it proposes, because restraint is what gives the remaining calls their credibility. And the output is an instrument for prioritization rather than a substitute for confirmation: it does not remove the final experiment, it tells you which one to run.

The part of a dataset previously written off becomes a pathway-placed, evidence-graded set of hypotheses that can be scrutinized, followed up, and published.

## Key terms

**Table 3. Key terms.**

| Term | What it means here |
|---|---|
| Annotation | Any label attached to a feature, including its area of biochemistry, and not by itself a claim about which compound it is |
| Authentic standard | A pure sample of a compound run the same way as the study samples, comparison against which settles an identity outright |
| Candidate | One proposed structure for a feature. Each feature receives several, ordered by the evidence behind them |
| Catalogued compound | A compound already characterized and present in the reference library, so a matching signal can be named by lookup |
| Confidence grade | The report's statement of how much evidence stands behind a call. It limits what the report claims rather than giving the probability that a structure is correct |
| Confirmation | The follow-up experiment that turns a proposed structure into a settled identity, normally by obtaining the compound and comparing against it |
| Deliverable | One self-contained report per feature, plus a searchable index across the full set |
| Evidence | The measurements bearing on whether a candidate is right, including how well it accounts for the observed mass and the fragmentation behaviour |
| Feature | One reproducible signal in the data, quantifiable across samples whether or not it has a name |
| Fragmentation | The pattern of pieces a molecule breaks into during measurement, carrying structural information beyond the mass alone |
| Identification | A claim about which compound a feature is, carrying a grade that says what kind of evidence stands behind it |
| Identification level | The community-standard grade for that evidence, from an exact mass alone up to a structure confirmed against an authentic standard |
| Library match | Naming a feature by finding it in the reference library. It is the most direct answer wherever the compound is catalogued, and it produces the named part of every deliverable |
| Metabolite | A small molecule produced or consumed by the biochemistry of a living system |
| Molecular formula | The elemental composition of a compound, which is what a report presents when the evidence supports composition but not shape |
| Observed mass | The mass measured for a feature, which constrains the structures chemically possible for it |
| Pathway and subpathway | The area of biochemistry a feature belongs to, from a broad class down to a narrower grouping. Every feature receives both, including those with no structure |
| Ranked shortlist | The ordered candidates for one feature. Working down the list is the intended use, because the correct structure is often present without being first |
| Reference library | The catalogue of already-characterized compounds an observed signal is matched against. Its size bounds how much of a dataset can be named by matching |
| Structure hypothesis | A specific proposed structure, offered as a lead rather than as a settled identity |
| Triage | Deciding which features are worth following up. Grading every feature against the same criteria is what makes that decision comparable across a study |
| Unnamed feature | A feature that has not been matched to a known compound, and so arrives with no structure attached |
| Untargeted metabolomics | Measurement designed to detect every molecule a method can, rather than a predefined list of them |
| Withheld structure | A deliberate decision not to propose a structure the evidence does not support. It is what gives the structures that are shown their weight |

## Availability

[REQUIRED: describe commercial availability, engagement model, and how a prospective customer initiates a project.]

[REQUIRED: confirm how this capability is positioned relative to existing Metabolon platform deliverables, and whether those deliverables may be named here.]

[REQUIRED: confirm trademark usage, legal disclaimer, and marketing approval before external release.]

Detailed methodology, the full validation analysis, and the accompanying limitations are documented in a separate scientific manuscript available under appropriate agreement. [REQUIRED: state the disclosure process and any confidentiality terms.]

---

*Prepared by Metabolon, Inc. [REQUIRED: document version, date, and approval reference.]*
