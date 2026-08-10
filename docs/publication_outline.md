# Publication outline

The white-paper writing guidance asks for an outline before another full draft. This
file is that outline, and it covers both publication surfaces at once so the two stay
consistent with each other:

- **White paper** — `docs/public_brief/` (customer-facing, methodology withheld).
- **Scientific manuscript** — `docs/paper/` (Nature initial-submission format).

Both draw their statistics from the same audited artifacts under `docs/paper/data/`
and the claim ledger at `docs/paper/supplement/claim_evidence.csv`.

## 1. Pre-writing capture

**The problem or unmet need.** A single untargeted LC-MS/MS study yields thousands of
reproducible signals, but only a fraction can be matched to an authentic chemical
standard, so the rest are delivered as unnamed features and set aside. That discarded
remainder is reproducible and quantifiable, and it frequently contains the compounds
that differentiate the study groups.

**Why existing approaches are insufficient.** Matching against a spectral library can
only return compounds the library already contains, so genuinely novel and
study-specific molecules are structurally invisible to it. Manual expert
interpretation followed by purchase or synthesis is authoritative but scales to a
handful of features per study, and neither route produces a comparable, explicit
statement of how much evidence stands behind each call.

**What the capability does.** For every unnamed feature it assembles candidate
structures that are not restricted to the reference library, orders them by how well
several independent lines of evidence agree, assigns biochemical pathway context, and
issues a graded evidence report — proposing a specific structure only where the
evidence supports one.

**What is technically differentiated.** Three things, in the order they matter:
candidate generation is not bounded by library membership; ranking draws on several
independent evidence axes rather than one similarity measure; and the reporting layer
is an explicit claim-control mechanism that withholds a structure when corroboration
is insufficient. The grade is therefore produced by the process rather than asserted
after it.

**What evidence supports the claims.** On 92 controls whose structures are known
independently, the correct structure was the single top-ranked proposal 47.8% of the
time, appeared within the top three 64.1% of the time, and within the top ten 73.9%
of the time. A leave-one-out test of the pathway assignment mechanism across 5,095
reference metabolites recovered the correct super-pathway 86.7% of the time. Across
3,897 delivered features, every feature received a report and pathway context, a
structure was proposed for 2,726, and a structure was withheld for 1,171.

**What users can now do that they could not do before.** They can act on the unnamed
fraction of their data instead of setting it aside: rank it, place it in a biochemical
pathway, and take a short, chemically specific shortlist into confirmation rather than
choosing a few features by intuition. Because the capability is computational, it is
applied uniformly to every feature, so the grades are comparable across a dataset and
across studies.

**The three most important takeaways.**

1. Every unnamed feature receives a documented, graded result. Nothing is silently
   dropped from the deliverable.
2. What the workflow declines to claim matters as much as what it proposes.
   Withholding is designed behaviour, and it is what makes the results that are shown
   worth acting on.
3. This is a prioritization instrument, not a replacement for authentic-standard
   confirmation. It makes confirmation cheaper by pointing it at the right molecules.

## 2. Audience and objective

**Primary audience.** Study leads and data users — in biopharma, biotech, and academic
research — who receive an untargeted metabolomics deliverable and have to decide what
to do about the unnamed features in it. They are scientifically literate but are not
mass-spectrometry method developers, and their decision is where to spend confirmation
effort.

**Secondary audience.** Analytical chemists, metabolomics core directors, and
computational biologists who will scrutinize how the confidence grades were derived
and how far they can be trusted before recommending the capability internally.

**Objective.** After reading this paper, the reader should understand *that the
unnamed fraction of an untargeted dataset can be turned into ranked, pathway-placed,
evidence-graded structure hypotheses*, believe *that the grades are earned rather than
asserted, because the workflow withholds a structure whenever the evidence does not
support one*, and feel confident that *they can use the output to decide which
features are worth the cost of confirmation, with the limits of each claim stated
plainly*.

## 3. Section jobs — white paper

| Section | What it needs to accomplish |
|---|---|
| Executive summary | Deliver the whole problem-solution-evidence-impact arc in miniature and close on the three takeaways, so a reader who stops here still leaves with the argument. |
| Introduction | Establish the unnamed-feature problem and price it in both directions: biologically, discarded features often carry the differentiating signal; operationally, the conventional route to naming one is too slow to apply at scale. |
| Background | Explain how identification confidence is communicated in the field today, and what the two existing routes — library matching and manual expert interpretation — each cover, so the gap is located in catalogue coverage and expert throughput rather than in the quality of either method. Library matching produces the named fraction of every deliverable and is the first step of this workflow, so it is described as the route being extended rather than one being displaced. Ends by naming the requirement the rest of the paper answers. |
| What the capability does | Describe the four stages so a non-specialist can follow the logic, with each paragraph ending on the consequence for the reader rather than on the mechanism. Carries the capability-to-value table. |
| How biochemical context is assigned | Explain the two paths by which a feature gets its pathway and subpathway, at a level a non-specialist can follow, and attach to each path the evidence that actually bounds it. A reader cannot judge a pathway call without knowing whether it was looked up or inferred. |
| Validation and evidence | Keep what was done, what was observed, what it means, and why it matters as four distinct moves, and state the limits of the evidence in the same breath as the result. |
| What you receive | Turn the deliverable into something concrete a reader can picture, and show that coverage and restraint arrive together. |
| How to read a result | Give the reader an operating manual: the key terms in plain language, then what each reporting outcome asserts and what the reasonable next step is. |
| Conclusion | Restate the three takeaways as impact, define where the capability fits against authentic-standard confirmation, and avoid any comparative claim the evidence cannot carry. |
| Availability | Carry the required commercial, legal, and contact placeholders through to marketing and legal review. |

Two things run through every section rather than living in one of them. **User value**
is the closing note of each substantive paragraph. **Limitations** are placed next to
the claims they qualify — the control set is modest and was assembled from available
material rather than sampled at random; a ranked list is not a confirmed answer; the
confidence output ranks and limits claims rather than stating a probability of
correctness; and pathway context orients rather than settles a classification.

## 4. Technical capability, and why it matters to the reader

| Technical capability | Why it matters to the target audience |
|---|---|
| A report for every unnamed feature | Nothing disappears between the instrument and the deliverable, so the unnamed fraction stays inside the biological interpretation instead of being written off. |
| Candidate generation not restricted to the reference library | Novel and study-specific molecules can be described rather than skipped, which is where the differentiating biology often sits. |
| Several independent evidence axes combined into one ordering | Reduces dependence on a single evidence source, so a candidate that looks good on one measure alone does not reach the top of the list unchallenged. |
| Graded reporting aligned to community identification levels | Increases transparency and supports more informed prioritization, because the output is interpretable without learning a bespoke vocabulary. |
| Deliberate withholding when evidence is insufficient | What is shown has earned the right to be shown, so confirmation effort is not spent chasing a structure the data never supported. |
| Pathway context on every feature, including those without a structure | Extends the biological usefulness of features that fall short of confirmation, because an unnamed feature can still be placed in a narrative. |
| Uniform automated application across the whole study | Enables consistent application at scale, so grades are comparable between the first feature and the three-thousandth, and between studies. |
| The supporting evidence shown alongside the call | A reviewer can see what the call rests on instead of trusting a single score, which is what makes the result defensible to a third party. |

## 5. Figure plan

Five figures in the white paper, six in the manuscript. Each entry completes the
sentence the guidance asks for, and the white-paper entries are listed in the order
they appear.

The last entry on each side is a screen capture rather than a drawn chart, so both
are PNG. They are regenerated by `tools/capture_report_screens.py` rather than by
the figure generators, because capture needs a headless browser that neither the
documentation dependencies nor CI install.

**White paper.**

- `fig_workflow_public.svg` — This figure demonstrates that *each stage adds a distinct
  kind of evidence and the final stage decides how much may responsibly be claimed*,
  which matters because *it shows the confidence grade is a product of the workflow
  rather than a label attached at the end*.
- `fig_pathway_public.svg` — This figure demonstrates that *biochemical context reaches a
  feature by one of two paths, one a lookup against a catalogued compound and the other
  an inference from chemically similar ones*, which matters because *the two carry very
  different weight, and a reader who cannot tell them apart cannot tell which pathway
  calls to lean on*.
- `fig_recovery_public.svg` — This figure demonstrates that *the correct structure is
  usually near the top of a short ranked list on controls with independently known
  structures*, which matters because *it converts an unbounded chemical space into a
  shortlist a chemist can actually work through*.
- `fig_report_public.png` — This figure demonstrates that *a delivered report states the
  proposed structure and the grade limiting it in the same place*, which matters because
  *it turns an abstract description of the deliverable into something the reader can
  picture, and shows the grade travelling with the call rather than sitting apart from
  it*.
- `fig_coverage_public.svg` — This figure demonstrates that *every delivered feature is
  reported and placed in a pathway while a specific structure appears only where the
  evidence supports one*, which matters because *it shows completeness and restraint
  are delivered together rather than traded against each other*.

**Manuscript.**

- `fig1_workflow.svg` — This figure demonstrates that *evidence generation, ranking,
  taxonomy assignment, and claim control are separable stages with distinct failure
  modes*, which matters because *each result in the paper can be attributed to a
  specific stage rather than to the pipeline as a whole*.
- `fig2_structure_accuracy.svg` — This figure demonstrates that *recovery improves
  substantially from rank 1 to top 10 and that catalog provenance does not confer
  higher accuracy*, which matters because *it justifies reporting ranked hypotheses and
  refuses the inference that an exact-source call is a confirmed one*.
- `fig3_confidence_reliability.svg` — This figure demonstrates that *the aggregate
  score is poorly calibrated under leave-one-out validation*, which matters because *it
  is the evidence for presenting the score as an ordering device rather than a
  probability of correctness*.
- `fig4_identification_distributions.svg` — This figure demonstrates that *the reporting
  rubric concentrates results in hypothesis-level categories and leaves the strongest
  tier empty*, which matters because *it shows the claim-control rules bind in practice
  rather than only in principle*.
- `fig5_pathway_evaluations.svg` — This figure demonstrates that *apparent pathway
  performance depends entirely on which evaluation context is used*, which matters
  because *quoting a single pathway accuracy number without its provenance would
  misstate the result by tens of percentage points*.
- `fig6_report_example.png` — This figure demonstrates that *the reporting level is
  enforced in the rendered page, with the public and internal renderings of one feature
  differing in how much of the candidate list each may expose*, which matters because
  *claim control was otherwise only ever described, and this shows it binding on the
  artifact a customer actually receives*.

## 6. Manuscript deltas

The manuscript keeps its IMRaD structure and Nature packaging; the guidance is applied
within it. The specific edits:

- **Abstract** — reshape into problem, approach, evidence, and consequence, and close
  on what the workflow is for rather than on what it computed.
- **Introduction** — end the contributions list on the consequence for a reader, not on
  the contribution itself.
- **Toolkit landscape and Discussion** — remove any implicit superiority; restrict
  comparative statements to defensible wording, since no head-to-head benchmark exists.
- **Results** — separate what was observed from what it means, so each subsection
  reports the measurement first and interprets it second.
- **Figure captions** — lead with the message from section 5, keeping the
  methodological detail that follows it.
- **Applications** — state what a user can now do that they could not before, each
  anchored to a specific result and its matching limitation.
- **Conclusion** — end on impact rather than on a restatement of the numbers.
- **Terminology** — use *recovery*, *accuracy*, *agreement*, and *consistency* for
  distinct things and say which is which; keep *identification* separate from
  *annotation* and *structure hypothesis*; keep identification-confidence tiers
  distinct from pathway-confidence tiers. These live in a Methods glossary rather than
  in a paragraph, together with the provenance vocabulary — reference library,
  retrieval, analogue search, de novo generation, catalog provenance, authentic
  standard — which is defined neutrally so *exact source* reads as a record of where a
  candidate came from rather than as a verdict on the route that produced it.
- **Coverage counts** — reconcile the pathway assignment counts in Results with
  `docs/paper/data/pathway/audit_pathway_summary.json`, which they had drifted from.

## 7. Constraints that shape the drafting

The white paper is machine-gated by `tools/verify_public_brief.py`, which is stricter
than a style guide. Three constraints change how sentences get written:

- **Numbers must resolve through the claim ledger.** Any numeral the scanner finds that
  does not trace to `claim_evidence.csv` fails the build, and there is no exemption for
  years or ordinary counts. Everything that is not an audited statistic is spelled out
  in words.
- **Parameters are not disclosable; the conceptual approach is.** Model families,
  descriptors, scoring parameters, thresholds, internal library numbers, pipeline
  versions, and infrastructure names are blocked by a scanner, and existing approaches
  are described as classes rather than named tools. What each stage does in concept may
  be stated, and for pathway assignment it is, because a customer cannot judge a pathway
  call without knowing whether the label was looked up or inferred. Nothing that would
  let a reader reproduce the workflow appears.
- **Overclaim language is blocked**, and budgets cap the document at 2,600 words, five
  figures, three tables, and ten pages. Pages bind before words do once a screenshot is
  in the document, so check `--check-pages` after any edit near a figure.

The manuscript is gated by `tools/verify_paper_package.py` and `tests/test_manuscript_docx.py`:
a 200-word abstract limit, a 75-character title, exactly six figures and one table,
ascending figure numbering, citation parity, and a fixed section order. No new
top-level sections may be introduced. Linked figures must be SVG or PNG and are
validated by format rather than by extension (`tools/figure_formats.py`).

## 8. Review passes

Recorded after drafting, following the four reviews the guidance specifies: argument,
audience, scientific rigor, and style. See `docs/publication_review.md`.
