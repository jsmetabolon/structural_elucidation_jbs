# Decisions on the mass axis

For whoever owns the inference pipeline, the manuscript, and the validation
budget. Each item states what is needed, the evidence, what happens if it is
deferred, and a recommendation.

**Status.** This began as three decisions that analysis could not settle. Decision
1 turned out to be answerable from the delivered data and is now closed apart from
a sanity check. Decision 2 keeps its shape but its recommendation changed once the
measurement was separated from the intervention. Decision 3 is still a real budget
decision, but it now has a cheaper option that did not exist when this was written.

Background and the full measurement are in the decision record at the end of
[docs/paper/REVISION_LOG.md](paper/REVISION_LOG.md) and in
`docs/paper/data/thresholds/audit_mass_policy.json`.

## What was found

The mass sub-score in every delivered report was produced by
`chemistry.conjugates.linear_mass_match_score` at a 1.0 Da default, not by the
`inference.predict.mass_match_score` the documentation described. The two shared
a name and differ 200-fold in tolerance.

Worse, the tolerance varies by which generator produced the candidate:

| Candidate source | n | Tolerance it was scored under |
|---|---:|---|
| `library_anchored_candidates_v1` | 13,127 | 1.0 Da |
| `peptide_first_denovo_v1` | 9,293 | 1.0 Da |
| `balanced` | 2,584 | 5.0 Da |
| `dehydro_v1`, `wgcna_seeded_v1`, `hydroxyl_v1`, `methylation_v1` | ~580 | 0.005 Da |

Each generator is consistent with itself. They disagree with each other. So
candidates are ranked against one another on scales that do not correspond, and
the confidence model averages the result with the other axes as though it were a
single quantity. `PLATT_A` and `PLATT_B` were fitted on that mixture.

## Decision 1: is a 1.0 Da tolerance intended? — answered

**This no longer needs the pipeline's history.** It was posed as a question for
whoever knew why the value was chosen. The history does not contain the answer and
the data does not need it.

The history is a dead end. Both `--mass-score-tolerance-da default=1.0` and the
conjugates default of 0.15 arrive in a single bulk refactor commit, so their origin
predates this repository.

The data settles the substance regardless of intent. Across the truth controls,
candidates that are correct have a **median mass error of 1.07 ppm**, and 81.4% of
them fall within 5 ppm. The measurement is accurate to roughly 1 to 2 ppm, so a
1.0 Da window — about 3,300 ppm at these masses — is some three thousand times
looser than the instrument warrants. **It is indefensible as a mass-accuracy
window whatever was intended**, and the evidence under decision 2 shows it is
actively counterproductive at that setting rather than merely permissive.

Two things follow that a tighter window alone does not fix.

- **17 of 97 candidates scored correct carry a conjugate-scale mass error**,
  concentrated in seven controls, five of them exactly one glucuronide (176.032 Da)
  out. A Tanimoto of 0.70 to the truth structure does not require the same
  molecular formula, so a candidate missing a conjugate passes as correct while
  being wrong by the conjugate's mass. Any window tight enough to be useful
  discards them.
- **One of those is not a conjugate relationship at all, and not what was first
  reported.** `X-21470` carries a candidate structurally *identical* to its truth
  annotation and still 233 Da out. An earlier version of this bullet concluded
  that the control's annotation or its observed mass must be wrong, and called it
  a truth-set integrity item. That was wrong, and it pointed at the wrong owner.

  The candidate mass, 466.1331 Da, is **exactly twice** the observed neutral mass
  of 233.0668. The compound is 5-androstenetriol 3,17-disulfate, `C19H30O9S2`: a
  disulfate, which readily forms a doubly charged `[M-2H]2-` anion in negative
  mode. The structure, the Tanimoto of 1.0 and the annotation are all correct.
  What is wrong is the conversion from m/z to neutral mass, which assumed a single
  charge. It is a **pipeline defect, not a curation task**, and it is now routed to
  handoff item 6.

So a single linear tolerance is the wrong instrument for two jobs. Mass accuracy
wants a ppm-scale window; conjugate relationships want the conjugate-aware
generation already proposed as fix 3 of handoff item 2.

**What remains for a person:** confirming that 1.0 Da was not a deliberate
scientific choice that this analysis has misread. That is a sanity check rather
than a blocker, and decisions 2 and 3 no longer wait on it.

## Decision 2: unify the scorers, rescore, and recalibrate?

**Owner: pipeline and manuscript jointly.**

Evidence that it may not be worth it. Sweeping the tolerance from 0.005 Da to
5.0 Da, a thousand-fold range:

| Tolerance | Rank-1 correct of 92 | Gained | Lost | Rank-1 strong-mass claims |
|---:|---:|---:|---:|---:|
| 0.005 Da | 41 | 0 | 3 | 2,057 |
| 0.05 Da | 42 | 0 | 2 | 2,756 |
| 0.5 Da | 42 | 0 | 2 | 3,068 |
| 1.0 Da | 42 | 0 | 2 | 3,085 |
| 5.0 Da | 43 | 0 | 1 | 3,105 |

Recovery moves by two controls across the entire range and no setting gains one.
What does move is what the reports claim: a thousand-report swing in strong-mass
calls.

**Correction to an earlier reading.** This section previously concluded that the
axis "relabels without discriminating". That took a result about the *gate* and
stated it about the *measurement*, and they turn out to differ. Splitting mass
error by known correctness across the truth controls:

| Window | Correct within | Incorrect within | Separation |
|---:|---:|---:|---:|
| 1 ppm | 48.5% | 37.3% | +11.2% |
| 2 ppm | 75.3% | 53.4% | +21.9% |
| 5 ppm | 81.4% | 59.5% | +21.9% |
| 20 ppm | 82.5% | 62.9% | +19.6% |
| ~3,300 ppm (the shipped 1.0 Da) | 82.5% | 90.1% | **−7.6%** |

The measurement discriminates, most sharply at 2 to 5 ppm. The intervention built
on it does not. And at the shipped tolerance the separation inverts: a window that
loose admits proportionally more incorrect candidates than correct ones, because
17 of 97 candidates scored correct carry a conjugate-scale error and sit outside
any window. So the axis at 1.0 Da is not neutral, it is very slightly worse than
neutral.

The cost of acting is a pipeline re-run, a Platt refit, and re-auditing every
number in both documents.

**Recommendation: unify the scorers onto one instrument-scale window.** Two
reasons, and the second is new. Ranking candidates scored on scales that differ
200-fold is wrong on its own terms, whatever the truth set can detect. And the
measurement does carry signal that the shipped tolerance spends: separation is
about 22 points at 2 to 5 ppm and negative at 1.0 Da. What is *not* supported is
keeping the current gate and tuning its threshold, which is what the sweep tested
and what it ruled out.

## Decision 3: fund truth-set expansion?

**Owner: whoever holds the validation budget.** This is the one that matters most.

Only **9 of the 598** rank-1 mass-mismatch reports have a known correct answer.
That is the entire evidence base for any claim about the other 589, including
handoff item 2's projection of roughly 523 recoverable structures.

The consequences are already visible. Decisions 1 and 2 rest on judgement because
measurement cannot reach them. Handoff item 2 cannot be validated before shipping.
And rank-1 recovery on 92 controls carries a twenty-point confidence interval, so
no future rerank change will be measurable either.

**Recommendation: treat handoff item 4 as a prerequisite rather than a parallel
improvement.** Every other item in the refinement handoff is gated behind it.
Expanding truth coverage of the mass-mismatch cohort specifically would be worth
more than expanding the truth set uniformly, since that cohort is where the
proposed changes act.

### A second option, which handoff item 4 does not consider

Item 4 proposes one route: acquire more authentic standards. That is a chemistry
budget and a long lead time. There is a second route that costs compute instead,
and the reason it is available is the thing the writing guidance already says
about the library — a large Level-1 library empowers prediction. It is also the
validation set.

**Metabolon already holds thousands of compounds with both a known structure and
an acquired spectrum: the reference library itself.** Running those spectra
through the pipeline with each compound suppressed from retrieval turns every one
of them into a held-out control. No new chemistry, no new acquisition.

There is precedent in this repository for exactly that reasoning. The pathway
audit runs leave-one-out over **5,095 reference metabolites**, hiding each
metabolite's own labels and asking the mechanism to recover them, and reports the
result as mechanism validation rather than as accuracy.

**Two honest limits, and neither is fatal.**

- **Catalogued compounds are not unknowns**, and this work already measures the
  gap: exact-source rank-1 recovery is 38.1% against 56.0% for predicted-source.
  So this cannot replace truth on real unknowns, and any absolute rate from it
  must be labelled as measured on catalogued chemistry.
- **It is not as cheap as the pathway version.** That one replays offline from the
  reference database because a neighbour vote over fingerprints can be recomputed
  from the database alone. Structural candidate generation cannot: it needs the
  spectra, the models and the inference stack. This is an EC2 run, not a laptop
  script.

**What it is good for is the question that blocks everything else.** Deciding
whether a reranking change helps needs paired comparisons, not a representative
sample. Thousands of them would resolve shifts that 92 controls cannot see at all,
and the same harness would then measure handoff item 2 before it ships rather
than after.

Costed in [docs/library_holdout_validation_scope.md](library_holdout_validation_scope.md).
Read the first section before funding anything: the whole proposition rests on how
many library standards carry a usable acquired spectrum, and that number is not
answerable from this repository.

## Settled since this was written

**The manuscript typeface.** It was open here as a decision needing a person, and
it has been taken: the manuscript is an internal document and is branded, body
face and heading colour together. Reversing it is no longer the five-line job
described here. `BRANDED` in
[tools/build_manuscript_docx.py](../tools/build_manuscript_docx.py) governs face
and colour together, and the release gate and the test read `BODY_FONT` from the
builder rather than naming a face, so a journal-plain submission is one flag.
