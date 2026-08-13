# Three decisions needed on the mass axis

For whoever owns the inference pipeline, the manuscript, and the validation
budget. Each item states what is needed, the evidence, what happens if it is
deferred, and a recommendation. None of the three can be settled by more
analysis of the current data; that is why they are here rather than in a ticket.

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

## Decision 1: is a 1.0 Da tolerance intended?

**Owner: whoever knows the pipeline's history.** This one is cheap and blocks the
other two.

At 300 Da, 1.0 Da is roughly 3,300 ppm. It admits essentially every candidate the
generator proposes, which is why rank-1 strong-mass claims sit at 3,085 of 3,897
under it. The profile carrying the value is named `strict_1da`, which reads as
deliberate, but one Dalton is not strict for high-resolution data.

The answer decides whether the sweep below is a diagnosis or a post-mortem.

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
calls. The axis relabels without discriminating.

The cost of acting is a pipeline re-run, a Platt refit, and re-auditing every
number in both documents.

**Recommendation: fix the incomparability, do not chase the tolerance.** Putting
all generators on one scale is defensible on its own terms, because ranking
candidates scored on different scales is wrong regardless of what the truth set
can detect. Tuning the tolerance is not, because the evidence says no value of it
changes correctness. If the re-run happens, do it to make the axis coherent, not
to make it tighter.

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

## Also outstanding, unrelated to the above

The manuscript body is currently set in Montserrat. It was Times New Roman until
the brand pass, and `docs/paper/README.md` still describes the file as a Nature
initial submission. If it is going to a journal, the body should revert while the
figures stay branded: about five lines in
[tools/build_manuscript_docx.py](../tools/build_manuscript_docx.py) plus two
assertions. Flagged at the time and never decided.
