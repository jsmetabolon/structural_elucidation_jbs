# Decisions on MSI sub-levels

For whoever owns the report layer, the client deliverable, and the validation
budget. Each item states what is needed, the evidence, what happens if it is
deferred, and a recommendation.

**Status.** This began as two questions — what the sub-level confidence bands
should be, and which score they apply to. Both turned out to be answerable from
the delivered data, and both answers are negative, so they are closed below
rather than routed to a person. What replaced them is a decision nobody had
written down: a month of finished sub-level work is sitting unmerged in
[PR #6](https://github.com/jsmetabolon/structural_elucidation_jbs/pull/6), and
its central validation claim does not survive holding out the reports that are
correct by construction.

One item moved while this was being written. Decision 7 was a question about an
unmerged branch; that branch reached `main` in `27e1f60`, so it is now a fix to
shipped code. Every measurement here has been recomputed against the post-merge
tables and none of the conclusions changed — see "What the 012b merge did to this
audit" below.

The measurement is [tools/audit_sublevel_bands.py](../tools/audit_sublevel_bands.py),
regenerating [docs/paper/data/identification/audit_sublevel_bands.json](paper/data/identification/audit_sublevel_bands.json).
The reasoning is in [docs/msi_sublevel_band_findings.md](msi_sublevel_band_findings.md).
The staged work is in [plans/msi_sublevel_bands.plan.md](../plans/msi_sublevel_bands.plan.md).

## What was found

Three things, on the delivered 3,897-report set.

**The tier cutoffs cannot all be reached.** The evidence score is
`sigmoid(A·raw + B)` with `A = 1.887084`, `B = -1.520154`, and `raw` is a
weighted average of 0..1 axes, so the score is bounded to `[0.1794, 0.5907]`.
`TIER_HIGH = 0.65` needs `raw = 1.13`. The delivered set agrees: there is no
`High` tier at all across 3,897 reports.

**The truth cohort has a trapdoor in it.** 255 compounds are truth-verified, but
44 of them are authentic-standard-confirmed Level 1 reports, which are correct
*by construction* — truth matching is what makes a report Level 1 — and which
carry a policy-pinned display score rather than model output. Any statistic that
leaves them in the denominator measures the pin. This is what separates the two
numbers in circulation: AUC 0.7865 with them, **0.378** without.

**Neither candidate score bands.** `evidence_confidence` is flat where the bands
would go: 62.4% of compounds at or above 0.98, and P(correct) moving only between
13.6% and 19.1% across its top seven deciles. `candidate_confidence` is inverted
on the model-scored population. Fine-letter median spread is 0.0139 and 0.0024
respectively, both inside the noise.

## Decision 1: what happens to PR #6 — needs a person

**Owner: report layer, jointly with whoever owns client delivery.** This is the
most urgent item here and the only one with a running cost.

[PR #6](https://github.com/jsmetabolon/structural_elucidation_jbs/pull/6) has
been open as a draft since 15 July. It already builds most of what the sub-level
work needs: the curated `SubLevel` surfaced on 2,484 reports, a derived `3A`,
client include/exclude columns in the master index, a rank-consistent
`candidate_confidence`, the 2A/2B deferral handoff, and a regenerated
3,897-report bundle. It is not speculative work — it is finished work.

Three of its parts stand on their own and one does not.

**Stands: the rank-consistency fix.** 58.8% of reports had at least one
rank/score inversion in the candidate widget, and the fix takes that to zero. A
widget that lists candidates 1 to 10 and scores them out of order is wrong on its
own terms, with no reference to truth.

**Stands: surfacing the curated `SubLevel`.** Expert-labeled, authoritative, and
it needs no model. 2,494 of 3,897 delivered compounds carry one on the
committed exec summary, 949 of them a fine letter.

**Stands, with a caveat: the un-deflation.** Re-anchoring rank 1 from the Platt
`overall` to the pre-Platt `raw_evidence` fixes a real display defect. Recovering
`raw_evidence` by inverting the Platt map puts its median at 0.6577, against the
0.657 PR #6 reports from its own bundle, so the two are the same quantity. On
that scale the tier cutoffs become reachable and the top band goes from 44
reports to 1,996 — the empty `2a`/`3a` problem is genuinely fixed.

**Does not stand: the validation attached to it.** PR #6 cites "still tracks
truth (AUC ~0.79)" and proposes `candidate_confidence` as the client
include/exclude score. But `sigmoid(A·raw + B)` is strictly increasing, so
re-anchoring is a rescale and cannot reorder anything. Measured rather than
argued, AUC on the model-scored population is **0.378 for both scores, the same
number to four decimal places**. The cited 0.79 is carried by the 44 pinned
reports. Verified separately: no non-Level-1 report exceeds the Platt ceiling
(maximum 0.5860 against 0.5907), so every one of them carries an unpinned score
and the equivalence is exact rather than approximate.

**If deferred:** the rank-inversion defect stays in the delivered reports, the
curated sub-level stays invisible, and the branch keeps drifting from `main`.
That cost is being paid now.

**Recommendation: split it.** Merge the rank-consistency fix, the curated
`SubLevel` surfacing, the un-deflation, and the "N of 10" note. Hold the
`Match conf` include/exclude column pending decision 2, and drop the AUC 0.79
claim from the description wherever it appears, since it describes a cohort
nobody filters on.

### Commit-level routing

Ten of the fifteen commits merge unchanged. Three need a wording or naming
amendment, one needs splitting, and none needs to be abandoned.

| Commit | Disposition |
|---|---|
| `2d13f6a` rank-consistent candidate-confidence helper + audit | merge |
| `01d3c85` rank-consistent candidate confidence in the widget | merge |
| `6dbcdf5` emit per-rank `candidate_confidence` into the sidecar | merge — closes a live gap, see below |
| `dcb6cf5` v2 candidate-score patcher | merge |
| `dd5c764` monotonicity, ties and patcher tests | merge |
| `56b9159` emit monotone `candidate_confidence` from the reranker | merge |
| `d025edd` reports-tree-only refresh | merge |
| `394c748` cover canonical `X_<id>.html`, preserve evidence on refresh | merge |
| `d79417c` rank-1 tie rate | merge |
| `6806176` real-bundle before/after | merge |
| `9dd75e8` handoff for the evidence the fine sub-levels need | merge |
| `e175b43` sub-level feasibility and defensibility audit | merge — conclusions corroborated independently |
| `efa280c` un-deflate the best match, "N of 10" | **amend** the justification: it cites AUC ~0.79 for `raw_evidence`, which is the whole-cohort figure |
| `4599439` surface curated `SubLevel` + derived `3A` | **split**: curated surfacing merges as is; the derived `3A` needs a name that is not `3A` |
| `9e6e359` index `Curated` + `Match conf` columns | **split**: `Curated` merges; `Match conf` is gated on decision 2 |

One merge is worth calling out because it fixes something already broken on
`main` rather than adding anything. `tools/build_exec_summary.py:242` reads
`candidate_confidence` out of the confidence sidecar, and nothing on `main`
writes it — `6dbcdf5` is the missing producer. Until it lands, that column
resolves to empty for any freshly enriched bundle, and the populated values in
the committed exec summary came from a bundle enriched off this branch.

## Decision 2: may `candidate_confidence` be a client-facing include/exclude score — needs a person

**Owner: product and commercial.**

A filter acts on the compounds a client is unsure about, which is the non-Level-1
population — exactly where AUC is 0.378. Below 0.50 means the ordering is
reversed, so a threshold filter on this score removes proportionally more correct
calls than incorrect ones. Keeping everything would outperform it.

The band precision makes the same point in the form a client would see. Under the
tier cutoffs, on the model-scored cohort, P(rank-1 correct) reads 6.1% in band
`b`, 9.9% in band `c`, and 23.5% below every cutoff.

**If deferred:** nothing ships that depends on it, and decision 1's larger part
can proceed without it.

**Recommendation: no, on current evidence.** Not "not yet, pending tuning" — a
threshold cannot repair a reversed ordering. This reopens only if decision 5 is
funded and produces a score that passes the gate in stage 6 of the plan.

## Decision 3: is the sub-level a filter or a display grade — needs a person

**Owner: product.**

PR #6 already assumes the former: its master-index columns are described as "for
client include/exclude". The screenshot that started this carries
`INCLUDE IN DATA SET` and `INCLUDE IN PATHWAYS` columns, which points the same
way. This asks for that to be confirmed rather than inferred, because the two
uses have different requirements and only one of them is blocked.

A **filter** needs a defensible cut point, which decision 2 says does not exist
yet — but a filter on the *curated* letter needs no model at all and could ship
today for the 2,494 compounds that carry one. A **display grade** needs an
ordered scale across all 3,897, which is what nothing currently provides.

**If deferred:** the sub-level work optimises for the wrong requirement.

**Recommendation: confirm it is the filter, and scope it to the curated letter.**
That is deliverable now and needs none of decisions 5 or 6.

## Decision 4: is the screenshot's scheme a real Metabolon specification — needs a person

**Owner: whoever produced the slide.**

The requested vocabulary — `2a`, `2b`, `2c`, `3a`, `3b`, `3c`, `5a`, with ranges
like `2a (0.9–1.0)` — does not match the curated data.
[data/known_X.csv](../data/known_X.csv) carries `2C` and `3A`–`3D` and nothing
else; there is no `2A`, no `2B` and no `5A` anywhere in it. The screenshot is
also watermarked "Fake Data".

There is a second, deeper mismatch. In Schymanski 2014 — reference [1] of the
manuscript, and the scheme `assign_msi_level` says it follows — the Level 2
sub-letters are *evidence-type* distinctions: `2a` is a probable structure by
reference library spectrum match, `2b` by diagnostic evidence. They say what kind
of evidence exists, not how much. Attaching a confidence range to them is a
category change, not a parameterisation, and this repo can compute neither
distinguishing input.

**If deferred:** any definition table written now is a reconstruction, and may
have to be redone.

**Recommendation: get the source before writing definitions.** If a real
specification exists it replaces stage 1 of the plan. If the slide is
illustrative, say so on the record, because the `2a (0.9–1.0)` framing is what
made these bands look like a parameter choice rather than a research problem.

## Decision 5: fund the score fix — needs a person

**Owner: pipeline, plus whoever holds the validation budget.**

Two pieces, and only the second is expensive. Re-fitting the Platt `(A, B)`
against the retention axis fixes reachability, but the un-deflation in decision 1
already achieves that for display purposes, so it is no longer urgent on its own.
Establishing a score that is *ordered* is the real work, and it needs a truth
cohort larger than 211 model-scored reports of which 23 are correct.

**Note the overlap, because it changes the price.** This is the same ask as
Decision 3 of [docs/mass_axis_decisions_needed.md](mass_axis_decisions_needed.md),
and the library-holdout route costed in
[docs/library_holdout_validation_scope.md](library_holdout_validation_scope.md)
would serve both. Funding them separately would buy the same cohort twice.

**If deferred:** decisions 2 and 3's display-grade option stay closed
indefinitely, and no future reranking change is measurable either.

**Recommendation: treat it as one budget line with the mass-axis item, and do not
re-fit Platt on its own** — a monotone rescale of an inverted score is still
inverted, which is precisely what decision 1 measured.

## Decision 6: fund the upstream evidence for genuine 2a/2b — needs a person

**Owner: chemistry and acquisition budget.**

Three inputs, none of which exist here: a reference MS² spectral library, an
observed retention time per X-compound (feature tables carry only m/z), and
candidate-versus-observed diagnostic fragment matching. These are what the
Schymanski sub-letters name, so this is the only route that produces `2a` and
`2b` as the cited scheme defines them rather than as renamed score bins.

**If deferred:** `2a`/`2b` remain unavailable at any threshold, which is a
statement worth being able to make to a client with a reason attached.

**Recommendation: gate on decision 4.** If the slide is illustrative this may not
be wanted at all.

## Decision 7: the matrix workbook's "Confidence Score" — needs a person

**Owner: whoever owns the client workbook.**

**Now a fix rather than a gate.** This was written while
`cursor/evidence-confidence-calculation-012b` was an unmerged branch; it landed
on `main` in `27e1f60` before the decision was taken, so both items below are
shipped code. Neither is hard to change and neither has reached a client yet as
far as this repository can tell, but the deferral cost is no longer zero.

It carries two things that should not have shipped as they stand.
`tools/build_all_metabolites_matrix_workbook.py` presents
`rank1_evidence_confidence` in a column headed **Confidence Score**, and that
column is 62.4% saturated above 0.98 and flat against truth above 0.9955. The
single real distinction it makes is at zero: the bottom two deciles sit at
exactly 0.0 and score 0% correct, everything above scores a flat 13.6–19.1%. It
is an evidence-availability flag presented as a confidence.

It also carries a hand-committed `docs/paper/data/msi_level_counts.csv` with no
generator behind it, so it cannot be refreshed or checked. It mixes case as
though `2c` and `2C` were different buckets, and lists `Exclude`, `Artifact` and
`Internal Standard` as if they were sub-levels.

**If deferred:** the workbook keeps presenting an evidence-availability flag as a
confidence score, and the count table stays unrefreshable. Both are on `main`
now, so deferring means shipping them rather than not shipping them — which
inverts the usual direction of caution on this page.

**Recommendation: rename the column on `main`.** Something that describes it, or
demote the float and expose the boolean beside it; update the workbook's
`Data_Dictionary` sheet, which currently reads "Canonical rank-1
evidence_confidence". Replace the count table with the generator in stage 2 of
the plan. Note the numbers above are recomputed against the post-merge tables and
still hold: the column's AUC improved to 0.5794 whole-cohort when 012b
canonicalised it, and it still bands into nothing ordered.

## Merge order

Three lines of work are outstanding and two of them touch the same files, so the
order matters. It cuts across decisions 1 and 7 rather than belonging to either.

**`cursor/evidence-confidence-calculation-012b` merged to `main` in `27e1f60`
while this was being written**, ahead of the order recommended here and without
the amendment decision 7 asks for. That is not reversible work and there is no
need to reverse it, but it changes what decision 7 is about: the matrix workbook
and its `Confidence Score` column are now on `main` rather than pending on a
branch, so the question is a fix to shipped code rather than a gate on a merge.

Two lines of work remain, and the order is unchanged for them.

**1. [PR #35](https://github.com/jsmetabolon/structural_elucidation_jbs/pull/35),
the audit, first.** Zero file overlap with PR #6, so it cannot conflict, and it
is the measurement PR #6 is judged against. Merging it first means the CI
regenerate-and-diff gate is in place before PR #6 lands.

**2. [PR #6](https://github.com/jsmetabolon/structural_elucidation_jbs/pull/6)
next, split per the commit routing above.** Now that 012b is on `main`, expect
the two overlaps it used to have with 012b to appear as conflicts against `main`
instead: `structural_elucidation/pipeline/rerank.py` (PR #6 adds 23 lines, 012b
changed 18 and removed 4) and `docs/pipeline_refinement_handoff.md`. Neither is
large, but PR #6 was branched before 012b existed, so resolve them deliberately
rather than taking either side wholesale.

### What the 012b merge did to this audit

It regenerated `docs/paper/data/exec_summary/exec_summary_by_compound.csv` and
`exec_summary_top10.csv`, which are this audit's inputs, and added a
`rank1_evidence_confidence_basis` column. So
`docs/paper/data/identification/audit_sublevel_bands.{json,md}` had to be
regenerated on top of the merge. CI catches this — `tools/audit_sublevel_bands.py`
is in the regenerate-and-diff block — so the failure mode is a red build rather
than a stale conclusion, and that is what happened.

**The conclusions themselves survived the refresh.** Every load-bearing result
reproduces unchanged: `TIER_HIGH` still unreachable, the re-anchor still fixing
reachability and not ordering with AUC 0.378 either side, derived `3A` still at
16.1% agreement, and no band scheme ordered for any of the three scores. The only
figure that moves is `evidence_confidence`, whose AUC goes from 0.4762 to 0.5794
whole-cohort and 0.6566 to 0.6005 model-scored, because that branch canonicalises
the column. Still no ordered bands, so no decision here changes. None of these
findings is an artifact of which exec-summary revision you read.

## Answered without a person

**Which score drives the bands.** Neither. Both failure modes are properties of
the scores rather than of any threshold, so no choice of cutoff repairs either.

**Whether the example cutoffs work.** No. `TIER_HIGH = 0.65` is above the Platt
ceiling, so `2a` and `3a` are empty at any cutoff above 0.5907. On the
re-anchored `raw_evidence` scale this specific problem goes away — the top band
holds 1,996 reports — but the ordering problem does not.

**Whether the fine letters can be derived.** No, and two independent audits
agree. This one measures a fine-letter median spread of 0.0139; PR #6's own
`audit_sublevel_evidence.py` reproduces the curated letter 9.5% of the time with
the fine ordering inverted.

**Whether the derived `3A` may keep its name.** No. The rule
(`msi_level == 3 and n_isomers >= 2`) fires on 2,129 reports, 85.5% of all Level
3, against 358 curated `3A`. Of the 1,245 it labels that a curator also labeled,
**16.1%** carry curated `3A`; the largest group inside it is 623 compounds the
curator called Level `5`. The underlying signal — multiple same-formula
candidates in the top 10 — is real and worth surfacing as isomeric ambiguity.
It is not the curator's `3A`, and sharing the letter guarantees it is read as
though it were. This one is a naming fix, not a research problem, and it belongs
with decision 1: rename it to what it measures — isomeric ambiguity — keep it in
a field of its own, and never render it where the curated letter goes.
