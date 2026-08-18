# MSI sub-level confidence bands: what the cutoffs can and cannot be

Answers the two questions that block the sub-level work — *what should the new
confidence bands be*, and *which score do they apply to* — against the delivered
3,897-report set. Every number here is produced by
[tools/audit_sublevel_bands.py](../tools/audit_sublevel_bands.py) and committed to
[docs/paper/data/identification/audit_sublevel_bands.json](paper/data/identification/audit_sublevel_bands.json);
the readable rendering is
[audit_sublevel_bands.md](paper/data/identification/audit_sublevel_bands.md).

## The short answer

**No cutoff on either candidate score produces a usable band scheme, so the
banding decision cannot be made yet.** This is not a case of picking the wrong
numbers. Both scores fail a precondition that no choice of cutoff can repair:

- `evidence_confidence` (the matrix **Confidence Score**) is **flat** where the
  bands would go. 62.4% of compounds sit at or above 0.98, and across the top
  seven deciles P(rank-1 correct) moves only between 13.6% and 19.1%. It
  separates *some evidence* from *no evidence at all* and nothing finer.
- `candidate_confidence` (the report overall) is **inverted** once the
  authentic-standard-confirmed reports are held out: AUC 0.378 against truth,
  where 0.50 means no information. Its lowest decile is its most accurate.

Two further findings constrain whatever is decided next, and both are
independent of which score is chosen. They are in sections 1 and 4.

## 1. `TIER_HIGH = 0.65` is unreachable, so `2a` and `3a` would be empty

The evidence score is `sigmoid(A * raw + B)` with `A = 1.887084`,
`B = -1.520154`, and `raw` is a weighted average of 0..1 axes, so `raw` cannot
exceed 1.0. That bounds the score to **[0.1794, 0.5907]** no matter how strong
the evidence is. Inverting the map at each tier cutoff:

| Cutoff | Value | Raw evidence required | Reachable |
|---|---:|---:|---|
| `TIER_HIGH` | 0.65 | 1.1336 | **no** |
| `TIER_MODERATE` | 0.48 | 0.7631 | yes |
| `TIER_LOW` | 0.38 | 0.5461 | yes |

A report can only land above 0.65 through a policy override — `CONFIRMED_OVERALL`
for an authentic standard, or the curator-identity promotion — and both of those
already imply Level 1. The delivered set agrees: the tier distribution in
[audit_presentation_summary.json](paper/data/presentation/audit_presentation_summary.json)
contains **no `High` tier at all** across 3,897 reports, only
`Moderate` 1,416, `Low` 1,344, `Tentative` 1,093 and `Confirmed` 44.

So the example bands in the request, applied to `candidate_confidence`:

| MSI level | ≥0.65 (`a`) | 0.48–0.65 (`b`) | 0.38–0.48 (`c`) | <0.38 |
|---|---:|---:|---:|---:|
| L1 | 44 | 0 | 0 | 0 |
| L2 | **0** | 172 | 21 | 0 |
| L3 | **0** | 1,246 | 1,062 | 181 |
| L4 | 0 | 0 | 785 | 380 |
| L5 | 0 | 0 | 3 | 3 |

`2a` and `3a` come out empty, and the entire `a` column is the 44 reports that
are already Level 1. This holds for `TIER_HIGH` at any value above 0.5907.

## 2. The bands are not ordered by correctness

P(rank-1 candidate correct) inside each band, over the 211 model-scored
truth-verified reports. A scheme that reads as confidence has to fall left to
right.

**`candidate_confidence`:**

| Cuts | `a` | `b` | `c` | below | Ordered |
|---|---:|---:|---:|---:|---|
| 0.65 / 0.48 / 0.38 (tiers) | n/a | 6.1% | 9.9% | 23.5% | **no** |
| 0.585 / 0.48 / 0.40 | 25.0% | 1.8% | 8.5% | 19.1% | **no** |
| 0.55 / 0.46 / 0.40 | 13.0% | 4.3% | 7.7% | 19.1% | **no** |

The ordering is backwards in every scheme: `c` beats `b`, and the reports below
every cutoff beat both. Shipping this would tell a reader that a `2c` assignment
is a weaker claim than a `2b` one when the opposite is measured.

The apparent AUC of 0.7865 quoted for this score elsewhere is an artifact of the
44 policy-pinned reports. They are correct *by construction* — truth matching is
what makes a report Level 1 — and their displayed score is a display value, not
model output. Held out, AUC is **0.378**.

**`evidence_confidence`:** the bands are degenerate rather than inverted. Because
65.7% of compounds clear 0.65, band `a` absorbs every correct call and `b`, `c`
and below all read 0.0%. Its AUC on model-scored reports is 0.6566, but the
deciles show where that comes from — deciles 1–2 sit at exactly 0.0 and score 0%
correct, while deciles 4–10 all sit above 0.9955 and score a flat 13.6–19.1%. It
is a usable binary flag for *no evidence at all*; it is not a graded score.

## 3. Neither score separates the curated sub-levels

Metabolon's curated `SubLevel` in `data/known_X.csv` is the only ground truth for
what a sub-level means here. Spread between the highest and lowest median across
the fine letters `2C`, `3A`, `3B`, `3C`:

| Score | Fine-letter median spread |
|---|---:|
| `candidate_confidence` | 0.0139 |
| `evidence_confidence` | 0.0024 |

Both are inside the noise. This reproduces, on a different cohort and by a
different method, the conclusion already reached in `docs/msi_sublevel_findings.md`
on branch `cursor/candidate-confidence-rework-f8a2` (not merged, so there is
nothing to link to here): a derived rubric reproduces the curated letter only
9.5% of the time and inverts the fine ordering at the top of Level 3.

## 4. `2a`/`2b` are not confidence bands in the scheme being cited

Worth settling before cutoffs are chosen at all, because it decides whether a
score is the right input. In Schymanski 2014 — reference [1] of the manuscript,
and the scheme
[`assign_msi_level`](../structural_elucidation/reports/confidence.py) says it
follows — the Level 2 sub-letters are **evidence-type** distinctions, not score
ranges:

- **2a** — probable structure by *reference library spectrum match*
- **2b** — probable structure by *diagnostic evidence*, where no standard or
  library spectrum is available but the evidence admits one structure

Metabolon's curated letters are the same kind of thing: they are expert
judgments encoded in the `Tentative information` free text, which is why they are
`2C` and `3A`–`3D` rather than numbers. Attaching `2a (0.9–1.0)` to a score band
makes a claim about *how much* evidence there is, where the scheme makes a claim
about *what kind*. A compound with a library MS² match is `2a` whether its score
is 0.95 or 0.55.

This matters practically because the repo cannot currently compute the
distinguishing evidence. There is no reference MS² spectral library and no
candidate-versus-observed diagnostic-fragment matching, so `2a` and `2b` are not
separable from anything on hand at any threshold.

## What is available instead

Three things are defensible today, in decreasing order of confidence:

1. **The curated `SubLevel`, surfaced not derived.** 2,494 of the 3,897
   delivered compounds carry an expert letter from `data/known_X.csv`, 949 of
   them a fine letter (`2C`, `3A`–`3D`). It is authoritative and needs no model.
   This also makes a live `msi_level_counts.csv` computable now — the copy on
   branch `cursor/evidence-confidence-calculation-012b` was hand-committed with
   no generator behind it.
2. **The coarse MSI level 1–5**, which is already assigned for 100% of reports.
   One caveat that should not be buried: on this truth-verified cohort it is not
   cleanly monotone either — L2 is 3.3% correct (n=30) against L3 at 13.3%
   (n=128). The L2 cell is one correct call out of thirty, so it is thin, but it
   is not evidence that the level ordering holds.
3. **A two-way split on `evidence_confidence`** at the zero boundary, which is
   the one distinction that score genuinely makes. It is a data-availability
   flag, and should be labeled as one rather than as confidence.

## What would unblock the original request

Bands become answerable once the score they are cut on is monotone with
correctness over its own output. Concretely, one of:

- **Recalibrate the Platt fit** so the score spans a usable range. The current
  `(A, B)` predate the retention axis and were fit in-sample on 92 controls; the
  ceiling at 0.5907 is a direct consequence. This is the smallest change and it
  fixes section 1, but on its own it will not fix section 2 — rescaling a score
  that is inverted leaves it inverted.
- **Establish a monotone per-candidate probability**, e.g. the isotonic
  `rerank_score_calibrated` that `candidate_confidence_column` on branch
  `cursor/candidate-confidence-rework-f8a2` prefers when present, validated
  out-of-sample on a truth cohort larger than 211.
- **Acquire the evidence the sub-letters actually name** — a reference MS²
  spectral library, observed retention time per X-compound, and diagnostic
  fragment matching. This is the only route that yields `2a`/`2b` as the cited
  scheme defines them rather than as a renamed score bin.

The ordered work for each of these is in
[plans/msi_sublevel_bands.plan.md](../plans/msi_sublevel_bands.plan.md).

## Reproducing

```
python3 tools/audit_sublevel_bands.py
```

Read-only, standard library only, and reproducible from committed inputs. It
reads `docs/paper/data/exec_summary/exec_summary_by_compound.csv` and the
constants in `structural_elucidation/reports/confidence.py`, and rewrites the two
artifacts under `docs/paper/data/identification/`.
