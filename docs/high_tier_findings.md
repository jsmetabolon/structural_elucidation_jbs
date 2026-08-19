# What "High" can be made to mean

Answers whether the identification-confidence `High` tier can be given a meaning
it does not currently have. `High` is empty: no report reaches it, across all
3,897. The question was whether to lower the cutoff, recalibrate underneath it,
or leave it alone.

Every number here is produced by [tools/audit_thresholds.py](../tools/audit_thresholds.py)
and committed to
[docs/paper/data/thresholds/audit_thresholds_summary.json](paper/data/thresholds/audit_thresholds_summary.json),
with the reachability arithmetic in
[tools/audit_sublevel_bands.py](../tools/audit_sublevel_bands.py). The staged work
is in [plans/msi_sublevel_bands.plan.md](../plans/msi_sublevel_bands.plan.md).

## The short answer

**Leave `High` alone, and leave it empty.** Not because moving it is hard, but
because the cohort that would have to justify a new cutoff cannot justify any
cutoff — its correct answers and its already-confirmed reports are the same 44
compounds.

Three findings, in the order they close the question.

## 1. No cutoff on today's score is both reachable and defensible

The evidence score is `sigmoid(A·raw + B)` with `A = 1.887084`, `B = -1.520154`,
and `raw` is a weighted average of 0..1 axes, so the score is bounded to
`[0.1794, 0.5907]`. That gives three regions and no landing zone:

| `TIER_HIGH` | Reachable | Defensible |
|---|---|---|
| above 0.5907 | no — the tier stays empty | n/a |
| 0.5233 to 0.5907 | yes | no — 52.4% accurate at 0.55 (n=21) against a 47.8% base rate, and **0.0% at 0.60** (n=9) |
| at or below 0.5233 | yes | no — exact catalog matches begin reading High, and they are 38.1% correct |

The middle row is the interesting one. It is the only region where a cutoff both
works mechanically and clears the exact-match guard in
`tests/test_confidence.py::test_exact_is_not_pinned_high`, and everything in it
is measured at a coin flip or worse.

## 2. Calibrating per source helps, and still is not enough

`source_stratification` has reported `confounded` for some time, with a docstring
prescribing the fix: pooled raw AUC 0.5914 sits *below* both strata (exact
0.6719, predicted 0.7784), because the two sources have opposite base rates
(0.381 and 0.560) and pooling them ranks provenance rather than correctness.

That prescription had never been carried out. It now is
(`per_source_calibration`), leave-one-out inside each stratum:

| Source | n | Base rate | Raw AUC | Pooled calibration | Per-source calibration |
|---|---:|---:|---:|---:|---:|
| exact | 42 | 38.1% | 0.672 | 0.517 | 0.522 |
| predicted | 50 | 56.0% | 0.778 | 0.675 | 0.739 |

Two results. Splitting the calibration **does** help — it beats the pooled fit in
both strata. And the **uncalibrated score beats both**, in both strata. The Platt
step is spending ordering rather than adding calibration, so if a boundary is
ever drawn it belongs on the raw evidence score, within a source, not on a
calibrated probability.

Leave-one-out is load-bearing rather than fastidious here: a single Platt fit is
monotone, so an in-sample "calibrated AUC" would restate the raw one exactly and
the comparison would be vacuous.

## 3. The cohort cannot answer the question at all

This is what actually closes it, and it invalidates the encouraging numbers in
section 2 as evidence for a tier.

**All 44 correct calls in the 92-control cohort are MSI Level 1, and the two sets
coincide exactly — zero correct calls fall outside Level 1.** That is structural,
not unlucky. The controls are Metabolon Level-1 compounds, which is what makes
them a truth set, and `assign_msi_level` returns Level 1 exactly when a Level-1
standard's top candidate matches it. "The pipeline got this right" and "this
report is already confirmed" are the same event on this cohort.

So every AUC computed on it — pooled, per-source, calibrated or raw — measures
how well a score separates the already-confirmed. A `High` tier is for reports
that have no standard, and this cohort contains no example of one being right.

The apparent operating point bears this out. On de novo candidates at raw >= 0.70
the precision is 87.5% (14/16, 95% CI [64.0%, 96.5%]) and the interval excludes
the 56.0% base rate — the only cutoff anywhere in this analysis that does. All 14
of those correct calls are Level 1 compounds, which are already `Confirmed`. The
bucket is not a new confidence tier; it is a rediscovery of the standards.

## What follows

**`High` stays as it is.** An empty tier that means "nothing reached this bar" is
honest. A populated one built on this cohort would mean "these look like the
compounds we already confirmed", which is not what a reader would take it to say.

**If `High` is needed for a client**, the tier that means authentic-standard
confirmation already exists and is populated: `Confirmed`, 44 reports. That is a
labelling decision, not a calibration project.

**The blocker is the truth cohort, not the score.** A tier for unconfirmed
reports needs truth data containing correct calls on unconfirmed reports. The
255-compound delivered truth set has 23 of them, and on that set the raw evidence
score runs backwards (AUC 0.378, see
[docs/msi_sublevel_band_findings.md](msi_sublevel_band_findings.md)). Expanding it
is decision 5 in
[docs/msi_sublevel_decisions_needed.md](msi_sublevel_decisions_needed.md), and the
same budget line as Decision 3 of
[docs/mass_axis_decisions_needed.md](mass_axis_decisions_needed.md).

**When the cohort grows, this reopens automatically.**
`test_the_committed_cohort_cannot_validate_a_tier_for_unconfirmed_reports`
asserts `correct_outside_level_1 == 0`. It is a trigger rather than an invariant:
the day it fails, the truth set has gained the evidence this question needs.

## Reproducing

```
python3 tools/audit_thresholds.py
```

Read-only, standard library only, reproducible from committed inputs. It reads
the structure controls and the per-compound exec summary and rewrites both
artifacts under `docs/paper/data/thresholds/`.

The calibration itself is also reproducible from committed data, which was not
previously true. `docs/paper/data/structure/audit_structure_controls.csv` carries
`raw_full` and `top1_correct_0p70` for all 92 controls, and fitting one against
the other reproduces the shipped `PLATT_A` / `PLATT_B` to five decimal places
(`tests/test_calibrate_confidence.py`). Re-fitting no longer waits on a report
bundle:

```
python3 tools/calibrate_confidence.py <reports-root> --per-source
```
