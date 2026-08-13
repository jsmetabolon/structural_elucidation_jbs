# Mass-policy audit

Read-only comparison of the mass-match window against a proposed one, over a delivered report bundle. Measures reordering only: the pool filter drops candidates before they reach a report, so a tighter tolerance's effect on what was available cannot be seen from here.

## Which formula produced the stored scores

| Rule | Candidates matching |
|---|---:|
| `predict.py (plateau 0.001 Da, tolerance 0.005 Da)` | 56 |
| `1 - err_da` | 22204 |
| `1 - err_da / 5` | 715 |
| `5 / (5 + err_da)` | 1987 |
| `1 / (1 + err_da)` | 167 |
| stored score exactly 0.0 (uninformative) | 5631 |
| unexplained | 885 |

Out of 31645 scored candidates.

## Which generator was scored on which scale

Inverting the score recovers the tolerance each candidate was scored under.

| Candidate source | n | Dominant tolerance | Share |
|---|---:|---|---:|
| `library_anchored_candidates_v1` | 13127 | 1 Da | 95% |
| `peptide_first_denovo_v1` | 9293 | 1 Da | 99% |
| `balanced` | 2584 | 5 Da | 68% |
| `dehydro_v1` | 218 | 0.005 Da | 100% |
| `wgcna_seeded_v1` | 125 | 0.005 Da | 100% |
| `hydroxyl_v1` | 124 | 0.005 Da | 98% |
| `methylation_v1` | 113 | 0.005 Da | 95% |
| `strict_1da` | 99 | other | 100% |
| `protected_peptide_permutation_baseline` | 92 | other | 52% |

Each generator is internally consistent and they disagree with each other, so the mass sub-score is not one quantity across the delivered set. Candidates scored on different scales are ranked against each other and averaged into one confidence number.

## Does the mass measurement carry signal

Whether the axis separates correct candidates from incorrect ones is a different question from whether a gate built on it helps, and the two answers differ.

Correct candidates have a median error of 1.066 ppm against 1.705 ppm for incorrect ones (97 and 676 candidates over the truth controls).

| Window | Correct within | Incorrect within | Separation |
|---:|---:|---:|---:|
| 1 ppm | 48.5% | 37.3% | +11.2% |
| 2 ppm | 75.3% | 53.4% | +21.9% |
| 5 ppm | 81.4% | 59.5% | +21.9% |
| 10 ppm | 81.4% | 61.2% | +20.2% |
| 20 ppm | 82.5% | 62.9% | +19.6% |
| 50 ppm | 82.5% | 69.1% | +13.4% |
| 1000 ppm | 82.5% | 90.1% | -7.6% |
| 3333 ppm | 82.5% | 90.1% | -7.6% |

Separation is widest at 2 ppm. At the shipped 1.0 Da, which is roughly 3,300 ppm at these masses, it inverts: a window that loose admits a higher share of incorrect candidates than correct ones, because the correct ones include a tail that is genuinely far out. So the axis is not merely uninformative at the shipped setting, it is very slightly worse than that.

### The tail, and what it says about the truth definition

17 of 97 candidates scored correct carry a mass error too large for their formula to match, concentrated in 7 controls rather than spread thinly. A Tanimoto of 0.70 to the truth structure does not require the same molecular formula, so a candidate missing a conjugate passes as correct while being wrong by the conjugate's mass.

Grouped by conjugate: {'unmatched': 12, 'glucuronide': 5}.

Worth separating from the rest: `X-21470` carries a candidate structurally identical to its truth annotation and still hundreds of Daltons out. That is a question about the control's annotation rather than about the scorer, and it should not be counted as a conjugate relationship.

## The screen

Candidates between the two windows: **1532** (95 of them at rank 1) out of 31645 with a mass error.

| Mass error | Candidates |
|---|---:|
| <=1 ppm | 8797 |
| <=2 ppm | 5393 |
| <=3 ppm | 1510 |
| <=5 ppm | 1418 |
| <=10 ppm | 1103 |
| <=20 ppm | 1012 |
| <=50 ppm | 1782 |
| <=100 ppm | 1850 |
| <=1000 ppm | 2693 |
| >1000 ppm | 6087 |

Bands are how the confidence model reads a mass error. Rank-1 candidates by band:

| Policy | strong | partial | fail | no mass |
|---|---:|---:|---:|---:|
| `current-absolute-da` | 2436 | 230 | 787 | 444 |
| `drafted-ppm` | 1424 | 985 | 1044 | 444 |

## Does the mass gate improve recovery

The gate from handoff item 2, applied to the whole candidate list under each policy, scored against the truth controls. Paired, because the same controls are scored twice and two independent intervals could not resolve a handful of flips.

| Policy | Depth | Correct before | Correct after | Gained | Lost | Exact p |
|---|---|---:|---:|---:|---:|---:|
| `current-absolute-da` | rank1 | 44 | 41 | 0 | 3 | 0.25 |
| `current-absolute-da` | rank3 | 59 | 58 | 3 | 4 | 1.0 |
| `drafted-ppm` | rank1 | 44 | 39 | 0 | 5 | 0.0625 |
| `drafted-ppm` | rank3 | 59 | 55 | 3 | 7 | 0.34375 |

### Does the tolerance matter at all

The delivered scorer is linear from zero with no plateau. Sweeping its tolerance across three orders of magnitude, against 44 controls correct at rank 1 before any gate:

| Tolerance | Rank-1 correct | Gained | Lost | Exact p | Rank-1 strong | partial | fail |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.005 Da | 41 | 0 | 3 | 0.25 | 2057 | 609 | 787 |
| 0.05 Da | 42 | 0 | 2 | 0.5 | 2756 | 196 | 501 |
| 0.5 Da | 42 | 0 | 2 | 0.5 | 3068 | 30 | 355 |
| 1 Da | 42 | 0 | 2 | 0.5 | 3085 | 20 | 348 |
| 5 Da | 43 | 0 | 1 | 1.0 | 3105 | 18 | 330 |

Recovery moves by 2 controls across the whole sweep, and no setting gains a single one. What does move is what the reports say: rank-1 strong-mass claims range from 2057 to 3105 across the same range. A spread of 2 across a 1000-fold change in tolerance means the mass axis is not separating correct candidates from incorrect ones at any setting, so no choice of window is defensible on this cohort. Note also that losses fall as the tolerance loosens: the gate does least harm where it fires least, which is what an intervention carrying no signal looks like.

### How far that generalises

Only **9 of the 598** rank-1 mass-mismatch reports have a known correct answer, so that handful carries the entire evidence base for any projection about the rest. The gate is aimed at a cohort the truth set barely covers, which bounds what these counts can settle either way.

## Scope

Bands describe how a measurement is read, not how often that reading is right.

Two limits are structural rather than incidental. Reordering a fixed list of ten cannot change top-10 recovery, so no top-10 figure appears above. And the pool filter in `structural_elucidation/utils/merge.py` drops candidates before they reach a report, so delivered data cannot show what a tighter tolerance would have discarded; answering that needs a pipeline re-run, not this bundle.

