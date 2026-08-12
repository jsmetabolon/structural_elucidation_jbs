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

### How far that generalises

Only **9 of the 598** rank-1 mass-mismatch reports have a known correct answer, so that handful carries the entire evidence base for any projection about the rest. The gate is aimed at a cohort the truth set barely covers, which bounds what these counts can settle either way.

## Scope

Bands describe how a measurement is read, not how often that reading is right.

Two limits are structural rather than incidental. Reordering a fixed list of ten cannot change top-10 recovery, so no top-10 figure appears above. And the pool filter in `structural_elucidation/utils/merge.py` drops candidates before they reach a report, so delivered data cannot show what a tighter tolerance would have discarded; answering that needs a pipeline re-run, not this bundle.

