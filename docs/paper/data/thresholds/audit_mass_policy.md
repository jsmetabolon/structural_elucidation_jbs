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

## Scope

Reordering a fixed list of ten cannot change top-10 recovery, so top-10 is not reported. Bands describe how a measurement is read, not how often that reading is right.

