# Mummichog X-compound pathway benchmark

- Evaluable controls (real known_X super label) with feature m/z: **881**
- Baseline (structure-based) audit controls: 909 (with m/z: 881)

## Arm 2 (mummichog MFN, variant = primary)

- Called: 109/881 evaluable (coverage 12.4%)
- Super accuracy on called: **0.5321** (n=109, 95% CI [0.4389, 0.6231])
- Head-to-head (n=109): mummichog 53.2% vs baseline 63.3%; discordant 17 vs 28; McNemar p = 0.135
- Known-structure diagnostics: true mass in MFN model for 9/78 (11.5%); when in model, matched by the feature 88.9% of the time

| vote share >= | n | coverage | accuracy |
|---|---|---|---|
| 0.0 | 109 | 12.4% | 53.2% |
| 0.2 | 109 | 12.4% | 53.2% |
| 0.4 | 104 | 11.8% | 54.8% |
| 0.5 | 104 | 11.8% | 54.8% |
| 0.6 | 95 | 10.8% | 56.8% |
| 0.8 | 86 | 9.8% | 59.3% |
| 0.9 | 85 | 9.7% | 58.8% |

## Arm 2 (mummichog MFN, variant = all)

- Called: 335/881 evaluable (coverage 38.0%)
- Super accuracy on called: **0.3791** (n=335, 95% CI [0.3288, 0.4322])
- Head-to-head (n=335): mummichog 37.9% vs baseline 58.2%; discordant 38 vs 106; McNemar p = 1.27e-08
- Known-structure diagnostics: true mass in MFN model for 9/78 (11.5%); when in model, matched by the feature 88.9% of the time

| vote share >= | n | coverage | accuracy |
|---|---|---|---|
| 0.0 | 335 | 38.0% | 37.9% |
| 0.2 | 335 | 38.0% | 37.9% |
| 0.4 | 324 | 36.8% | 38.6% |
| 0.5 | 318 | 36.1% | 38.0% |
| 0.6 | 295 | 33.5% | 39.0% |
| 0.8 | 253 | 28.7% | 36.8% |
| 0.9 | 241 | 27.4% | 38.2% |

## Arm 3 (mass-only vote vs Metabolon reference DB)

- Skipped: metabolite synonym SQLite asset not available in this checkout

