# Pathway abstention precision-coverage analysis

Cohort: 909 evaluable controls from `docs/paper/data/pathway/audit_pathway_accuracy.csv`.

| policy | n called | coverage | raw acc | harmonized acc |
|---|---|---|---|---|
| all_assignments | 909 | 100.0% | 55.3% | 63.1% |
| drop_low_tier | 546 | 60.1% | 67.6% | 74.0% |
| high_tier_only | 380 | 41.8% | 70.3% | 77.1% |
| exact_status_only | 265 | 29.1% | 63.0% | 72.1% |
| similarity_ge_0.5 | 754 | 83.0% | 60.5% | 68.6% |
| similarity_ge_0.7 | 449 | 49.4% | 67.9% | 76.6% |

Recommended operating point (harmonized accuracy >= 70%, maximize coverage): **drop_low_tier**.
