# Pathway/subpathway prediction audit

- Controls predicted: **3897**
- Status: {'strong_similarity_vote': 1995, 'exact_reference_match': 1459, 'weak_similarity_backfill': 443}
- Tier: {'high': 1863, 'low': 1418, 'medium': 616}
- best_reference_similarity p10/p50/p90: {'p10': 0.375, 'p50': 0.733, 'p90': 1.0}

## Reference-taxonomy consistency (not independent prediction accuracy)
- Exact-structure controls (InChIKey found in reference DB): **1466**
- Overall super/sub consistency: **100.0% / 100.0%**
- exact_reference_lookup: super 100.0% (n=1459), sub 100.0% (n=1459)
- neighbor_consensus: super 100.0% (n=7), sub 100.0% (n=7)
- Direct exact-reference assignments copy the database label and are circular as accuracy tests.
- Caveat: the reference DB and known_X annotations agree only 62.5% (n=267) even on identical structures, so any known_X-based score is capped near that.

## Agreement with known_X annotation (independent labels, not the predictor's taxonomy)
- Evaluable (super): **909** controls
- Super-pathway agreement with known_X: 55.3% (understated: known_X uses a different Metabolon annotation than the reference DB the predictor targets)
- Sub-pathway agreement: 23.5% (given correct super: 42.5%; n=909)

### Super accuracy by confidence tier
| tier | n | accuracy |
|------|---|----------|
| high | 380 | 70.3% |
| medium | 166 | 61.5% |
| low | 363 | 36.9% |

### Super accuracy by best-similarity bin
| bin | n | accuracy |
|-----|---|----------|
| [0.0,0.3) | 39 | 25.6% |
| [0.3,0.5) | 116 | 31.9% |
| [0.5,0.7) | 305 | 49.5% |
| [0.7,1.0) | 449 | 67.9% |

### Top true-super confusions (true -> predicted counts)
- **Lipid**: {'Lipid': 159, 'Peptide': 59, 'Xenobiotics': 35, 'Amino Acid': 11, 'Carbohydrate': 5}
- **Peptide**: {'Peptide': 135, 'Lipid': 56, 'Amino Acid': 47, 'Xenobiotics': 8, 'Nucleotide': 6}
- **Xenobiotics**: {'Xenobiotics': 119, 'Lipid': 26, 'Amino Acid': 23, 'Peptide': 10, 'Carbohydrate': 4}
- **Amino Acid**: {'Amino Acid': 71, 'Xenobiotics': 28, 'Peptide': 24, 'Lipid': 10, 'Carbohydrate': 4}
- **Cofactors and Vitamins**: {'Xenobiotics': 5, 'Energy': 3, 'Cofactors and Vitamins': 3, 'Nucleotide': 2, 'Peptide': 2}
- **Nucleotide**: {'Nucleotide': 7, 'Xenobiotics': 6, 'Lipid': 1, 'Peptide': 1, 'Carbohydrate': 1}
- **Carbohydrate**: {'Carbohydrate': 8, 'Lipid': 3, 'Xenobiotics': 3, 'Nucleotide': 2}
- **Energy**: {'Energy': 1}

### Predicted supers not in truth vocabulary (possible drift)
- ['partially characterized molecules']
