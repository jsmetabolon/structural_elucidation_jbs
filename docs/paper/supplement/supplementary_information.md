# Supplementary Information

## S1. Internal toolkit-benchmarking design (prospective)

The earlier internal study design compares MS2Query, SIRIUS, Spec2Mol, and MSNovelist with the legacy workflow using CASMI categories: adduct, elemental formula, structure class, and 2D structure. Library 1 contains 2,091 known serum/plasma metabolites spanning 9 super-pathways and 135 subpathways. Library 2 contains 120 internally confirmed NMD/X-compounds plus 500 CASMI 2022 compounds. No comparative toolkit results were available for this manuscript; the design is prospective and its counts are not pooled with the 3,897-feature deployed evaluation.

## S2. Cohorts and confidence vocabularies

| Term | Definition | Source / status |
|---|---|---|
| Delivered set | 3,897 unique X-compounds; first library wins duplicates | report bundle |
| Structure truth set | 92 controls with nonzero `truth_tanimoto` | candidate sidecars |
| Reference-structure cohort | 1,466 rank-1 structures present in reference DB | 1,459 direct lookups |
| `known_X` cohort | 909 controls with independent super-/subpathway labels | taxonomy-discordant |
| P0 | true unknowns with `has_p0_reference=true` | authoritative P0 flags file unavailable |
| P1 | Metabolon Level 1 plus Level 2/2C | authoritative P1 truth directory unavailable |
| Pathway-confidence tier | high / medium / low neighbor-vote confidence | not structure confidence |
| Identification-confidence tier | Confirmed / High / Moderate / Low / Tentative | exploratory Platt score |
| MSI reporting level | evidence-derived L1–L5 rubric | custom; not externally validated |

## S3. Structural recovery and threshold sensitivity

**Table S1. Cumulative structural recovery on n=92 truth controls (Wilson 95% confidence intervals).**

| Truth threshold | Rank 1 | Top 3 | Top 10 |
|---|---:|---:|---:|
| Tanimoto ≥ 0.50 | 59.8% (49.6–69.2%) | 72.8% (63.0–80.9%) | 77.2% (67.6–84.6%) |
| Tanimoto ≥ 0.70 | 47.8% (37.9–57.9%) | 64.1% (53.9–73.2%) | 73.9% (64.1–81.8%) |
| Tanimoto ≥ 0.85 | 40.2% (30.8–50.4%) | 56.5% (46.3–66.2%) | 65.2% (55.0–74.2%) |
| Exact 2D | 38.0% (28.8–48.2%) | 55.4% (45.3–65.2%) | 64.1% (53.9–73.2%) |

**Table S2. Rank-1 recovery by provenance.**

| Source | n | Accuracy (95% CI) | Median truth Tanimoto (bootstrap 95% CI) |
|---|---:|---:|---:|
| exact | 42 | 38.1% (25.0–53.2%) | 0.506 (0.251–0.707) |
| predicted | 50 | 56.0% (42.3–68.8%) | 0.811 (0.565–1.000) |

### S3.1 Cross-validated calibration

Full model leave-one-out Brier score: **0.2546**; calibration slope 0.056, intercept -0.082. Without the pathway axis, Brier score is **0.2528**. These results indicate weak discrimination/calibration and motivate treating the score as exploratory.

![Figure S1. Independent pathway agreement by pathway-confidence tier.](../figures/figS1_pathway_tier.svg)

## S4. Pathway evaluation details

**Table S3. Independent-label agreement by assignment status.**

| Status | Super-pathway (95% CI) | Subpathway (95% CI) |
|---|---:|---:|
| exact_reference_match | 63.0% (57.1–68.6%) | 35.1% (29.6–41.0%) |
| strong_similarity_vote | 54.5% (50.4–58.6%) | 20.1% (17.0–23.7%) |
| weak_similarity_backfill | 36.1% (26.6–46.9%) | 9.6% (5.0–17.9%) |

Leave-one-out uses k=10, similarity floor=0.3, and similarity exponent=2.0. It evaluates the neighbor-vote mechanism only, not the complete production predictor.

![Figure S2. Neighbor-vote leave-one-out accuracy by similarity.](../figures/figS2_loo_similarity.svg)

![Figure S3. Independent-label agreement by assignment status.](../figures/figS3_pathway_status.svg)

**Table S4. Predicted vs independent `known_X` super-pathway distribution (n=909).**

| Super-pathway | Predicted | Independent label |
|---|---:|---:|
| Lipid | 257 | 276 |
| Peptide | 231 | 257 |
| Xenobiotics | 204 | 184 |
| Amino Acid | 153 | 141 |
| Cofactors and Vitamins | 9 | 18 |
| Carbohydrate | 26 | 16 |
| Nucleotide | 21 | 16 |
| Energy | 5 | 1 |
| Partially Characterized Molecules | 3 | 0 |

**Table S5. Independent-label super-pathway confusion for the four most common classes.**

| Independent label | Lipid | Peptide | Xenobiotics | Amino Acid | Other |
|---|---:|---:|---:|---:|---:|
| Lipid | 159 | 59 | 35 | 11 | 5 |
| Peptide | 56 | 135 | 8 | 47 | 6 |
| Xenobiotics | 26 | 10 | 119 | 23 | 4 |
| Amino Acid | 10 | 24 | 28 | 71 | 4 |

## S5. Mass-mismatch audit

There are **598** rank-1 mass mismatches among 3897 reports. **590** are in the MSI-L4 deep-audit subset. A lower-ranked mass-consistent candidate exists for **523/598** overall and **519/590** in the L4 subset. These are projected reranking opportunities, not demonstrated correct structures.

![Figure S4. Rank-1 mass-mismatch cohorts.](../figures/figS4_mass_mismatch.svg)

## S6. Per-library delivered cohort

| Library | Ion mode | n | Percentage |
|---|---|---:|---:|
| 209 | negative | 2,389 | 61.3% |
| 305 | negative | 202 | 5.2% |
| 400 | positive | 875 | 22.5% |
| 402 | positive | 431 | 11.1% |

## S7. Required provenance before submission

- Complete author list, corresponding-author email, CRediT roles, funding, and conflict disclosures.
- Source-study matrix/species metadata and ethics/consent/IRB determination.
- Spec2Mol encoder/decoder checkpoint hashes, training snapshot, random seeds, software/CUDA versions, and EC2 hardware.
- Production neighbor-builder parameters and override rules.
- Authoritative P0 flags and P1 truth files for cohort-specific rank-1/top-10 results.
- Immutable repository release/DOI and controlled-access data instructions.

Full-set identification-confidence median: 0.4322; MSI coverage: 100.0%.
