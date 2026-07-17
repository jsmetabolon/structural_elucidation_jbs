# MSI sublevel feasibility + defensibility findings

Basis for how sublevels are assigned. Produced by
[tools/audit_sublevel_evidence.py](../tools/audit_sublevel_evidence.py) on the
delivered 3,897-report set joined with the curated `known_X` `SubLevel` (2,331
labeled) and the new `se_report_data_bundle.tgz` evidence.

## 1. The curated scheme already exists (and is the ground truth)
`known_X.SubLevel` carries Metabolon's curated sublevels: `2C` (163), `3A/3B/3C/3D`
(356/253/173/3), `4`/`4?`, `5`, `1-KNOWN`/`Released`. These are expert judgments of
structural specificity encoded in the `Tentative information` free text. They are
**surfaced directly** in the reports (authoritative).

## 2. Our evidence axes do NOT separate the curated sublevels
Per-sublevel medians of the available evidence overlap and are non-monotone:

| SubLevel | n | struct | mass | raw_evid | ms2 | iso | ms2_present |
|---|---:|---:|---:|---:|---:|---:|---:|
| 2C | 163 | 0.659 | 1.0 | 0.639 | 0.121 | 0.5 | 0.67 |
| 3A | 356 | 0.671 | 1.0 | 0.627 | 0.154 | 0.5 | 0.63 |
| 3B | 253 | 0.731 | 1.0 | 0.636 | 0.122 | 0.5 | 0.89 |
| 3C | 173 | 0.618 | 1.0 | 0.609 | 0.140 | 0.5 | 0.78 |
| 4  | 102 | 0.660 | 1.0 | 0.610 | 0.023 | 0.5 | 0.76 |
| 5  | 1281 | 0.700 | 1.0 | 0.613 | 0.176 | 0.5 | 0.64 |

`iso` is a constant 0.5 (default; not populated), `mass` is saturated, `ms2` is
uniformly weak (observed fragment matching was never computed). Structure
similarity is non-monotone across sublevels.

## 3. A derived rubric cannot reproduce or correctly order the fine sublevels
A proposed evidence-based rubric ([tools/audit_sublevel_evidence.py](../tools/audit_sublevel_evidence.py) `derive_sublevel`) reproduces the curated letter only **9.5% exactly** (44% within one order). Its *ordering* only partly tracks truth (overall tau 0.59, driven by the coarse 3-series > 4 > 5 split), and the fine top-of-L3 is **inverted**: derived `2C` median truth 0.452 < derived `3B` 0.731. So the letters 2C/3B/3C/3D cannot be assigned defensibly from current evidence.

## 4. What IS defensible
- **Surface the curated `SubLevel`** wherever `known_X` has it (2,331 compounds) - expert-labeled, authoritative.
- **Coarse MSI level (1-5)** for all reports - already validated and monotone with truth.
- **`candidate_confidence` (`raw_evidence`)** as the continuous include/exclude score - validated to track correctness (AUC 0.79; correct median 0.78 vs incorrect 0.61).
- **`3A` positional-isomer** is the one fine category with a real, semantic signal we can compute (multiple same-formula top candidates), distinct from a raw-evidence bin.

Not defensible from current evidence: deriving `2C`, `3B`, `3C`, `3D`, `2A`, `2B` for uncurated unknowns.

## 5. Evidence required to derive the fine sublevels (upstream handoff)
- **Observed RT** per X-compound (feature tables carry only m/z) -> "Consistent RT".
- A **reference MS² spectral library** -> "Library MS²" (2a, and 2C vs 3-series).
- **Candidate-vs-observed diagnostic-fragment matching** computed from the observed spectra -> 2b/3b typing and homologue evidence.

## 6. Recommended delivery
1. Surface curated `SubLevel` (authoritative) in reports + master index; add a client include/exclude filter on level/sublevel + `candidate_confidence`.
2. For uncurated unknowns, present the validated MSI level + `candidate_confidence`, and the semantic `3A` positional-isomer flag; do **not** fabricate fine 2C/3B/3C/3D letters.
3. Ship the handoff (S5) for the evidence needed to derive the full rubric.
