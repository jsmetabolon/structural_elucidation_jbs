# Would mummichog be a better way to assign pathways for X-compounds?

**Verdict: No — not for per-compound pathway assignment.** Mummichog-style
m/z-network inference can only produce a super-pathway call for 12–38% of
evaluable X-compounds, and even on that self-selected (easier) subset it is
less accurate than the production structure-based predictor. The dominant
failure is metabolic-model coverage: only ~12% of curated-known X-compound
masses exist in mummichog's human MFN model at all, and the model is
effectively blind to Xenobiotics and Peptides — the two largest X-compound
classes after Lipids. Mummichog remains a reasonable *complement* for
study-level pathway-activity screening, which is the task it was designed for,
but that is a different deliverable from per-compound Metabolon
super/sub-pathway labels.

## Question and setup

The production pipeline assigns Metabolon super/sub-pathways to X-compounds
from their *predicted structures*
(`structural_elucidation/pipeline/enrichment/pathway.py`): exact
canonical-SMILES lookup in the Metabolon metabolite synonym DB, else a
Morgan-fingerprint neighbor vote. Its audited accuracy
(`docs/paper/data/pathway/audit_pathway_summary.json`) is **55.3%**
super-pathway agreement with independent `known_X` labels (n=909 evaluable;
70.3% for high-confidence tier), against a ~**62.5%** ceiling caused by
taxonomy disagreement between the reference DB and `known_X` annotations.

[Mummichog](https://pypi.org/project/mummichog/) (v2.7.0) instead maps raw
m/z features onto a KEGG-derived metabolic network (human MFN: 119 pathways,
3,557 compounds, 2,124 mass-searchable) using adduct rules and 10 ppm
tolerance, requiring no structure identification. The question: would that be
a better way to label X-compounds with pathways?

## Benchmark design

Implemented in `tools/mummichog_pathway_benchmark.py`; all arms are scored
against the same ground truth (`data/known_X.csv` `SuperPathway`, sentinel
labels excluded) with the same normalization and Wilson CIs as
`tools/audit_pathway_predictions.py`.

- **Arm 1 — baseline (structure-based)**: the committed per-control audit
  (`docs/paper/data/pathway/audit_pathway_accuracy.csv`), not recomputed.
- **Arm 2 — mummichog**: each X-compound's parent m/z (from
  `data/feature_rows/lib{209,305,400,402}_*_feature_rows.csv`, first-lib-wins,
  ionization-mode aware) is matched against the MFN model using mummichog's
  own adduct rules (`mummichog.ions`), wanted-adduct lists and 10 ppm
  tolerance (`mummichog.config`). Candidate compounds vote a Metabolon
  super-pathway through a hand-reviewed 119-row crosswalk
  (`tools/data/mummichog_pathway_crosswalk.csv`; 25 rows flagged ambiguous),
  each candidate contributing one vote split across its mapped supers.
  Two variants:
  - `primary`: matches restricted to mummichog's primary ions
    (M+H/M+Na/M-H2O+H; M-H/M-2H/M-H2O-H) — the faithful analog of its
    `force_primary_ion` default given one quantified parent ion per compound.
  - `all`: every adduct in the wanted list for the mode (permissive).
- **Arm 3 — mass-only ablation (gated, skipped)**: the same vote against the
  Metabolon synonym DB itself, isolating "mass-only vs structure" from "KEGG
  coverage". The SQLite asset is not in this checkout (`fetch_assets.sh`
  needs EC2/local access); the subcommand is implemented and will run once
  the asset is fetched.

Data coverage: 2,473 X-compound controls have a parent m/z; **881** of the
909 baseline-evaluable controls have both m/z and a real truth label, and all
881 are in the baseline audit — so head-to-head comparisons are exact on
shared compounds.

## Results

| Metric | Baseline (structure) | Arm 2 `primary` | Arm 2 `all` |
|---|---|---|---|
| Coverage of evaluable set | 100% (909/909) | **12.4%** (109/881) | **38.0%** (335/881) |
| Super accuracy on called | 55.3% (n=909) | 53.2% (n=109, CI 43.9–62.3) | 37.9% (n=335, CI 32.9–43.2) |
| Baseline accuracy on the same called subset | — | 63.3% | 58.2% |
| Discordant pairs (mummichog-only vs baseline-only correct) | — | 17 vs 28 | 38 vs 106 |
| McNemar p (two-sided) | — | 0.135 | 1.27e-08 |
| Correct calls as share of all 881 evaluable | 55.3% | 6.6% (58) | 14.4% (127) |

Notes:

- The compounds mummichog *can* call skew canonical/endogenous, which is
  exactly where the baseline is also strongest (63.3% on the primary-variant
  subset vs 55.3% overall) — mummichog loses the head-to-head even on its own
  best terrain, significantly so in the higher-coverage variant.
- Raising the vote-share threshold does not rescue accuracy: the best
  operating point of the `primary` variant is 59.3% at ≥0.8 vote share with
  9.8% coverage — still below the baseline on the same rows.

### Per-super accuracy (Arm 2 `all`, n=335)

| True super | n | Mummichog accuracy |
|---|---|---|
| Lipid | 78 | 78.2% |
| Amino Acid | 74 | 71.6% |
| Carbohydrate | 8 | 50.0% |
| Nucleotide | 11 | 18.2% |
| Cofactors and Vitamins | 11 | 9.1% |
| Xenobiotics | 86 | **5.8%** |
| Peptide | 66 | **0.0%** |

The evaluable X-compound population is dominated by Lipid (276), Peptide
(257), Xenobiotics (184), and Amino Acid (141) — i.e. roughly half of all
X-compounds belong to classes (Peptide, Xenobiotics) where a KEGG-style
metabolic model has essentially no vocabulary. True peptides are
systematically mis-called as Lipid (52/66): their masses collide with
glycerophospholipid/fatty-acyl model compounds.

### Why it fails: model coverage, not matcher wiring

For the 78 evaluable controls with curated known structures, only **9 (11.5%)**
have their true neutral mass anywhere in the MFN model at ±10 ppm. When the
true mass *is* in the model, the feature matcher finds it **88.9%** of the
time — the adduct/tolerance implementation is validated; the limitation is
that X-compounds are, almost by definition, molecules outside canonical
metabolic models (uncatalogued xenobiotics, food components, conjugates,
partial peptides, modified lipids).

A secondary failure mode is mass-only ambiguity: a single m/z at 10 ppm
matches a median of 2 (up to 26) model compounds spanning multiple supers,
whereas the structure-based predictor votes over MS2-verified candidate
structures.

## What mummichog would still be good for

Mummichog's actual design target is **study-level pathway activity**: given a
full feature table with per-feature statistics, it tests whether significant
features concentrate in pathways, tolerating per-feature ambiguity because it
aggregates over many features. That use is complementary to (not a
replacement for) per-compound labels, and could be offered as a study-level
screen on `untargeted_workflow` exports (`feature_meta_analysis.csv` +
statistics). It was out of scope here because neither repo contains a study
feature table with statistics (requires an Oracle-backed export).

## Recommendation

1. **Keep the structure-based predictor** as the mechanism for per-X-compound
   super/sub-pathway assignment. It has full coverage, higher accuracy on
   every shared subset, and produces sub-pathways, which mummichog cannot
   (KEGG pathway granularity does not map onto Metabolon sub-pathways).
2. **Do not adopt mummichog for per-compound assignment.** Its usable signal
   (Amino Acid / Lipid calls at ~72–78%) is already covered better by the
   baseline on the same compounds.
3. Optional follow-ups, in decreasing order of value:
   - Run **Arm 3** (mass-only vote against the Metabolon synonym DB) once the
     SQLite asset is fetched, to quantify how much of the baseline's accuracy
     is attributable to structure information vs mass alone within a
     Metabolon-native reference. If mass-only against the in-house DB scores
     close to the baseline, cheap mass votes could pre-fill pathways before
     structure prediction completes.
   - Offer study-level mummichog enrichment as a separate deliverable on
     pipeline exports (pathway activity per study, not per-compound labels).

## Reproduction

```bash
pip install mummichog==2.7.0 pandas rdkit
python tools/mummichog_pathway_benchmark.py inventory
python tools/mummichog_pathway_benchmark.py arm2
python tools/mummichog_pathway_benchmark.py arm3   # skips unless synonym DB fetched
python tools/mummichog_pathway_benchmark.py score
```

Artifacts:

- `docs/paper/data/pathway/mummichog_arm2_predictions.csv` — per-compound
  candidates, votes, and known-structure diagnostics (both variants).
- `docs/paper/data/pathway/mummichog_eval_summary.json` / `.md` — all metrics
  reported above (coverage, accuracy, per-super, confusion, vote-share
  curves, head-to-head, McNemar).
- `tools/data/mummichog_pathway_crosswalk.csv` — reviewed MFN→Metabolon
  super crosswalk (`ambiguous`/`notes` columns mark judgment calls).
