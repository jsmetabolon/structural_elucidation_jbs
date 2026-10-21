# Two exports: the top-ten summary, and feces pathway composition

Both tools are read-only over a delivered report bundle. Neither changes any
published number. The generated CSVs are **not committed**, because both depend
on inputs that do not live in this repository; the commands to produce them are
below.

## 1. Top ten candidates per X-compound

`tools/build_exec_summary.py` emits, for every X-compound, its top ten proposed
structures with scores and an indication of which of the ten matched the truth.

```bash
python3 tools/build_exec_summary.py \
  --reports-root /path/to/reports \
  --out-dir docs/paper/data/exec_summary
```

Outputs `exec_summary_top10.csv` (one row per candidate),
`exec_summary_by_compound.csv` (one row per compound) and `exec_summary.json`.

### What "matched the truth" can mean, per cohort

The request was for a match indication covering known unknowns and Level 2C
identified metabolites. Those two cohorts do not have the same kind of answer
behind them, and conflating them would overstate the second.

| Cohort | N | Reference available | What the match column means |
|---|---:|---|---|
| `truth_control` | 92 | Blinded reference structure, per-candidate `truth_tanimoto` | Structural match at Tanimoto >= 0.70, the published definition |
| `level_2c` | 163 | Tentative **name**; a formula for only 26 | Formula agreement, where a formula exists; blank for the other 137 |
| `other` | 3,642 | None | Blank, never `FALSE` |

The Level 2C limitation is worth stating plainly, because it is a property of
the data rather than of the tool:

- **No 2C compound has a reference structure.** Structural truth exists for
  exactly 92 compounds, and all of them are essentially Metabolon Level 1
  (65 `KNOWN`, 25 `Released`). Zero 2C compounds carry a `truth_tanimoto`.
- **No 2C compound has a single named candidate in its top ten**, so the
  tentative name cannot be matched against a candidate name.
- **Only 26 of the 163 record a molecular formula**, so formula agreement is
  all that can be checked, and only for those 26.
- The annotation itself is tentative. Where the delivered report and `known_X`
  both name a 2C compound they sometimes disagree: `X-26082` is
  `N-acetylglycine` in the report and `N-formylalanine` in `known_X`, two
  isomers of `C4H7NO3`.

A `truth_basis` column records which of these applies on every row, and
unassessable compounds are left blank rather than false so they never enter a
match rate.

### Measured on the delivered bundle

Run against the 3,897-report bundle, the export reproduced the published
recovery figures exactly, which is the check that it is reading the same truth
the manuscript does:

- 3,897 compounds, 34,296 candidate rows
- Truth controls: 44/92 at rank 1 (**47.8%**), **64.1%** top 3, **73.9%** top 10
- Level 2C: 23 of the 26 assessable compounds have a top-ten candidate whose
  formula agrees with the annotation
- 3,642 compounds not assessable on either basis

### Feature names and the workbook

`feature_name` is the name the report's naming policy allows. For Naming v2
bundles (`output/reports_ppm10_full_naming_v2`) that is: a bare name at Level 1
(authentic standard), `Probable <name>` at Level 2, `Possible <name>` at Level 3
(curated `known_X` hypothesis), `<formula>; structure unresolved` at Level 4 and
`Exact mass <m>; unidentified` at Level 5. One sheet-level rule sits on top: a
compound whose rank-1 structure is already a catalogued global-library
metabolite, and which the report therefore left at a formula or exact-mass
placeholder, is written `Isomer of <library metabolite>` (an X-compound cannot
*be* a metabolite the library already names, so the relationship is the
informative label). The report's own label is kept in `report_display_name`,
and every excluded library metabolite is listed in `excluded_library_reference`.
Names from a compound's own `known_X` row and authentic-standard identities are
never re-qualified.

`tools/build_exec_summary_workbook.py --exec-dir <out-dir>` turns the CSVs into
`executive_summary_top10.xlsx` (top-10 sheet, by-compound sheet, a
`naming_conventions` legend and a `summary` sheet). `scripts/generate_ppm10_full_reports.sh`
runs it after the CSV export. The reviewed copy lives at
`docs/paper/data/exec_summary/executive_summary_top10.xlsx`.

## 2. Updated feces pathway composition

`tools/build_feces_pathway_update.py` recomputes the chemical-pathway breakdown
of human feces with X-compound pathway assignments added to the named
metabolites the panel already detects.

```bash
python3 tools/build_feces_pathway_update.py \
  --global-workbook "20260130 Microbiome vs. Global Detection (9).xlsx" \
  --detection-history "20250728 hFeces X-compound (Unknown) Detection History.xlsx" \
  --reports-root /path/to/reports \
  --stratify
```

Detection is read from **GDP** (Global Discovery Panel), column
`GDP Bucket - hFecal` on sheet `GLOBAL 20260129`, not from the narrower MFP
columns beside it.

### Known pathways where they exist, predicted only where they do not

An X-compound's pathway is taken from the curated `SuperPathway` in
`data/known_X.csv` when one is recorded there, and from the pipeline's
prediction otherwise. Of the 704 feces-detected compounds, **344 are curated and
360 predicted**. Where both exist they agree 54.1% of the time, close to the
52.2% super-pathway agreement measured elsewhere, so preferring the curated
label changes a meaningful number of rows.

A fifth of that column is not taxonomy. `Unknown` (155 rows), `DO NOT INCLUDE`
(17), `Ask DTS` (11) and `Can add this` (12) are review notes, the same defect
the manuscript already records for the subpathway labels. A curated value is
only accepted if it appears in the Global panel's own vocabulary; anything else
falls back to the prediction. Without that filter the results table grows
columns called `Unknown` and `Ask DTS` and reports them as chemistry.

Pass `--no-curated-pathways` to predict everything and see the difference.

### Which "known" set to compare against

Two defensible denominators, and they give different pictures:

- **Everything the panel detects in feces** — 2,265 of the 5,189 metabolites on
  the Global panel. This is the full chemical inventory feces is known to
  contain. Written to `docs/paper/data/feces_pathway/`.
- **The curated fecal panel** — the 944 metabolites on the Fecal Inclusion List,
  of which 913 are GDP-detected. This is what the microbiome panel is actually
  meant to carry. Pass `--inclusion-list`; written to
  `docs/paper/data/feces_pathway_panel/`.

The inclusion list is far more xenobiotic- and peptide-poor than the full panel,
so the X-compounds move it much more. All 944 join to Global by `CHEMICAL_ID`
with a pathway, so nothing is lost in the join.

### Two ways to count, and they differ

Baseline is the 2,265 named metabolites the Global panel detects in feces at all
(bucket other than `Never`, out of 5,189). Added are the 704 X-compounds that
appear in the feces detection history and have a delivered report with a pathway
annotation; one of the 705, `X-25434`, has no report. Reported averages always
include both sides.

How much each compound contributes is a real choice, not a detail:

- **By compound.** Every compound counts once, whether it turns up in every
  sample or one. This describes the *inventory*: how many distinct things of
  each class the matrix can contain.
- **By detection frequency** (`--weight-by-detection`). Each compound
  contributes the fraction of studies it is seen in. For a pathway the reported
  average is the **expected number of metabolites detected**:
  `Σ detection_rate` over knowns and pathway-annotated X-compounds in that
  pathway. Five metabolites each seen 20% of the time contribute 1.0. That sum
  is the average. Weighted exports do **not** carry an `updated_pct` share
  column; compositional shares remain only on the by-compound counting mode.

The two sources use different denominators — Global fecal detection is over 82
SSIDs and the X-compound history over 57 sets — so rates are comparable between
them while raw detection counts are not. The weighting uses rates for that
reason.

Both are written out: `feces_pathway_update.csv` /
`feces_pathway_update_subpathway.csv` and
`feces_pathway_update_weighted.csv` /
`feces_pathway_update_weighted_subpathway.csv`. The weighted CSVs carry
`avg_metabolites_detected_known`, `avg_metabolites_detected_x`, and
`avg_metabolites_detected` (known + X) at both super-pathway and sub-pathway
level, each rounded to the nearest integer. `feces_pathway_update_weighted.xlsx`
puts the same results in one workbook, with `Pathway` and `Subpathway` tabs.

### All four views (compositional shares)

Previous to updated, in percentage points, across both denominators and both
counting rules. These are shares of detection mass (or of compound counts), not
the average metabolite detections:

| Chemical pathway | Panel, by compound | Panel, weighted | Global, by compound | Global, weighted |
|---|---:|---:|---:|---:|
| Amino Acid | 22.23 → 21.95 | 29.50 → 27.89 | 15.76 → 17.14 | 22.10 → 22.96 |
| Lipid | 19.50 → 20.78 | 23.40 → 22.08 | 34.70 → 31.80 | 35.40 → 31.60 |
| Xenobiotics | 33.73 → 30.06 | 21.34 → 20.91 | 26.23 → 26.00 | 16.16 → 17.19 |
| **Peptide** | **2.41 → 11.26** | **1.71 → 11.35** | **3.44 → 8.02** | **3.85 → 9.22** |
| Secondary metabolism \* | 10.51 → 5.94 | 9.95 → 5.90 | 5.21 → 3.97 | 5.71 → 4.29 |
| Nucleotide | 2.41 → 2.72 | 4.11 → 3.85 | 4.33 → 4.04 | 5.85 → 5.26 |
| Cofactors and Vitamins | 2.63 → 2.60 | 3.36 → 3.35 | 2.83 → 2.76 | 3.67 → 3.58 |
| Carbohydrate | 4.60 → 3.40 | 4.13 → 2.93 | 4.15 → 3.60 | 3.92 → 3.24 |
| Energy | 0.88 → 0.68 | 1.45 → 1.13 | 0.79 → 0.71 | 1.06 → 0.96 |
| Cofactors, Prosthetic Groups \* | 0.55 → 0.31 | 0.64 → 0.38 | 0.75 → 0.57 | 0.92 → 0.69 |
| Hormone metabolism \* | 0.44 → 0.25 | 0.40 → 0.24 | 0.18 → 0.13 | 0.19 → 0.15 |
| Partially Characterized Molecules | 0.11 → 0.06 | 0.00 → 0.00 | 1.63 → 1.25 | 1.16 → 0.87 |

**Peptide rises in all four**, between 4.6 and 9.6 points, and it is the only
class that does. Every other pathway moves by less than about 4.6 points and
most by under 1.5. The choice of denominator and counting rule changes the size
of the peptide result substantially — it is largest against the curated panel,
which carries only 22 named peptides — but not its direction or its rank.

Nothing else here is stable enough to report as a finding. Xenobiotics falls 3.7
points against the panel and 0.2 against Global; Amino Acid rises against Global
and falls against the panel. Those are artefacts of which denominator is chosen,
not statements about feces.

### Measured result, by detection frequency

The frequency-weighted view answers "how many metabolites of each class are
detected in a typical sample": each cell is the sum of detection rates for
knowns and pathway-annotated X-compounds in that pathway (Global denominator):

| Chemical pathway | Known avg | +X avg | Avg metabolites detected |
|---|---:|---:|---:|
| Lipid | 385 | 73 | 458 |
| Amino Acid | 240 | 92 | 333 |
| Xenobiotics | 176 | 73 | 249 |
| Peptide | 42 | 92 | **134** |
| Nucleotide | 64 | 13 | 76 |
| Secondary metabolism | 62 | 0 | 62 \* |
| Cofactors and Vitamins | 40 | 12 | 52 |
| Carbohydrate | 43 | 4 | 47 |
| Energy | 12 | 2 | 14 |
| Partially Characterized Molecules | 13 | 0 | 13 \* |
| Cofactors, Prosthetic Groups, Electron Carriers | 10 | 0 | 10 \* |
| Hormone metabolism | 2 | 0 | 2 \* |

Peptide's average rises from 42 known detections to 134 once X-compounds
are included — more than a threefold increase in expected detections. Weighting
changes the **baseline share under by-compound counting** more than it changes
the absolute averages: xenobiotics fall from 26.23% to 16.16% of compound-count
share and amino acids rise from 15.76% to 22.10% before a single X-compound is
added when that inventory is re-expressed as detection rates, because the panel
lists many xenobiotics that each turn up rarely and fewer amino acids that turn
up almost always. Counting compounds and summing detection rates answer different
questions; the xenobiotic-heavy inventory is the one that does not survive the
switch.

The Peptide finding does not depend on which counting rule is chosen: its
absolute average rises under weighting, and its by-compound share rises +4.57
points.

### Measured result, by detection frequency (sub-pathway)

The same sum-of-rates average, at sub-pathway resolution. Written to
the `Subpathway` tab of `feces_pathway_update_weighted.xlsx` and to
`feces_pathway_update_weighted_subpathway.csv` (223 sub-pathways on the Global
denominator). Top fifteen by expected metabolites detected:

| Sub-pathway | Known avg | +X avg | Avg metabolites detected |
|---|---:|---:|---:|
| Food Component/Plant | 81 | 37 | 118 |
| Dipeptide | 27 | 29 | 57 |
| Polypeptide | 3 | 46 | 49 |
| Chemical | 13 | 15 | 28 |
| Tryptophan Metabolism | 14 | 13 | 27 |
| Fatty Acid, Dicarboxylate | 13 | 13 | 26 |
| Secondary Bile Acid Metabolism | 20 | 5 | 24 |
| Steroid | 23 | 0 | 23 \* |
| Benzoate Metabolism | 11 | 11 | 22 |
| Bile acids | 21 | 0 | 21 \* |
| Glutamate, Arginine, and Proline Metabolism | 21 | 0 | 21 \* |
| Lysine Metabolism | 16 | 4 | 20 |
| Fatty acid, hydroxy | 19 | 0 | 19 \* |
| Fatty Acid, Monohydroxy | 9 | 8 | 17 |
| Aspartate family (OAA derived) | 16 | 0 | 16 \* |

The Peptide super-pathway rise shows up here as **Dipeptide** (27 → 57)
and especially **Polypeptide** (3 → 49): almost all of the polypeptide
average is from pathway-annotated X-compounds.

### Measured result, by compound

| Chemical pathway | Known | Baseline | +X | Updated | Delta |
|---|---:|---:|---:|---:|---:|
| Lipid | 786 | 34.70% | 158 | 31.80% | −2.91 |
| Xenobiotics | 594 | 26.23% | 178 | 26.00% | −0.22 |
| Amino Acid | 357 | 15.76% | 152 | 17.14% | +1.38 |
| Peptide | 78 | 3.44% | 160 | 8.02% | **+4.57** |
| Nucleotide | 98 | 4.33% | 22 | 4.04% | −0.28 |
| Secondary metabolism | 118 | 5.21% | 0 | 3.97% | −1.24 \* |
| Carbohydrate | 94 | 4.15% | 13 | 3.60% | −0.55 |
| Cofactors and Vitamins | 64 | 2.83% | 18 | 2.76% | −0.06 |
| Partially Characterized Molecules | 37 | 1.63% | 0 | 1.25% | −0.39 \* |
| Energy | 18 | 0.79% | 3 | 0.71% | −0.09 |
| Cofactors, Prosthetic Groups, Electron Carriers | 17 | 0.75% | 0 | 0.57% | −0.18 \* |
| Hormone metabolism | 4 | 0.18% | 0 | 0.13% | −0.04 \* |

\* The pipeline never predicts this class, so its share can only fall here. The
movement reflects a growing denominator, not anything about feces.

The one substantive shift is **Peptide, from 3.4% to 8.0%**.
Everything else moves by about a point or less, and most of the negative moves
are dilution rather than displacement.

Note that shares are compositional and must total 100%, so adding 704 compounds
to 2,265 grows the denominator by 31% and any class that does not gain at least
that proportionally shows a smaller share while losing nothing. Lipid is the
clearest case: it takes 182 additions, the second most of any class, and still
falls 2.10 points because Peptide grew faster.

### How much to trust it

These are predicted pathways. Measured super-pathway agreement against
independent labels is 52.2% against a ceiling of 62.5%, so roughly half the
added assignments are expected to be wrong at the class level. That is why
`--stratify` exists.

The Peptide finding survives that test: it is present at a similar size in every
confidence tier, rather than resting on the weak ones.

| Tier | N | Peptide delta |
|---|---:|---:|
| high | 266 | +1.69 |
| medium | 144 | +1.33 |
| low | 294 | +1.36 |

Restricting **both** sides to compounds detected in at least half of samples
(the panel's `Often` and `Very Often` buckets, and X-compounds at >= 50%)
sharpens the same picture: Peptide 4.28% to 8.92%, Lipid 36.31% to 33.43%.
Applying that filter to only one side inflates the shift, so the tool applies it
to both.

Two further caveats hold regardless of tier. Compounds are counted once each
regardless of detection frequency, with `--weight-by-detection` as the
alternative that reports expected metabolites detected (`Σ` detection rates for
knowns plus pathway-annotated X-compounds); and the sheet spells
`Secondary metabolism` two ways, which the tool collapses so one class does not
report as two.

Double counting was checked and is negligible. The Global panel contains no
X-compounds by name, and of the 704 added only one, `X-24643` annotated as
`sulfate`, collides by name with a feces-detected named metabolite.

## Regenerating after a fresh VM

Neither export survives a rebuilt machine, and neither set of inputs is in the
repository:

- The report bundle. Re-export per
  [docs/mass_policy_bundle_handoff.md](mass_policy_bundle_handoff.md).
- `20260130 Microbiome vs. Global Detection (9).xlsx` and
  `20250728 hFeces X-compound (Unknown) Detection History.xlsx`, which are
  uploads rather than repository files.

`data/known_X.csv` and `data/known_X_all_available_structures_simple.csv` are in
the repository, so cohort assignment needs no external input.
