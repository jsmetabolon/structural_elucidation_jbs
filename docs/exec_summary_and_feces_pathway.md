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

### Measured result

Baseline is the 2,265 named metabolites the Global panel detects in feces at all
(bucket other than `Never`, out of 5,189). Added are the 704 X-compounds that
appear in the feces detection history and have a delivered report; one of the
705, `X-25434`, has no report.

| Chemical pathway | Identified | Baseline | +X | Updated | Delta |
|---|---:|---:|---:|---:|---:|
| Lipid | 786 | 34.70% | 182 | 32.60% | −2.10 |
| Xenobiotics | 594 | 26.23% | 191 | 26.44% | +0.21 |
| Amino Acid | 357 | 15.76% | 142 | 16.81% | +1.05 |
| Peptide | 78 | 3.44% | 134 | 7.14% | **+3.70** |
| Secondary metabolism | 118 | 5.21% | 0 | 3.97% | −1.24 \* |
| Nucleotide | 98 | 4.33% | 17 | 3.87% | −0.45 |
| Carbohydrate | 94 | 4.15% | 20 | 3.84% | −0.31 |
| Cofactors and Vitamins | 64 | 2.83% | 11 | 2.53% | −0.30 |
| Partially Characterized Molecules | 37 | 1.63% | 2 | 1.31% | −0.32 |
| Energy | 18 | 0.79% | 5 | 0.77% | −0.02 |
| Cofactors, Prosthetic Groups, Electron Carriers | 17 | 0.75% | 0 | 0.57% | −0.18 \* |
| Hormone metabolism | 4 | 0.18% | 0 | 0.13% | −0.04 \* |

\* The pipeline never predicts this class, so its share can only fall here. The
movement reflects a growing denominator, not anything about feces.

The one substantive shift is **Peptide, roughly doubling from 3.4% to 7.1%**.
Everything else moves by about a point or less, and most of the negative moves
are dilution rather than displacement.

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
alternative; and the sheet spells `Secondary metabolism` two ways, which the
tool collapses so one class does not report as two.

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
