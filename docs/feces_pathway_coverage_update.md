# Chemical-pathway coverage in feces, with X-compounds included

The panel's picture of what feces contains is built from metabolites it can
name. Every X-compound is a reproducibly detected feature it cannot, so all of
them are missing from that picture by construction. This is what the composition
looks like once the pipeline's pathway assignments for those features are added.

Regenerate with:

```bash
python3 tools/build_feces_pathway_update.py \
    --global-workbook "20260130 Microbiome vs. Global Detection (9).xlsx" \
    --detection-history "20250728 hFeces X-compound (Unknown) Detection History.xlsx" \
    --reports-root /path/to/reports --stratify
```

## What counts as detected

Detection is read from the **GDP** column (`GDP Bucket - hFecal`, Global
Discovery Panel) of the `GLOBAL 20260129` sheet, not the MFP columns beside it,
which describe the narrower microbiome-focused panel. A metabolite counts as
detected when that bucket is anything other than `Never`: 2,265 of the 5,189 on
the panel.

On the X-compound side, 705 appear in the human-feces detection history, with
detection rates from 1.75% to 98.25%. 704 of those have a delivered report and
so carry a pathway assignment; `X-25434` does not. Pathway is the rank-1
candidate's `predicted_pathway`.

## Coverage by pathway

| Chemical pathway | Identified | Previous | X added | Updated | Change |
|---|---:|---:|---:|---:|---:|
| Lipid | 786 | 34.70% | 182 | 32.60% | −2.10 |
| Xenobiotics | 594 | 26.23% | 191 | 26.44% | +0.21 |
| Amino Acid | 357 | 15.76% | 142 | 16.81% | +1.05 |
| Peptide | 78 | 3.44% | 134 | 7.14% | **+3.70** |
| Secondary metabolism † | 118 | 5.21% | 0 | 3.97% | −1.24 |
| Nucleotide | 98 | 4.33% | 17 | 3.87% | −0.45 |
| Carbohydrate | 94 | 4.15% | 20 | 3.84% | −0.31 |
| Cofactors and Vitamins | 64 | 2.83% | 11 | 2.53% | −0.30 |
| Partially Characterized Molecules | 37 | 1.63% | 2 | 1.31% | −0.32 |
| Energy | 18 | 0.79% | 5 | 0.77% | −0.02 |
| Cofactors, Prosthetic Groups, Electron Carriers † | 17 | 0.75% | 0 | 0.57% | −0.18 |
| Hormone metabolism † | 4 | 0.18% | 0 | 0.13% | −0.04 |
| **Total** | **2,265** | | **704** | **2,969** | |

† The pipeline never predicts this class, so its share falls only because the
denominator grew. Read those three rows as arithmetic, not as a finding about
feces.

## Reading it

**Peptide is the result.** It roughly doubles, from 3.44% to 7.14%, and it is
the only movement larger than the noise the rest of the table sits in. 134 of
the 704 additions are peptides against a named baseline of only 78. That is
consistent with the panel being weakest exactly where the pipeline's peptide
track is strongest, and it survives every cut below.

**Lipid falls without losing anything.** It drops 2.10 points on 182 additions,
which is the largest absolute gain of any class. Its share falls because other
classes grew faster. Nothing about the lipid inventory changed.

**Everything else moves less than about one point**, which is well inside what
the pathway predictor's own error can account for.

## How much to trust it

These are predicted pathways. Super-pathway agreement against independent labels
was measured at 52.2% against a ceiling of 62.5%
(`docs/paper/data/pathway/audit_pathway_summary.json`), so roughly half of the
704 assignments are expected to be wrong at the class level. The table shows
direction and rough size; it is not a measurement of fecal composition.

Three checks make the Peptide result more robust than that caveat alone implies:

- **It holds at every confidence tier.** Splitting the 704 by the pipeline's own
  pathway-confidence tier moves Peptide +1.69 points on the 266 high-tier
  assignments, +1.33 on 144 medium and +1.36 on 294 low. If it were an artefact
  of weak predictions, the high tier would not show it most strongly.
- **It holds among frequently detected compounds.** Restricting both sides to
  things seen in at least half of sample sets, Peptide goes 4.28% to 8.92%.
- **It is not double counting.** The Global panel contains no X-compounds by
  name. Of the 704, 105 are separately identified at Metabolon Level 1 or 2C,
  but only one of those (`X-24643`, annotated `sulfate`) matches a named Global
  metabolite that is also detected in feces.

## What is not settled

- **Whether a predicted class should count the same as a measured one.** The
  table gives each compound one vote regardless of provenance. Weighting the
  predictions down, or reporting them as a separate stratum rather than merged,
  is a presentation decision nobody has made yet.
  `--weight-by-detection` reports the frequency-weighted alternative.
- **The Peptide baseline may be understated for reasons unrelated to this work.**
  78 named peptides out of 2,265 is low for feces, and part of the doubling may
  be the panel's coverage rather than the pipeline's discovery.
