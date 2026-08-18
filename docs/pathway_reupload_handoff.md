# Handoff: the three small files needed to redo the feces pathway sums

Self-contained task for whoever can reach the host holding the rendered
X-compound reports. Assume no other context. You are only exporting data.

**Do not send the report bundle.** An earlier version of this work moved all
11,691 sidecar files, roughly 190 MB, over a presigned S3 link. That is far more
than this needs. Of those files, the pathway work reads exactly four columns
from the first row of one file per report. Extracting them locally turns a
190 MB transfer into about 250 KB, and removes S3 from the process entirely.

## Goal

Let the cloud agent recompute the chemical-pathway composition of human feces on
an expected-count basis: each metabolite contributes its detection frequency
rather than a whole unit, so five metabolites each seen 20% of the time sum to
1.0 rather than 5. Both named metabolites and X-compounds carrying our pathway
annotations must be included.

## What is needed

Three files, about 7.5 MB in total. All go through the normal chat upload; none
needs S3, a tarball, or a presigned URL.

| # | File | Approx size | Where it comes from |
|---|---|---:|---|
| 1 | `x_compound_pathways.csv` | ~250 KB | extracted from the report host, script below |
| 2 | `20260130 Microbiome vs. Global Detection (9).xlsx` | 7.2 MB | already exists, re-upload as is |
| 3 | `20250728 hFeces X-compound (Unknown) Detection History.xlsx` | 34 KB | already exists, re-upload as is |

Items 2 and 3 were uploaded before and were lost when the cloud VM was rebuilt.
They are unchanged; no regeneration is needed, just send the same files.

### Why each one

- **(1)** supplies each X-compound's predicted pathway. Nothing else in the
  repository has it for all 3,897 compounds; the committed
  `docs/paper/data/pathway/audit_pathway_accuracy.csv` covers only the 909 with
  independent labels.
- **(2)** supplies the named metabolites: their pathway and their fecal
  detection rate, read from sheet `GLOBAL 20260129`, columns
  `Orig SuperPathway`, `Orig SubPathway` and `hFecal HD4` (GDP, not MFP).
- **(3)** supplies each X-compound's fecal detection rate, from sheet
  `unk_feces_detection_rates`.

Nothing else is required. The top-ten candidate export is already committed at
`docs/paper/data/exec_summary/` for all 3,897 compounds and needs no new data.

## Step 1: extract the pathway table

Run on the host holding the reports. Pure standard library, no dependencies.

```bash
REPORTS="/path/to/reports"   # the directory containing X-01911/, X-12410/, ...
OUT="/tmp/x_compound_pathways.csv"

python3 - "$REPORTS" "$OUT" <<'PY'
import csv, sys, pathlib
csv.field_size_limit(10**9)
root, out = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2])
cols = ["predicted_pathway", "predicted_subpathway",
        "pathway_confidence_tier", "pathway_assignment_status"]
rows, missing = [], 0
for d in sorted(p for p in root.iterdir() if p.is_dir() and p.name.upper().startswith("X-")):
    f = d / "data" / f"{d.name.lower()}_top_candidates.csv"
    if not f.exists():
        missing += 1
        continue
    with f.open(newline="") as fh:
        first = next(csv.DictReader(fh), None)
    if first is None:
        missing += 1
        continue
    rows.append({"x_compound_id": d.name,
                 **{c: (first.get(c) or "").strip() for c in cols}})
with out.open("w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["x_compound_id"] + cols, lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
print(f"wrote {len(rows)} rows to {out}; {missing} reports skipped")
PY

ls -lh "$OUT" && head -3 "$OUT"
```

The first row of `*_top_candidates.csv` is the rank-1 candidate, and the pathway
columns are per report rather than per candidate, so only that row is read.

## Step 2: check it before sending

- The row count is **3,897**. A materially smaller number means the report root
  is wrong or partial; a `missing` count above zero should be explained.
- `head -3` shows real pathway names, not blanks. Expect values such as
  `Xenobiotics`, `Lipid`, `Amino Acid`, `Peptide`.
- The file is a few hundred kilobytes. Megabytes means the wrong columns were
  captured.

Sanity figures from the last run, useful as a comparison: among the 704
X-compounds detected in feces, the rank-1 pathway split was Xenobiotics 191,
Lipid 182, Amino Acid 142, Peptide 134, then Carbohydrate 20, Nucleotide 17,
Cofactors and Vitamins 11, Energy 5, Partially Characterized Molecules 2.

## Step 3: send

Upload all three files directly. If any is rejected for size, S3 remains
available as a fallback, per
[docs/mass_policy_bundle_handoff.md](mass_policy_bundle_handoff.md); the cloud
VM has no AWS CLI or credentials, so it can only fetch a presigned HTTPS URL
with `curl`, and the link must return bytes rather than an HTML interstitial.

Report back the row count from step 1 and the `missing` count.

## Known problem with a previously supplied file

`Feces (Unspecified Type) - Human.xlsx` arrived corrupt and was never usable. It
is not needed for this task, but if it is ever sent again it needs re-exporting
rather than re-sending:

- 4,206,592 bytes, beginning with a valid local file header (`PK\x03\x04`), which
  is why `file` reports "Microsoft Excel 2007+".
- No end-of-central-directory record (`PK\x05\x06` absent) and no central
  directory (`PK\x01\x02` absent).
- Only 8 local headers in the whole file, the last at offset 355,685.

That is a truncated or partially written upload rather than a format problem.
Verify any replacement with `unzip -t` before sending.

## If you are exporting from macOS

Archive utilities there add AppleDouble `._*` sidecars, which broke an earlier
transfer. Set `COPYFILE_DISABLE=1` before any `tar`, and confirm with
`tar -tzf <archive> | grep -c '/\._'` returning 0. This only matters if you fall
back to a tarball; the plain CSV upload above is unaffected.

## Unresolved, for a person

- **Which known set should the denominator use?** Two variants are committed and
  they tell different stories. Using everything the Global panel detects in
  feces puts the known-side expected count at 1,087.79 and makes Lipid largest.
  Using the MFP panel as the known set puts it at 526.45 and makes Amino Acid
  largest. The X-compound contribution, 361.37, is the same either way. The
  cloud agent has no basis for choosing.
- **Should low-confidence pathway predictions be included?** Measured
  super-pathway agreement is 52.2% against a 62.5% ceiling, and 294 of the 704
  feces-detected assignments are the pipeline's own `low` tier. Every figure
  above includes them.
