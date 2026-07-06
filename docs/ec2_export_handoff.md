# Handoff: export X-compound rerank inputs from EC2 for cloud HTML generation

Self-contained task for an agent running in the Cursor window that **can reach the
EC2 host** `ec2-user@10.12.5.28`. Assume you have no other context. Do not
generate any HTML yourself; your only job is to export inputs and deliver them.

## Goal

Export the pipeline result files that the HTML report generator needs, bundle
them in a fixed layout, and deliver the bundle to the requesting (cloud) agent
via a URL it can `curl`. The cloud agent will do all report generation.

## Background (enough to act)

- This is the Metabolon **X-compound structural elucidation** pipeline. "X
  compounds" are unknowns with ids like `X-25003`; each has ranked SMILES
  predictions in per-library rerank CSVs.
- A cloud agent already built a batch HTML report generator
  (`structural_elucidation.reports.batch` + `generate_all_reports.sh`) but it
  runs on an isolated VM that **cannot SSH into** `10.12.5.28`. You can. So you
  export the data; it generates.

## Inputs to export

1. **Required — per-lib rerank predictions** (small, a few MB each):
   `<OUT_ROOT>/rerank_lib{209,305,400,402}/unified_rerank_all_top10.csv`
   where `OUT_ROOT` is typically `models/pipeline_equalization` (full pipeline)
   or `output/pipeline_equalization`. Each CSV has a `control_id` column (the X
   id), `predicted_smiles`, `rank`, and score columns.
2. **Optional — metabolite synonyms DB** (can be multiple GB):
   a SQLite file usually named like `metabolite_synonyms_merged_with_mtx_giga.sqlite`.
   Enables library matches, synonyms, pathways, and the known-metabolite
   similarity network. Reports still build without it (degraded, smaller).
3. **Optional — observed MS2 lookups** (for MS-fit scoring):
   `data/observed_ms2_lookups/lib{lib}_{ion}_observed_ms2_lookup.csv`.

## Hard constraints

- The cloud VM's SSH to `10.12.5.28` is reset pre-auth (server-side IP
  allowlist); its egress IP `54.159.155.113` is not whitelisted. Do **not** try
  to "fix" cloud SSH — deliver data out-of-band instead.
- The cloud VM **does** have outbound HTTPS, so a presigned S3 URL (or any
  HTTPS-reachable location) works for delivery.
- Data is Metabolon-internal: prefer your **own** S3 bucket, not public hosts.
- The bundle's internal paths must be **repo-root-relative** so the cloud agent
  can `tar -C <repo> -xzf` cleanly (see layout below).

## Required bundle layout (inside the .tgz)

```
output/pipeline_equalization/rerank_lib209/unified_rerank_all_top10.csv
output/pipeline_equalization/rerank_lib305/unified_rerank_all_top10.csv
output/pipeline_equalization/rerank_lib400/unified_rerank_all_top10.csv
output/pipeline_equalization/rerank_lib402/unified_rerank_all_top10.csv
vendor/metabolite_synonyms.sqlite                # only if including the DB
```

## Decisions (defaults; change if needed)

- **Include the DB?** Default yes (full reports). For a fast first pass, omit it
  (skip the DB block) — the cloud agent can fold it in later.
- **Transfer = S3 presigned URL**, 24h expiry. Alternatives: any HTTPS file
  host; or, for CSVs only (never the multi-GB DB), commit them to the existing
  branch `cursor/batch-x-compound-html-reports-f8a2` if this host has git auth.

## Steps

If your shell is **on the EC2 host**, run as-is. If your shell is a **local
machine with SSH to EC2**, run these on EC2 (e.g. `ssh ec2-user@10.12.5.28`)
or adapt paths via `scp`.

```bash
set -euo pipefail
WS="${WS:-/data/josh/spec2mol/workspaces/metabolon}"   # adjust to your workspace
cd "$WS"

# 1) Locate rerank CSVs (works regardless of models/ vs output/ OUT_ROOT)
mapfile -t PREDS < <(find . -type f -path '*/rerank_lib*/unified_rerank_all_top10.csv' 2>/dev/null | sort)
echo "found ${#PREDS[@]} prediction CSVs:"; printf '  %s\n' "${PREDS[@]}"
test "${#PREDS[@]}" -gt 0 || { echo "No rerank CSVs found under $WS — fix WS and retry"; exit 1; }

# 2) Stage into the required repo-relative layout (normalizes to output/pipeline_equalization)
STAGE="$(mktemp -d)/se_inputs"; mkdir -p "$STAGE/output/pipeline_equalization" "$STAGE/vendor"
for p in "${PREDS[@]}"; do
  libdir="$(basename "$(dirname "$p")")"          # rerank_libNNN
  mkdir -p "$STAGE/output/pipeline_equalization/$libdir"
  cp "$p" "$STAGE/output/pipeline_equalization/$libdir/"
done

# 3) OPTIONAL: include the metabolite synonyms DB (skip for a fast first pass)
DB="$(find . -iname 'metabolite_synonyms*giga*.sqlite' 2>/dev/null | head -1)"
[ -z "$DB" ] && DB="$(find . -iname 'metabolite_synonyms*.sqlite' 2>/dev/null | head -1)"
[ -n "$DB" ] && cp "$DB" "$STAGE/vendor/metabolite_synonyms.sqlite" && echo "included DB: $DB"

# 4) Per-lib + total unique X-compound counts (report these back)
python3 - "$STAGE" <<'PY'
import sys, glob, os, pandas as pd
root = sys.argv[1]
tot = set()
for f in sorted(glob.glob(os.path.join(root, "output/pipeline_equalization/rerank_lib*/unified_rerank_all_top10.csv"))):
    df = pd.read_csv(f, usecols=lambda c: c in ("control_id", "x_compound_id"))
    col = "control_id" if "control_id" in df.columns else "x_compound_id"
    ids = set(df[col].dropna().astype(str)) - {"", "nan", "none"}
    tot |= ids
    print(f"  {os.path.basename(os.path.dirname(f))}: {len(ids)} unique compounds")
print(f"TOTAL unique X compounds: {len(tot)}")
PY

# 5) Bundle + integrity info (report size + sha256 back)
tar -C "$STAGE" -czf /tmp/se_report_inputs.tgz .
ls -lh /tmp/se_report_inputs.tgz
sha256sum /tmp/se_report_inputs.tgz
tar -tzf /tmp/se_report_inputs.tgz | sed 's#^\./##' | grep -E 'unified_rerank_all_top10.csv|metabolite_synonyms.sqlite'

# 6) Upload + presign (adjust bucket/region)
aws s3 cp /tmp/se_report_inputs.tgz s3://YOUR_BUCKET/tmp/se_report_inputs.tgz
aws s3 presign s3://YOUR_BUCKET/tmp/se_report_inputs.tgz --expires-in 86400
```

## What to return to the requesting (cloud) agent

- The **presigned HTTPS URL** from step 6.
- The **sha256** and **size** from step 5.
- The **per-lib + total unique compound counts** from step 4.
- Whether the **DB was included** (full vs degraded run).

## Verification checklist (before handing off)

- `find` returned the expected number of lib CSVs (ideally all of 209/305/400/402).
- Step 5's `tar -tzf | grep` lists each `rerank_lib*/unified_rerank_all_top10.csv`
  (and `vendor/metabolite_synonyms.sqlite` if DB included).
- Total unique compound count is non-zero and plausible.
- Presigned URL is freshly generated (24h expiry) and points at the uploaded object.

## Unresolved / adjust if needed

- Confirm `WS` (workspace path on EC2). If `find` returns nothing, fix it.
- Confirm the DB filename/path if the `find` heuristic misses it.
- Confirm an S3 bucket you can write to (or pick another HTTPS delivery method).
- If duplicate `rerank_lib*` dirs exist under multiple roots, prefer the
  canonical `OUT_ROOT` (latest/full run) to avoid mixing versions.
