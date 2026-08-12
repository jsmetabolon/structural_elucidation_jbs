# Handoff: export the report sidecars needed for the mass-policy audit

Self-contained task for whoever can reach the host holding the rendered X-compound
reports. Assume no other context. You are only exporting data; the cloud agent runs
the audit.

## Goal

Deliver three small sidecar files per report to a cloud agent over HTTPS, so it can
measure whether changing the mass-match window would move any result. See
`plans/mass_policy_impact_measurement.plan.md` for what the audit does with them.

## What is needed, and what is not

For each report directory `X-*/data/`, exactly three files, where `<slug>` is the
report id lowercased (`X-25984` gives `x-25984`):

| File | Supplies |
|---|---|
| `<slug>_candidate_score_summary.csv` | per-candidate `abs_mass_error_da`, `mass_match_score`, `predicted_smiles`, `predicted_formula` |
| `<slug>_top_candidates.csv` | per-candidate `rank`, `truth_tanimoto` — the truth labels recovery is scored against |
| `<slug>_confidence_inputs.json` | `observed_neutral_mass` for the ppm conversion, plus `msi_level` and `source` |

**Do not include the HTML, images, or the synonyms DB.** They are large and the
audit does not read them. The three files above are a few kilobytes per report, so
the whole set should compress to well under a gigabyte.

The audit is read-only and never writes back to the bundle.

## Do the five-report smoke bundle first

This costs a minute and saves a round trip. Column names are the one thing the
cloud agent cannot verify without real data, and if an assumption is wrong it is
better to find out on five reports than on four thousand.

```bash
set -euo pipefail
REPORTS="/path/to/reports"          # the directory containing X-25984/, X-21434/, ...
STAGE="$(mktemp -d)/mass_policy_smoke"

mkdir -p "$STAGE"
for d in $(ls -d "$REPORTS"/X-*/ | head -5); do
  id="$(basename "$d")"; slug="$(printf '%s' "$id" | tr 'A-Z' 'a-z')"
  dst="$STAGE/reports/$id/data"; mkdir -p "$dst"
  for f in "${slug}_candidate_score_summary.csv" \
           "${slug}_top_candidates.csv" \
           "${slug}_confidence_inputs.json"; do
    [ -f "$d/data/$f" ] && cp "$d/data/$f" "$dst/" || echo "MISSING: $id/$f"
  done
done

tar -C "$STAGE" -czf /tmp/mass_policy_smoke.tgz .
ls -lh /tmp/mass_policy_smoke.tgz && sha256sum /tmp/mass_policy_smoke.tgz
```

Upload and presign that, send the link, and wait for confirmation before the full
export.

## The full bundle

```bash
set -euo pipefail
REPORTS="/path/to/reports"
STAGE="$(mktemp -d)/mass_policy_bundle"

mkdir -p "$STAGE"
n=0; missing=0
for d in "$REPORTS"/X-*/; do
  id="$(basename "$d")"; slug="$(printf '%s' "$id" | tr 'A-Z' 'a-z')"
  [ -d "$d/data" ] || continue
  dst="$STAGE/reports/$id/data"; mkdir -p "$dst"
  for f in "${slug}_candidate_score_summary.csv" \
           "${slug}_top_candidates.csv" \
           "${slug}_confidence_inputs.json"; do
    if [ -f "$d/data/$f" ]; then cp "$d/data/$f" "$dst/"; else missing=$((missing+1)); fi
  done
  n=$((n+1))
done
echo "staged $n reports, $missing missing files"

du -sh "$STAGE"
tar -C "$STAGE" -czf /tmp/mass_policy_bundle.tgz .
ls -lh /tmp/mass_policy_bundle.tgz && sha256sum /tmp/mass_policy_bundle.tgz
```

## Establishing the GET link

The cloud VM has unrestricted outbound HTTPS, so any host it can `curl` works. The
data is Metabolon-internal, so prefer a bucket you control over a public file host,
and keep the expiry short.

```bash
aws s3 cp /tmp/mass_policy_bundle.tgz s3://YOUR_BUCKET/tmp/mass_policy_bundle.tgz
aws s3 presign s3://YOUR_BUCKET/tmp/mass_policy_bundle.tgz --expires-in 86400
```

If `aws` is unavailable, any HTTPS endpoint serving a plain `GET` is fine. Avoid
anything that returns an HTML interstitial rather than the bytes, which is the
usual failure with consumer file-sharing links.

## What to send back

- The presigned HTTPS URL.
- The `sha256` and size, so the download can be verified before it is trusted.
- The report count `n` and the `missing` count from the staging loop. A large
  `missing` number means the naming assumption is wrong and should be resolved
  before the audit runs.

## Verification before handing off

- `tar -tzf /tmp/mass_policy_bundle.tgz | head` shows paths of the form
  `./reports/X-*/data/x-*_top_candidates.csv`.
- `tar -tzf /tmp/mass_policy_bundle.tgz | wc -l` is roughly three times the report
  count.
- The report count is plausible. The delivered set audited to date is 3,897.
- The archive contains no `.html` and no `.sqlite`:
  `tar -tzf /tmp/mass_policy_bundle.tgz | grep -cE '\.(html|sqlite)$'` should be 0.

## Notes

- Paths inside the archive are relative to a `reports/` root, so the cloud agent
  extracts with `tar -xzf` into a scratch directory and points `--reports-root` at
  the `reports/` directory inside it.
- If the sidecars are absent because reports were generated before the enrichment
  step, `tools/enrich_report_evidence.py` is what writes
  `<slug>_confidence_inputs.json`. Without it the ppm conversion has no observed
  neutral mass and the audit degrades to counting band changes on the stored score
  only.
