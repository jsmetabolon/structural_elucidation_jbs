#!/usr/bin/env bash
set -euo pipefail

# Generate self-contained HTML reports for every X compound found in the rerank
# prediction CSVs, then build a master index. Thin wrapper around
# `python -m structural_elucidation.reports.batch`.
#
# Typical use (inputs already present locally):
#   OUT_ROOT=output/pipeline_equalization JOBS=8 ./generate_all_reports.sh
#
# Quick scope check without generating anything:
#   DRY_RUN=1 OUT_ROOT=output/pipeline_equalization ./generate_all_reports.sh
#
# Env vars:
#   PY            python to use (default: .venv/bin/python if present, else python3)
#   OUT_ROOT      root holding rerank_lib*/<PRED_NAME> (default: output/pipeline_equalization)
#   RUN_LIBS      comma libs to include (default: 209,305,400,402)
#   PRED_NAME     per-lib predictions filename (default: unified_rerank_all_top10.csv)
#   REPORTS_OUT   master output dir (default: output/reports)
#   JOBS          parallel compound builds (default: 1)
#   DRY_RUN       1 = enumerate only (default: unset)
#   OVERWRITE     1 = rebuild existing reports (default: unset, i.e. resume)
#   LIMIT         only build first N compounds (default: unset)
#   OBSERVED_LOOKUP_DIR  dir of lib{N}_{ion}_observed_ms2_lookup.csv for MS scoring (optional)

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$REPO_ROOT"

PY="${PY:-}"
if [[ -z "$PY" ]]; then
  if [[ -x "$REPO_ROOT/.venv/bin/python" ]]; then
    PY="$REPO_ROOT/.venv/bin/python"
  else
    PY="python3"
  fi
fi

OUT_ROOT="${OUT_ROOT:-output/pipeline_equalization}"
RUN_LIBS="${RUN_LIBS:-209,305,400,402}"
PRED_NAME="${PRED_NAME:-unified_rerank_all_top10.csv}"
REPORTS_OUT="${REPORTS_OUT:-output/reports}"
JOBS="${JOBS:-1}"

args=(
  -m structural_elucidation.reports.batch
  --equalization-root "$OUT_ROOT"
  --libs "$RUN_LIBS"
  --predictions-name "$PRED_NAME"
  --out-dir "$REPORTS_OUT"
  --jobs "$JOBS"
)

[[ "${DRY_RUN:-0}" == "1" ]] && args+=(--dry-run)
[[ "${OVERWRITE:-0}" == "1" ]] && args+=(--overwrite)
[[ -n "${LIMIT:-}" ]] && args+=(--limit "$LIMIT")
[[ -n "${OBSERVED_LOOKUP_DIR:-}" ]] && args+=(--observed-lookup-dir "$OBSERVED_LOOKUP_DIR")

echo "[generate_all_reports] PY=$PY"
echo "[generate_all_reports] OUT_ROOT=$OUT_ROOT RUN_LIBS=$RUN_LIBS JOBS=$JOBS REPORTS_OUT=$REPORTS_OUT"
echo "+ $PY ${args[*]}"
exec "$PY" "${args[@]}"
