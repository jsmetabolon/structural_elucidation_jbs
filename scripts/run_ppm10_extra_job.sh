#!/usr/bin/env bash
# Score extra compounds at a 10 ppm library window and add only new <=10 ppm
# compounds to the main summary.
#
# Usage: scripts/run_ppm10_extra_job.sh <inputs_dir> <extra_out_root>
#
# <inputs_dir> holds lib*_{feature_rows,neutral_mass_hypotheses,gap_manifest}.csv
# and selection.json (run_libs). Reports already in MAIN_ROOT are not rebuilt.
# EXTRA_START_AT=verifier skips candidate generation and scores
# <inputs_dir>/lib*_candidate_metrics_all.csv instead.
set -uo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

inputs_dir="${1:?inputs dir}"
extra_root="${2:?extra output root}"
[[ "$inputs_dir" = /* ]] || inputs_dir="$REPO_ROOT/$inputs_dir"
MAIN_ROOT="${MAIN_ROOT:-output/pipeline_equalization_adduct_derive}"
: "${PY:?set PY}"
export PY
export PYTHONPATH="${PYTHONPATH:-$REPO_ROOT}"

run_libs="$("$PY" -c 'import json,sys; print(json.load(open(sys.argv[1]))["run_libs"])' "$inputs_dir/selection.json")"
echo "RUN_LIBS=$run_libs"
if [[ -z "$run_libs" ]]; then
  echo "EXIT=0"
  echo "nothing to score"
  echo GATHER_EXIT=0
  exit 0
fi

OUT_ROOT="$extra_root" \
  FEATURE_ROWS_DIR="$inputs_dir" \
  CANDIDATE_INPUTS_DIR="$inputs_dir" \
  EMPTY_BASE_POOL=1 \
  START_AT="${EXTRA_START_AT:-candidates}" \
  CANDIDATES_EXTRA_ARGS="--derive-neutral-mass-from-adduct --ppm-tolerance 10" \
  VERIFY_EXTRA_ARGS=--ignore-acceptance \
  JOBS=4 \
  RUN_LIBS="$run_libs" \
  ./run_inference.sh
status=$?
echo "EXIT=$status"

present_libs=""
IFS=',' read -r -a lib_items <<< "$run_libs"
for lib in "${lib_items[@]}"; do
  if [[ -s "$extra_root/rerank_lib$lib/unified_rerank_all_top10.csv" ]]; then
    present_libs="${present_libs:+$present_libs,}$lib"
  fi
done
if [[ -z "$present_libs" ]]; then
  echo "no rerank output to merge"
  echo GATHER_EXIT=$status
  exit "$status"
fi

"$PY" tools/refill_top10_ppm.py \
  --equalization-root "$extra_root" \
  --libs "$present_libs" \
  --top-k 10 \
  --ppm-ceiling 10 \
  ${MIN_SPECTRUM_SCORE:+--min-spectrum-score "$MIN_SPECTRUM_SCORE"} \
  --summary-json "$extra_root/ppm10_refill_summary.json" || { echo GATHER_EXIT=$?; exit 1; }
merge_out="$("$PY" tools/merge_ppm10_recovery.py \
  --main-root "$MAIN_ROOT" \
  --recovery-root "$extra_root" \
  --libs "$present_libs")" || { echo GATHER_EXIT=$?; exit 1; }
printf '%s\n' "$merge_out"
new="$(printf '%s\n' "$merge_out" | sed -n 's/.*"new_compounds": \([0-9][0-9]*\).*/\1/p' | tail -1)"
if [[ "${new:-0}" == "0" ]]; then
  echo "no new compounds at or under 10 ppm"
  echo GATHER_EXIT=0
  exit 0
fi

OVERWRITE=0 \
  OUT_ROOT="$MAIN_ROOT" \
  PRED_NAME=unified_rerank_all_top10_ppm10.csv \
  REPORTS_OUT="$MAIN_ROOT/reports" \
  EXEC_OUT="$MAIN_ROOT/exec_summary" \
  CANDIDATE_TAB_N=10 \
  RUN_LIBS=209,305,400,402 \
  JOBS=4 \
  ./scripts/generate_ppm10_full_reports.sh
echo GATHER_EXIT=$?
