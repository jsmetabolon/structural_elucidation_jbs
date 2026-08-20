#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY="${PY:-python3}"
if [[ -x "$REPO_ROOT/.venv/bin/python" ]]; then
  PY="$REPO_ROOT/.venv/bin/python"
elif [[ -x "$REPO_ROOT/../.venv/bin/python" ]]; then
  PY="$REPO_ROOT/../.venv/bin/python"
fi

OUT_ROOT="${OUT_ROOT:-$REPO_ROOT/output/pipeline_equalization}"
RUN_LIBS="${RUN_LIBS:-209,305,400,402}"
START_AT="${START_AT:-verifier}"
# JOBS>1 runs independent library scopes concurrently (per-scope logs under
# OUT_ROOT). Per-stage worker counts are divided across concurrent scopes
# unless explicitly provided.
JOBS="${JOBS:-1}"
TOTAL_CORES="$(nproc 2>/dev/null || echo 4)"
default_stage_workers=$(( TOTAL_CORES / (JOBS > 0 ? JOBS : 1) ))
(( default_stage_workers < 1 )) && default_stage_workers=1
export SE_VERIFIER_WORKERS="${SE_VERIFIER_WORKERS:-$default_stage_workers}"
export SE_RERANK_JOBS="${SE_RERANK_JOBS:-$default_stage_workers}"
MERGED_IONS_CSV="${MERGED_IONS_CSV:-$REPO_ROOT/data/ions_unknown_merged_clean.csv}"
FEATURE_TABLE="${FEATURE_TABLE:-$MERGED_IONS_CSV}"
ADDUCT_LADDER="${ADDUCT_LADDER:-$REPO_ROOT/data/chrolib_adduct_ladder.csv}"
OBSERVED_209="${OBSERVED_209:-$REPO_ROOT/data/observed_ms2_lookups/per_control_all_x.csv}"

scope_ion() {
  case "$1" in
    209|305) echo negative ;;
    400|402) echo positive ;;
    *) echo "Unsupported lib: $1" >&2; exit 2 ;;
  esac
}

current_integrated() {
  case "$1" in
    209) echo "$REPO_ROOT/data/integrated_best/lib_209_negative.csv" ;;
    305) echo "$REPO_ROOT/data/integrated_best/lib_305_negative.csv" ;;
    400) echo "$REPO_ROOT/data/integrated_best/lib_400_positive.csv" ;;
    402) echo "$REPO_ROOT/data/integrated_best/lib_402_positive.csv" ;;
  esac
}

base_pool_for() {
  local lib="$1"
  local ion
  ion=$(scope_ion "$lib")
  echo "$REPO_ROOT/data/base_pools/lib${lib}_${ion}.csv"
}

TIMINGS_CSV="${TIMINGS_CSV:-$OUT_ROOT/stage_timings.csv}"

record_timing() {
  local stage="$1" lib="$2" seconds="$3" rc="$4"
  if [[ ! -f "$TIMINGS_CSV" ]]; then
    echo "timestamp_utc,stage,lib_id,seconds,exit_code" > "$TIMINGS_CSV"
  fi
  echo "$(date -u +%Y-%m-%dT%H:%M:%SZ),$stage,$lib,$seconds,$rc" >> "$TIMINGS_CSV"
}

run_stage() {
  local stage="$1" lib="$2"
  shift 2
  echo
  echo "+ $*"
  local start_ns end_ns rc=0
  start_ns=$(date +%s%N)
  if "$@"; then rc=0; else rc=$?; fi
  end_ns=$(date +%s%N)
  local seconds
  seconds=$(awk -v s="$start_ns" -v e="$end_ns" 'BEGIN {printf "%.2f", (e - s) / 1e9}')
  echo "[timing] stage=$stage lib=$lib seconds=$seconds exit=$rc"
  record_timing "$stage" "$lib" "$seconds" "$rc"
  return "$rc"
}

append_existing() {
  local joined="$1"
  local path="$2"
  if [[ -s "$path" ]]; then
    if [[ -n "$joined" ]]; then
      printf "%s,%s" "$joined" "$path"
    else
      printf "%s" "$path"
    fi
  else
    printf "%s" "$joined"
  fi
}

run_scope() {
  local lib="$1"
  local ion
  ion=$(scope_ion "$lib")
  local manifest="$REPO_ROOT/data/feature_rows/lib${lib}_${ion}_gap_manifest.csv"
  local feature_rows="$REPO_ROOT/data/feature_rows/lib${lib}_${ion}_feature_rows.csv"
  local mass_hyp="$REPO_ROOT/data/feature_rows/lib${lib}_${ion}_neutral_mass_hypotheses.csv"
  local observed_lookup="$REPO_ROOT/data/observed_ms2_lookups/lib${lib}_${ion}_observed_ms2_lookup.csv"
  local base_pool
  base_pool=$(base_pool_for "$lib")

  if [[ ! -s "$manifest" ]]; then
    echo "[run_inference] missing manifest for lib${lib}_${ion}; skipping"
    return 0
  fi

  local inputs=""
  local la_path="$OUT_ROOT/library_anchored_lib${lib}/runs/library_anchored_v1/track_control_blind/lib_${lib}_${ion}/selfies_profiles/strict_1da/candidate_metrics_all.csv"

  if [[ "$START_AT" != "verifier" ]]; then
    run_stage candidates "$lib" "$PY" -m structural_elucidation.pipeline.candidates \
      --output-root "$OUT_ROOT/library_anchored_lib${lib}" \
      --diagnosis-csv "$manifest" \
      --feature-rows "$feature_rows" \
      --neutral-mass-hypotheses "$mass_hyp" \
      --lib-id "$lib" --ion-mode "$ion" \
      --controls "$(awk 'NR>1{print $1}' "$manifest" | paste -sd, -)"
  fi
  inputs=$(append_existing "$inputs" "$la_path")

  if [[ -z "$inputs" ]]; then
    echo "[run_inference] no candidate inputs for lib${lib}_${ion}; skipping verifier/pool/rerank"
    return 0
  fi

  run_stage verify "$lib" "$PY" -m structural_elucidation.pipeline.verify \
    --output-root "$OUT_ROOT/verifier_lib${lib}" \
    --observed-lookup "$observed_lookup" \
    --all-x-observed-lookup "$observed_lookup" \
    --current-integrated-best "$(current_integrated "$lib")" \
    --input-candidate-metrics "$inputs" \
    --lib-id "$lib" --ion-mode "$ion" --include-all-x

  run_stage pool "$lib" "$PY" -m structural_elucidation.pipeline.pool \
    --output-root "$OUT_ROOT/unified_pool_lib${lib}" \
    --scored-candidates "$OUT_ROOT/verifier_lib${lib}/scored_candidates.csv" \
    --base-unified-pool-csv "$base_pool" \
    --current-integrated-best "$(current_integrated "$lib")" \
    --lib-id "$lib" --ion-mode "$ion" \
    --include-all-x

  run_stage rerank "$lib" "$PY" -m structural_elucidation.pipeline.rerank \
    --unified-pool "$OUT_ROOT/unified_pool_lib${lib}/unified_candidates_per_control.csv" \
    --non-s0-pool "$OUT_ROOT/unified_pool_lib${lib}/unified_candidates_per_control_non_s0.csv" \
    --output-root "$OUT_ROOT/rerank_lib${lib}" \
    --lib-id "$lib" --ion-mode "$ion" \
    --model-selection-mode auto
}

run_scopes_parallel() {
  local pending=("$@")
  local overall=0
  while (( ${#pending[@]} )); do
    local batch=("${pending[@]:0:JOBS}")
    pending=("${pending[@]:JOBS}")
    local -A pid_of=()
    local lib
    for lib in "${batch[@]}"; do
      echo "[run_inference] starting scope lib${lib} (log: $OUT_ROOT/lib${lib}_scope.log)"
      run_scope "$lib" > "$OUT_ROOT/lib${lib}_scope.log" 2>&1 &
      pid_of[$lib]=$!
    done
    for lib in "${batch[@]}"; do
      if wait "${pid_of[$lib]}"; then
        echo "[run_inference] scope lib${lib} completed"
      else
        local rc=$?
        overall=$rc
        echo "[run_inference] scope lib${lib} FAILED rc=$rc (see $OUT_ROOT/lib${lib}_scope.log)" >&2
      fi
    done
  done
  return "$overall"
}

mkdir -p "$OUT_ROOT"
IFS=',' read -r -a libs <<< "$RUN_LIBS"
clean_libs=()
for lib in "${libs[@]}"; do
  lib=$(echo "$lib" | xargs)
  [[ -z "$lib" ]] && continue
  clean_libs+=("$lib")
done

if (( JOBS > 1 && ${#clean_libs[@]} > 1 )); then
  run_scopes_parallel "${clean_libs[@]}"
else
  for lib in "${clean_libs[@]}"; do
    run_scope "$lib"
  done
fi

echo "[run_inference] completed $(date -u +%Y-%m-%dT%H:%M:%SZ)"
