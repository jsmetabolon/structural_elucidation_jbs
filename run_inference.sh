#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PYTHONPATH="$REPO_ROOT/src:${PYTHONPATH:-}"
PY="${PY:-python3}"
if [[ -x "$REPO_ROOT/.venv/bin/python" ]]; then
  PY="$REPO_ROOT/.venv/bin/python"
elif [[ -x "$REPO_ROOT/../.venv/bin/python" ]]; then
  PY="$REPO_ROOT/../.venv/bin/python"
fi

OUT_ROOT="${OUT_ROOT:-$REPO_ROOT/output/pipeline_equalization_v21}"
RUN_LIBS="${RUN_LIBS:-209,305,400,402}"
START_AT="${START_AT:-verifier}"
MERGED_IONS_CSV="${MERGED_IONS_CSV:-$REPO_ROOT/data/ions_unknown_merged_v21_clean.csv}"
FEATURE_TABLE="${FEATURE_TABLE:-$MERGED_IONS_CSV}"
ADDUCT_LADDER="${ADDUCT_LADDER:-$REPO_ROOT/data/chrolib_adduct_ladder_v21.csv}"
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

run_cmd() {
  echo
  echo "+ $*"
  "$@"
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
    run_cmd "$PY" "$REPO_ROOT/src/build_library_anchored_candidates_v1.py" \
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

  run_cmd "$PY" "$REPO_ROOT/src/score_candidates_with_verifier_v1.py" \
    --output-root "$OUT_ROOT/verifier_lib${lib}" \
    --observed-lookup "$observed_lookup" \
    --all-x-observed-lookup "$observed_lookup" \
    --current-integrated-best "$(current_integrated "$lib")" \
    --input-candidate-metrics "$inputs" \
    --lib-id "$lib" --ion-mode "$ion" --include-all-x

  run_cmd "$PY" "$REPO_ROOT/src/build_unified_candidate_pool_v16.py" \
    --output-root "$OUT_ROOT/unified_pool_lib${lib}" \
    --scored-candidates "$OUT_ROOT/verifier_lib${lib}/scored_candidates.csv" \
    --base-unified-pool-csv "$base_pool" \
    --current-integrated-best "$(current_integrated "$lib")" \
    --lib-id "$lib" --ion-mode "$ion" \
    --include-all-x

  run_cmd "$PY" "$REPO_ROOT/src/train_unified_rerank_v1.py" \
    --unified-pool "$OUT_ROOT/unified_pool_lib${lib}/unified_candidates_per_control.csv" \
    --non-s0-pool "$OUT_ROOT/unified_pool_lib${lib}/unified_candidates_per_control_non_s0.csv" \
    --output-root "$OUT_ROOT/rerank_lib${lib}" \
    --lib-id "$lib" --ion-mode "$ion" \
    --model-selection-mode auto
}

mkdir -p "$OUT_ROOT"
IFS=',' read -r -a libs <<< "$RUN_LIBS"
for lib in "${libs[@]}"; do
  lib=$(echo "$lib" | xargs)
  [[ -z "$lib" ]] && continue
  run_scope "$lib"
done

echo "[run_inference] completed $(date -u +%Y-%m-%dT%H:%M:%SZ)"
