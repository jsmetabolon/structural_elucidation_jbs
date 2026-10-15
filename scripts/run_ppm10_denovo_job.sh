#!/usr/bin/env bash
# Scoped SELFIES de novo for compounds with no library candidate at 10 ppm.
#
# Decodes only the compounds in output/ppm10_remaining/nohit10_ids.txt that
# have an encoder feature row, adds the June bundle's decoded structures for
# the same compounds, rescores mass against adduct-derived neutral masses, then
# scores, refills at 10 ppm, and merges new compounds through
# scripts/run_ppm10_extra_job.sh.
set -uo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"
: "${PY:?set PY}"
: "${WS:?set WS (workspace holding generate_selfies_denovo.py)}"
export PY
export PYTHONPATH="${PYTHONPATH:-$REPO_ROOT}"

BUNDLE="${BUNDLE:-$WS/models/merged_structuremasst_bundle_excl_x_unresolved_2026-04-17_metabolon}"
ENCODER="${ENCODER:-$WS/models/merged_main_excl_x_unresolved_2026-04-17/best_encoder.pt}"
NUM_SAMPLES="${NUM_SAMPLES:-60}"
IDS="${IDS:-output/ppm10_remaining/nohit10_ids.txt}"
SELECTION="${SELECTION:-output/ppm10_recovery_inputs/selection.json}"
WORK="${WORK:-output/ppm10_denovo}"
INPUTS="$WORK/merged_structuremasst_bundle_ppm10/selfies_profiles"

rm -rf "$WORK"
"$PY" tools/prepare_denovo_ppm10_inputs.py slice \
  --bundle-dir "$BUNDLE" --ids "$IDS" --selection "$SELECTION" --out-dir "$WORK/features" \
  || { echo "slice failed"; exit 1; }

for features in "$WORK"/features/lib*_features.csv; do
  [[ -e "$features" ]] || continue
  name="$(basename "$features" _features.csv)"
  lib="${name#lib}"; lib="${lib%%_*}"
  ion="${name##*_}"
  if [[ "$ion" == "positive" ]]; then pair=pos_400_402; else pair=neg_209_305; fi
  out="$REPO_ROOT/$WORK/decode/lib_${lib}_${ion}/selfies_denovo"
  echo "[denovo] decoding lib $lib ($ion)"
  (cd "$WS" && env -u PYTHONPATH "$PY" generate_selfies_denovo.py \
    --model_ckpt "models/selfies_decoder_${pair}/selfies_decoder_${pair}.pt" \
    --vocab_json "models/selfies_decoder_${pair}/selfies_vocab_${pair}.json" \
    --feature_table "$REPO_ROOT/$features" \
    --encoder_ckpt "$ENCODER" --encoder_kind merged \
    --lib_id "$lib" --control_lib_id "$lib" --ion_mode "$ion" \
    --output_dir "$out" --artifact_prefix "merged_${lib}" \
    --num_samples "$NUM_SAMPLES" --device cpu \
    --minmass 50.0 --maxmass 950.0 --parent_mz_min 50.0 --parent_mz_max 950.0 \
    --disable_ri_conditioning --skip_control_eval) \
    || { echo "decode failed for lib $lib"; exit 1; }
done

"$PY" tools/prepare_denovo_ppm10_inputs.py inputs \
  --candidate-roots "$WORK/decode" "$BUNDLE" \
  --ids "$IDS" --selection "$SELECTION" \
  --feature-dir "${FEATURE_DIR:-output/ppm10_recovery_inputs}" \
  --out-dir "$INPUTS" \
  || { echo "inputs failed"; exit 1; }

DENOVO_ROOT="${DENOVO_ROOT:-output/pipeline_equalization_adduct_derive_denovo}"
EXTRA_START_AT=verifier MIN_SPECTRUM_SCORE="${MIN_SPECTRUM_SCORE:-0}" \
  ./scripts/run_ppm10_extra_job.sh "$INPUTS" "$DENOVO_ROOT"
