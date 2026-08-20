#!/usr/bin/env bash
set -euo pipefail

# Profile one library scope end-to-end: serial baseline vs optimized vs warm
# score-cache rerun, with output-parity verification between all arms.
#
# Prerequisites: required data staged (./fetch_get_assets.sh --verify passes).
# Expect multi-minute wall clock per arm on real data.
#
# Usage:
#   LIB=209 ./tools/profile_baseline.sh
# Optional env:
#   LIB           library scope to profile (default 209)
#   PROFILE_ROOT  output root (default output/profile_lib<LIB>)
#   SKIP_BASELINE=1  skip the serial baseline arm (e.g. on re-runs)
#
# After single-scope profiling, measure the production shape with:
#   OUT_ROOT=output/full_parallel JOBS=4 START_AT=full ./run_inference.sh
#
# For statistical profiles of one stage, py-spy works well, e.g.:
#   py-spy record -o verify.svg -- .venv/bin/python -m structural_elucidation.pipeline.verify --workers 1 ...

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

LIB="${LIB:-209}"
PROFILE_ROOT="${PROFILE_ROOT:-$REPO_ROOT/output/profile_lib${LIB}}"
NPROC="$(nproc 2>/dev/null || echo 4)"

PY="${PY:-python3}"
if [[ -x "$REPO_ROOT/.venv/bin/python" ]]; then
  PY="$REPO_ROOT/.venv/bin/python"
fi

echo "== 1/5: verify staged assets =="
./fetch_get_assets.sh --verify

echo
echo "== 2/5: validate NCE modelfiles load in this environment =="
SE_PROFILE_LIB="$LIB" "$PY" - <<'PYEOF'
import os
lib = int(os.environ["SE_PROFILE_LIB"])
from structural_elucidation.nce.runtime import get_predicted_spectrum_matcher
matcher = get_predicted_spectrum_matcher(lib)
print(f"matcher OK for lib {lib}: {type(matcher).__name__} "
      f"(model bins: {len(getattr(matcher, 'model_dict', {}))})")
PYEOF

run_arm() {
  local name="$1"
  shift
  echo
  echo "== arm: $name =="
  local out="$PROFILE_ROOT/$name"
  mkdir -p "$out"
  ( export "$@"; \
    OUT_ROOT="$out" RUN_LIBS="$LIB" START_AT=full JOBS=1 ./run_inference.sh )
  echo "-- $name stage timings:"
  cat "$out/stage_timings.csv" || true
}

if [[ "${SKIP_BASELINE:-0}" != "1" ]]; then
  echo
  echo "== 3/5: serial baseline (workers=1, score cache disabled) =="
  run_arm serial_baseline \
    SE_VERIFIER_WORKERS=1 SE_RERANK_JOBS=1 SE_CANDIDATES_WORKERS=1 SE_VERIFIER_NO_CACHE=1
fi

echo
echo "== 4/5: optimized cold (workers=$NPROC, fresh score cache) =="
rm -rf "$REPO_ROOT/cache/verifier_scores"
run_arm optimized_cold \
  SE_VERIFIER_WORKERS="$NPROC" SE_RERANK_JOBS="$NPROC" SE_CANDIDATES_WORKERS="$NPROC" SE_VERIFIER_NO_CACHE=0

echo
echo "== 4b: optimized warm rerun (score cache reused) =="
run_arm optimized_warm \
  SE_VERIFIER_WORKERS="$NPROC" SE_RERANK_JOBS="$NPROC" SE_CANDIDATES_WORKERS="$NPROC" SE_VERIFIER_NO_CACHE=0

echo
echo "== 5/5: output parity checks =="
if [[ -d "$PROFILE_ROOT/serial_baseline" ]]; then
  "$PY" tools/compare_run_outputs.py "$PROFILE_ROOT/serial_baseline" "$PROFILE_ROOT/optimized_cold"
fi
"$PY" tools/compare_run_outputs.py "$PROFILE_ROOT/optimized_cold" "$PROFILE_ROOT/optimized_warm"

echo
echo "== timing summary =="
for arm in serial_baseline optimized_cold optimized_warm; do
  f="$PROFILE_ROOT/$arm/stage_timings.csv"
  [[ -f "$f" ]] || continue
  total=$(awk -F, 'NR>1 {s+=$4} END {printf "%.1f", s}' "$f")
  echo "$arm: total ${total}s ($f)"
done
echo "[profile_baseline] done"
