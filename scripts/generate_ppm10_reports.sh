#!/usr/bin/env bash
set -euo pipefail
# Full ppm-gated report rebuild + enrich/patch + exec summary + ppm audit.
# Expected runtime: tens of minutes at JOBS=8.

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

export PY="${PY:-/Users/jsimpson/miniconda3/bin/python}"
export PYTHONPATH="${PYTHONPATH:-$REPO_ROOT}"
export PRED_NAME="${PRED_NAME:-unified_rerank_all_top10_ppm10.csv}"
export OUT_ROOT="${OUT_ROOT:-output/pipeline_equalization}"
export REPORTS_OUT="${REPORTS_OUT:-output/reports_ppm10}"
export OVERWRITE="${OVERWRITE:-1}"
export JOBS="${JOBS:-8}"
export CANDIDATE_TAB_N="${CANDIDATE_TAB_N:-10}"

echo "[ppm10] generating reports into $REPORTS_OUT"
./generate_all_reports.sh

BUNDLE="${BUNDLE:-$REPO_ROOT/data}"
KNOWN_X="${KNOWN_X:-$REPO_ROOT/data/known_X.csv}"

echo "[ppm10] enrich"
"$PY" tools/enrich_report_evidence.py "$REPORTS_OUT" "$BUNDLE" "$KNOWN_X" || true
echo "[ppm10] patch confidence"
"$PY" tools/patch_report_confidence.py "$REPORTS_OUT"
echo "[ppm10] patch candidate scores"
"$PY" tools/patch_candidate_scores.py "$REPORTS_OUT"
echo "[ppm10] patch consensus crowding"
"$PY" tools/patch_consensus_crowding.py "$REPORTS_OUT"
echo "[ppm10] rebuild master index"
"$PY" tools/rebuild_master_index.py "$REPORTS_OUT"

EXEC_OUT="${EXEC_OUT:-$REPO_ROOT/docs/paper/data/exec_summary_ppm10}"
mkdir -p "$EXEC_OUT"
echo "[ppm10] exec summary -> $EXEC_OUT"
"$PY" tools/build_exec_summary.py --reports-root "$REPORTS_OUT" --out-dir "$EXEC_OUT"

THRESH="${THRESH:-$REPO_ROOT/docs/paper/data/thresholds}"
echo "[ppm10] ppm audit (gated) -> ${THRESH}/audit_ppm_qc_gated.*"
"$PY" tools/audit_ppm_qc.py \
  --by-compound "$EXEC_OUT/exec_summary_by_compound.csv" \
  --top10 "$EXEC_OUT/exec_summary_top10.csv" \
  --json-out "$THRESH/audit_ppm_qc_gated.json" \
  --md-out "$THRESH/audit_ppm_qc_gated.md" \
  --csv-out "$THRESH/audit_ppm_qc_gated.csv" \
  --filtered-candidates-out "$THRESH/audit_ppm_qc_gated_filtered_candidates.csv" \
  --filtered-by-compound-out "$THRESH/audit_ppm_qc_gated_filtered_by_compound.csv"

echo "[ppm10] done"
