#!/usr/bin/env bash
set -euo pipefail
# Rebuild ppm-gated reports + exec summary from the full peptide_corrected pools
# (3897 delivered set). Refill must already have written:
#   output/pipeline_equalization_peptide_corrected/rerank_lib*/unified_rerank_all_top10_ppm10.csv
#
# Expected runtime: ~1.5–2.5 h at JOBS=8 for ~3200 compounds with ≥1 ≤10 ppm candidate.
# Compounds with zero ≤10 ppm survivors are omitted (~700).

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

export PY="${PY:-/Users/jsimpson/miniconda3/bin/python}"
export PYTHONPATH="${PYTHONPATH:-$REPO_ROOT}"
export OUT_ROOT="${OUT_ROOT:-output/pipeline_equalization_peptide_corrected}"
export PRED_NAME="${PRED_NAME:-unified_rerank_all_top10_ppm10.csv}"
export REPORTS_OUT="${REPORTS_OUT:-output/reports_ppm10_full}"
export OVERWRITE="${OVERWRITE:-1}"
export JOBS="${JOBS:-8}"
export CANDIDATE_TAB_N="${CANDIDATE_TAB_N:-10}"

echo "[ppm10-full] generating reports into $REPORTS_OUT"
./generate_all_reports.sh

BUNDLE="${BUNDLE:-$REPO_ROOT/data}"
KNOWN_X="${KNOWN_X:-$REPO_ROOT/data/known_X.csv}"

if [[ "${NAMING_POLICY_VERSION:-v2}" == "v1" ]]; then
  echo "[ppm10-full] enrich"
  "$PY" tools/enrich_report_evidence.py "$REPORTS_OUT" "$BUNDLE" "$KNOWN_X" || true
  echo "[ppm10-full] patch confidence"
  "$PY" tools/patch_report_confidence.py "$REPORTS_OUT"
  echo "[ppm10-full] patch candidate scores"
  "$PY" tools/patch_candidate_scores.py "$REPORTS_OUT"
  echo "[ppm10-full] patch consensus crowding"
  "$PY" tools/patch_consensus_crowding.py "$REPORTS_OUT"
fi
echo "[ppm10-full] rebuild master index"
"$PY" tools/rebuild_master_index.py "$REPORTS_OUT"

EXEC_OUT="${EXEC_OUT:-$REPO_ROOT/docs/paper/data/exec_summary_ppm10}"
mkdir -p "$EXEC_OUT"
echo "[ppm10-full] exec summary -> $EXEC_OUT"
"$PY" tools/build_exec_summary.py --reports-root "$REPORTS_OUT" --out-dir "$EXEC_OUT"
echo "[ppm10-full] exec summary workbook -> $EXEC_OUT/executive_summary_top10.xlsx"
"$PY" tools/build_exec_summary_workbook.py --exec-dir "$EXEC_OUT"

THRESH="${THRESH:-$REPO_ROOT/docs/paper/data/thresholds}"
mkdir -p "$THRESH"
echo "[ppm10-full] ppm audit (gated) -> ${THRESH}/audit_ppm_qc_gated.*"
"$PY" tools/audit_ppm_qc.py \
  --by-compound "$EXEC_OUT/exec_summary_by_compound.csv" \
  --top10 "$EXEC_OUT/exec_summary_top10.csv" \
  --json-out "$THRESH/audit_ppm_qc_gated.json" \
  --md-out "$THRESH/audit_ppm_qc_gated.md" \
  --csv-out "$THRESH/audit_ppm_qc_gated.csv" \
  --filtered-candidates-out "$THRESH/audit_ppm_qc_gated_filtered_candidates.csv" \
  --filtered-by-compound-out "$THRESH/audit_ppm_qc_gated_filtered_by_compound.csv"

echo "[ppm10-full] done"
