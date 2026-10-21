#!/usr/bin/env bash
# Handoff: gated adduct-derive run on EC2. Do not run locally (hours; NCE/GPU).
# Skip analog: not consumed by run_inference.sh (eval assembly only).
# Skip run_full_pipeline_ec2.sh and START_AT=verifier.
#
# Prerequisites on this laptop:
#   python3 tools/inventory_adduct_ladder.py --limit-per-file 100000
#   (2026-09-02: 378/10568 hypotheses move >5 mDa; derive gate PASS)
#
# On a machine that can SSH to EC2:
set -euo pipefail

HOST="${EC2_HOST:-ec2-user@10.12.5.28}"
REMOTE_ROOT="${EC2_JBS:-/data/josh/spec2mol/workspaces/metabolon/spec2mol_v21_inference}"
EC2_PY="${EC2_PY:-/data/josh/spec2mol/workspaces/metabolon/.venv/bin/python}"
LOCAL_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT_ROOT_NAME="${OUT_ROOT_NAME:-output/pipeline_equalization_adduct_derive}"

usage() {
  cat <<'EOF'
Usage:
  ./run_adduct_derive_ec2.sh sync     # rsync this repo to EC2 (excludes output/, .venv, large models)
  ./run_adduct_derive_ec2.sh start    # tmux: START_AT=candidates + derive flag, fresh OUT_ROOT
  ./run_adduct_derive_ec2.sh recover  # tmux: 10 ppm window for compounds that already have inputs
  ./run_adduct_derive_ec2.sh library10 # tmux: 10 ppm re-search of scored compounds still over 10 ppm
  ./run_adduct_derive_ec2.sh denovo   # tmux: SELFIES de novo for compounds with no library hit at 10 ppm
  ./run_adduct_derive_ec2.sh compare  # after job: McNemar vs models/pipeline_equalization baseline
  ./run_adduct_derive_ec2.sh analog   # only if eval assembly still needs analog pools (not inference)

Do not change closest_neutral_mass until compare says rank-1/rank-3 moved.
EOF
}

WS="${EC2_WS:-/data/josh/spec2mol/workspaces/metabolon}"
CLEAN="${EC2_CLEAN:-$WS/models/pipeline_equalization_v21_clean}"

sync_repo() {
  ssh "$HOST" "mkdir -p '$REMOTE_ROOT/data/adducts' '$REMOTE_ROOT/data/feature_rows' '$REMOTE_ROOT/data/integrated_best' '$REMOTE_ROOT/tools' '$REMOTE_ROOT/tests'"
  rsync -az \
    --exclude '.venv/' \
    --exclude 'output/' \
    --exclude 'logs/' \
    --exclude '.git/' \
    --exclude 'models/encoder/' \
    --exclude 'models/decoder_weights/' \
    --exclude 'models/embedding_banks/' \
    --exclude 'vendor/modelfiles/' \
    --exclude 'vendor/metabolite_synonyms.sqlite' \
    "$LOCAL_ROOT/structural_elucidation/" "$HOST:$REMOTE_ROOT/structural_elucidation/"
  rsync -az "$LOCAL_ROOT/data/adducts/" "$HOST:$REMOTE_ROOT/data/adducts/"
  rsync -az "$LOCAL_ROOT/data/feature_rows/" "$HOST:$REMOTE_ROOT/data/feature_rows/"
  rsync -az "$LOCAL_ROOT/data/integrated_best/" "$HOST:$REMOTE_ROOT/data/integrated_best/"
  rsync -az "$LOCAL_ROOT/run_inference.sh" "$LOCAL_ROOT/run_adduct_derive_ec2.sh" \
    "$LOCAL_ROOT/generate_all_reports.sh" "$HOST:$REMOTE_ROOT/"
  rsync -az "$LOCAL_ROOT/scripts/generate_ppm10_full_reports.sh" \
    "$LOCAL_ROOT/scripts/run_ppm10_extra_job.sh" \
    "$LOCAL_ROOT/scripts/run_ppm10_denovo_job.sh" "$HOST:$REMOTE_ROOT/scripts/"
  rsync -az \
    "$LOCAL_ROOT/tools/inventory_adduct_ladder.py" \
    "$LOCAL_ROOT/tools/compare_adduct_derive_rerank.py" \
    "$LOCAL_ROOT/tools/refill_top10_ppm.py" \
    "$LOCAL_ROOT/tools/build_exec_summary.py" \
    "$LOCAL_ROOT/tools/build_exec_summary_workbook.py" \
    "$LOCAL_ROOT/tools/reported_identification.py" \
    "$LOCAL_ROOT/tools/audit_ppm_qc.py" \
    "$LOCAL_ROOT/tools/rebuild_master_index.py" \
    "$LOCAL_ROOT/tools/prepare_ppm10_recovery_inputs.py" \
    "$LOCAL_ROOT/tools/select_ppm10_remaining.py" \
    "$LOCAL_ROOT/tools/prepare_denovo_ppm10_inputs.py" \
    "$LOCAL_ROOT/tools/merge_ppm10_recovery.py" \
    "$HOST:$REMOTE_ROOT/tools/"
  rsync -az "$LOCAL_ROOT/tests/test_adduct_ladder.py" "$HOST:$REMOTE_ROOT/tests/"
  echo "synced code + adducts + feature_rows + integrated_best to $HOST:$REMOTE_ROOT"
  link_inputs
}

link_inputs() {
  ssh "$HOST" bash -s <<EOF
set -euo pipefail
INF='$REMOTE_ROOT'
WS='$WS'
CLEAN='$CLEAN'
mkdir -p "\$INF/data/observed_ms2_lookups" "\$INF/data/base_pools" "\$INF/vendor/nce/elucidation/table_builders" "\$INF/models/x_known_curation_v1"
# pool/verify --include-all-x read this catalog relative to the repo root. Without it
# the pool keeps only the integrated-best controls (~80 compounds instead of ~2,500).
ln -sfn "\$WS/models/x_known_curation_v1/known_x_evidence_catalog.csv" "\$INF/models/x_known_curation_v1/known_x_evidence_catalog.csv"
test -s "\$INF/models/x_known_curation_v1/known_x_evidence_catalog.csv"
ln -sfn "\$WS/known_X.csv" "\$INF/known_X.csv"
ln -sfn "\$WS/known_X.csv" "\$INF/data/known_X.csv"
ln -sfn "\$WS/known_X_all_available_structures_simple.csv" "\$INF/known_X_all_available_structures_simple.csv"
ln -sfn "\$WS/known_X_all_available_structures_simple.csv" "\$INF/data/known_X_all_available_structures_simple.csv"
ln -sfn "\$WS/models/original_nce_current_data_v1_all_libs/nce_full_smiles_candidates.csv" "\$INF/data/nce_full_smiles_candidates.csv"
ln -sfn /data/josh/spec2mol/shared/references/metabolite_synonyms.sqlite "\$INF/vendor/metabolite_synonyms.sqlite"
ln -sfn "\$WS/SMILES_DB/hmdb_minimal.csv" "\$INF/vendor/nce/elucidation/table_builders/hmdb_minimal.csv"
if [[ -f "\$WS/SMILES_DB/HMDB_36_classyfire_21_annotations.csv" ]]; then
  ln -sfn "\$WS/SMILES_DB/HMDB_36_classyfire_21_annotations.csv" "\$INF/vendor/nce/elucidation/table_builders/HMDB_36_classyfire_21_annotations.csv"
fi
if [[ -f "\$WS/SMILES_DB/hmdb_spreadsheet.csv" ]]; then
  ln -sfn "\$WS/SMILES_DB/hmdb_spreadsheet.csv" "\$INF/vendor/nce/elucidation/table_builders/hmdb_spreadsheet.csv"
fi
for spec in 209:negative 305:negative 400:positive 402:positive; do
  lib="\${spec%%:*}"
  ion="\${spec##*:}"
  ln -sfn "\$CLEAN/lib\${lib}_\${ion}_observed_ms2_lookup.csv" "\$INF/data/observed_ms2_lookups/lib\${lib}_\${ion}_observed_ms2_lookup.csv"
  ln -sfn "\$CLEAN/unified_pool_lib\${lib}/unified_candidates_per_control.csv" "\$INF/data/base_pools/lib\${lib}_\${ion}.csv"
done
test -s "\$INF/data/feature_rows/lib209_negative_gap_manifest.csv"
test -s "\$INF/data/adducts/adduct_ladder.csv"
test -s "\$INF/vendor/nce/elucidation/table_builders/hmdb_minimal.csv"
echo "linked EC2 inputs under \$INF"
EOF
}

start_derive() {
  # Stored neutral mass. Adduct-derived mass is only the fallback job for
  # compounds with no ion-supported candidate within 10 ppm.
  ssh "$HOST" "cat > '$REMOTE_ROOT/run_adduct_derive_job.sh'" <<EOF
#!/usr/bin/env bash
set -uo pipefail
cd '$REMOTE_ROOT'
export PY='$EC2_PY'
export PYTHONPATH='$REMOTE_ROOT'
test -x "\$PY" || { echo "missing \$PY"; exit 2; }
test -f data/adducts/adduct_ladder.csv
export OUT_ROOT='$OUT_ROOT_NAME'
export START_AT=candidates
export RUN_LIBS=209,305,400,402
export VERIFY_EXTRA_ARGS=--ignore-acceptance
export JOBS=4
./run_inference.sh
status=\$?
echo EXIT=\$status
if [[ \$status -eq 0 ]]; then
  "\$PY" tools/refill_top10_ppm.py \\
    --equalization-root '$OUT_ROOT_NAME' \\
    --libs 209,305,400,402 \\
    --top-k 10 \\
    --ppm-ceiling 10 \\
    --summary-json '$OUT_ROOT_NAME/ppm10_refill_summary.json'
  OUT_ROOT='$OUT_ROOT_NAME' \\
    PRED_NAME=unified_rerank_all_top10_ppm10.csv \\
    REPORTS_OUT='$OUT_ROOT_NAME/reports' \\
    EXEC_OUT='$OUT_ROOT_NAME/exec_summary' \\
    CANDIDATE_TAB_N=10 \\
    RUN_LIBS=209,305,400,402 \\
    JOBS=4 \\
    ./scripts/generate_ppm10_full_reports.sh
  echo GATHER_EXIT=\$?
fi
exec bash
EOF
  ssh "$HOST" bash -s <<EOF
set -euo pipefail
cd '$REMOTE_ROOT'
chmod +x run_adduct_derive_job.sh
test -x '$EC2_PY' || { echo "missing $EC2_PY"; exit 2; }
test -f data/adducts/adduct_ladder.csv
test -s models/x_known_curation_v1/known_x_evidence_catalog.csv || { echo "missing all-X catalog; run link first" >&2; exit 2; }
test -s data/known_X.csv || { echo "missing data/known_X.csv; run link first" >&2; exit 2; }
test -s vendor/metabolite_synonyms.sqlite || { echo "missing vendor/metabolite_synonyms.sqlite; run link first" >&2; exit 2; }
if [[ -d $OUT_ROOT_NAME ]]; then
  if find $OUT_ROOT_NAME -name 'unified_rerank_all_top10.csv' | grep -q .; then
    echo "refusing to overwrite non-empty $OUT_ROOT_NAME" >&2
    exit 2
  fi
  rm -rf $OUT_ROOT_NAME
fi
if tmux has-session -t adduct_derive 2>/dev/null; then
  tmux kill-session -t adduct_derive
fi
tmux new-session -d -s adduct_derive '$REMOTE_ROOT/run_adduct_derive_job.sh'
echo "started tmux session adduct_derive. attach: ssh $HOST -t tmux attach -t adduct_derive"
EOF
}

# Launch scripts/run_ppm10_extra_job.sh in tmux after a selector writes its inputs.
# Args: <session> <inputs_dir> <extra_root> <selector command...>
launch_ppm10_extra() {
  local session="$1" inputs_dir="$2" extra_root="$3"
  shift 3
  local selector="$*"
  ssh "$HOST" "cat > '$REMOTE_ROOT/run_${session}_job.sh'" <<EOF
#!/usr/bin/env bash
set -uo pipefail
cd '$REMOTE_ROOT'
export PY='$EC2_PY'
export PYTHONPATH='$REMOTE_ROOT'
test -x "\$PY" || { echo "missing \$PY"; exit 2; }
"\$PY" $selector || { echo "selector failed"; exec bash; }
MAIN_ROOT='$OUT_ROOT_NAME' ./scripts/run_ppm10_extra_job.sh '$inputs_dir' '$extra_root'
exec bash
EOF
  ssh "$HOST" bash -s <<EOF
set -euo pipefail
cd '$REMOTE_ROOT'
chmod +x run_${session}_job.sh scripts/run_ppm10_extra_job.sh
test -s models/x_known_curation_v1/known_x_evidence_catalog.csv
test -d '$WS/models/pipeline_equalization_v1_recovered_ms2'
if [[ -d $extra_root ]] && find $extra_root -name 'unified_rerank_all_top10.csv' | grep -q .; then
  echo "refusing to overwrite non-empty $extra_root" >&2
  exit 2
fi
rm -rf $extra_root
if tmux has-session -t $session 2>/dev/null; then
  tmux kill-session -t $session
fi
tmux new-session -d -s $session '$REMOTE_ROOT/run_${session}_job.sh'
echo "started tmux session $session. attach: ssh $HOST -t tmux attach -t $session"
EOF
}

start_recover() {
  launch_ppm10_extra ppm10_recover output/ppm10_recovery_inputs \
    output/pipeline_equalization_adduct_derive_recover \
    tools/prepare_ppm10_recovery_inputs.py \
    --catalog models/x_known_curation_v1/known_x_evidence_catalog.csv \
    --scored-root "$OUT_ROOT_NAME" \
    --feature-dir data/feature_rows \
    --recovered-dir "$WS/models/pipeline_equalization_v1_recovered_ms2" \
    --out-dir output/ppm10_recovery_inputs
}

start_library10() {
  launch_ppm10_extra ppm10_library10 output/ppm10_remaining/library10_inputs \
    output/pipeline_equalization_adduct_derive_library10 \
    tools/select_ppm10_remaining.py \
    --catalog models/x_known_curation_v1/known_x_evidence_catalog.csv \
    --main-root "$OUT_ROOT_NAME" \
    --recovery-root output/pipeline_equalization_adduct_derive_recover \
    --recovery-selection output/ppm10_recovery_inputs/selection.json \
    --summary "$OUT_ROOT_NAME/exec_summary/exec_summary_by_compound.csv" \
    --feature-dir data/feature_rows \
    --recovered-dir "$WS/models/pipeline_equalization_v1_recovered_ms2" \
    --out-dir output/ppm10_remaining
}

start_denovo() {
  ssh "$HOST" bash -s <<EOF
set -euo pipefail
cd '$REMOTE_ROOT'
chmod +x scripts/run_ppm10_denovo_job.sh scripts/run_ppm10_extra_job.sh
test -s output/ppm10_remaining/nohit10_ids.txt || { echo "run library10 first" >&2; exit 2; }
test -s output/ppm10_recovery_inputs/selection.json
test -s '$WS/generate_selfies_denovo.py'
extra=output/pipeline_equalization_adduct_derive_denovo
if [[ -d \$extra ]] && find \$extra -name 'unified_rerank_all_top10.csv' | grep -q .; then
  echo "refusing to overwrite non-empty \$extra" >&2
  exit 2
fi
rm -rf \$extra
if tmux has-session -t ppm10_denovo 2>/dev/null; then
  tmux kill-session -t ppm10_denovo
fi
tmux new-session -d -s ppm10_denovo "cd '$REMOTE_ROOT' && PY='$EC2_PY' WS='$WS' MAIN_ROOT='$OUT_ROOT_NAME' ./scripts/run_ppm10_denovo_job.sh; exec bash"
echo "started tmux session ppm10_denovo. attach: ssh $HOST -t tmux attach -t ppm10_denovo"
EOF
}

start_analog() {
  echo "Analog is NOT on run_inference.sh. Only run this if evaluation.assembly still reads analog pools."
  ssh "$HOST" bash -s <<EOF
set -euo pipefail
cd '$REMOTE_ROOT'
tmux has-session -t analog_ladder 2>/dev/null && { echo 'tmux session analog_ladder already exists'; exit 2; }
tmux new-session -d -s analog_ladder "cd '$REMOTE_ROOT' && export PY='$EC2_PY' PYTHONPATH='$REMOTE_ROOT' && \\
  \$PY -m structural_elucidation.evaluation.analog \\
    --output-root models/denovo_x_domain_analog_ladder_v1; echo EXIT=\$?; exec bash"
echo "started tmux analog_ladder"
EOF
}

run_compare() {
  ssh "$HOST" bash -s <<EOF
set -euo pipefail
cd '$REMOTE_ROOT'
export PYTHONPATH='$REMOTE_ROOT'
'$EC2_PY' tools/compare_adduct_derive_rerank.py \\
  --baseline /data/josh/spec2mol/workspaces/metabolon/models/pipeline_equalization_v21_clean \\
  --derive $OUT_ROOT_NAME \\
  --json-out $OUT_ROOT_NAME/mcnemar_vs_v21_clean.json
EOF
}

cmd="${1:-}"
case "$cmd" in
  sync) sync_repo ;;
  start) start_derive ;;
  recover) start_recover ;;
  library10) start_library10 ;;
  denovo) start_denovo ;;
  analog) start_analog ;;
  compare) run_compare ;;
  link) link_inputs ;;
  *) usage; exit 2 ;;
esac
