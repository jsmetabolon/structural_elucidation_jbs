#!/usr/bin/env bash
# Cohort the unpublished 551 and score a 20-compound same-formula isomer pilot.
# Runs on ec2-user@10.12.5.28. Does not touch the merged delivery workbook.
set -euo pipefail

HOST="${HOST:-ec2-user@10.12.5.28}"
REPO="${REPO:-/data/josh/spec2mol/workspaces/metabolon/spec2mol_v21_inference}"
PY="${PY:-/data/josh/spec2mol/workspaces/metabolon/.venv/bin/python}"
STORED_ROOT="${STORED_ROOT:-/data/josh/spec2mol/workspaces/metabolon/models/pipeline_equalization_v21}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
REMOTE="$REPO/output/pipeline_equalization_spectrum_pilot"
SSH=(ssh -o ConnectTimeout=20 -o ConnectionAttempts=1 -o StrictHostKeyChecking=accept-new)

echo "[wait] looking for $HOST"
ok=0
for attempt in $(seq 1 8); do
  if "${SSH[@]}" "$HOST" 'echo up'; then
    ok=1
    break
  fi
  echo "[wait] attempt $attempt failed"
  sleep 10
done
if [[ "$ok" -ne 1 ]]; then
  echo "[wait] GIVE_UP host unreachable"
  exit 2
fi

rsync -a -e "ssh -o ConnectTimeout=20 -o ConnectionAttempts=1" \
  "$ROOT/tools/spectrum_pilot.py" \
  "$ROOT/tools/refill_top10_ppm.py" \
  "$HOST:$REPO/tools/"

"${SSH[@]}" "$HOST" "mkdir -p '$REMOTE/inputs'"
rsync -a -e "ssh -o ConnectTimeout=20 -o ConnectionAttempts=1" \
  "$ROOT/output/pipeline_equalization_merged_delivery/exec_summary_by_compound.csv" \
  "$ROOT/output/pipeline_equalization_stored_ppm10/fallback_ids.txt" \
  "$ROOT/output/_merge_inputs/selection.json" \
  "$ROOT/models/pipeline_equalization/lib209_negative_observed_ms2_lookup.csv" \
  "$ROOT/models/pipeline_equalization/lib305_negative_observed_ms2_lookup.csv" \
  "$ROOT/models/pipeline_equalization/lib400_positive_observed_ms2_lookup.csv" \
  "$ROOT/models/pipeline_equalization/lib402_positive_observed_ms2_lookup.csv" \
  "$HOST:$REMOTE/inputs/"

# Cohort is a full-pool scan. Pilot scoring follows in the same tmux session
# and stops after the 20-compound pilot so the scores can be read first.
"${SSH[@]}" "$HOST" "tmux has-session -t spectrum_pilot 2>/dev/null || tmux new-session -d -s spectrum_pilot \
  'cd \"$REPO\" && export PYTHONPATH=\"$REPO\" && \
   $PY tools/spectrum_pilot.py cohort \
     --by-compound \"$REMOTE/inputs/exec_summary_by_compound.csv\" \
     --lookup-dir \"$REMOTE/inputs\" \
     --stored-root \"$STORED_ROOT\" \
     --fallback-root \"$REPO/output/pipeline_equalization_adduct_fallback\" \
     --fallback-ids \"$REMOTE/inputs/fallback_ids.txt\" \
     --fallback-selection \"$REMOTE/inputs/selection.json\" \
     --out \"$REMOTE/cohort.csv\" && \
   $PY tools/spectrum_pilot.py pilot \
     --cohort \"$REMOTE/cohort.csv\" \
     --out \"$REMOTE/pilot_scores.csv\" \
     --n 20 --min-peaks 3 --isomer-limit 50; \
   echo PILOT_EXIT:\$?'"
echo "[wait] TMUX_STARTED spectrum_pilot"
