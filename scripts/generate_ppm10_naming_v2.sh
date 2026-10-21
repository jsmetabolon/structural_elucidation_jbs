#!/usr/bin/env bash
set -euo pipefail

# Versioned report delivery. Runtime is roughly 1.5–2.5 hours at JOBS=8.
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"
export NAMING_POLICY_VERSION=v2
export REPORTS_OUT=output/reports_ppm10_full_naming_v2
export EXEC_OUT="$REPO_ROOT/$REPORTS_OUT/exec_summary"
export THRESH="$REPO_ROOT/$REPORTS_OUT/thresholds"
export OVERWRITE=1
./scripts/generate_ppm10_full_reports.sh
"${PY:-/Users/jsimpson/miniconda3/bin/python}" tools/build_naming_v2_release.py \
  --old-root output/reports_ppm10_full --new-root "$REPORTS_OUT"
