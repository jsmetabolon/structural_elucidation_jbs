#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"

PY=${PY:-"$ROOT/.venv/bin/python"}
RUN_ID=${RUN_ID:-$(date -u +%Y%m%d_%H%M%S)}
LOG_DIR=${LOG_DIR:-"$ROOT/logs/full_pipeline_v21_clean_${RUN_ID}"}
OUT_ROOT=${OUT_ROOT:-models/pipeline_equalization_v21_clean}
ADDUCT_OUT=${ADDUCT_OUT:-models/upstream_adduct_evidence_v21_clean}
PEPTIDE_OUT=${PEPTIDE_OUT:-models/peptide_first_denovo_v21_clean}
RUN_LIBS=${RUN_LIBS:-209,305,400,402}
ENABLE_PATHWAY_ENRICHMENT=${ENABLE_PATHWAY_ENRICHMENT:-1}
PATHWAY_REFERENCE_DB=${PATHWAY_REFERENCE_DB:-"Metabolite Synonym Database/metabolite_synonyms_merged_with_mtx_giga.sqlite"}

KNOWN_CSV=${KNOWN_CSV:-merged_encoder/data/ions_w_ri_v3_clean.csv}
KNOWN_OUT=${KNOWN_OUT:-merged_encoder/data/ions_merged_v21_clean.csv}
UNKNOWN_OUT=${UNKNOWN_OUT:-merged_encoder/data/ions_unknown_merged_v21_clean.csv}

mkdir -p "$LOG_DIR"
exec > >(tee -a "$LOG_DIR/run.log") 2>&1

echo "[v21-clean] run_id=$RUN_ID"
echo "[v21-clean] root=$ROOT"
echo "[v21-clean] log_dir=$LOG_DIR"

require_file() {
  local path="$1"
  if [[ ! -s "$path" ]]; then
    echo "[v21-clean] missing required file: $path" >&2
    exit 2
  fi
}

require_file "$PY"
require_file build_recovered_ms2_from_chrolib_v3.py
require_file run_full_merged_pipeline_metabolon.sh
require_file run_pipeline_equalization_ec2_v3.sh
if [[ "$ENABLE_PATHWAY_ENRICHMENT" == "1" ]]; then
  require_file build_pathway_enrichment_v1.py
  require_file "$PATHWAY_REFERENCE_DB"
fi
require_file "$KNOWN_CSV"
require_file "$UNKNOWN_OUT"
require_file ions_with_ri.csv
require_file SMILES_DB/comp_id_to_smiles_209.csv
require_file SMILES_DB/fasstrecords_prevalence.csv
require_file SMILES_DB/fasstrecords_prevalence_pos.csv
require_file decoder/models/model.pt
require_file decoder/models/selfies_vocab_expanded.json
require_file models/observed_ms2_lookup_v8/per_control_all_x.csv

for lib in 209 305 400 402; do
  ion=positive
  [[ "$lib" == "209" || "$lib" == "305" ]] && ion=negative
  require_file "models/unified_pool_v20_lib${lib}_${ion}/unified_candidates_per_control.csv"
done
require_file models/integrated_denovo_audit_v15/lib_209_negative/integrated_per_control_best.csv
require_file models/integrated_denovo_audit_v19b/lib_305_negative/integrated_per_control_best.csv
require_file models/integrated_denovo_audit_v19b/lib_400_positive/integrated_per_control_best.csv
require_file models/integrated_denovo_audit_v19b/lib_402_positive/integrated_per_control_best.csv

if [[ -d "$OUT_ROOT" ]] && [[ -n "$(find "$OUT_ROOT" -mindepth 1 -maxdepth 2 -type f 2>/dev/null | head -n 1)" ]]; then
  echo "[v21-clean] refusing to overwrite non-empty output root: $OUT_ROOT" >&2
  echo "[v21-clean] set OUT_ROOT to a fresh directory to rerun." >&2
  exit 2
fi

if [[ "${PREFLIGHT_ONLY:-0}" == "1" ]]; then
  echo "[v21-clean] preflight complete"
  exit 0
fi

echo "[v21-clean] phase0: regenerate known and unknown ions"
"$PY" build_recovered_ms2_from_chrolib_v3.py \
  --chrolib-csv "$KNOWN_CSV" \
  --known-out "$KNOWN_OUT" \
  --unknown-out "$UNKNOWN_OUT" \
  --out-raw merged_encoder/data/ions_full_v21_clean_raw.csv \
  --out-trimmed merged_encoder/data/ions_full_v21_clean_trimmed.csv \
  --audit-out merged_encoder/data/ions_full_v21_clean_audit.csv \
  --summary-json merged_encoder/data/ions_full_v21_clean_summary.json \
  --chrolib-source-version v21_clean \
  --force

"$PY" - <<'PY'
import pandas as pd
from pathlib import Path

known = Path("merged_encoder/data/ions_merged_v21_clean.csv")
unknown = Path("merged_encoder/data/ions_unknown_merged_v21_clean.csv")
for path in [known, unknown]:
    rows = sum(1 for _ in path.open()) - 1
    print(f"[v21-clean] {path}: {rows} data rows")
PY

echo "[v21-clean] rebuild references/ pickles from v21_clean CSVs"
"$PY" - <<PY
from pathlib import Path
import pandas as pd

run_id = "${RUN_ID}"
refs = Path("references")
refs.mkdir(parents=True, exist_ok=True)

drop_cols = ["lib_id", "ion_mode", "ion_mode_label", "dataset_split"]

def backup(path: Path) -> None:
    if path.exists() or path.is_symlink():
        target = path.with_name(f"{path.name}.{run_id}.bak")
        if target.exists():
            raise FileExistsError(target)
        path.rename(target)
        print(f"[v21-clean] backed up {path} -> {target}")

def write_by_lib(csv_path: str, template: str, label: str) -> None:
    ions = pd.read_csv(csv_path, low_memory=False)
    if "lib_id" not in ions.columns:
        raise KeyError(f"{csv_path} missing lib_id")
    ions["lib_id"] = pd.to_numeric(ions["lib_id"], errors="raise").astype(int)
    for lib_id in [209, 305, 400, 402]:
        subset = ions.loc[ions["lib_id"].eq(lib_id)].copy()
        if subset.empty:
            raise ValueError(f"{csv_path} has zero rows for lib {lib_id}")
        subset = subset.drop(columns=drop_cols, errors="ignore")
        path = refs / template.format(lib_id=lib_id)
        backup(path)
        subset.to_pickle(path)
        n_compounds = subset["comp_id"].nunique() if "comp_id" in subset.columns else 0
        print(f"[v21-clean] {label} {path}: {len(subset)} rows, {n_compounds} compounds")

write_by_lib("merged_encoder/data/ions_merged_v21_clean.csv", "ions_all_{lib_id}.pkl", "known")
write_by_lib("merged_encoder/data/ions_unknown_merged_v21_clean.csv", "ions_unknown_or_x_{lib_id}_new.pkl", "unknown")
PY

echo "[v21-clean] phase1: full merged metabolon pipeline"
bash run_full_merged_pipeline_metabolon.sh

require_file models/merged_main_excl_x_unresolved_2026-04-17/best_encoder.pt
if [[ ! -d embeddings_latent_from_smiles_runtime_metabolon ]]; then
  echo "[v21-clean] missing embedding bank: embeddings_latent_from_smiles_runtime_metabolon" >&2
  exit 2
fi

echo "[v21-clean] phase2: pipeline equalization"
OUT_ROOT="$OUT_ROOT" \
MERGED_IONS_CSV="$UNKNOWN_OUT" \
FEATURE_TABLE="$UNKNOWN_OUT" \
RUN_LIBS="$RUN_LIBS" \
START_AT=full \
ADDUCT_OUT="$ADDUCT_OUT" \
PEPTIDE_OUT="$PEPTIDE_OUT" \
bash run_pipeline_equalization_ec2_v3.sh

echo "[v21-clean] post-run verification"
for lib in 209 305 400 402; do
  path="$OUT_ROOT/rerank_lib${lib}/unified_rerank_all_top10.csv"
  require_file "$path"
  wc -l "$path"
done

if [[ "$ENABLE_PATHWAY_ENRICHMENT" == "1" ]]; then
  echo "[v21-clean] phase3: pathway enrichment"
  "$PY" build_pathway_enrichment_v1.py \
    --equalization-root "$OUT_ROOT" \
    --reference-db "$PATHWAY_REFERENCE_DB" \
    --libs "$RUN_LIBS"

  for lib in 209 305 400 402; do
    path="$OUT_ROOT/rerank_lib${lib}/unified_rerank_all_top10_with_pathways.csv"
    require_file "$path"
    wc -l "$path"
  done
else
  echo "[v21-clean] phase3: pathway enrichment disabled"
fi

find "$OUT_ROOT" -maxdepth 3 -type f \
  \( -path "*/rerank_lib*/unified_rerank_all_top10.csv" \
  -o -path "*/rerank_lib*/unified_rerank_all_top10_with_pathways.csv" \
  -o -path "*/pathway_enrichment_summary.json" \
  -o -path "*/verifier_lib*/scored_candidates.csv" \
  -o -path "*/unified_pool_lib*/unified_candidates_per_control.csv" \) \
  -print | sort > "$LOG_DIR/final_outputs.txt"

echo "[v21-clean] completed $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "[v21-clean] final outputs listed in $LOG_DIR/final_outputs.txt"
