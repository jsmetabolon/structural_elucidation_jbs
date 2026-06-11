#!/usr/bin/env bash
set -euo pipefail

# Fetch large model assets for structural_elucidation_jbs.
# Usage:
#   ./fetch_assets.sh --source ec2 --host user@ec2-host
#   ./fetch_assets.sh --source local --parent /path/to/spec2mol_merged_fullrun_metabolon_only

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOURCE="${SOURCE:-ec2}"
EC2_HOST="${EC2_HOST:-}"
EC2_WORKSPACE="${EC2_WORKSPACE:-/data/josh/spec2mol/workspaces/metabolon}"
PARENT_REPO="${PARENT_REPO:-$(cd "$REPO_ROOT/.." && pwd)}"

usage() {
  cat <<EOF
Usage: $0 [--source ec2|local] [--host user@host] [--parent /path/to/parent/repo]

Downloads:
  - models/encoder/best_encoder.pt
  - models/decoder_weights/*.pt
  - models/embedding_banks/embeddings_clean/
  - vendor/modelfiles/
  - vendor/metabolite_synonyms.sqlite
  - data/wgcna_cache/ (optional)
  - data/ions_unknown_merged_clean.csv
  - data/chrolib_adduct_ladder.csv
  - data/base_pools/
  - data/nce_full_smiles_candidates.csv
  - data/smiles_db/hmdb_minimal.csv
  - data/smiles_db/name_dict_runtime_combined.csv
  - data/smiles_db/comp_id_to_smiles_expanded_950.csv
  - vendor/nce/elucidation/table_builders/*.csv (HMDB tables)
  - data/smiles_db/library_*_fragment_database_element0.csv
  - data/smiles_db/cid_to_parent_mz_*.pkl
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --source) SOURCE="$2"; shift 2 ;;
    --host) EC2_HOST="$2"; shift 2 ;;
    --parent) PARENT_REPO="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown arg: $1" >&2; usage; exit 2 ;;
  esac
done

mkdir -p "$REPO_ROOT/models/encoder" "$REPO_ROOT/models/decoder_weights" \
  "$REPO_ROOT/models/embedding_banks" "$REPO_ROOT/vendor/modelfiles" \
  "$REPO_ROOT/vendor" "$REPO_ROOT/data/wgcna_cache" \
  "$REPO_ROOT/data/base_pools" "$REPO_ROOT/data/observed_ms2_lookups" \
  "$REPO_ROOT/vendor/nce/elucidation/table_builders" \
  "$REPO_ROOT/data/smiles_db"

rsync_one() {
  local src="$1"
  local dst="$2"
  echo "+ rsync $src -> $dst"
  rsync -avz --progress "$src" "$dst"
}

first_existing_glob() {
  local pattern="$1"
  shopt -s nullglob
  local matches=( $pattern )
  shopt -u nullglob
  if (( ${#matches[@]} > 0 )); then
    printf '%s\n' "${matches[0]}"
    return 0
  fi
  return 1
}

rsync_encoder_from_parent() {
  local src=""
  if [[ -n "${PARENT_ENCODER:-}" && -e "$PARENT_ENCODER" ]]; then
    src="$PARENT_ENCODER"
  elif [[ -e "$PARENT_REPO/models/merged_main_clean/best_encoder.pt" ]]; then
    src="$PARENT_REPO/models/merged_main_clean/best_encoder.pt"
  elif [[ -e "$PARENT_REPO/models/x400_positive_v2_clean/encoder/best_encoder.pt" ]]; then
    src="$PARENT_REPO/models/x400_positive_v2_clean/encoder/best_encoder.pt"
  else
    src=$(first_existing_glob "$PARENT_REPO/models/merged_main_*_clean/best_encoder.pt" || true)
  fi
  [[ -z "$src" ]] && return 0
  rsync_one "$src" "$REPO_ROOT/models/encoder/" || \
    rsync_one "$src" "$REPO_ROOT/models/encoder/best_encoder.pt"
}

rsync_embeddings_from_parent() {
  local src=""
  if [[ -n "${PARENT_EMBEDDINGS:-}" && -e "$PARENT_EMBEDDINGS" ]]; then
    src="$PARENT_EMBEDDINGS"
  elif [[ -d "$PARENT_REPO/merged_encoder/data/embeddings_clean" ]]; then
    src="$PARENT_REPO/merged_encoder/data/embeddings_clean"
  else
    src=$(first_existing_glob "$PARENT_REPO/merged_encoder/data/embeddings_*_clean" || true)
  fi
  [[ -z "$src" ]] && return 0
  rsync_one "$src/" "$REPO_ROOT/models/embedding_banks/embeddings_clean/"
}

rsync_ions_from_parent() {
  local src=""
  if [[ -n "${PARENT_IONS:-}" && -e "$PARENT_IONS" ]]; then
    src="$PARENT_IONS"
  elif [[ -e "$PARENT_REPO/merged_encoder/data/ions_unknown_merged_clean.csv" ]]; then
    src="$PARENT_REPO/merged_encoder/data/ions_unknown_merged_clean.csv"
  else
    src=$(first_existing_glob "$PARENT_REPO/merged_encoder/data/ions_unknown_merged*_clean.csv" || true)
  fi
  if [[ -z "$src" && -e "$PARENT_REPO/merged_encoder/data/ions_w_ri_v3_clean.csv" ]]; then
    src="$PARENT_REPO/merged_encoder/data/ions_w_ri_v3_clean.csv"
  fi
  [[ -z "$src" ]] && return 0
  rsync_one "$src" "$REPO_ROOT/data/ions_unknown_merged_clean.csv"
}

rsync_chrolib_from_parent() {
  local src=""
  if [[ -n "${PARENT_CHROLIB:-}" && -e "$PARENT_CHROLIB" ]]; then
    src="$PARENT_CHROLIB"
  elif [[ -e "$PARENT_REPO/data/chrolib_adduct_ladder.csv" ]]; then
    src="$PARENT_REPO/data/chrolib_adduct_ladder.csv"
  else
    src=$(first_existing_glob "$PARENT_REPO/data/chrolib_adduct_ladder*.csv" || true)
  fi
  [[ -z "$src" ]] && return 0
  rsync_one "$src" "$REPO_ROOT/data/chrolib_adduct_ladder.csv"
}

rsync_encoder_from_ec2() {
  local src="${EC2_ENCODER:-$EC2_HOST:$EC2_WORKSPACE/models/merged_main_clean/best_encoder.pt}"
  rsync_one "$src" "$REPO_ROOT/models/encoder/" || \
    rsync_one "$EC2_HOST:$EC2_WORKSPACE/models/x400_positive_v2_clean/encoder/best_encoder.pt" "$REPO_ROOT/models/encoder/best_encoder.pt" || \
    rsync_one "$EC2_HOST:$EC2_WORKSPACE/models/merged_main_*_clean/best_encoder.pt" "$REPO_ROOT/models/encoder/" || true
}

rsync_embeddings_from_ec2() {
  rsync_one "${EC2_EMBEDDINGS:-$EC2_HOST:$EC2_WORKSPACE/merged_encoder/data/embeddings_clean/}" "$REPO_ROOT/models/embedding_banks/embeddings_clean/" || \
    rsync_one "$EC2_HOST:$EC2_WORKSPACE/merged_encoder/data/embeddings_*_clean/" "$REPO_ROOT/models/embedding_banks/embeddings_clean/" || true
}

rsync_ions_from_ec2() {
  rsync_one "${EC2_IONS:-$EC2_HOST:$EC2_WORKSPACE/merged_encoder/data/ions_unknown_merged_clean.csv}" "$REPO_ROOT/data/ions_unknown_merged_clean.csv" || \
    rsync_one "$EC2_HOST:$EC2_WORKSPACE/merged_encoder/data/ions_unknown_merged*_clean.csv" "$REPO_ROOT/data/ions_unknown_merged_clean.csv" || true
}

rsync_chrolib_from_ec2() {
  rsync_one "${EC2_CHROLIB:-$EC2_HOST:$EC2_WORKSPACE/data/chrolib_adduct_ladder.csv}" "$REPO_ROOT/data/chrolib_adduct_ladder.csv" || \
    rsync_one "$EC2_HOST:$EC2_WORKSPACE/data/chrolib_adduct_ladder*.csv" "$REPO_ROOT/data/chrolib_adduct_ladder.csv" || true
}

if [[ "$SOURCE" == "local" ]]; then
  rsync_encoder_from_parent
  rsync_one "$PARENT_REPO/decoder/models/model.pt" "$REPO_ROOT/models/decoder_weights/"
  rsync_one "$PARENT_REPO/decoder/models/model_mass_neg_209_305_best.pt" "$REPO_ROOT/models/decoder_weights/" || true
  rsync_embeddings_from_parent
  rsync_one "$HOME/dev/modelfiles/" "$REPO_ROOT/vendor/modelfiles/"
  rsync_one "$PARENT_REPO/SMILES_DB/library_209_fragment_database_element0.csv" "$REPO_ROOT/data/smiles_db/" || true
  rsync_one "$PARENT_REPO/SMILES_DB/library_400_fragment_database_element0.csv" "$REPO_ROOT/data/smiles_db/" || true
  rsync_one "$PARENT_REPO/SMILES_DB/cid_to_parent_mz_runtime_combined_neg.pkl" "$REPO_ROOT/data/smiles_db/" || true
  rsync_one "$PARENT_REPO/SMILES_DB/cid_to_parent_mz_runtime_combined_pos.pkl" "$REPO_ROOT/data/smiles_db/" || true
  rsync_one "$PARENT_REPO/../spec2mol/Metabolite Synonym Database/metabolite_synonyms_merged_with_mtx_giga.sqlite" "$REPO_ROOT/vendor/metabolite_synonyms.sqlite" || \
    rsync_one "$HOME/dev/spec2mol/Metabolite Synonym Database/metabolite_synonyms_merged_with_mtx_giga.sqlite" "$REPO_ROOT/vendor/metabolite_synonyms.sqlite" || true
  rsync_one "$PARENT_REPO/models/wgcna_inputs_v1/cache/" "$REPO_ROOT/data/wgcna_cache/" || true
  rsync_ions_from_parent
  rsync_chrolib_from_parent
  rsync_one "$PARENT_REPO/data/base_pools/" "$REPO_ROOT/data/base_pools/" || true
  rsync_one "$PARENT_REPO/models/original_nce_current_data_v1_all_libs_20260427/nce_full_smiles_candidates.csv" "$REPO_ROOT/data/" || true
  rsync_one "$PARENT_REPO/SMILES_DB/hmdb_minimal.csv" "$REPO_ROOT/data/smiles_db/" || true
  rsync_one "$PARENT_REPO/SMILES_DB/name_dict_runtime_combined.csv" "$REPO_ROOT/data/smiles_db/" || true
  rsync_one "$PARENT_REPO/SMILES_DB/comp_id_to_smiles_expanded_950.csv" "$REPO_ROOT/data/smiles_db/" || true
  rsync_one "$PARENT_REPO/SMILES_DB/hmdb_minimal.csv" "$REPO_ROOT/vendor/nce/elucidation/table_builders/" || true
  rsync_one "$PARENT_REPO/SMILES_DB/HMDB_36_classyfire_21_annotations.csv" "$REPO_ROOT/vendor/nce/elucidation/table_builders/" || \
    rsync_one "$HOME/dev/NewCompoundElucidation/elucidation/table_builders/HMDB_36_classyfire_21_annotations.csv" "$REPO_ROOT/vendor/nce/elucidation/table_builders/" || true
  rsync_one "$PARENT_REPO/SMILES_DB/hmdb_spreadsheet.csv" "$REPO_ROOT/vendor/nce/elucidation/table_builders/" || true
  rsync_one "$PARENT_REPO/data/observed_ms2_lookups/" "$REPO_ROOT/data/observed_ms2_lookups/" || true
elif [[ "$SOURCE" == "ec2" ]]; then
  if [[ -z "$EC2_HOST" ]]; then
    echo "Set EC2_HOST=user@host for --source ec2" >&2
    exit 2
  fi
  rsync_encoder_from_ec2
  rsync_one "$EC2_HOST:$EC2_WORKSPACE/decoder/models/model.pt" "$REPO_ROOT/models/decoder_weights/"
  rsync_one "$EC2_HOST:$EC2_WORKSPACE/decoder/models/model_mass_neg_209_305_best.pt" "$REPO_ROOT/models/decoder_weights/" || true
  rsync_embeddings_from_ec2
  rsync_one "$EC2_HOST:$EC2_WORKSPACE/../shared/modelfiles/" "$REPO_ROOT/vendor/modelfiles/" || \
    rsync_one "$EC2_HOST:~/dev/modelfiles/" "$REPO_ROOT/vendor/modelfiles/"
  rsync_one "$EC2_HOST:$EC2_WORKSPACE/SMILES_DB/library_209_fragment_database_element0.csv" "$REPO_ROOT/data/smiles_db/" || true
  rsync_one "$EC2_HOST:$EC2_WORKSPACE/SMILES_DB/library_400_fragment_database_element0.csv" "$REPO_ROOT/data/smiles_db/" || true
  rsync_one "$EC2_HOST:$EC2_WORKSPACE/../spec2mol/Metabolite Synonym Database/metabolite_synonyms_merged_with_mtx_giga.sqlite" "$REPO_ROOT/vendor/metabolite_synonyms.sqlite" || true
  rsync_one "$EC2_HOST:$EC2_WORKSPACE/models/wgcna_inputs_v1/cache/" "$REPO_ROOT/data/wgcna_cache/" || true
  rsync_ions_from_ec2
  rsync_chrolib_from_ec2
  rsync_one "$EC2_HOST:$EC2_WORKSPACE/data/base_pools/" "$REPO_ROOT/data/base_pools/" || true
  rsync_one "$EC2_HOST:$EC2_WORKSPACE/models/original_nce_current_data_v1_all_libs_20260427/nce_full_smiles_candidates.csv" "$REPO_ROOT/data/" || true
  rsync_one "$EC2_HOST:$EC2_WORKSPACE/SMILES_DB/hmdb_minimal.csv" "$REPO_ROOT/data/smiles_db/" || true
  rsync_one "$EC2_HOST:$EC2_WORKSPACE/SMILES_DB/name_dict_runtime_combined.csv" "$REPO_ROOT/data/smiles_db/" || true
  rsync_one "$EC2_HOST:$EC2_WORKSPACE/SMILES_DB/comp_id_to_smiles_expanded_950.csv" "$REPO_ROOT/data/smiles_db/" || true
  rsync_one "$EC2_HOST:$EC2_WORKSPACE/SMILES_DB/hmdb_minimal.csv" "$REPO_ROOT/vendor/nce/elucidation/table_builders/" || true
  rsync_one "$EC2_HOST:$EC2_WORKSPACE/SMILES_DB/HMDB_36_classyfire_21_annotations.csv" "$REPO_ROOT/vendor/nce/elucidation/table_builders/" || true
  rsync_one "$EC2_HOST:$EC2_WORKSPACE/SMILES_DB/hmdb_spreadsheet.csv" "$REPO_ROOT/vendor/nce/elucidation/table_builders/" || true
  rsync_one "$EC2_HOST:$EC2_WORKSPACE/data/observed_ms2_lookups/" "$REPO_ROOT/data/observed_ms2_lookups/" || true
else
  echo "Unsupported SOURCE=$SOURCE" >&2
  exit 2
fi

echo "[fetch_assets] done $(date -u +%Y-%m-%dT%H:%M:%SZ)"
