#!/usr/bin/env bash
set -euo pipefail

# Stage pipeline data assets from plain HTTP(S) GET links.
#
# This is the GET-link counterpart of fetch_assets.sh (which needs rsync/SSH).
# Each asset has a stable NAME; provide its URL either as an environment
# variable SE_GET_<NAME> or as a NAME=URL line in a links file (default:
# ./get_links.env, which is gitignored because presigned URLs are secrets).
#
# Usage:
#   ./fetch_get_assets.sh --verify              # report present/missing, no downloads
#   ./fetch_get_assets.sh                       # download missing assets with known URLs
#   ./fetch_get_assets.sh --force               # redownload even if present
#   ./fetch_get_assets.sh --links-file my.env   # use an alternate links file
#
# Asset names (REQUIRED for the core candidates->verify->pool->rerank path):
#   OBS_LOOKUP_209, OBS_LOOKUP_305, OBS_LOOKUP_400, OBS_LOOKUP_402,
#   OBS_PER_CONTROL_ALL_X, BASE_POOL_209, BASE_POOL_305, BASE_POOL_400,
#   BASE_POOL_402, MODELFILES_TGZ, NCE_FULL_SMILES, HMDB_MINIMAL,
#   HMDB_CLASSYFIRE, HMDB_SPREADSHEET, NAME_DICT, COMP_ID_TO_SMILES_950,
#   FRAG_DB_209, FRAG_DB_400, CID_TO_PARENT_MZ_NEG, CID_TO_PARENT_MZ_POS
# Optional (enrichment/reports): METABOLITE_SQLITE, WGCNA_CACHE_TGZ
# Optional (torch inference lane): ENCODER_PT, DECODER_PT, DECODER_NEG_PT,
#   EMBEDDINGS_TGZ, IONS_MERGED, CHROLIB_LADDER

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LINKS_FILE="$REPO_ROOT/get_links.env"
VERIFY_ONLY=0
FORCE=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --verify) VERIFY_ONLY=1; shift ;;
    --force) FORCE=1; shift ;;
    --links-file) LINKS_FILE="$2"; shift 2 ;;
    -h|--help) grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) echo "Unknown arg: $1" >&2; exit 2 ;;
  esac
done

if [[ -f "$LINKS_FILE" ]]; then
  # shellcheck disable=SC1090
  set -a; source "$LINKS_FILE"; set +a
fi

# NAME|KIND|REQUIRED|DEST[|EXTRA_DEST]
# KIND: file (plain download) or tgz (download + extract into DEST directory).
ASSETS=(
  "OBS_LOOKUP_209|file|required|data/observed_ms2_lookups/lib209_negative_observed_ms2_lookup.csv"
  "OBS_LOOKUP_305|file|required|data/observed_ms2_lookups/lib305_negative_observed_ms2_lookup.csv"
  "OBS_LOOKUP_400|file|required|data/observed_ms2_lookups/lib400_positive_observed_ms2_lookup.csv"
  "OBS_LOOKUP_402|file|required|data/observed_ms2_lookups/lib402_positive_observed_ms2_lookup.csv"
  "OBS_PER_CONTROL_ALL_X|file|required|data/observed_ms2_lookups/per_control_all_x.csv"
  "BASE_POOL_209|file|required|data/base_pools/lib209_negative.csv"
  "BASE_POOL_305|file|required|data/base_pools/lib305_negative.csv"
  "BASE_POOL_400|file|required|data/base_pools/lib400_positive.csv"
  "BASE_POOL_402|file|required|data/base_pools/lib402_positive.csv"
  "MODELFILES_TGZ|tgz|required|vendor/modelfiles"
  "NCE_FULL_SMILES|file|required|data/nce_full_smiles_candidates.csv"
  "HMDB_MINIMAL|file|required|vendor/nce/elucidation/table_builders/hmdb_minimal.csv|data/smiles_db/hmdb_minimal.csv"
  "HMDB_CLASSYFIRE|file|required|vendor/nce/elucidation/table_builders/HMDB_36_classyfire_21_annotations.csv"
  "HMDB_SPREADSHEET|file|required|vendor/nce/elucidation/table_builders/hmdb_spreadsheet.csv"
  "NAME_DICT|file|required|data/smiles_db/name_dict_runtime_combined.csv"
  "COMP_ID_TO_SMILES_950|file|required|data/smiles_db/comp_id_to_smiles_expanded_950.csv"
  "FRAG_DB_209|file|required|data/smiles_db/library_209_fragment_database_element0.csv"
  "FRAG_DB_400|file|required|data/smiles_db/library_400_fragment_database_element0.csv"
  "CID_TO_PARENT_MZ_NEG|file|required|data/smiles_db/cid_to_parent_mz_runtime_combined_neg.pkl"
  "CID_TO_PARENT_MZ_POS|file|required|data/smiles_db/cid_to_parent_mz_runtime_combined_pos.pkl"
  "METABOLITE_SQLITE|file|optional|vendor/metabolite_synonyms.sqlite"
  "WGCNA_CACHE_TGZ|tgz|optional|data/wgcna_cache"
  "ENCODER_PT|file|optional|models/encoder/best_encoder.pt"
  "DECODER_PT|file|optional|models/decoder_weights/model.pt"
  "DECODER_NEG_PT|file|optional|models/decoder_weights/model_mass_neg_209_305_best.pt"
  "EMBEDDINGS_TGZ|tgz|optional|models/embedding_banks/embeddings_clean"
  "IONS_MERGED|file|optional|data/ions_unknown_merged_clean.csv"
  "CHROLIB_LADDER|file|optional|data/chrolib_adduct_ladder.csv"
)

asset_url() {
  local name="$1"
  local var="SE_GET_${name}"
  printf '%s' "${!var:-${!name:-}}"
}

present() {
  local kind="$1" dest="$2"
  if [[ "$kind" == "tgz" ]]; then
    [[ -d "$REPO_ROOT/$dest" ]] && [[ -n "$(ls -A "$REPO_ROOT/$dest" 2>/dev/null)" ]]
  else
    [[ -s "$REPO_ROOT/$dest" ]]
  fi
}

download() {
  local url="$1" out="$2"
  mkdir -p "$(dirname "$out")"
  curl -fL --retry 4 --retry-delay 4 -o "$out.part" "$url"
  mv "$out.part" "$out"
}

missing_required=0
downloaded=0
for row in "${ASSETS[@]}"; do
  IFS='|' read -r name kind requirement dest extra_dest <<<"$row"
  url="$(asset_url "$name")"
  if present "$kind" "$dest" && [[ "$FORCE" -eq 0 ]]; then
    echo "[present ] $name -> $dest"
    if [[ -n "${extra_dest:-}" && ! -s "$REPO_ROOT/$extra_dest" ]]; then
      mkdir -p "$(dirname "$REPO_ROOT/$extra_dest")"
      cp "$REPO_ROOT/$dest" "$REPO_ROOT/$extra_dest"
      echo "[copied  ] $name -> $extra_dest"
    fi
    continue
  fi
  if [[ "$VERIFY_ONLY" -eq 1 || -z "$url" ]]; then
    if [[ "$requirement" == "required" ]]; then
      echo "[MISSING ] $name ($requirement) -> $dest ${url:+(url known)}"
      missing_required=$((missing_required + 1))
    else
      echo "[missing ] $name ($requirement) -> $dest ${url:+(url known)}"
    fi
    continue
  fi
  echo "[fetch   ] $name -> $dest"
  if [[ "$kind" == "tgz" ]]; then
    tmp="$(mktemp --suffix=.tgz)"
    curl -fL --retry 4 --retry-delay 4 -o "$tmp" "$url"
    mkdir -p "$REPO_ROOT/$dest"
    tar -xzf "$tmp" -C "$REPO_ROOT/$dest" --strip-components=0
    rm -f "$tmp"
  else
    download "$url" "$REPO_ROOT/$dest"
    if [[ -n "${extra_dest:-}" ]]; then
      mkdir -p "$(dirname "$REPO_ROOT/$extra_dest")"
      cp "$REPO_ROOT/$dest" "$REPO_ROOT/$extra_dest"
    fi
  fi
  downloaded=$((downloaded + 1))
done

echo
echo "[fetch_get_assets] downloaded=$downloaded missing_required=$missing_required"
if [[ "$missing_required" -gt 0 ]]; then
  echo "[fetch_get_assets] provide URLs via SE_GET_<NAME> env vars or $LINKS_FILE (NAME=URL per line)"
  exit 1
fi
echo "[fetch_get_assets] all required assets staged"
