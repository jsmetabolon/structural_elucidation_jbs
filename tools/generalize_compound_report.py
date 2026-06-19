#!/usr/bin/env python3
"""Post-process copied X-25003 report builder into generic build_compound_report.py.

Rewrites compound-specific identifiers (X-25003, lactoyl anchors) into
parameterized ``control_id`` / ``target_name`` fields and injects a generic
``build_report`` CLI. Run once after copying the legacy report from parent repo.
"""

from __future__ import annotations

import re
from pathlib import Path

# ── Target file ───────────────────────────────────────────────────────────────

TARGET = Path(__file__).resolve().parents[1] / "src" / "build_compound_report.py"

BUILD_REPORT = '''
def control_slug(control_id: str) -> str:
    return re.sub(r"[^a-z0-9]+", "-", control_id.lower()).strip("-") or "compound"


def predictions_to_candidate_pool(predictions: pd.DataFrame, control_id: str) -> pd.DataFrame:
    work = predictions.copy()
    id_col = "control_id" if "control_id" in work.columns else "x_compound_id"
    work = work[work[id_col].astype(str).eq(control_id)].copy()
    if work.empty:
        return work
    smiles_col = next(
        (c for c in ("predicted_smiles", "smiles", "candidate_smiles", "canonical_2d_smiles") if c in work.columns),
        None,
    )
    if smiles_col is None:
        return pd.DataFrame()
    rank_col = "rank" if "rank" in work.columns else "_rank" if "_rank" in work.columns else None
    if rank_col:
        work = work.sort_values(rank_col, kind="mergesort")
        work["rank"] = pd.to_numeric(work[rank_col], errors="coerce")
    else:
        work["rank"] = range(1, len(work) + 1)
    work["control_id"] = control_id
    work["predicted_smiles"] = work[smiles_col].astype(str)
    if "predicted_formula" not in work.columns and "candidate_formula" in work.columns:
        work["predicted_formula"] = work["candidate_formula"]
    if "candidate_name" not in work.columns and "base_name" in work.columns:
        work["candidate_name"] = work["base_name"]
    score_col = next(
        (c for c in ("unified_final_score", "rerank_score_calibrated", "rerank_score", "score_rank") if c in work.columns),
        None,
    )
    if score_col:
        work["candidate_report_score"] = pd.to_numeric(work[score_col], errors="coerce")
    return work


def lookup_known_x_row(known_x: pd.DataFrame, control_id: str) -> dict[str, Any]:
    unk = known_x.get("Unknown", pd.Series("", index=known_x.index)).astype(str)
    rows = known_x[unk.eq(control_id)]
    return rows.iloc[0].to_dict() if not rows.empty else {}


def resolve_chemical_id(db_path: Path, known_x_row: dict[str, Any], target_name: str) -> str:
    if known_x_row.get("KNOWN"):
        return str(known_x_row["KNOWN"]).strip()
    if not db_path.exists():
        return ""
    df = query_df(
        db_path,
        "SELECT chemical_id FROM metabolites WHERE lower(canonical_name) = lower(?) LIMIT 1",
        (target_name,),
    )
    if not df.empty:
        return str(df.iloc[0]["chemical_id"])
    return ""


def load_real_correlation_rows(
    control_id: str,
    chemical_id: str,
    target_name: str,
) -> tuple[list[dict[str, Any]], str | None]:
    anchor_names = {control_id.lower(), target_name.lower()}
    if target_name:
        anchor_names.add(target_name.strip().lower())
    for cohort, spec in wgcna.COHORT_INPUTS.items():
        modules = wgcna.load_module_assignments(cohort)
        if modules.empty:
            continue
        matches = modules[
            modules["plot_name"].astype(str).str.strip().str.lower().isin(anchor_names)
            | modules["plot_name"].astype(str).eq(control_id)
        ]
        if chemical_id:
            matches = matches | modules["chem_id"].astype(str).eq(str(chemical_id))
        if not matches.empty:
            break
    blocker = (
        f"Real sample-abundance correlation rows for {control_id} are not available in the configured WGCNA exports."
    )
    return [], blocker


def build_placeholder_correlation_rows(
    extended_unknown_rows: list[dict[str, Any]],
    *,
    exclude_control_id: str,
    limit: int = 5,
) -> list[dict[str, Any]]:
    proxy_rows: list[dict[str, Any]] = []
    for row in extended_unknown_rows:
        x_id = str(row.get("x_compound_id", "") or "")
        if not x_id or x_id == exclude_control_id:
            continue
        try:
            proxy_score = float(row.get("similarity_to_target", 0.0) or 0.0)
        except (TypeError, ValueError):
            proxy_score = 0.0
        proxy_rows.append(
            {
                "id": x_id,
                "compound": str(row.get("known_name", "") or x_id),
                "correlation_r": f"{proxy_score:.3f}",
                "pathway": "Verified X structure proxy",
            }
        )
        if len(proxy_rows) >= limit:
            break
    return proxy_rows


def render_guides_html(control_id: str, guide_targets: list[tuple[str, str]]) -> str:
    guide_blocks = "".join(
        f"<section class='guide-card' id='guide-{slug}'><h2>{html.escape(label)}</h2><p>{html.escape(label)} guide will go here.</p></section>"
        for slug, label in guide_targets
    ) or "<section class='guide-card'><h2>No guides yet</h2><p>Guide content will go here.</p></section>"
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>{html.escape(control_id)} Guide Placeholders</title>
  <style>
    :root {{
      --navy:#001D37; --blue:#28496E; --paper:#FFFFFF; --wash:#F1FAFB; --ink:#15304B; --muted:#5C6B7A; --line:rgba(40,73,110,.18);
      --shadow:0 24px 58px rgba(0,29,55,.10);
    }}
    * {{ box-sizing:border-box; }}
    body {{ margin:0; font-family:"Montserrat","Avenir Next","Segoe UI",Arial,sans-serif; color:var(--ink); background:linear-gradient(180deg,#fff 0,#f7fbfc 100%); }}
    main {{ max-width:1080px; margin:0 auto; padding:32px 24px 56px; }}
    h1 {{ margin:0 0 8px; color:var(--navy); }}
    .lede {{ color:var(--muted); margin:0 0 22px; }}
    .guide-grid {{ display:grid; gap:14px; }}
    .guide-card {{ background:var(--paper); border:1px dashed var(--line); border-radius:10px; padding:18px; box-shadow:var(--shadow); scroll-margin-top:18px; }}
    .guide-card h2 {{ margin:0 0 6px; font-size:20px; color:var(--navy); }}
    .guide-card p {{ margin:0; color:var(--muted); }}
  </style>
</head>
<body>
  <main>
    <h1>{html.escape(control_id)} guides</h1>
    <p class="lede">Placeholder guide cards for dashboard info buttons.</p>
    <div class="guide-grid">{guide_blocks}</div>
  </main>
</body>
</html>"""


def build_report(args: argparse.Namespace) -> dict[str, Any]:
    repo_root = REPO_ROOT
    control_id = str(args.control_id).strip()
    slug = control_slug(control_id)
    out_dir = resolve(args.out_dir, repo_root)
    assets = out_dir / "assets"
    data_dir = out_dir / "data"
    out_dir.mkdir(parents=True, exist_ok=True)
    assets.mkdir(parents=True, exist_ok=True)
    data_dir.mkdir(parents=True, exist_ok=True)

    db_path = resolve(args.db, repo_root)
    predictions_path = resolve(args.predictions, repo_root)
    predictions = pd.read_csv(predictions_path, low_memory=False)
    if args.candidate_pool and Path(resolve(args.candidate_pool, repo_root)).exists():
        candidate_pool = pd.read_csv(resolve(args.candidate_pool, repo_root), low_memory=False)
        top = candidate_pool[candidate_pool.get("control_id", pd.Series("", index=candidate_pool.index)).astype(str).eq(control_id)].copy()
    else:
        top = predictions_to_candidate_pool(predictions, control_id)
    if top.empty:
        raise SystemExit(f"No prediction rows found for control_id={control_id} in {predictions_path}")

    known_x_path = resolve(args.known_x, repo_root)
    known_x = pd.read_csv(known_x_path, low_memory=False)
    known_x_all_path = resolve(args.known_x_all, repo_root)
    known_x_all = pd.read_csv(known_x_all_path, low_memory=False) if known_x_all_path.exists() else pd.DataFrame()

    peak_csv = resolve(args.peak_csv, repo_root) if args.peak_csv else None
    peaks = pd.read_csv(peak_csv, low_memory=False) if peak_csv and peak_csv.exists() else pd.DataFrame()

    rank1 = top.sort_values("rank", kind="mergesort").iloc[0]
    smiles_col = "predicted_smiles" if "predicted_smiles" in rank1.index else "smiles"
    target_smiles = str(args.target_smiles or rank1.get(smiles_col) or "").strip()
    if not target_smiles:
        raise SystemExit(f"Could not resolve target SMILES for {control_id}")

    target_name = str(args.target_name or rank1.get("candidate_name") or rank1.get("base_name") or control_id).strip()
    kx_row = lookup_known_x_row(known_x, control_id)
    chemical_id = str(args.chemical_id or resolve_chemical_id(db_path, kx_row, target_name) or "").strip()

    target = {
        "control_id": control_id,
        "name": target_name,
        "smiles": target_smiles,
        "formula": formula_from_smiles(target_smiles),
        "exact_mass": exact_mass_from_smiles(target_smiles),
        "inchikey": inchikey_from_smiles(target_smiles),
    }
    target_structure_svg = render_smiles_svg(target_smiles)
    asset_base = f"{slug}_structure"
    if target_structure_svg:
        (assets / f"{asset_base}.svg").write_text(target_structure_svg, encoding="utf-8")
        structure_panel_html = (
            "<div class='structure-svg-frame'>"
            f"{target_structure_svg}"
            "<p class='structure-note'>2D structure rendered from SMILES.</p>"
            "</div>"
        )
    else:
        structure_image(target_smiles, assets / f"{asset_base}.png", f"{control_id}: {target_name}")
        structure_panel_html = f"<img class='structure' src='assets/{asset_base}.png' alt='{html.escape(target_name)} structure'>"

    top = top.sort_values("rank", kind="mergesort")
    if not peaks.empty:
        feature_col = next((c for c in ("feature_key", "FEATURE_KEY") if c in peaks.columns), "")
        peak_mask = pd.Series(False, index=peaks.index)
        if feature_col:
            peak_mask = peak_mask | peaks[feature_col].astype(str).str.contains(re.escape(control_id), regex=True, na=False)
        if peak_mask.any():
            peaks = peaks[peak_mask].copy()

    top.head(int(args.candidate_tab_n)).copy().to_csv(data_dir / f"{slug}_top_candidates.csv", index=False)
    candidate_score_input = score_smiles_spectrum_alignment(
        top.head(int(args.candidate_tab_n)).copy(),
        peaks,
        tolerance_da=float(args.spectrum_tolerance_da),
    )
    candidate_score_df = compact_candidate_scores(candidate_score_input)
    candidate_score_df.to_csv(data_dir / f"{slug}_candidate_score_summary.csv", index=False)
    top10_candidates = candidate_score_df.head(10).copy()
    candidate_consensus_maps, candidate_consensus_summaries, default_consensus_rank = build_candidate_consensus_maps(top10_candidates)
    (data_dir / f"{slug}_candidate_consensus_maps.json").write_text(
        json.dumps(
            {
                "maps": {str(k): v for k, v in candidate_consensus_maps.items()},
                "summaries": {str(k): v for k, v in candidate_consensus_summaries.items()},
                "default_rank": default_consensus_rank,
            },
            indent=2,
        )
    )
    consensus_rows, consensus_svg = consensus_atom_support(target_smiles, top.head(int(args.consensus_top_n)).copy())
    write_csv(data_dir / f"{slug}_consensus_atom_support.csv", consensus_rows)

    pd.DataFrame([kx_row]).to_csv(data_dir / f"{slug}_known_x_row.csv", index=False)

    db_detail = {}
    syn_df = pd.DataFrame()
    if db_path.exists() and chemical_id:
        db_detail_df = query_df(
            db_path,
            """
            SELECT chemical_id, canonical_name, smiles, inchikeys, super_pathway, pathway,
                   metabolon_super, metabolon_sub, resolution_bucket, resolution_status,
                   resolution_method, resolution_sources, resolution_query,
                   resolution_match_name, resolution_score, resolution_candidate_count
            FROM metabolites
            WHERE chemical_id = ?
            """,
            (chemical_id,),
        )
        db_detail = db_detail_df.iloc[0].to_dict() if not db_detail_df.empty else {}
        db_detail_df.to_csv(data_dir / "metabolite_db_detail.csv", index=False)
        syn_df = query_df(
            db_path,
            """
            SELECT entry_value AS synonym, normalized_value AS normalized, sources
            FROM synonym_index
            WHERE chemical_id = ? AND entry_type = 'synonym'
            ORDER BY lower(entry_value), entry_value
            """,
            (chemical_id,),
        )
        syn_df.to_csv(data_dir / "metabolite_db_synonyms.csv", index=False)

    correlation_rows, correlation_blocker = load_real_correlation_rows(control_id, chemical_id, target_name)

    all_structured = pd.DataFrame()
    if db_path.exists():
        all_structured = query_df(
            db_path,
            """
            SELECT chemical_id, canonical_name, smiles, inchikeys, super_pathway, pathway,
                   metabolon_super, metabolon_sub, resolution_status
            FROM metabolites
            WHERE smiles IS NOT NULL AND trim(smiles) <> ''
            """,
        )

    extended = build_extended_similarity_network(
        target_smiles=target_smiles,
        all_structured=all_structured,
        known_x_all=known_x_all,
        max_known_render=int(args.extended_known_render_limit),
        known_min_similarity=float(args.extended_known_min_similarity),
        x_min_similarity=float(args.extended_x_min_similarity),
        target_edge_min=float(args.extended_target_edge_min),
        pair_edge_min=float(args.extended_pair_edge_min),
        max_edges_per_node=int(args.extended_max_edges_per_node),
        control_id=control_id,
        target_name=target_name,
    )
    if not extended["known_all"].empty:
        extended["known_all"].to_csv(data_dir / f"all_known_similarity_to_{slug}.csv", index=False)
    if not extended["x_all"].empty:
        extended["x_all"].to_csv(data_dir / f"all_unknown_similarity_to_{slug}.csv", index=False)
    write_csv(data_dir / "extended_similarity_network_nodes.csv", extended["render_nodes"])
    write_csv(data_dir / "extended_similarity_network_edges.csv", extended["render_edges"])
    (data_dir / "extended_similarity_network.json").write_text(
        json.dumps({"nodes": extended["render_nodes"], "edges": extended["render_edges"]}, indent=2)
    )

    rows: list[dict[str, Any]] = []
    _, target_fp = mol_fp(target_smiles)
    for _, row in all_structured.iterrows():
        smiles = str(row.get("smiles") or "").split("|")[0].strip()
        mol, fp = mol_fp(smiles)
        if mol is None or fp is None or target_fp is None:
            continue
        sim = float(TanimotoSimilarity(target_fp, fp))
        if sim < float(args.min_similarity):
            continue
        item = row.to_dict()
        item["smiles"] = smiles
        item["similarity_to_target"] = sim
        rows.append(item)
    sim_df = pd.DataFrame(rows)
    if not sim_df.empty:
        sim_df = sim_df.sort_values("similarity_to_target", ascending=False, kind="mergesort").head(int(args.max_neighbors))
        sim_df["confidence"] = sim_df.apply(confidence_bucket, axis=1)
    neighbors = sim_df.to_dict(orient="records") if not sim_df.empty else []
    write_csv(data_dir / "similarity_neighbors.csv", neighbors)

    structural_possibility_rows: list[dict[str, Any]] = []
    structural_possibility_tree_svg = ""
    if args.hypotheses and Path(resolve(args.hypotheses, repo_root)).exists():
        structural_possibility_rows = x25003_structural_possibility_summary(
            candidate_scores=candidate_score_df,
            neighbors=neighbors,
            target_smiles=target_smiles,
            target_formula=target["formula"],
        )
        structural_possibility_tree_svg = x25003_structural_possibility_tree_svg(structural_possibility_rows)
        write_csv(data_dir / f"{slug}_structural_possibility_summary.csv", structural_possibility_rows)

    for row in neighbors[:12]:
        safe = "".join(ch if ch.isalnum() else "-" for ch in str(row.get("canonical_name", "")).lower()).strip("-")[:64]
        structure_image(str(row["smiles"]), assets / f"neighbor-{row['chemical_id']}-{safe}.png", str(row.get("canonical_name", "")))

    nodes = extended["render_nodes"]
    edges = extended["render_edges"]
    write_csv(data_dir / "similarity_network_nodes.csv", nodes)
    write_csv(data_dir / "similarity_network_edges.csv", edges)
    (data_dir / "similarity_network.json").write_text(json.dumps({"nodes": nodes, "edges": edges}, indent=2))

    parent_mz = float(pd.to_numeric(rank1.get("parent_mz"), errors="coerce") or 0.0)
    evidence = {
        "top10_best_tanimoto": float(pd.to_numeric(rank1.get("score_rank"), errors="coerce") or 0.0),
        "full_pool_oracle_best_tanimoto": float(pd.to_numeric(rank1.get("unified_final_score"), errors="coerce") or 0.0),
        "parent_mz": parent_mz,
        "mass_error_da": float(target["exact_mass"] - parent_mz) if parent_mz else 0.0,
        "lactoyl_hypothesis_score": "",
        "diagnostic_bin_score": "",
        "matched_diagnostic_bins_json": "",
    }
    write_csv(data_dir / f"{slug}_evidence_summary.csv", [{"field": k, "value": v} for k, v in evidence.items()])
    (data_dir / f"{slug}_evidence_summary.json").write_text(
        json.dumps({"target": target, "evidence": evidence, "db_detail": db_detail}, indent=2)
    )

    extended_x_all = extended["x_all"]
    if not extended_x_all.empty and "x_compound_id" in extended_x_all.columns:
        extended_unknown_rows = extended_x_all[extended_x_all["x_compound_id"].astype(str) != control_id].head(20).fillna("").to_dict(orient="records")
    else:
        extended_unknown_rows = []
    if not correlation_rows:
        correlation_rows = build_placeholder_correlation_rows(extended_unknown_rows, exclude_control_id=control_id, limit=5)
        correlation_blocker = (
            "Proxy correlation rows use nearest verified X-structure neighbors from the expanded similarity network."
        )
    write_csv(data_dir / f"{slug}_top_correlating_metabolites.csv", correlation_rows)
    (data_dir / f"{slug}_top_correlating_metabolites_status.json").write_text(
        json.dumps({"n_rows": len(correlation_rows), "blocker": correlation_blocker}, indent=2)
    )
    extended_known_rows = extended["known_all"].head(20).fillna("").to_dict(orient="records") if not extended["known_all"].empty else []
    known_render = extended.get("known_render", pd.DataFrame())
    x_render = extended.get("x_render", pd.DataFrame())

    render_kwargs = dict(
        out_dir=out_dir,
        target=target,
        evidence=evidence,
        db_detail=db_detail,
        synonyms=syn_df.fillna("").to_dict(orient="records"),
        neighbors=neighbors,
        nodes=nodes,
        edges=edges,
        candidate_scores=candidate_score_df.fillna("").to_dict(orient="records"),
        structure_panel_html=structure_panel_html,
        consensus_svg=consensus_svg,
        candidate_consensus_maps=candidate_consensus_maps,
        candidate_consensus_summaries=candidate_consensus_summaries,
        default_consensus_rank=default_consensus_rank,
        structural_possibility_tree_svg=structural_possibility_tree_svg,
        structural_possibility_rows=structural_possibility_rows,
        feature_total_count=int(len(peaks)),
        extended_network_svg_html=extended_network_svg(extended["render_nodes"], extended["render_edges"]),
        extended_unknown_rows=extended_unknown_rows,
        extended_known_rows=extended_known_rows,
        extended_stats={
            "n_known_all": int(len(extended["known_all"])),
            "n_x_all": int(len(extended["x_all"])),
            "n_known_render": int(len(known_render)),
            "n_x_render": int(len(x_render)),
        },
        correlation_rows=correlation_rows,
        correlation_blocker=correlation_blocker,
    )
    public_html = render_html(**render_kwargs, variant="public")
    internal_html = render_html(**render_kwargs, variant="internal")
    public_guide_targets = [
        (slugify_label(label), label)
        for label in [
            "Top Candidate Compound",
            "Parameters Contributing to Identification Confidence",
            "Top 10 potential structures",
            "Top 10 Metabolon Library Matches",
            "Fragment Scores",
            "Mass Match Scores",
            "MS Fit Scores",
            "Top Correlating Metabolites",
            "Similarity Network",
        ]
    ]
    guides_html = render_guides_html(control_id, public_guide_targets)
    (out_dir / "public.html").write_text(public_html)
    (out_dir / "guides.html").write_text(guides_html)
    (out_dir / "internal.html").write_text(internal_html)
    (out_dir / "index.html").write_text(public_html)
    compound_html = out_dir / f"{control_id.replace('-', '_')}.html"
    compound_html.write_text(public_html)
    payload = {
        "out_dir": str(out_dir),
        "html": str(out_dir / "index.html"),
        "compound_html": str(compound_html),
        "public_html": str(out_dir / "public.html"),
        "guides_html": str(out_dir / "guides.html"),
        "internal_html": str(out_dir / "internal.html"),
        "n_neighbors": int(len(neighbors)),
        "target": target,
        "db_chemical_id": chemical_id,
        "control_id": control_id,
    }
    (out_dir / "manifest.json").write_text(json.dumps(payload, indent=2))
    print(json.dumps(payload, indent=2, sort_keys=True))
    return payload
'''

MAIN_AND_HELPERS = BUILD_REPORT  # helpers + build_report only

BUILD_REPORT = MAIN_AND_HELPERS.replace(
    "\n\ndef main() -> None:",
    "",
).rsplit("\nif __name__", 1)[0]

MAIN_BLOCK = '''
def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--control-id", required=True, help="X compound control id, e.g. X-25003")
    parser.add_argument("--predictions", type=Path, required=True, help="Rerank or candidate pool CSV")
    parser.add_argument("--out-dir", type=Path, default=None)
    parser.add_argument("--db", type=Path, default=DEFAULT_DB)
    parser.add_argument("--known-x", type=Path, default=DEFAULT_KNOWN_X)
    parser.add_argument("--known-x-all", type=Path, default=DEFAULT_KNOWN_X_ALL)
    parser.add_argument("--candidate-pool", type=Path, default=None)
    parser.add_argument("--hypotheses", type=Path, default=None)
    parser.add_argument("--peak-csv", type=Path, default=None)
    parser.add_argument("--target-smiles", type=str, default="")
    parser.add_argument("--target-name", type=str, default="")
    parser.add_argument("--chemical-id", type=str, default="")
    parser.add_argument("--spectrum-tolerance-da", type=float, default=1.0)
    parser.add_argument("--min-similarity", type=float, default=0.28)
    parser.add_argument("--max-neighbors", type=int, default=10)
    parser.add_argument("--candidate-tab-n", type=int, default=20)
    parser.add_argument("--consensus-top-n", type=int, default=20)
    parser.add_argument("--extended-known-render-limit", type=int, default=120)
    parser.add_argument("--extended-known-min-similarity", type=float, default=0.30)
    parser.add_argument("--extended-x-min-similarity", type=float, default=0.0)
    parser.add_argument("--extended-target-edge-min", type=float, default=0.32)
    parser.add_argument("--extended-pair-edge-min", type=float, default=0.42)
    parser.add_argument("--extended-max-edges-per-node", type=int, default=6)
    args = parser.parse_args()
    if args.out_dir is None:
        args.out_dir = Path("output/reports") / args.control_id
    build_report(args)


if __name__ == "__main__":
    main()
'''


def main() -> None:
    """Apply regex substitutions to generalize the legacy X-25003 report builder."""
    text = TARGET.read_text(encoding="utf-8")

    text = text.replace(
        '"""Build a high-level biological exploration report for X-25003."""',
        '"""Build a self-contained HTML exploration report for any X compound."""',
    )
    text = text.replace("from lactoyl_amino_acid_v1_utils import THEORETICAL_N_LACTOYL_SMILES\n", "")
    text = re.sub(
        r"DEFAULT_OUT_DIR = Path\([^\)]+\)\nDEFAULT_DB = Path\([^\)]+\)\n(?:DEFAULT_[^\n]+\n)+",
        """REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT_DIR = REPO_ROOT / "output/reports"
DEFAULT_DB = REPO_ROOT / "vendor/metabolite_synonyms.sqlite"
DEFAULT_KNOWN_X = REPO_ROOT / "data/known_X.csv"
DEFAULT_KNOWN_X_ALL = REPO_ROOT / "data/known_X_all_available_structures_simple.csv"
""",
        text,
        count=1,
    )
    text = text.replace("similarity_to_x25003", "similarity_to_target")
    text = text.replace("load_x25003_real_correlation_rows", "load_real_correlation_rows_LEGACY")
    text = text.replace("build_x25003_placeholder_correlation_rows", "build_placeholder_correlation_rows_LEGACY")
    text = text.replace("render_x25003_guides_html", "render_guides_html_LEGACY")

    text = re.sub(
        r"def build_extended_similarity_network\(\n    \*,\n    target_smiles: str,\n    all_structured: pd\.DataFrame,\n    known_x_all: pd\.DataFrame,\n    max_known_render: int,\n    known_min_similarity: float,\n    x_min_similarity: float,\n    target_edge_min: float,\n    pair_edge_min: float,\n    max_edges_per_node: int,\n\) -> dict\[str, Any\]:",
        """def build_extended_similarity_network(
    *,
    target_smiles: str,
    all_structured: pd.DataFrame,
    known_x_all: pd.DataFrame,
    max_known_render: int,
    known_min_similarity: float,
    x_min_similarity: float,
    target_edge_min: float,
    pair_edge_min: float,
    max_edges_per_node: int,
    control_id: str = "X-00000",
    target_name: str = "target",
) -> dict[str, Any]:""",
        text,
        count=1,
    )
    text = text.replace(
        '            "id": "X-25003",\n            "label": "X-25003 · N-lactoylalanine",\n            "display_label": "X-25003",',
        '            "id": control_id,\n            "label": f"{control_id} · {target_name}",\n            "display_label": control_id,',
    )
    text = text.replace(
        'x_render = x_all[\n        x_all["similarity_to_target"] >= x_min_similarity\n    ].copy()\n    if not x_render.empty:\n        x_render = x_render[x_render["x_compound_id"].astype(str) != "X-25003"].copy()',
        'x_render = x_all[x_all["similarity_to_target"] >= x_min_similarity].copy()\n    if not x_render.empty:\n        x_render = x_render[x_render["x_compound_id"].astype(str) != control_id].copy()',
    )

    text = re.sub(
        r"def build_report\(args: argparse\.Namespace\) -> dict\[str, Any\]:.*?^def main\(\) -> None:",
        BUILD_REPORT + "\n\ndef main() -> None:",
        text,
        count=1,
        flags=re.MULTILINE | re.DOTALL,
    )
    text = re.sub(
        r"def main\(\) -> None:.*?^if __name__ == \"__main__\":\n    main\(\)\s*$",
        MAIN_BLOCK.strip(),
        text,
        count=1,
        flags=re.MULTILINE | re.DOTALL,
    )

    TARGET.write_text(text, encoding="utf-8")
    print(f"Patched {TARGET}")


if __name__ == "__main__":
    main()
