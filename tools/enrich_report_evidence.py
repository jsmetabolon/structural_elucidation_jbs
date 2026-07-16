#!/usr/bin/env python3
"""Fold bundle evidence into each report as confidence inputs + display extras.

Writes ``<slug>_confidence_inputs.json`` into every report's ``data/`` dir with
the normalized confidence ``fields`` plus display extras used by the patcher:
  - observed parent m/z  <- feature_rows + known_X m/z fallback  (gap 4)
  - observed MS2 present  <- observed_ms2_lookups                 (gap 3)
  - pipeline MS2 fit      <- <slug>_top_candidates.csv
  - Metabolon Level       <- known_X.csv (authentic-standard)     (gap 6)
  - display name/synonyms <- metabolite_synonyms.sqlite / systematic descriptor (gap 2)
  - mass-consistent alt   <- <slug>_candidate_score_summary.csv   (gap 5)
  - co-abundance neighbors<- wgcna parent_neighbors (placeholder)  (gap 1)

Pure standard library. Downstream reads a single self-contained sidecar.

Usage:
    python3 tools/enrich_report_evidence.py <reports_root> <bundle_dir> [known_X.csv]
"""
from __future__ import annotations

import csv
import json
import sqlite3
import sys
from pathlib import Path

csv.field_size_limit(10 ** 7)
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "structural_elucidation" / "reports"))
import confidence as C  # noqa: E402

LEVEL_JUNK = {"", "nan", "artifact", "exclude", "4?"}
TRUTH_CORRECT = 0.70


def _num(v):
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return None if f != f else f


def _read(path: Path) -> list[dict]:
    if not path.exists():
        return []
    with path.open(newline="") as h:
        return list(csv.DictReader(h))


def load_observed_mz(bundle: Path, known_x: Path) -> dict[str, dict]:
    out: dict[str, dict] = {}
    for fr in sorted((bundle / "feature_rows").glob("*_feature_rows.csv")):
        for row in _read(fr):
            cid = str(row.get("control_id", "")).strip()
            mz = _num(row.get("parent_mz"))
            if cid and cid not in out and mz is not None:
                out[cid] = {"parent_mz": mz, "neutral_mass": _num(row.get("neutral_mass_fallback"))}
    # known_X m/z fallback (raises coverage to ~97%)
    for row in _read(known_x):
        cid = str(row.get("Unknown", "")).strip()
        mz = _num(row.get("m/z"))
        if cid and cid not in out and mz is not None:
            out[cid] = {"parent_mz": mz, "neutral_mass": None}
    return out


def load_ms2(bundle: Path) -> dict[str, dict]:
    out: dict[str, dict] = {}
    path = bundle / "observed_ms2_lookups" / "per_control_all_x.csv"
    for row in _read(path):
        cid = str(row.get("control_id", "")).strip()
        if cid:
            out[cid] = {
                "present": str(row.get("ms2_present", "")).strip().lower() in ("true", "1"),
                "n_pairs": _num(row.get("n_ms2_pairs")) or 0,
                "n_ion_bins": _num(row.get("n_observed_ion_bins")) or 0,
                "n_loss_bins": _num(row.get("n_observed_loss_bins")) or 0,
                "parent_mz": _num(row.get("parent_mz")),
                "neutral_mass": _num(row.get("neutral_mass")),
            }
    return out


def load_levels(known_x: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    for row in _read(known_x):
        cid = str(row.get("Unknown", "")).strip()
        lvl = str(row.get("Level", "")).strip()
        if cid and cid not in out and lvl.lower() not in LEVEL_JUNK:
            out[cid] = lvl
    return out


def load_db_names(bundle: Path) -> dict[str, dict]:
    """InChIKey skeleton -> {name, synonyms}."""
    db = bundle / "metabolite_synonyms.sqlite"
    out: dict[str, dict] = {}
    if not db.exists():
        return out
    con = sqlite3.connect(str(db))
    for ik, name, shorts, chem in con.execute(
        "SELECT inchikeys, canonical_name, short_names, chemical_names FROM metabolites"
    ):
        if not ik or not name:
            continue
        syns: list[str] = []
        for src in (shorts, chem):
            if src:
                syns += [s.strip() for s in str(src).replace(";", ",").split(",") if s.strip()]
        for part in str(ik).replace(";", ",").split(","):
            skel = part.strip().split("-")[0]
            if skel and skel not in out:
                out[skel] = {"name": str(name).strip(), "synonyms": syns[:8]}
    con.close()
    return out


def load_wgcna(bundle: Path) -> dict[str, list]:
    out: dict[str, list] = {}
    for f in sorted((bundle / "wgcna").glob("lib*_parent_neighbors.csv")):
        for row in _read(f):
            cid = str(row.get("control_id", "")).strip()
            if not cid:
                continue
            out.setdefault(cid, [])
            if len(out[cid]) < 6:
                out[cid].append({
                    "name": str(row.get("parent_neighbor_name", "")).strip(),
                    "id": str(row.get("parent_neighbor_id", "")).strip(),
                    "score": _num(row.get("parent_neighbor_score")),
                    "source": str(row.get("parent_neighbor_source", "")).strip(),
                })
    return out


def _rank1(rows: list[dict]) -> dict:
    best, br = {}, None
    for r in rows:
        k = _num(r.get("rank"))
        k = 1e9 if k is None else k
        if br is None or k < br:
            br, best = k, r
    return best


def candidate_evidence_rows(css: list[dict], tc: list[dict], feature_overall: float | None) -> list[dict]:
    """Per-rank evidence + rank-consistent confidence for the top-10 widget.

    The featured per-candidate number is ``candidate_confidence``
    (``confidence.candidate_confidences``): monotone non-increasing with rank and
    anchored so rank 1 equals the report overall confidence. The class-level
    ``evidence_confidence`` + ``mass_match_score`` (from ``top_candidates.csv``,
    matched by SMILES, falling back to rank-ordered row order) are retained as
    supporting context. ``tied`` flags genuinely unresolved same-formula isomers.
    """
    ev_by_smiles: dict[str, dict] = {}
    for row in tc:
        s = str(row.get("predicted_smiles", "")).strip()
        if s and s not in ev_by_smiles:
            ev_by_smiles[s] = {
                "evidence_confidence": _num(row.get("evidence_confidence")),
                "mass_match_score": _num(row.get("mass_match_score")),
            }
    top = css[:10]
    cc = C.candidate_confidences(top, feature_overall)
    out = []
    for i, row in enumerate(top):
        s = str(row.get("predicted_smiles", "")).strip()
        ev = ev_by_smiles.get(s)
        if ev is None and i < len(tc):
            ev = {
                "evidence_confidence": _num(tc[i].get("evidence_confidence")),
                "mass_match_score": _num(tc[i].get("mass_match_score")),
            }
        ev = ev or {}
        rank = _num(row.get("rank"))
        out.append({
            "rank": int(rank) if rank is not None else i + 1,
            "candidate_confidence": cc[i]["candidate_confidence"],
            "candidate_confidence_tied": cc[i]["tied"],
            "evidence_confidence": ev.get("evidence_confidence"),
            "mass_match_score": ev.get("mass_match_score"),
        })
    return out


def enrich(report_dir: Path, obs_mz, ms2, levels, db_names, wgcna) -> str:
    data = report_dir / "data"
    if not data.is_dir():
        return "no_data"
    slug = report_dir.name.lower()
    cid = report_dir.name
    css = _read(data / f"{slug}_candidate_score_summary.csv")
    pw = _read(data / f"{slug}_pathway_prediction.csv")
    if not css and not pw:
        return "no_sidecars"
    tc = _read(data / f"{slug}_top_candidates.csv")

    mz_rec = obs_mz.get(cid, {})
    ms2_rec = ms2.get(cid, {})
    observed_mz = mz_rec.get("parent_mz") if mz_rec.get("parent_mz") is not None else ms2_rec.get("parent_mz")
    ms2_present = bool(ms2_rec.get("present"))
    r1 = _rank1(tc) if tc else {}
    ms2_cands = [_num(r1.get(k)) for k in ("predicted_spectrum_score_v2", "forward_spectrum_score", "ms2_match_score")]
    ms2_cands = [c for c in ms2_cands if c is not None]
    ms2_score = max(ms2_cands) if ms2_cands else None
    truth = _num(r1.get("truth_tanimoto"))
    truth_match = truth is not None and truth >= TRUTH_CORRECT

    neigh = sum(1 for row in _read(data / "similarity_neighbors.csv")
                if (_num(row.get("similarity_to_target")) or 0) >= 0.50)
    correlations_real = True
    status = data / f"{slug}_top_correlating_metabolites_status.json"
    if status.exists():
        try:
            correlations_real = not json.loads(status.read_text()).get("blocker")
        except (ValueError, OSError):
            correlations_real = True

    fields = C.fields_from_parts(
        css[0] if css else {}, css[1] if len(css) > 1 else None, pw[0] if pw else {},
        observed_mz_available=observed_mz is not None,
        ms2_spectrum_present=ms2_present, ms2_score=ms2_score,
        library_neighbors_ge_050=neigh, correlations_real=correlations_real,
    )
    # display + confidence extras
    fields["observed_mz"] = observed_mz
    fields["observed_neutral_mass"] = mz_rec.get("neutral_mass") or ms2_rec.get("neutral_mass")
    fields["ms2_n_pairs"] = ms2_rec.get("n_pairs") or 0
    fields["ms2_detail"] = {
        "present": ms2_present,
        "n_pairs": ms2_rec.get("n_pairs") or 0,
        "n_ion_bins": ms2_rec.get("n_ion_bins") or 0,
        "n_loss_bins": ms2_rec.get("n_loss_bins") or 0,
        "forward": _num(r1.get("forward_spectrum_score")),
        "predicted_v2": _num(r1.get("predicted_spectrum_score_v2")),
        "match": _num(r1.get("ms2_match_score")),
    }
    fields["metabolon_level"] = levels.get(cid, "")
    fields["truth_match"] = truth_match
    fields["mass_alt"] = C.mass_consistent_alt(css)
    # display name / synonyms (gap 2)
    ev = data / f"{slug}_evidence_summary.json"
    target = {}
    if ev.exists():
        try:
            target = json.loads(ev.read_text()).get("target", {})
        except (ValueError, OSError):
            target = {}
    ik = str(target.get("inchikey", "")).strip()
    formula = str(target.get("formula", "")).strip()
    dbn = db_names.get(ik.split("-")[0]) if ik else None
    if dbn:
        fields["display_name"] = dbn["name"]
        fields["synonyms"] = dbn["synonyms"]
    else:
        fields["display_name"] = f"De novo structure ({formula})" if formula else cid
        fields["synonyms"] = []
    fields["name_resolved"] = bool(dbn)
    # co-abundance neighbors (gap 1, placeholder)
    fields["correlations"] = wgcna.get(cid, [])
    # derive + store the MSI identification level (also derivable at render time)
    conf = C.compute_confidence(fields)
    fields["msi_level"] = conf.msi_level
    fields["msi_label"] = conf.msi_label
    fields["msi_rationale"] = conf.msi_rationale
    # Per-candidate confidence is anchored to the report's evidence match strength
    # (raw_evidence), which is fuller-range than the compressed Platt overall.
    fields["candidate_evidence"] = candidate_evidence_rows(css, tc, conf.raw_evidence)
    (data / f"{slug}_confidence_inputs.json").write_text(json.dumps(fields, indent=2))
    return "enriched"


def main() -> None:
    root = Path(sys.argv[1] if len(sys.argv) > 1 else "_audit_pull/reports")
    bundle = Path(sys.argv[2] if len(sys.argv) > 2 else "_bundle")
    known_x = Path(sys.argv[3] if len(sys.argv) > 3 else "data/known_X.csv")
    obs_mz = load_observed_mz(bundle, known_x)
    ms2 = load_ms2(bundle)
    levels = load_levels(known_x)
    db_names = load_db_names(bundle)
    wgcna = load_wgcna(bundle)
    print(f"[enrich] m/z:{len(obs_mz)} ms2:{len(ms2)} levels:{len(levels)} "
          f"db_names:{len(db_names)} wgcna:{len(wgcna)}", file=sys.stderr)
    counts: dict[str, int] = {}
    for d in sorted(p for p in root.iterdir() if p.is_dir()):
        st = enrich(d, obs_mz, ms2, levels, db_names, wgcna)
        counts[st] = counts.get(st, 0) + 1
    print(f"[enrich] {counts}", file=sys.stderr)


if __name__ == "__main__":
    main()
