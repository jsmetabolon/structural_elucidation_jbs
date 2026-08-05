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
import re
import sqlite3
import sys
from pathlib import Path
from typing import Any

csv.field_size_limit(10 ** 7)
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "structural_elucidation" / "reports"))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "structural_elucidation" / "nce"))
import confidence as C  # noqa: E402
import tentative_naming as TN  # noqa: E402
import tentative_naming_ext as TN_EXT  # noqa: E402

try:
    import evidence as NCE_EVIDENCE  # noqa: E402
except ImportError:
    NCE_EVIDENCE = None

LEVEL_JUNK = {"", "nan", "artifact", "exclude", "4?"}
# Feature keys look like ..._mz100.0404_ri1248_X-23917. Mirrors
# structural_elucidation.utils.domain.ri_from_feature_key, which this tool cannot
# import without pulling in pandas.
FEATURE_KEY_RI = re.compile(r"_ri(\d+)_")
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
            ion_bins: list[Any] = []
            loss_bins: list[Any] = []
            for key, target in (
                ("observed_ion_bins_json", ion_bins),
                ("observed_loss_bins_json", loss_bins),
            ):
                raw = str(row.get(key, "") or "").strip()
                if raw:
                    try:
                        parsed = json.loads(raw)
                        if isinstance(parsed, list):
                            target.extend(parsed)
                    except (ValueError, TypeError):
                        pass
            out[cid] = {
                "present": str(row.get("ms2_present", "")).strip().lower() in ("true", "1"),
                "n_pairs": _num(row.get("n_ms2_pairs")) or 0,
                "n_ion_bins": _num(row.get("n_observed_ion_bins")) or 0,
                "n_loss_bins": _num(row.get("n_observed_loss_bins")) or 0,
                "parent_mz": _num(row.get("parent_mz")),
                "neutral_mass": _num(row.get("neutral_mass")),
                "observed_ion_bins": ion_bins,
                "observed_loss_bins": loss_bins,
                "lib_id": _num(row.get("lib_id")) or 209,
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


def candidate_evidence_rows(css: list[dict], tc: list[dict]) -> list[dict]:
    """Per-rank evidence for the top-10 widget.

    Joins the pipeline's ``evidence_confidence`` + ``mass_match_score`` (from
    ``top_candidates.csv``) onto the displayed candidate rows so the patcher can
    surface them alongside the raw rerank score. Matches by SMILES, falling back
    to row order (both files are rank-ordered).
    """
    ev_by_smiles: dict[str, dict] = {}
    for row in tc:
        s = str(row.get("predicted_smiles", "")).strip()
        if s and s not in ev_by_smiles:
            ev_by_smiles[s] = {
                "evidence_confidence": _num(row.get("evidence_confidence")),
                "mass_match_score": _num(row.get("mass_match_score")),
            }
    out = []
    for i, row in enumerate(css[:10]):
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
            "evidence_confidence": ev.get("evidence_confidence"),
            "mass_match_score": ev.get("mass_match_score"),
        })
    return out


def fragment_candidate_fits(tc: list[dict], ms2_rec: dict) -> list[dict]:
    """Per-candidate matched ion/loss bin counts (non-raw; public-safe)."""

    if not tc or not ms2_rec.get("present"):
        return []
    ion_bins = {str(item) for item in (ms2_rec.get("observed_ion_bins") or []) if str(item) not in {"", "0"}}
    loss_bins = {str(item) for item in (ms2_rec.get("observed_loss_bins") or []) if str(item) not in {"", "0"}}
    if not ion_bins and not loss_bins:
        return []
    lib_id = int(ms2_rec.get("lib_id") or 209)
    out: list[dict] = []
    for row in tc[:3]:
        rank = _num(row.get("rank"))
        smiles = str(row.get("predicted_smiles", "")).strip()
        matched_ions = None
        matched_losses = None
        if smiles and NCE_EVIDENCE is not None:
            pred_ions, pred_losses = NCE_EVIDENCE.candidate_fragment_bins(smiles, lib_id)
            if ion_bins:
                matched_ions = len({str(item) for item in pred_ions} & ion_bins)
            if loss_bins:
                matched_losses = len({str(item) for item in pred_losses} & loss_bins)
        out.append(
            {
                "rank": int(rank) if rank is not None else len(out) + 1,
                "matched_ion_bins": matched_ions,
                "matched_loss_bins": matched_losses,
            }
        )
    return out


def retention_inputs(rank1: dict, candidate_summary: dict) -> tuple[Any, float | None, float | None]:
    """Recover the retention evidence for the top candidate.

    The observed index is not a column in any sidecar; it is embedded in the
    feature key, which the reranker carries through. The predicted index needs
    RDKit, so it is computed here rather than in the confidence model.
    """
    lib_id = rank1.get("lib_id") or candidate_summary.get("lib_id")
    observed = _ri_from_feature_key(rank1.get("feature_key") or candidate_summary.get("feature_key"))
    smiles = str(rank1.get("predicted_smiles") or candidate_summary.get("predicted_smiles") or "").strip()
    proxy = None
    if smiles and NCE_EVIDENCE is not None:
        value = NCE_EVIDENCE.ri_proxy_from_smiles(smiles)
        proxy = float(value) if value == value else None  # NaN check without importing math
    return lib_id, observed, proxy


def _ri_from_feature_key(value: Any) -> float | None:
    match = FEATURE_KEY_RI.search(str(value or ""))
    if not match:
        return None
    ri = float(match.group(1))
    return ri if ri > 0 else None


def enrich(report_dir: Path, obs_mz, ms2, levels, db_names, wgcna, ledger=None) -> str:
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

    lib_id, observed_ri, retention_proxy = retention_inputs(r1, css[0] if css else {})
    fields = C.fields_from_parts(
        css[0] if css else {}, css[1] if len(css) > 1 else None, pw[0] if pw else {},
        observed_mz_available=observed_mz is not None,
        ms2_spectrum_present=ms2_present, ms2_score=ms2_score,
        library_neighbors_ge_050=neigh, correlations_real=correlations_real,
        lib_id=lib_id,
        observed_retention_index=observed_ri,
        candidate_retention_proxy=retention_proxy,
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
    fields["fragment_candidate_fits"] = fragment_candidate_fits(tc, ms2_rec)
    fields["metabolon_level"] = levels.get(cid, "")
    fields["truth_match"] = truth_match
    fields["mass_alt"] = C.mass_consistent_alt(css)
    fields["candidate_evidence"] = candidate_evidence_rows(css, tc)
    # Reviewed verdict: drives the headline, the MSI level and the identity metric.
    if ledger:
        C.apply_curation_entry(fields, ledger.get(cid))
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
    smiles = str(target.get("smiles", "")).strip()
    fields["candidate_exact_mass"] = target.get("exact_mass")
    dbn = db_names.get(ik.split("-")[0]) if ik else None
    name_resolved = bool(dbn)
    tentative = TN.resolve(
        smiles,
        formula,
        r1,
        db_names,
        name_resolved=name_resolved,
        bundle_dir=report_dir,
        extended_ctx=TN_EXT.default_extended_context(reports_root=report_dir.parent),
        inchikey=ik,
        control_id=cid,
    )
    if dbn:
        resolved_name = dbn["name"]
        fields["synonyms"] = dbn["synonyms"]
    else:
        resolved_name = tentative.label if tentative else "De novo structure"
        fields["synonyms"] = []
    fields["smiles"] = smiles
    fields["name_resolved"] = name_resolved
    fields.update(TN.descriptor_fields(tentative))
    # co-abundance neighbors (gap 1, placeholder)
    fields["correlations"] = wgcna.get(cid, [])
    # derive + store the MSI identification level (also derivable at render time)
    conf = C.compute_confidence(fields)
    fields["msi_level"] = conf.msi_level
    fields["msi_label"] = conf.msi_label
    fields["msi_rationale"] = conf.msi_rationale
    # The level decides what the compound may be called. Levels 4-5 withhold the
    # structure, so they cannot carry a structural name however well it resolved.
    fields["resolved_name"] = resolved_name
    fields["display_name"] = TN.display_name_for_level(
        resolved_name,
        msi_level=conf.msi_level,
        formula=formula,
        observed_mz=fields.get("observed_mz"),
        name_is_reference=name_resolved,
    )
    if conf.msi_level >= 4:
        fields["synonyms"] = []
    (data / f"{slug}_confidence_inputs.json").write_text(json.dumps(fields, indent=2))
    return "enriched"


DEFAULT_CURATION_LEDGER = (
    Path(__file__).resolve().parent.parent / "data" / "curation_ledger.csv"
)


def main() -> None:
    root = Path(sys.argv[1] if len(sys.argv) > 1 else "_audit_pull/reports")
    bundle = Path(sys.argv[2] if len(sys.argv) > 2 else "_bundle")
    known_x = Path(sys.argv[3] if len(sys.argv) > 3 else "data/known_X.csv")
    ledger_path = Path(sys.argv[4]) if len(sys.argv) > 4 else DEFAULT_CURATION_LEDGER
    obs_mz = load_observed_mz(bundle, known_x)
    ms2 = load_ms2(bundle)
    levels = load_levels(known_x)
    db_names = load_db_names(bundle)
    wgcna = load_wgcna(bundle)
    ledger = C.load_curation_ledger(ledger_path)
    print(f"[enrich] m/z:{len(obs_mz)} ms2:{len(ms2)} levels:{len(levels)} "
          f"db_names:{len(db_names)} wgcna:{len(wgcna)} reviewed_verdicts:{len(ledger)}",
          file=sys.stderr)
    counts: dict[str, int] = {}
    for d in sorted(p for p in root.iterdir() if p.is_dir()):
        st = enrich(d, obs_mz, ms2, levels, db_names, wgcna, ledger=ledger)
        counts[st] = counts.get(st, 0) + 1
    print(f"[enrich] {counts}", file=sys.stderr)


if __name__ == "__main__":
    main()
