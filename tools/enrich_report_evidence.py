#!/usr/bin/env python3
"""Fold the report_data_bundle evidence into each report as confidence inputs.

Writes ``<slug>_confidence_inputs.json`` into every report's ``data/`` dir with
the normalized confidence ``fields`` plus display extras, combining:
  - observed parent m/z  <- feature_rows/*.csv (fixes the ~85% missing-m/z gap)
  - observed MS2 present  <- observed_ms2_lookups/per_control_all_x.csv
  - pipeline MS2 fit      <- <slug>_top_candidates.csv (predicted/forward spectrum)
  - structure / mass / pathway / orthogonal <- existing sidecars

Pure standard library (no RDKit/pandas). Downstream (confidence patcher, audit,
calibration) then reads a single self-contained sidecar.

Usage:
    python3 tools/enrich_report_evidence.py <reports_root> <bundle_dir>
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

csv.field_size_limit(10 ** 7)
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "structural_elucidation" / "reports"))
import confidence as C  # noqa: E402


def _num(v):
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return None if f != f else f


def load_observed_mz(bundle: Path) -> dict[str, dict]:
    """control_id -> {parent_mz, neutral_mass, ion_mode} from feature rows."""
    out: dict[str, dict] = {}
    for fr in sorted((bundle / "feature_rows").glob("*_feature_rows.csv")):
        with fr.open(newline="") as h:
            for row in csv.DictReader(h):
                cid = str(row.get("control_id", "")).strip()
                if not cid or cid in out:
                    continue
                mz = _num(row.get("parent_mz"))
                if mz is None:
                    continue
                out[cid] = {
                    "parent_mz": mz,
                    "neutral_mass": _num(row.get("neutral_mass_fallback")),
                    "ion_mode": str(row.get("ion_mode", "")).strip(),
                }
    return out


def load_ms2(bundle: Path) -> dict[str, dict]:
    """control_id -> {present, n_pairs, parent_mz, neutral_mass} from MS2 lookup."""
    out: dict[str, dict] = {}
    path = bundle / "observed_ms2_lookups" / "per_control_all_x.csv"
    if not path.exists():
        return out
    with path.open(newline="") as h:
        for row in csv.DictReader(h):
            cid = str(row.get("control_id", "")).strip()
            if not cid:
                continue
            present = str(row.get("ms2_present", "")).strip().lower() in ("true", "1")
            out[cid] = {
                "present": present,
                "n_pairs": _num(row.get("n_ms2_pairs")) or 0,
                "parent_mz": _num(row.get("parent_mz")),
                "neutral_mass": _num(row.get("neutral_mass")),
            }
    return out


def _read_rows(path: Path) -> list[dict]:
    if not path.exists():
        return []
    with path.open(newline="") as h:
        return list(csv.DictReader(h))


def _rank1(rows: list[dict]) -> dict:
    best, br = {}, None
    for r in rows:
        k = _num(r.get("rank"))
        k = 1e9 if k is None else k
        if br is None or k < br:
            br, best = k, r
    return best


def pipeline_ms2_score(tc_rows: list[dict]) -> float | None:
    if not tc_rows:
        return None
    row = _rank1(tc_rows)
    cands = [
        _num(row.get("predicted_spectrum_score_v2")),
        _num(row.get("forward_spectrum_score")),
        _num(row.get("ms2_match_score")),
    ]
    cands = [c for c in cands if c is not None]
    return max(cands) if cands else None


def enrich(report_dir: Path, obs_mz: dict, ms2: dict) -> str:
    data = report_dir / "data"
    if not data.is_dir():
        return "no_data"
    slug = report_dir.name.lower()
    cid = report_dir.name
    css = _read_rows(data / f"{slug}_candidate_score_summary.csv")
    pw = _read_rows(data / f"{slug}_pathway_prediction.csv")
    if not css and not pw:
        return "no_sidecars"
    tc = _read_rows(data / f"{slug}_top_candidates.csv")

    mz_rec = obs_mz.get(cid, {})
    ms2_rec = ms2.get(cid, {})
    observed_mz = mz_rec.get("parent_mz") if mz_rec.get("parent_mz") is not None else ms2_rec.get("parent_mz")
    ms2_present = bool(ms2_rec.get("present"))
    ms2_score = pipeline_ms2_score(tc)

    neigh = 0
    for row in _read_rows(data / "similarity_neighbors.csv"):
        s = _num(row.get("similarity_to_target"))
        if s is not None and s >= 0.50:
            neigh += 1
    correlations_real = True
    status = data / f"{slug}_top_correlating_metabolites_status.json"
    if status.exists():
        try:
            correlations_real = not json.loads(status.read_text()).get("blocker")
        except (ValueError, OSError):
            correlations_real = True

    fields = C.fields_from_parts(
        css[0] if css else {},
        css[1] if len(css) > 1 else None,
        pw[0] if pw else {},
        observed_mz_available=observed_mz is not None,
        ms2_spectrum_present=ms2_present,
        ms2_score=ms2_score,
        library_neighbors_ge_050=neigh,
        correlations_real=correlations_real,
    )
    # display extras (patcher reads these for the Mass / MS Fit sections)
    fields["observed_mz"] = observed_mz
    fields["observed_neutral_mass"] = mz_rec.get("neutral_mass") or ms2_rec.get("neutral_mass")
    fields["ms2_n_pairs"] = ms2_rec.get("n_pairs") or 0
    (data / f"{slug}_confidence_inputs.json").write_text(json.dumps(fields, indent=2))
    return "enriched"


def main() -> None:
    root = Path(sys.argv[1] if len(sys.argv) > 1 else "_audit_pull/reports")
    bundle = Path(sys.argv[2] if len(sys.argv) > 2 else "_bundle")
    obs_mz = load_observed_mz(bundle)
    ms2 = load_ms2(bundle)
    print(f"[enrich] observed m/z map: {len(obs_mz)} | MS2 lookup: {len(ms2)}", file=sys.stderr)
    counts: dict[str, int] = {}
    for d in sorted(p for p in root.iterdir() if p.is_dir()):
        st = enrich(d, obs_mz, ms2)
        counts[st] = counts.get(st, 0) + 1
    print(f"[enrich] {counts}", file=sys.stderr)


if __name__ == "__main__":
    main()
