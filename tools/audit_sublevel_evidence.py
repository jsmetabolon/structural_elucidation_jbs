#!/usr/bin/env python3
"""Feasibility + validation audit for MSI sublevel assignment.

Answers three questions on the labeled cohort (controls with a curated ``SubLevel``
in ``known_X``):
  1. Separability - do our available evidence axes distinguish the curated
     sublevels (2C, 3A-3D, 4, 5) at all?  (per-sublevel feature medians)
  2. Agreement - how well does a proposed *derived* sublevel rubric reproduce the
     curated ``SubLevel``?  (confusion matrix + adjacent-order agreement)
  3. Defensibility - does the derived sublevel *ordering* track correctness, i.e.
     do higher sublevels have higher truth Tanimoto?  (monotonicity vs truth) -
     this is what justifies client include/exclude by sublevel.

Reads the delivered reports tree (per-report ``data/`` sidecars) plus the curated
``known_X``/known_x_row sublevel. Pure standard library.

Usage:
    python3 tools/audit_sublevel_evidence.py <reports_root> [--json out.json]
"""
# ruff: noqa: E501  (dense print/format lines)
from __future__ import annotations

import argparse
import csv
import json
import statistics as st
import sys
from pathlib import Path

csv.field_size_limit(10 ** 7)
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "structural_elucidation" / "reports"))
import confidence as C  # noqa: E402

# Confidence ordering (lower index = higher identification confidence).
SUBLEVEL_ORDER = ["1", "2A", "2B", "2C", "3A", "3B", "3C", "3D", "4", "5"]
ORDER_RANK = {s: i for i, s in enumerate(SUBLEVEL_ORDER)}


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


def _rank1(rows: list[dict]) -> dict:
    best, br = {}, None
    for r in rows:
        k = _num(r.get("rank"))
        k = 1e9 if k is None else k
        if br is None or k < br:
            br, best = k, r
    return best


def report_features(report_dir: Path) -> dict | None:
    data = report_dir / "data"
    if not data.is_dir():
        return None
    slug = report_dir.name.lower()
    kx = _read(data / f"{slug}_known_x_row.csv")
    curated = str(kx[0].get("SubLevel", "")).strip() if kx else ""
    ci_path = data / f"{slug}_confidence_inputs.json"
    tc = _read(data / f"{slug}_top_candidates.csv")
    if not ci_path.exists() or not tc:
        return None
    try:
        fields = json.loads(ci_path.read_text())
    except (ValueError, OSError):
        return None
    r1 = _rank1(tc)
    # same-formula multiplicity among the top candidates (positional-isomer signal)
    r1_formula = str(r1.get("predicted_formula", "") or "").strip()
    n_isomers = sum(1 for r in tc[:10] if str(r.get("predicted_formula", "") or "").strip() == r1_formula and r1_formula)
    conf = C.compute_confidence(fields)
    return {
        "control_id": report_dir.name,
        "curated": curated,
        "struct": _num(fields.get("structure_similarity")) or 0.0,
        "mass": _num(fields.get("mass_match_score")),
        "raw_evidence": conf.raw_evidence,
        "ms2_present": bool(fields.get("ms2_spectrum_present")),
        "ms2_score": _num(fields.get("ms2_score")) or 0.0,
        "observed": bool(fields.get("observed_mz_available")),
        "iso_score": _num(r1.get("isotope_formula_score")),
        "iso_like": str(r1.get("feature_is_likely_isotopologue", "")).strip().lower() in ("true", "1"),
        "homolog": bool(str(r1.get("conjugate_name", "")).strip())
        or str(r1.get("is_multi_conjugate", "")).strip().lower() in ("true", "1"),
        "n_isomers": n_isomers,
        "truth": _num(r1.get("truth_tanimoto")),
        "source": str(fields.get("source", "")).strip(),
    }


def derive_sublevel(f: dict) -> str:
    """Proposed evidence-based derived sublevel (aligned to the rubric intent).

    Uses only signals we actually have: accurate mass, isotope pattern, structure
    similarity, same-formula multiplicity (positional-isomer signal), observed-MS2
    presence, and conjugate/homologue hints. 2A/2B are intentionally not derived
    (require spectral-library MS2 + RT, unavailable). L1 is curated-only.
    """
    mass = f["mass"]
    mass_ok = mass is not None and mass >= C.MASS_STRONG
    mass_bad = mass is not None and mass < C.MASS_WEAK
    iso_ok = (f["iso_score"] is None) or (f["iso_score"] >= 0.5)
    struct = f["struct"]
    ms2 = f["ms2_present"] and f["ms2_score"] >= C.MS2_PRESENT

    if mass_bad or not f["observed"]:
        return "4" if f["observed"] else "5"
    if struct >= C.STRUCT_STRONG and mass_ok and iso_ok and (ms2 or f["homolog"]):
        return "2C"
    if mass_ok and iso_ok and struct >= C.STRUCT_MODERATE and f["n_isomers"] >= 2:
        return "3A"
    if mass_ok and struct >= 0.55 and f["ms2_present"]:
        return "3B"
    if mass_ok and struct >= C.STRUCT_WEAK:
        return "3C"
    if f["homolog"]:
        return "3D"
    return "4"


def _median(xs):
    xs = [x for x in xs if x is not None]
    return round(st.median(xs), 3) if xs else None


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("reports_root", type=Path, nargs="?", default=Path("_audit_pull/reports"))
    ap.add_argument("--json", type=Path, default=None)
    args = ap.parse_args()

    rows = []
    for d in sorted(p for p in args.reports_root.iterdir() if p.is_dir() and p.name.upper().startswith("X-")):
        rec = report_features(d)
        if rec:
            rec["derived"] = derive_sublevel(rec)
            rows.append(rec)

    curated_rows = [r for r in rows if r["curated"] in ORDER_RANK and r["curated"] not in ("1",)]
    feats = ["struct", "mass", "raw_evidence", "ms2_score", "iso_score", "n_isomers"]
    dist = {}
    for sl in ["2C", "3A", "3B", "3C", "3D", "4", "5"]:
        sub = [r for r in curated_rows if r["curated"] == sl]
        if not sub:
            continue
        dist[sl] = {"n": len(sub), **{k: _median([r[k] for r in sub]) for k in feats},
                    "ms2_present_rate": round(sum(r["ms2_present"] for r in sub) / len(sub), 2)}

    # agreement: derived vs curated (exact + within-one-order)
    comparable = [r for r in curated_rows if r["curated"] in ("2C", "3A", "3B", "3C", "3D", "4", "5")]
    exact = sum(r["derived"] == r["curated"] for r in comparable)
    within1 = sum(abs(ORDER_RANK[r["derived"]] - ORDER_RANK[r["curated"]]) <= 1 for r in comparable)
    confusion: dict = {}
    for r in comparable:
        confusion.setdefault(r["curated"], {}).setdefault(r["derived"], 0)
        confusion[r["curated"]][r["derived"]] += 1

    # defensibility: does derived sublevel ORDER track correctness (truth Tanimoto)?
    truth = [r for r in rows if r["truth"] is not None]
    by_derived_truth = {}
    for sl in ["2C", "3A", "3B", "3C", "3D", "4", "5"]:
        sub = [r["truth"] for r in truth if r["derived"] == sl]
        if sub:
            by_derived_truth[sl] = {"n": len(sub), "median_truth": round(st.median(sub), 3),
                                    "frac_correct_ge_0.7": round(sum(t >= 0.7 for t in sub) / len(sub), 3)}
    # Spearman-ish: correlation of derived order rank vs truth (lower order rank should = higher truth)
    pairs = [(ORDER_RANK[r["derived"]], r["truth"]) for r in truth if r["derived"] in ORDER_RANK]
    mono = None
    if len(pairs) > 10:
        conc = disc = 0
        for i in range(len(pairs)):
            for j in range(i + 1, len(pairs)):
                a, b = pairs[i], pairs[j]
                if a[0] == b[0] or a[1] == b[1]:
                    continue
                # expect: smaller order rank -> larger truth
                agree = (a[0] < b[0]) == (a[1] > b[1])
                conc += agree
                disc += not agree
        mono = round((conc - disc) / (conc + disc), 3) if (conc + disc) else None

    out = {
        "n_reports": len(rows),
        "n_curated": len(curated_rows),
        "curated_feature_medians": dist,
        "agreement": {
            "n": len(comparable),
            "exact": round(exact / len(comparable), 3) if comparable else None,
            "within_one_order": round(within1 / len(comparable), 3) if comparable else None,
            "confusion": confusion,
        },
        "defensibility_truth": {
            "by_derived": by_derived_truth,
            "order_vs_truth_tau": mono,
        },
    }
    if args.json:
        args.json.write_text(json.dumps(out, indent=2))

    print(f"[sublevel-audit] reports={out['n_reports']} curated={out['n_curated']}")
    print("  curated feature medians (do axes separate sublevels?):")
    for sl, d in dist.items():
        print(f"    {sl:3s} n={d['n']:4d} struct={d['struct']} mass={d['mass']} raw={d['raw_evidence']} "
              f"ms2={d['ms2_score']} iso={d['iso_score']} n_isomers={d['n_isomers']} ms2_present={d['ms2_present_rate']}")
    ag = out["agreement"]
    print(f"  derived vs curated: exact={ag['exact']} within-one-order={ag['within_one_order']} (n={ag['n']})")
    dt = out["defensibility_truth"]
    print(f"  derived order vs truth Tanimoto tau={dt['order_vs_truth_tau']} (higher sublevel -> more correct?)")
    for sl, d in dt["by_derived"].items():
        print(f"    derived {sl:3s} n={d['n']:4d} median_truth={d['median_truth']} correct>=0.7={d['frac_correct_ge_0.7']}")


if __name__ == "__main__":
    main()
