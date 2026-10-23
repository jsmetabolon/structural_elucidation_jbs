#!/usr/bin/env python3
"""Benchmark mummichog-style m/z pathway assignment for X-compounds.

Investigates whether mummichog (m/z + metabolic-network based pathway
inference, no structure required) would assign Metabolon super-pathways to
X-compounds better than the production structure-based predictor
(``structural_elucidation/pipeline/enrichment/pathway.py``, audited by
``tools/audit_pathway_predictions.py``).

Benchmark arms
--------------
- Arm 1 (baseline): the committed structure-based audit,
  ``docs/paper/data/pathway/audit_pathway_accuracy.csv`` (per-control rank-1
  predictions vs known_X truth). Not recomputed here.
- Arm 2 (mummichog): each X-compound's parent m/z is matched against the
  mummichog v2.7 human MFN metabolic model using mummichog's own adduct rules
  and 10 ppm tolerance; candidate compounds vote a Metabolon super-pathway via
  a hand-reviewed crosswalk (``tools/data/mummichog_pathway_crosswalk.csv``).
  Two variants: ``primary`` (primary ions only, the faithful analog of
  mummichog's force_primary_ion for single-feature evidence) and ``all``
  (every adduct in mummichog's wanted list for the ionization mode).
- Arm 3 (mass-only ablation, gated): the same neutral mass matched against the
  Metabolon metabolite synonym reference DB (the label source the baseline
  votes over). Requires the fetched SQLite asset; skipped when absent.

Ground truth and scoring conventions (evaluable set, label normalization,
Wilson CIs) mirror ``tools/audit_pathway_predictions.py``.

Usage::

    python tools/mummichog_pathway_benchmark.py inventory
    python tools/mummichog_pathway_benchmark.py arm2 \
        --out docs/paper/data/pathway/mummichog_arm2_predictions.csv
    python tools/mummichog_pathway_benchmark.py arm3 \
        --out docs/paper/data/pathway/mummichog_arm3_predictions.csv
    python tools/mummichog_pathway_benchmark.py score \
        --arm2 docs/paper/data/pathway/mummichog_arm2_predictions.csv \
        --out-json docs/paper/data/pathway/mummichog_eval_summary.json \
        --out-md docs/paper/data/pathway/mummichog_eval_summary.md
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_LIBS = ("209", "305", "400", "402")
PPM_DEFAULT = 10.0
FEATURE_DIR = REPO_ROOT / "data" / "feature_rows"
KNOWN_X = REPO_ROOT / "data" / "known_X.csv"
KNOWN_STRUCTURES = REPO_ROOT / "data" / "known_X_all_available_structures_simple.csv"
BASELINE_AUDIT = REPO_ROOT / "docs" / "paper" / "data" / "pathway" / "audit_pathway_accuracy.csv"
CROSSWALK = REPO_ROOT / "tools" / "data" / "mummichog_pathway_crosswalk.csv"
REFERENCE_DB_CANDIDATES = (
    REPO_ROOT / "Metabolite Synonym Database" / "metabolite_synonyms_merged_with_mtx_giga.sqlite",
    REPO_ROOT / "vendor" / "metabolite_synonyms.sqlite",
)
# mummichog ionization-mode keys per library ionization ('pos_default' is the
# mummichog CLI default for positive data).
MODE_BY_ION_MODE = {"positive": "pos_default", "negative": "negative"}
VARIANTS = ("primary", "all")


def _load_audit_module():
    """Load tools/audit_pathway_predictions.py for its shared helpers."""

    path = REPO_ROOT / "tools" / "audit_pathway_predictions.py"
    spec = importlib.util.spec_from_file_location("audit_pathway_predictions", path)
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    return module


_AUDIT = _load_audit_module()
norm = _AUDIT.norm
clean = _AUDIT.clean
wilson = _AUDIT.wilson
EXCLUDE_TRUTH = _AUDIT.EXCLUDE_TRUTH


# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------

def load_crosswalk(path: Path = CROSSWALK) -> dict[str, str]:
    """MFN pathway name -> Metabolon super-pathway. Raises on duplicates."""

    mapping: dict[str, str] = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            name = row["mfn_pathway"].strip()
            super_ = row["metabolon_super"].strip()
            if not name or not super_:
                raise ValueError(f"crosswalk row missing data: {row}")
            if name in mapping:
                raise ValueError(f"duplicate crosswalk pathway: {name}")
            mapping[name] = super_
    return mapping


def load_features(feature_dir: Path = FEATURE_DIR, libs: tuple[str, ...] = DEFAULT_LIBS) -> pd.DataFrame:
    """One feature row per control: first-lib-wins in audit lib order."""

    frames = []
    for lib in libs:
        for path in sorted(feature_dir.glob(f"lib{lib}_*_feature_rows.csv")):
            df = pd.read_csv(path, low_memory=False)
            frames.append(df[["control_id", "lib_id", "ion_mode", "parent_mz", "neutral_mass_fallback"]])
    if not frames:
        raise FileNotFoundError(f"no feature_rows CSVs under {feature_dir}")
    feat = pd.concat(frames, ignore_index=True)
    feat = feat.drop_duplicates(subset=["control_id"], keep="first").reset_index(drop=True)
    return feat


def load_truth(known_x_path: Path = KNOWN_X) -> pd.DataFrame:
    """control_id -> (true_super, true_sub), first occurrence wins (audit convention)."""

    kx = pd.read_csv(known_x_path, low_memory=False)
    out: dict[str, dict[str, str]] = {}
    for _, r in kx.iterrows():
        cid = clean(r.get("Unknown"))
        if not cid or cid in out:
            continue
        out[cid] = {"true_super": clean(r.get("SuperPathway")), "true_sub": clean(r.get("Draft Subpathway"))}
    return pd.DataFrame([{"control_id": k, **v} for k, v in out.items()])


def evaluable_truth(truth: pd.DataFrame) -> pd.DataFrame:
    """Rows whose super label is real (non-empty, non-sentinel)."""

    sup = truth["true_super"].map(norm)
    return truth[sup.ne("") & ~sup.isin(EXCLUDE_TRUTH)].copy()


def load_known_structures(path: Path = KNOWN_STRUCTURES) -> dict[str, str]:
    """control_id -> SMILES of the curated known identity."""

    df = pd.read_csv(path, low_memory=False)
    out: dict[str, str] = {}
    for _, r in df.iterrows():
        cid = clean(r.get("x_compound_id"))
        smi = clean(r.get("smiles"))
        if cid and smi and cid not in out:
            out[cid] = smi
    return out


# ---------------------------------------------------------------------------
# Mummichog model access (Arm 2)
# ---------------------------------------------------------------------------

def _mummichog_imports():
    from mummichog.config import MASS_RANGE, primary_ions, wanted_adduct_list
    from mummichog.ions import compute_isotopes_adducts
    from mummichog.JSON_metabolicModels import metabolicModels

    return metabolicModels, compute_isotopes_adducts, MASS_RANGE, primary_ions, wanted_adduct_list


IonIndex = list[list[tuple[str, str, float]]]


def build_ion_index(mode: str, restrict_primary: bool) -> tuple[IonIndex, dict[str, list[str]], dict[str, float]]:
    """Integer-mass-binned (compound, ion, theoretical m/z) index for one mode.

    Mirrors mummichog's ``IndexedTable.__build_cpdindex__``: adducts computed
    with mummichog's own rules, restricted to its wanted list for the mode
    (optionally to its primary ions), binned by int(mass) within MASS_RANGE.
    Returns (index, cpd2pathways, cpd_mw).
    """

    metabolicModels, compute_isotopes_adducts, MASS_RANGE, primary_ions, wanted_adduct_list = _mummichog_imports()
    model = metabolicModels["human_model_mfn"]
    wanted = set(wanted_adduct_list[mode])
    if restrict_primary:
        wanted &= set(primary_ions)
    index: IonIndex = [[] for _ in range(MASS_RANGE[1] + 2)]
    cpd_mw: dict[str, float] = {}
    for cpd_id, cpd in model["Compounds"].items():
        if not cpd.get("mw"):
            continue
        cpd_mw[cpd_id] = float(cpd["mw"])
        adducts = compute_isotopes_adducts(cpd, mode)
        for ion, mass in adducts.items():
            if ion in wanted and MASS_RANGE[0] < mass < MASS_RANGE[1]:
                index[int(mass)].append((cpd_id, ion, mass))
    return index, model["cpd2pathways"], cpd_mw


def match_mz(mz: float, index: IonIndex, ppm: float = PPM_DEFAULT) -> list[tuple[str, str, float]]:
    """Candidate (compound, ion, theoretical m/z) tuples within ppm of mz.

    Mirrors mummichog's ``__match_mz_ion__`` (integer bins floor-1..floor+1,
    flat ppm tolerance).
    """

    floor = int(mz)
    tol = 1e-6 * ppm * mz
    matched = []
    for ii in (floor - 1, floor, floor + 1):
        if 0 <= ii < len(index):
            for entry in index[ii]:
                if abs(entry[2] - mz) < tol:
                    matched.append(entry)
    return matched


def vote_super(
    candidate_cpds: list[str],
    cpd2pathways: dict[str, list[str]],
    crosswalk: dict[str, str],
) -> dict[str, Any]:
    """Super-pathway vote over candidate compounds.

    Each candidate compound with >=1 mapped pathway contributes weight 1,
    split equally across the distinct Metabolon supers its pathways map to.
    Deterministic tie-break: alphabetical; ties are flagged.
    """

    scores: dict[str, float] = defaultdict(float)
    n_voting = 0
    pathway_names: Counter[str] = Counter()
    for cpd in candidate_cpds:
        pathways = cpd2pathways.get(cpd) or []
        supers = sorted({crosswalk[p] for p in pathways if p in crosswalk})
        if not supers:
            continue
        n_voting += 1
        for s in supers:
            scores[s] += 1.0 / len(supers)
        pathway_names.update(pathways)
    if not scores:
        return {
            "predicted_super": "",
            "vote_share": float("nan"),
            "n_voting": 0,
            "tie": False,
            "top_pathways": "",
        }
    total = sum(scores.values())
    best = max(scores.values())
    winners = sorted(s for s, v in scores.items() if v == best)
    top_pathways = "; ".join(f"{name} ({count})" for name, count in pathway_names.most_common(3))
    return {
        "predicted_super": winners[0],
        "vote_share": round(best / total, 4),
        "n_voting": n_voting,
        "tie": len(winners) > 1,
        "top_pathways": top_pathways,
    }


def _true_mass_diagnostics(
    smiles: str,
    candidate_cpds: list[str],
    cpd_mw: dict[str, float],
    ppm: float,
) -> dict[str, Any]:
    """Is the curated-known structure's neutral mass represented in the model /
    among the matched candidates? Diagnoses coverage vs voting failures."""

    from rdkit import Chem
    from rdkit.Chem.Descriptors import ExactMolWt

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"true_neutral_mass": float("nan"), "model_has_true_mass": None, "candidates_include_true_mass": None}
    true_mass = ExactMolWt(mol)
    tol = 1e-6 * ppm * true_mass
    model_has = any(abs(mw - true_mass) < tol for mw in cpd_mw.values())
    cand_has = any(abs(cpd_mw.get(c, -1e9) - true_mass) < tol for c in candidate_cpds)
    return {
        "true_neutral_mass": round(true_mass, 5),
        "model_has_true_mass": model_has,
        "candidates_include_true_mass": cand_has,
    }


def run_arm2(out_path: Path, ppm: float = PPM_DEFAULT) -> pd.DataFrame:
    """Mummichog MFN mass-match + crosswalk vote for every X-compound feature."""

    crosswalk = load_crosswalk()
    features = load_features()
    known_smiles = load_known_structures()

    indexes: dict[tuple[str, bool], Any] = {}
    rows: list[dict[str, Any]] = []
    for _, feat in features.iterrows():
        mode = MODE_BY_ION_MODE.get(str(feat["ion_mode"]))
        if mode is None:
            continue
        for variant in VARIANTS:
            key = (mode, variant == "primary")
            if key not in indexes:
                indexes[key] = build_ion_index(mode, restrict_primary=(variant == "primary"))
            index, cpd2pathways, cpd_mw = indexes[key]
            matched = match_mz(float(feat["parent_mz"]), index, ppm)
            candidate_cpds = sorted({m[0] for m in matched})
            vote = vote_super(candidate_cpds, cpd2pathways, crosswalk)
            row: dict[str, Any] = {
                "control_id": feat["control_id"],
                "lib": feat["lib_id"],
                "ion_mode": feat["ion_mode"],
                "mz": feat["parent_mz"],
                "variant": variant,
                "n_candidate_compounds": len(candidate_cpds),
                "matched_ions": "; ".join(sorted({m[1] for m in matched})),
                "matched_compounds": "; ".join(candidate_cpds[:25]),
                **vote,
            }
            smiles = known_smiles.get(str(feat["control_id"]), "")
            if smiles:
                row.update(_true_mass_diagnostics(smiles, candidate_cpds, cpd_mw, ppm))
            rows.append(row)
    pred = pd.DataFrame(rows)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    pred.to_csv(out_path, index=False)
    print(f"[arm2] wrote {len(pred)} rows ({features.shape[0]} controls x {len(VARIANTS)} variants) -> {out_path}")
    return pred


# ---------------------------------------------------------------------------
# Arm 3: mass-only vote against the Metabolon synonym reference DB (gated)
# ---------------------------------------------------------------------------

def find_reference_db() -> Path | None:
    for candidate in REFERENCE_DB_CANDIDATES:
        if candidate.exists():
            return candidate
    return None


def run_arm3(out_path: Path, ppm: float = PPM_DEFAULT) -> pd.DataFrame | None:
    """Neutral-mass +/- ppm vote against the Metabolon synonym DB labels."""

    db_path = find_reference_db()
    if db_path is None:
        print("[arm3] SKIPPED: metabolite synonym SQLite not present "
              f"(looked for {[str(p) for p in REFERENCE_DB_CANDIDATES]}). "
              "Fetch it with ./fetch_assets.sh and re-run.")
        return None

    import sqlite3

    from rdkit import Chem
    from rdkit.Chem.Descriptors import ExactMolWt

    con = sqlite3.connect(str(db_path))
    ref = pd.read_sql_query(
        "SELECT smiles, metabolon_super, metabolon_sub FROM metabolites "
        "WHERE TRIM(COALESCE(smiles,'')) <> '' AND TRIM(COALESCE(metabolon_super,'')) <> ''",
        con,
    )
    con.close()

    entries: list[tuple[float, str, str]] = []  # (mass, super, sub)
    seen: set[str] = set()
    for _, r in ref.iterrows():
        smi = str(r["smiles"]).strip()
        if smi in seen:
            continue
        seen.add(smi)
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        entries.append((ExactMolWt(mol), clean(r["metabolon_super"]), clean(r["metabolon_sub"])))
    entries.sort(key=lambda e: e[0])
    masses = [e[0] for e in entries]

    import bisect

    features = load_features()
    rows: list[dict[str, Any]] = []
    for _, feat in features.iterrows():
        neutral = float(feat["neutral_mass_fallback"])
        tol = 1e-6 * ppm * neutral
        lo = bisect.bisect_left(masses, neutral - tol)
        hi = bisect.bisect_right(masses, neutral + tol)
        hits = entries[lo:hi]
        super_scores: dict[str, float] = defaultdict(float)
        for _, sup, _sub in hits:
            super_scores[sup] += 1.0
        if super_scores:
            total = sum(super_scores.values())
            best = max(super_scores.values())
            winners = sorted(s for s, v in super_scores.items() if v == best)
            pred_super = winners[0]
            share = round(best / total, 4)
            sub_scores: Counter[str] = Counter(
                _sub for _, sup, _sub in hits if sup == pred_super and norm(_sub)
            )
            pred_sub = sub_scores.most_common(1)[0][0] if sub_scores else ""
            tie = len(winners) > 1
        else:
            pred_super, pred_sub, share, tie = "", "", float("nan"), False
        rows.append(
            {
                "control_id": feat["control_id"],
                "lib": feat["lib_id"],
                "ion_mode": feat["ion_mode"],
                "neutral_mass": neutral,
                "n_candidates": len(hits),
                "predicted_super": pred_super,
                "predicted_sub": pred_sub,
                "vote_share": share,
                "tie": tie,
            }
        )
    pred = pd.DataFrame(rows)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    pred.to_csv(out_path, index=False)
    print(f"[arm3] wrote {len(pred)} rows -> {out_path}")
    return pred


# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------

def mcnemar_exact(b: int, c: int) -> float:
    """Two-sided exact McNemar p-value on discordant pair counts (binomial)."""

    n = b + c
    if n == 0:
        return 1.0
    k = min(b, c)
    p = sum(math.comb(n, i) for i in range(0, k + 1)) / 2.0 ** n
    return min(1.0, 2.0 * p)


def _accuracy_block(df: pd.DataFrame) -> dict[str, Any]:
    return wilson(int(df["correct_super"].sum()), len(df))


def score_variant(pred: pd.DataFrame, truth: pd.DataFrame, baseline: pd.DataFrame) -> dict[str, Any]:
    """Coverage/accuracy metrics for one Arm-2 variant, plus baseline head-to-head."""

    ev_truth = evaluable_truth(truth)
    joined = pred.merge(ev_truth, on="control_id", how="inner").copy()
    joined["has_prediction"] = joined["predicted_super"].map(norm).ne("")
    called = joined[joined["has_prediction"]].copy()
    called["correct_super"] = [
        norm(p) == norm(t) for p, t in zip(called["predicted_super"], called["true_super"])
    ]

    result: dict[str, Any] = {
        "n_controls_with_mz": int(pred["control_id"].nunique()),
        "n_evaluable_with_mz": int(len(joined)),
        "n_with_any_candidate": int((joined["n_candidate_compounds"] > 0).sum()),
        "n_called": int(len(called)),
        "coverage_called_of_evaluable": round(len(called) / len(joined), 4) if len(joined) else 0.0,
        "super_accuracy_on_called": _accuracy_block(called) if len(called) else {"n": 0},
        "accuracy_by_true_super": {
            str(k): _accuracy_block(g) for k, g in called.groupby("true_super") if len(g) >= 5
        },
        "candidate_count_pctiles": {
            p: float(joined["n_candidate_compounds"].quantile(q))
            for p, q in (("p10", 0.1), ("p50", 0.5), ("p90", 0.9))
        },
    }

    # Confusion (top true supers) and predicted-vs-truth distribution.
    confusion: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    for _, r in called.iterrows():
        confusion[clean(r["true_super"])][clean(r["predicted_super"])] += 1
    top_true = [t for t, _ in Counter(called["true_super"].map(clean)).most_common(8)]
    result["confusion_top_true_super"] = {
        t: dict(sorted(confusion[t].items(), key=lambda kv: kv[1], reverse=True)[:5]) for t in top_true if t
    }
    predicted_counts = Counter(called["predicted_super"].map(clean)).most_common()
    truth_counts = Counter(called["true_super"].map(clean)).most_common()
    result["bias"] = {
        "predicted_super_distribution": {k: int(v) for k, v in predicted_counts},
        "truth_super_distribution": {k: int(v) for k, v in truth_counts},
    }

    # Accuracy/coverage trade-off across vote-share thresholds.
    curve = []
    for thr in (0.0, 0.2, 0.4, 0.5, 0.6, 0.8, 0.9):
        sel = called[called["vote_share"] >= thr]
        if len(sel) == 0:
            continue
        curve.append(
            {
                "vote_share_min": thr,
                "n": int(len(sel)),
                "coverage_of_evaluable": round(len(sel) / len(joined), 4),
                "accuracy": round(float(sel["correct_super"].mean()), 4),
            }
        )
    result["vote_share_curve"] = curve

    # Head-to-head vs the structure-based baseline on the joint called set.
    base = baseline[["control_id", "predicted_pathway", "true_super", "correct_super", "tier"]].copy()
    base = base.rename(
        columns={
            "predicted_pathway": "baseline_super",
            "correct_super": "baseline_correct",
            "true_super": "baseline_truth",
        }
    )
    joint = called.merge(base, on="control_id", how="inner")
    if len(joint):
        joint["baseline_correct"] = joint["baseline_correct"].astype(bool)
        b = int((joint["correct_super"] & ~joint["baseline_correct"]).sum())  # mummichog only
        c = int((~joint["correct_super"] & joint["baseline_correct"]).sum())  # baseline only
        both = int((joint["correct_super"] & joint["baseline_correct"]).sum())
        neither = int((~joint["correct_super"] & ~joint["baseline_correct"]).sum())
        result["head_to_head_vs_baseline"] = {
            "n_joint": int(len(joint)),
            "mummichog_accuracy": round(float(joint["correct_super"].mean()), 4),
            "baseline_accuracy": round(float(joint["baseline_correct"].mean()), 4),
            "both_correct": both,
            "mummichog_only_correct": b,
            "baseline_only_correct": c,
            "neither_correct": neither,
            "mcnemar_p_two_sided": float(f"{mcnemar_exact(b, c):.3g}"),
            "baseline_accuracy_by_tier_on_joint": {
                str(t): wilson(int(g["baseline_correct"].sum()), len(g)) for t, g in joint.groupby("tier")
            },
            "mummichog_accuracy_by_baseline_tier_on_joint": {
                str(t): wilson(int(g["correct_super"].sum()), len(g)) for t, g in joint.groupby("tier")
            },
        }

    # Mechanistic diagnostics on the curated-known subset.
    if "model_has_true_mass" in joined.columns:
        diag = joined[joined["model_has_true_mass"].notna()]
        if len(diag):
            in_model = diag[diag["model_has_true_mass"] == True]  # noqa: E712
            result["known_structure_diagnostics"] = {
                "n_known_structures_evaluated": int(len(diag)),
                "true_mass_in_model": int(len(in_model)),
                "true_mass_in_model_frac": round(len(in_model) / len(diag), 4),
                "true_mass_among_candidates_when_in_model": (
                    round(float((in_model["candidates_include_true_mass"] == True).mean()), 4)  # noqa: E712
                    if len(in_model)
                    else None
                ),
            }
    return result


def score_arm3(pred: pd.DataFrame, truth: pd.DataFrame) -> dict[str, Any]:
    """Super/sub accuracy for the mass-only reference-DB vote."""

    ev_truth = evaluable_truth(truth)
    joined = pred.merge(ev_truth, on="control_id", how="inner").copy()
    called = joined[joined["predicted_super"].map(norm).ne("")].copy()
    called["correct_super"] = [
        norm(p) == norm(t) for p, t in zip(called["predicted_super"], called["true_super"])
    ]
    sub_ev = called[called["true_sub"].map(norm).ne("")].copy()
    sub_ev["correct_sub"] = [
        norm(p) == norm(t) for p, t in zip(sub_ev["predicted_sub"], sub_ev["true_sub"])
    ]
    return {
        "n_evaluable_with_mz": int(len(joined)),
        "n_called": int(len(called)),
        "coverage_called_of_evaluable": round(len(called) / len(joined), 4) if len(joined) else 0.0,
        "super_accuracy_on_called": _accuracy_block(called) if len(called) else {"n": 0},
        "sub_accuracy_on_called": wilson(int(sub_ev["correct_sub"].sum()), len(sub_ev)) if len(sub_ev) else {"n": 0},
    }


def run_inventory() -> dict[str, Any]:
    """Data-coverage report: which evaluable controls have m/z inputs."""

    features = load_features()
    truth = load_truth()
    ev = evaluable_truth(truth)
    baseline = pd.read_csv(BASELINE_AUDIT)
    feat_controls = set(features["control_id"])
    ev_controls = set(ev["control_id"])
    base_controls = set(baseline["control_id"])
    inv = {
        "n_controls_with_feature_mz": len(feat_controls),
        "ion_mode_counts": {k: int(v) for k, v in features["ion_mode"].value_counts().items()},
        "n_truth_controls": int(len(truth)),
        "n_evaluable_super_truth": len(ev_controls),
        "n_evaluable_with_mz": len(ev_controls & feat_controls),
        "n_baseline_audit_controls": len(base_controls),
        "n_baseline_with_mz": len(base_controls & feat_controls),
        "n_joint_baseline_evaluable_mz": len(base_controls & ev_controls & feat_controls),
    }
    print(json.dumps(inv, indent=2))
    return inv


def run_score(arm2_path: Path, arm3_path: Path | None, out_json: Path, out_md: Path) -> dict[str, Any]:
    truth = load_truth()
    baseline = pd.read_csv(BASELINE_AUDIT)
    arm2 = pd.read_csv(arm2_path)

    summary: dict[str, Any] = {
        "inventory": run_inventory(),
        "baseline_reference": {
            "source": str(BASELINE_AUDIT.relative_to(REPO_ROOT)),
            "note": "structure-based predictor audit (Arm 1); see audit_pathway_summary.json",
        },
        "arm2_mummichog": {},
    }
    for variant in VARIANTS:
        sub = arm2[arm2["variant"] == variant].copy()
        summary["arm2_mummichog"][variant] = score_variant(sub, truth, baseline)

    if arm3_path is not None and arm3_path.exists():
        arm3 = pd.read_csv(arm3_path)
        summary["arm3_mass_only_reference_db"] = score_arm3(arm3, truth)
    else:
        summary["arm3_mass_only_reference_db"] = {
            "skipped": True,
            "reason": "metabolite synonym SQLite asset not available in this checkout",
        }

    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(summary, indent=2, default=str))
    out_md.write_text(render_markdown(summary))
    print(f"[score] wrote {out_json} and {out_md}")
    return summary


def render_markdown(summary: dict[str, Any]) -> str:
    """Compact human-readable summary of the benchmark."""

    inv = summary["inventory"]
    lines = [
        "# Mummichog X-compound pathway benchmark",
        "",
        f"- Evaluable controls (real known_X super label) with feature m/z: **{inv['n_evaluable_with_mz']}**",
        f"- Baseline (structure-based) audit controls: {inv['n_baseline_audit_controls']}"
        f" (with m/z: {inv['n_baseline_with_mz']})",
        "",
    ]
    for variant, s in summary["arm2_mummichog"].items():
        acc = s["super_accuracy_on_called"]
        lines += [
            f"## Arm 2 (mummichog MFN, variant = {variant})",
            "",
            f"- Called: {s['n_called']}/{s['n_evaluable_with_mz']} evaluable"
            f" (coverage {s['coverage_called_of_evaluable']:.1%})",
            f"- Super accuracy on called: **{acc.get('accuracy')}** (n={acc.get('n')}, 95% CI {acc.get('ci95')})",
        ]
        h2h = s.get("head_to_head_vs_baseline")
        if h2h:
            lines += [
                f"- Head-to-head (n={h2h['n_joint']}): mummichog {h2h['mummichog_accuracy']:.1%}"
                f" vs baseline {h2h['baseline_accuracy']:.1%};"
                f" discordant {h2h['mummichog_only_correct']} vs {h2h['baseline_only_correct']};"
                f" McNemar p = {h2h['mcnemar_p_two_sided']}",
            ]
        diag = s.get("known_structure_diagnostics")
        if diag:
            lines += [
                f"- Known-structure diagnostics: true mass in MFN model for"
                f" {diag['true_mass_in_model']}/{diag['n_known_structures_evaluated']}"
                f" ({diag['true_mass_in_model_frac']:.1%}); when in model, matched by the feature"
                f" {diag['true_mass_among_candidates_when_in_model']:.1%} of the time"
                if diag.get("true_mass_among_candidates_when_in_model") is not None
                else f"- Known-structure diagnostics: true mass in MFN model for"
                f" {diag['true_mass_in_model']}/{diag['n_known_structures_evaluated']}",
            ]
        lines += [""]
        curve = s.get("vote_share_curve") or []
        if curve:
            lines += ["| vote share >= | n | coverage | accuracy |", "|---|---|---|---|"]
            for row in curve:
                coverage = f"{row['coverage_of_evaluable']:.1%}"
                accuracy = f"{row['accuracy']:.1%}"
                lines.append(f"| {row['vote_share_min']} | {row['n']} | {coverage} | {accuracy} |")
            lines += [""]
    arm3 = summary.get("arm3_mass_only_reference_db", {})
    lines += ["## Arm 3 (mass-only vote vs Metabolon reference DB)", ""]
    if arm3.get("skipped"):
        lines += [f"- Skipped: {arm3.get('reason')}", ""]
    else:
        acc = arm3["super_accuracy_on_called"]
        sub = arm3["sub_accuracy_on_called"]
        lines += [
            f"- Called: {arm3['n_called']}/{arm3['n_evaluable_with_mz']}"
            f" (coverage {arm3['coverage_called_of_evaluable']:.1%})",
            f"- Super accuracy on called: **{acc.get('accuracy')}** (n={acc.get('n')})",
            f"- Sub accuracy on called: {sub.get('accuracy')} (n={sub.get('n')})",
            "",
        ]
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)

    sub.add_parser("inventory", help="report data coverage for the benchmark")

    p2 = sub.add_parser("arm2", help="run mummichog MFN mass-match + pathway vote")
    p2.add_argument("--out", type=Path, default=REPO_ROOT / "docs/paper/data/pathway/mummichog_arm2_predictions.csv")
    p2.add_argument("--ppm", type=float, default=PPM_DEFAULT)

    p3 = sub.add_parser("arm3", help="run mass-only vote against the Metabolon synonym DB (gated)")
    p3.add_argument("--out", type=Path, default=REPO_ROOT / "docs/paper/data/pathway/mummichog_arm3_predictions.csv")
    p3.add_argument("--ppm", type=float, default=PPM_DEFAULT)

    ps = sub.add_parser("score", help="score predictions and write summary JSON/MD")
    ps.add_argument("--arm2", type=Path, default=REPO_ROOT / "docs/paper/data/pathway/mummichog_arm2_predictions.csv")
    ps.add_argument("--arm3", type=Path, default=REPO_ROOT / "docs/paper/data/pathway/mummichog_arm3_predictions.csv")
    ps.add_argument("--out-json", type=Path, default=REPO_ROOT / "docs/paper/data/pathway/mummichog_eval_summary.json")
    ps.add_argument("--out-md", type=Path, default=REPO_ROOT / "docs/paper/data/pathway/mummichog_eval_summary.md")

    args = parser.parse_args()
    if args.command == "inventory":
        run_inventory()
    elif args.command == "arm2":
        run_arm2(args.out, args.ppm)
    elif args.command == "arm3":
        run_arm3(args.out, args.ppm)
    elif args.command == "score":
        run_score(args.arm2, args.arm3 if args.arm3.exists() else None, args.out_json, args.out_md)


if __name__ == "__main__":
    main()
