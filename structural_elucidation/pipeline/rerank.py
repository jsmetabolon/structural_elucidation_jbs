#!/usr/bin/env python3
"""Train and apply a compact unified cross-source reranker.

Fits a gradient-boosted or logistic classifier on pooled candidate features,
validates with leave-one-control-out CV (G3 gate), calibrates per truth-blind
source with isotonic regression, and emits top-10 ranked outputs.

Pipeline context:
    Consumes ``unified_candidates_per_control.csv`` from ``pool``. Trains on
    X-control labels (Tanimoto >= 0.7). Writes scored pool, OOF/full top-10,
    and generalization diagnostics.

Key I/O:
    Inputs: unified pool CSV, known_X truth structures, optional non-S0 pool.
    Outputs: ``unified_candidates_scored.csv``, ``unified_rerank_top10.csv``,
    ``unified_rerank_all_top10.csv``, isotonic calibrators, G3 reports.
"""

from __future__ import annotations

import argparse
import json
import pickle
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.isotonic import IsotonicRegression
from sklearn.linear_model import LogisticRegressionCV

from structural_elucidation.chemistry.diagnostics import TRUTH_BLIND_COMBO_SOURCES
from structural_elucidation.utils.denovo import canonicalize_smiles
from structural_elucidation.evaluation.failure_modes import best_tanimoto_for_smiles, load_x_controls
from structural_elucidation.evaluation.rerank_assembly import eval_reranked, write_denovo_report
from structural_elucidation.paths import REPO_ROOT, resolve
from structural_elucidation.utils.domain import dataframe_to_markdown, write_csv, write_json

# ── Reranker feature columns and defaults ────────────────────────────
FEATURE_COLUMNS: list[str] = [
    "predicted_spectrum_score_v2",
    "mass_match_score",
    "formula_prior_score",
    "evidence_confidence",
    "ms2_match_score",
    "nce_pipeline_score",
    "forward_spectrum_score",
    "rerank_score",
    "generation_score",
    "novelty_score",
    "min_source_rank_score",
    "n_source_rows_score",
    "is_peptide",
    "is_library_hmdb",
    "is_library_mlims",
    "is_nonpep_gated",
    "is_multi_conjugate",
    "is_nce_direct",
    "is_biotransformation",
    "isotope_formula_score",
    "isotope_ratio_score",
    "carbon_count_agreement",
]

WGCNA_FEATURE_COLUMNS: list[str] = [
    "wgcna_pathway_match",
    "wgcna_consensus_strength",
]

CLASS_MATCH_FEATURE_COLUMNS: list[str] = [
    "class_match_score",
]

FEATURE_DEFAULTS: dict[str, float] = {
    "isotope_formula_score": 0.5,
    "isotope_ratio_score": 0.5,
    "carbon_count_agreement": 0.5,
    "class_match_score": 0.5,
}


# ── Feature engineering ──────────────────────────────────────────────
def numeric(df: pd.DataFrame, col: str, default: float = 0.0) -> pd.Series:
    """Coerce column to float with fillna."""
    if col not in df.columns:
        return pd.Series(default, index=df.index, dtype=float)
    return pd.to_numeric(df[col], errors="coerce").fillna(default).astype(float)


def load_wgcna_consensus(
    path: Path | None = None, *, lib_id: int | None = None, ion_mode: str | None = None
) -> pd.DataFrame:
    """Load per-control WGCNA consensus pathway labels for soft-prior features."""
    if path is None:
        if lib_id is None or ion_mode is None:
            return pd.DataFrame()
        path = (
            Path("models")
            / "wgcna_evidence_v1"
            / f"lib_{int(lib_id)}_{str(ion_mode)}"
            / "wgcna_evidence_per_control.csv"
        )
    path = Path(path)
    if not path.exists():
        return pd.DataFrame()
    try:
        raw = pd.read_csv(path, low_memory=False)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()
    if raw.empty or "control_id" not in raw.columns:
        return pd.DataFrame()
    consensus = raw[raw.get("cohort", pd.Series("", index=raw.index)).astype(str).eq("consensus")].copy()
    if consensus.empty:
        return pd.DataFrame()
    keep = ["control_id", "consensus_super_pathway", "consensus_super_strength"]
    for col in keep:
        if col not in consensus.columns:
            consensus[col] = ""
    return consensus[keep].drop_duplicates("control_id", keep="first")


def _infer_lib_ion(df: pd.DataFrame) -> tuple[int | None, str | None]:
    """Infer lib_id and ion_mode from the first non-null pool row."""
    lib_id = None
    ion_mode = None
    if "lib_id" in df.columns and not df.empty:
        values = pd.to_numeric(df["lib_id"], errors="coerce").dropna()
        if not values.empty:
            lib_id = int(values.iloc[0])
    if "ion_mode" in df.columns and not df.empty:
        values = df["ion_mode"].dropna().astype(str)
        if not values.empty:
            ion_mode = str(values.iloc[0])
    return lib_id, ion_mode


def add_wgcna_features(df: pd.DataFrame, wgcna_consensus: pd.DataFrame | None = None) -> pd.DataFrame:
    """Add pathway-match and consensus-strength features from WGCNA evidence."""
    out = df.copy()
    if wgcna_consensus is None:
        lib_id, ion_mode = _infer_lib_ion(out)
        wgcna_consensus = load_wgcna_consensus(lib_id=lib_id, ion_mode=ion_mode)
    if wgcna_consensus is None or wgcna_consensus.empty or "control_id" not in out.columns:
        out["wgcna_pathway_match"] = 0.5
        out["wgcna_consensus_strength"] = 0.0
        return out
    merged = out.merge(wgcna_consensus, on="control_id", how="left")
    predicted = (
        merged.get("predicted_super_pathway", pd.Series("", index=merged.index))
        .fillna("")
        .astype(str)
        .str.strip()
        .str.lower()
    )
    consensus = (
        merged.get("consensus_super_pathway", pd.Series("", index=merged.index))
        .fillna("")
        .astype(str)
        .str.strip()
        .str.lower()
    )
    has_consensus = consensus.str.len().gt(0)
    has_prediction = predicted.str.len().gt(0)
    merged["wgcna_pathway_match"] = np.where(
        ~has_consensus, 0.5, np.where(has_prediction & predicted.eq(consensus), 1.0, 0.0)
    )  # 0.5 when no cohort consensus avoids penalizing sparse WGCNA coverage
    merged["wgcna_consensus_strength"] = (
        pd.to_numeric(merged.get("consensus_super_strength"), errors="coerce").fillna(0.0).clip(0.0, 1.0)
    )
    return merged.drop(
        columns=[col for col in ["consensus_super_pathway", "consensus_super_strength"] if col in merged.columns]
    )


def add_features(df: pd.DataFrame, wgcna_consensus: pd.DataFrame | None = None) -> pd.DataFrame:
    """Fill missing feature columns and derive rank/source-density scores."""
    out = df.copy()
    out["min_source_rank_score"] = 1.0 / numeric(out, "min_source_rank", default=999999).clip(lower=1.0)
    out["n_source_rows_score"] = np.log1p(numeric(out, "n_source_rows", default=1.0)).clip(lower=0.0, upper=5.0) / 5.0
    for col in FEATURE_COLUMNS:
        if col not in out.columns:
            out[col] = FEATURE_DEFAULTS.get(col, 0.0)
        default = FEATURE_DEFAULTS.get(col, 0.0)
        out[col] = numeric(out, col, default=default).replace([np.inf, -np.inf], default).fillna(default)
    return out


def scope_training_pool(
    pool: pd.DataFrame,
    controls: list[Any],
    *,
    lib_id: int,
    ion_mode: str,
    control_scope: str = "x_controls",
) -> pd.DataFrame:
    """Restrict training rows to the requested lib/ion and control scope."""
    out = pool.copy()
    truth_ids = {str(control.control_id) for control in controls}
    if "ion_mode" in out.columns:
        out = out[out["ion_mode"].astype(str).eq(str(ion_mode))].copy()
    if str(control_scope) == "x_controls" and "control_id" in out.columns:
        out = out[out["control_id"].astype(str).isin(truth_ids)].copy()
    if "lib_id" in out.columns:
        lib_values = pd.to_numeric(out["lib_id"], errors="coerce")
        if str(control_scope) == "x_controls":
            out = out[lib_values.eq(int(lib_id)) | lib_values.isna()].copy()
        else:
            out = out[lib_values.eq(int(lib_id))].copy()
    return out


def add_labels(pool: pd.DataFrame, controls: list[Any]) -> pd.DataFrame:
    """Attach truth Tanimoto and binary label (>= 0.7) per candidate row."""
    truth = {str(control.control_id): list(control.smiles) for control in controls}
    out = pool.copy()
    tans: list[float] = []
    labels: list[int] = []
    for _, row in out.iterrows():
        canonical = canonicalize_smiles(row.get("predicted_smiles", "")) or str(
            row.get("canonical_2d_smiles", "") or ""
        )
        _exact, tan = best_tanimoto_for_smiles(canonical, truth.get(str(row.get("control_id", "")), []))
        tans.append(float(tan))
        labels.append(int(float(tan) >= 0.7))
    out["truth_tanimoto"] = tans
    out["label_ge_0p7"] = labels
    return out


# ── Model factories and prediction ───────────────────────────────────
def make_model() -> GradientBoostingClassifier:
    """Return a shallow regularized GBM tuned for small X-control training sets."""
    return GradientBoostingClassifier(
        n_estimators=80,
        learning_rate=0.03,
        max_depth=2,
        min_samples_leaf=20,
        subsample=0.7,
        random_state=17,
    )


def make_model_b() -> LogisticRegressionCV:
    """Return L2 logistic regression with ROC-AUC CV as the linear baseline."""
    return LogisticRegressionCV(Cs=10, cv=5, penalty="l2", scoring="roc_auc", max_iter=2000, random_state=17)


def _predict_proba(model: Any, frame: pd.DataFrame) -> np.ndarray:
    """Predict positive-class probability, falling back to point predictions."""
    if hasattr(model, "predict_proba"):
        return model.predict_proba(frame[FEATURE_COLUMNS])[:, 1]
    return np.asarray(model.predict(frame[FEATURE_COLUMNS]), dtype=float)


# ── LOOCV evaluation and model selection ─────────────────────────────
def loocv_scores(pool: pd.DataFrame, model_factory=make_model) -> pd.Series:
    """Leave-one-control-out predicted probabilities; heuristic fallback when no positives."""
    scores = pd.Series(0.0, index=pool.index, dtype=float)
    controls = sorted(pool["control_id"].astype(str).unique())
    for control_id in controls:
        train = pool[~pool["control_id"].astype(str).eq(control_id)].copy()
        test = pool[pool["control_id"].astype(str).eq(control_id)].copy()
        if pd.to_numeric(test.get("label_ge_0p7"), errors="coerce").fillna(0).sum() <= 0:
            scores.loc[test.index] = heuristic_score(test)
            continue
        if train.empty or test.empty or train["label_ge_0p7"].nunique() < 2:
            scores.loc[test.index] = heuristic_score(test)
            continue
        model = model_factory()
        model.fit(train[FEATURE_COLUMNS], train["label_ge_0p7"])
        scores.loc[test.index] = _predict_proba(model, test)
    return scores


def rank1_hits_from_scores(pool: pd.DataFrame, score_col: str) -> int:
    """Count controls where the top-scored candidate has label_ge_0p7."""
    if pool.empty:
        return 0
    ranked = pool.sort_values(["control_id", score_col], ascending=[True, False], kind="mergesort").drop_duplicates(
        "control_id", keep="first"
    )
    return int(pd.to_numeric(ranked.get("label_ge_0p7"), errors="coerce").fillna(0).ge(1).sum())


def top10_hits_from_scores(pool: pd.DataFrame, score_col: str) -> int:
    """Count controls with a label_ge_0p7 hit in the top 10 by *score_col*."""
    if pool.empty:
        return 0
    ranked = pool.sort_values(["control_id", score_col], ascending=[True, False], kind="mergesort").copy()
    ranked["_score_rank"] = ranked.groupby("control_id", dropna=False).cumcount() + 1
    return int(ranked[ranked["_score_rank"].le(10)].groupby("control_id")["label_ge_0p7"].max().fillna(0).ge(1).sum())


def evaluate_model_option(pool: pd.DataFrame, name: str, model_factory) -> dict[str, Any]:
    """Fit one model variant and return OOF/full rank-1 and top-10 metrics."""
    work = pool.copy()
    work[f"{name}_oof_score"] = loocv_scores(work, model_factory=model_factory)
    if work["label_ge_0p7"].nunique() >= 2:
        model = model_factory()
        model.fit(work[FEATURE_COLUMNS], work["label_ge_0p7"])
        work[f"{name}_full_score"] = _predict_proba(model, work)
    else:
        model = None
        work[f"{name}_full_score"] = heuristic_score(work)
    full_rank1 = rank1_hits_from_scores(work, f"{name}_full_score")
    oof_rank1 = rank1_hits_from_scores(work, f"{name}_oof_score")
    return {
        "name": name,
        "model": model,
        "scores": work[[f"{name}_oof_score", f"{name}_full_score"]],
        "full_rank1": int(full_rank1),
        "oof_rank1": int(oof_rank1),
        "rank1_delta": int(abs(full_rank1 - oof_rank1)),
        "full_top10": int(top10_hits_from_scores(work, f"{name}_full_score")),
        "oof_top10": int(top10_hits_from_scores(work, f"{name}_oof_score")),
    }


def select_regularized_model(pool: pd.DataFrame) -> dict[str, Any]:
    """Pick the model with smallest full-vs-OOF rank-1 gap (legacy alias)."""
    options = [
        evaluate_model_option(pool, "option_a_regularized_gbm", make_model),
        evaluate_model_option(pool, "option_b_logistic_cv", make_model_b),
    ]
    selected = sorted(options, key=lambda item: (item["rank1_delta"], -item["oof_rank1"], item["name"]))[0]
    selected["all_options"] = options
    return selected


def select_model(pool: pd.DataFrame, mode: str = "auto") -> dict[str, Any]:
    """Select GBM, logistic, or auto-pick by LOOCV rank-1 stability."""
    options = [
        evaluate_model_option(pool, "option_a_regularized_gbm", make_model),
        evaluate_model_option(pool, "option_b_logistic_cv", make_model_b),
    ]
    if str(mode) == "gbm":
        selected = options[0]
    elif str(mode) == "logistic":
        selected = options[1]
    else:
        selected = sorted(options, key=lambda item: (item["rank1_delta"], -item["oof_rank1"], item["name"]))[0]
    selected["all_options"] = options
    return selected


def g3_gate(full_rank1: int, oof_rank1: int, min_floor_delta: int = 3) -> dict[str, Any]:
    """G3 validates rerank generalization against the rerank full-fit baseline.

    Pass requires rank-1 delta <= 2 and OOF rank-1 within *min_floor_delta*
    of the full-fit rank-1 (default: rerank baseline minus 3).

    Returns:
        Dict with ``g3_pass``, rank-1 counts, delta, and floor policy metadata.
    """
    delta = int(abs(int(full_rank1) - int(oof_rank1)))
    min_oof = int(full_rank1) - int(min_floor_delta)
    return {
        "g3_pass": bool(delta <= 2 and int(oof_rank1) >= int(min_oof)),
        "full_fit_rank1": int(full_rank1),
        "oof_rank1": int(oof_rank1),
        "rank1_delta_abs": int(delta),
        "min_oof_rank1": int(min_oof),
        "g3_floor_policy": "rerank_baseline_minus_3",
    }


# ── Per-source isotonic calibration ──────────────────────────────────
def source_mask(df: pd.DataFrame, source: str) -> pd.Series:
    """Boolean mask for rows whose source label contains *source*."""
    labels = df.get("source_label", df.get("source_labels", pd.Series("", index=df.index))).fillna("").astype(str)
    return labels.str.contains(source, regex=False)


def fit_isotonic_per_source(pool: pd.DataFrame, score_col: str = "unified_oof_score") -> dict[str, Any]:
    """Fit isotonic calibrators per truth-blind enumerator source."""
    calibrators: dict[str, Any] = {}
    labels = pd.to_numeric(pool.get("label_ge_0p7"), errors="coerce").fillna(0).astype(int)
    scores = pd.to_numeric(pool.get(score_col), errors="coerce").fillna(0.0)
    for source in sorted(TRUTH_BLIND_COMBO_SOURCES):
        mask = source_mask(pool, source)
        x = scores[mask].to_numpy(dtype=float)
        y = labels[mask].to_numpy(dtype=int)
        if len(x) >= 2 and len(set(y.tolist())) >= 2:
            iso = IsotonicRegression(out_of_bounds="clip")
            iso.fit(x, y)
            calibrators[source] = {"model": iso, "status": "fit", "n": int(len(x))}
        else:
            calibrators[source] = {"model": None, "status": "raw_fallback", "n": int(len(x))}
    return calibrators


def apply_isotonic_per_source(
    df: pd.DataFrame, calibrators: dict[str, Any], score_col: str = "unified_final_score"
) -> pd.Series:
    """Apply per-source isotonic maps; uncalibrated rows keep raw scores."""
    raw = pd.to_numeric(df.get(score_col), errors="coerce").fillna(0.0)
    out = raw.copy()
    for source, payload in calibrators.items():
        mask = source_mask(df, source)
        model = payload.get("model") if isinstance(payload, dict) else None
        if model is not None and mask.any():
            out.loc[mask] = model.predict(raw.loc[mask].to_numpy(dtype=float))
    return out.clip(0.0, 1.0)


def oof_drop_attribution(
    pool: pd.DataFrame, full_col: str = "unified_final_score", oof_col: str = "unified_oof_score"
) -> pd.DataFrame:
    """List controls where full-fit rank-1 hits but LOOCV rank-1 misses."""
    rows: list[dict[str, Any]] = []
    work = pool.copy()
    for control_id, group in work.groupby("control_id", dropna=False):
        full_top = group.sort_values(full_col, ascending=False, kind="mergesort").head(1)
        oof_top = group.sort_values(oof_col, ascending=False, kind="mergesort").head(1)
        if full_top.empty or oof_top.empty:
            continue
        full_hit = int(pd.to_numeric(full_top.get("label_ge_0p7"), errors="coerce").fillna(0).iloc[0] >= 1)
        oof_hit = int(pd.to_numeric(oof_top.get("label_ge_0p7"), errors="coerce").fillna(0).iloc[0] >= 1)
        if full_hit and not oof_hit:
            top = full_top.iloc[0]
            feature_values = {col: float(top.get(col, 0.0)) for col in FEATURE_COLUMNS}
            leaking = sorted(feature_values.items(), key=lambda item: abs(item[1]), reverse=True)[:3]
            rows.append(
                {
                    "control_id": control_id,
                    "full_fit_rank1": full_hit,
                    "oof_rank1": oof_hit,
                    "dominant_source": top.get("source_label", top.get("source_labels", "")),
                    "leaking_features": json.dumps([name for name, _value in leaking]),
                }
            )
    return pd.DataFrame(rows)


def heuristic_score(df: pd.DataFrame) -> pd.Series:
    """Weighted fallback score when supervised training is ill-conditioned."""
    return (
        0.34 * numeric(df, "predicted_spectrum_score_v2")
        + 0.22 * numeric(df, "mass_match_score")
        + 0.14 * numeric(df, "formula_prior_score")
        + 0.12 * numeric(df, "evidence_confidence")
        + 0.08 * numeric(df, "is_library_mlims")
        + 0.06 * numeric(df, "is_peptide")
        + 0.04 * numeric(df, "min_source_rank_score")
    ).clip(0.0, 1.0)


def rank_frame(scored: pd.DataFrame, score_col: str, top_k: int, policy: str = "unified_rerank_v1") -> pd.DataFrame:
    """Sort by rerank score and assign per-control ranks up to *top_k*."""
    out = scored.copy()
    out["rerank_policy"] = policy
    out["rerank_score"] = pd.to_numeric(out[score_col], errors="coerce").fillna(0.0)
    out = out.sort_values(
        [
            "control_id",
            "rerank_score",
            "predicted_spectrum_score_v2",
            "mass_match_score",
            "formula_prior_score",
            "min_source_rank",
            "predicted_smiles",
        ],
        ascending=[True, False, False, False, False, True, True],
        kind="mergesort",
    )
    out["feature_key"] = out.get("feature_key", pd.Series("", index=out.index)).fillna("").astype(str)
    missing_key = out["feature_key"].eq("")
    out.loc[missing_key, "feature_key"] = out.loc[missing_key, "control_id"].map(
        lambda value: f"unified_rerank_{value}"
    )
    out["rank"] = out.groupby("control_id", dropna=False).cumcount() + 1
    return out[out["rank"].le(int(top_k))].reset_index(drop=True)


# ── CLI entry point ──────────────────────────────────────────────────
def run(args: argparse.Namespace) -> dict[str, Any]:
    """Train/select reranker, score pool, emit top-10 and evaluation artifacts."""
    global FEATURE_COLUMNS
    repo_root = REPO_ROOT
    output_root = resolve(args.output_root, repo_root)
    output_root.mkdir(parents=True, exist_ok=True)
    wgcna_soft_prior = bool(getattr(args, "wgcna_soft_prior", False))
    if wgcna_soft_prior:
        FEATURE_COLUMNS = list(dict.fromkeys(list(FEATURE_COLUMNS) + WGCNA_FEATURE_COLUMNS))
    if bool(getattr(args, "class_match_score_feature", False)):
        FEATURE_COLUMNS = list(dict.fromkeys(list(FEATURE_COLUMNS) + CLASS_MATCH_FEATURE_COLUMNS))
    if bool(getattr(args, "drop_novelty_score", False)):
        FEATURE_COLUMNS = [col for col in FEATURE_COLUMNS if col != "novelty_score"]
    pool = pd.read_csv(resolve(args.unified_pool, repo_root), low_memory=False)
    controls = load_x_controls(
        known_x=repo_root / "known_X.csv",
        known_x_all=repo_root / "known_X_all_available_structures_simple.csv",
        lib_id=int(args.lib_id),
        ion_mode=str(args.ion_mode),
    )
    wgcna_consensus = load_wgcna_consensus(lib_id=int(args.lib_id), ion_mode=str(args.ion_mode))
    pool = scope_training_pool(
        pool,
        controls,
        lib_id=int(args.lib_id),
        ion_mode=str(args.ion_mode),
        control_scope=str(getattr(args, "training_control_scope", "x_controls")),
    )
    if wgcna_soft_prior:
        pool = add_wgcna_features(pool, wgcna_consensus)
    pool = add_labels(add_features(pool), controls)
    selection_mode = str(getattr(args, "model_selection_mode", "auto"))
    if not bool(getattr(args, "v8_regularized", True)) and selection_mode == "auto":
        selection_mode = "gbm"
    selection = select_model(pool, selection_mode)
    selected_name = str(selection["name"])
    pool["unified_oof_score"] = selection["scores"][f"{selected_name}_oof_score"].to_numpy(dtype=float)
    pool["unified_final_score"] = selection["scores"][f"{selected_name}_full_score"].to_numpy(dtype=float)
    pool = add_wgcna_features(pool, wgcna_consensus)
    model = selection.get("model")
    if model is not None and hasattr(model, "feature_importances_"):
        importances = pd.DataFrame({"feature": FEATURE_COLUMNS, "importance": model.feature_importances_}).sort_values(
            "importance", ascending=False
        )
    elif model is not None and hasattr(model, "coef_"):
        coefs = np.abs(np.asarray(model.coef_).reshape(-1))
        importances = pd.DataFrame({"feature": FEATURE_COLUMNS, "importance": coefs}).sort_values(
            "importance", ascending=False
        )
    else:
        importances = pd.DataFrame({"feature": FEATURE_COLUMNS, "importance": 0.0})
    g3 = g3_gate(int(selection["full_rank1"]), int(selection["oof_rank1"]))
    calibrators = fit_isotonic_per_source(pool, "unified_oof_score")
    with open(output_root / "isotonic_per_source.pkl", "wb") as handle:
        pickle.dump(calibrators, handle)
    drop_attr = oof_drop_attribution(pool, "unified_final_score", "unified_oof_score")
    write_csv(output_root / "oof_drop_attribution.csv", drop_attr)
    option_payloads = selection.get("all_options", [selection])
    selection_rows = []
    for payload in option_payloads:
        row = {k: v for k, v in payload.items() if k not in {"model", "scores", "all_options"}}
        row["selected"] = bool(str(row.get("name", "")) == selected_name)
        selection_rows.append(row)
    write_csv(output_root / "model_selection.csv", pd.DataFrame(selection_rows))
    write_csv(output_root / "unified_candidates_scored.csv", pool)
    write_csv(output_root / "rerank_predictions_s0.csv", pool)
    write_csv(output_root / "feature_importances.csv", importances)
    non_s0_predictions = pd.DataFrame()
    non_s0_path = resolve(args.non_s0_pool, repo_root) if getattr(args, "non_s0_pool", None) else None
    if non_s0_path is not None and non_s0_path.exists():
        try:
            non_s0_predictions = pd.read_csv(non_s0_path, low_memory=False)
        except pd.errors.EmptyDataError:
            non_s0_predictions = pd.DataFrame()
        if not non_s0_predictions.empty:
            if wgcna_soft_prior:
                non_s0_predictions = add_wgcna_features(non_s0_predictions, wgcna_consensus)
            non_s0_predictions = add_features(non_s0_predictions)
            if model is not None and pool["label_ge_0p7"].nunique() >= 2:
                non_s0_predictions["unified_final_score"] = _predict_proba(model, non_s0_predictions)
            else:
                non_s0_predictions["unified_final_score"] = heuristic_score(non_s0_predictions)
            non_s0_predictions = add_wgcna_features(non_s0_predictions, wgcna_consensus)
            non_s0_predictions["rerank_score_calibrated"] = apply_isotonic_per_source(
                non_s0_predictions, calibrators, "unified_final_score"
            )
            non_s0_predictions["rerank_score"] = non_s0_predictions["unified_final_score"]
            write_csv(output_root / "non_s0_candidates_scored.csv", non_s0_predictions)
            write_csv(output_root / "rerank_predictions_non_s0.csv", non_s0_predictions)
    if non_s0_predictions.empty:
        write_csv(output_root / "non_s0_candidates_scored.csv", pd.DataFrame())
        write_csv(output_root / "rerank_predictions_non_s0.csv", pd.DataFrame())
    oof_top10 = rank_frame(pool, "unified_oof_score", int(args.top_k), policy="unified_rerank_v1_oof")
    final_top10 = rank_frame(pool, "unified_final_score", int(args.top_k), policy="unified_rerank_v1")
    non_s0_top10 = (
        rank_frame(non_s0_predictions, "unified_final_score", int(args.top_k), policy="unified_rerank_v1_non_s0")
        if not non_s0_predictions.empty
        else pd.DataFrame(columns=final_top10.columns)
    )
    all_top10 = (
        pd.concat([final_top10, non_s0_top10], ignore_index=True, sort=False)
        if not non_s0_top10.empty
        else final_top10.copy()
    )
    write_csv(output_root / "unified_rerank_oof_top10.csv", oof_top10)
    write_csv(output_root / "unified_rerank_top10.csv", final_top10)
    write_csv(output_root / "unified_rerank_non_s0_top10.csv", non_s0_top10)
    write_csv(output_root / "unified_rerank_all_top10.csv", all_top10)
    candidate_path = output_root / "unified_candidates_scored.csv"
    denovo_path = output_root / "unified_rerank_top10.csv"
    write_denovo_report(final_top10, denovo_path)
    metrics, per_control = eval_reranked(
        candidate_df=pool,
        ranked_df=final_top10,
        controls=controls,
        lib_id=int(args.lib_id),
        ion_mode=str(args.ion_mode),
        profile="strict_1da",
        policy="unified_rerank_v1",
        top_k=int(args.top_k),
        candidate_path=candidate_path,
        denovo_path=denovo_path,
    )
    oof_metrics, oof_per_control = eval_reranked(
        candidate_df=pool,
        ranked_df=oof_top10,
        controls=controls,
        lib_id=int(args.lib_id),
        ion_mode=str(args.ion_mode),
        profile="strict_1da_oof",
        policy="unified_rerank_v1_oof",
        top_k=int(args.top_k),
        candidate_path=candidate_path,
        denovo_path=output_root / "unified_rerank_oof_top10.csv",
    )
    metrics_df = pd.DataFrame(metrics + oof_metrics)
    per_df = pd.DataFrame(per_control)
    oof_per_df = pd.DataFrame(oof_per_control)
    write_csv(output_root / "unified_rerank_metrics.csv", metrics_df)
    write_csv(output_root / "unified_rerank_per_control.csv", per_df)
    write_csv(output_root / "unified_rerank_oof_per_control.csv", oof_per_df)
    top10_hits = (
        int(pd.to_numeric(per_df.get("top10_best_tanimoto"), errors="coerce").fillna(0.0).ge(0.7).sum())
        if not per_df.empty
        else 0
    )
    rank1_hits = (
        int(pd.to_numeric(per_df.get("rank1_best_tanimoto"), errors="coerce").fillna(0.0).ge(0.7).sum())
        if not per_df.empty
        else 0
    )
    oof_top10_hits = (
        int(pd.to_numeric(oof_per_df.get("top10_best_tanimoto"), errors="coerce").fillna(0.0).ge(0.7).sum())
        if not oof_per_df.empty
        else 0
    )
    oof_rank1_hits = (
        int(pd.to_numeric(oof_per_df.get("rank1_best_tanimoto"), errors="coerce").fillna(0.0).ge(0.7).sum())
        if not oof_per_df.empty
        else 0
    )
    summary = {
        "n_candidates": int(pool.shape[0]),
        "n_positive_labels": int(pool["label_ge_0p7"].sum()),
        "selected_model": selected_name,
        "model_selection_mode": selection_mode,
        "training_control_scope": str(getattr(args, "training_control_scope", "x_controls")),
        "wgcna_soft_prior": bool(wgcna_soft_prior),
        "class_match_score_feature": bool(getattr(args, "class_match_score_feature", False)),
        "drop_novelty_score": bool(getattr(args, "drop_novelty_score", False)),
        "final_top10_ge_0p7": top10_hits,
        "final_rank1_ge_0p7": rank1_hits,
        "oof_top10_ge_0p7": oof_top10_hits,
        "oof_rank1_ge_0p7": oof_rank1_hits,
        "g3_pass": bool(g3["g3_pass"]),
        "g3_rank1_delta_abs": int(g3["rank1_delta_abs"]),
        "g3_full_fit_rank1": int(g3["full_fit_rank1"]),
        "g3_oof_rank1": int(g3["oof_rank1"]),
        "min_oof_rank1": int(g3["min_oof_rank1"]),
        "g3_floor_policy": str(g3["g3_floor_policy"]),
        "n_non_s0_predictions": int(non_s0_predictions.shape[0]) if not non_s0_predictions.empty else 0,
    }
    write_json(output_root / "unified_rerank_summary.json", summary)
    write_json(output_root / "rerank_generalization_summary.json", summary)
    if not bool(g3["g3_pass"]):
        (output_root / "g3_failure_report.md").write_text(
            "\n".join(
                [
                    "# G3 Failure Report",
                    "",
                    f"Selected model: `{selected_name}`",
                    f"Full-fit rank1: `{g3['full_fit_rank1']}`",
                    f"OOF rank1: `{g3['oof_rank1']}`",
                    f"Absolute rank1 delta: `{g3['rank1_delta_abs']}`",
                    f"Minimum OOF rank1: `{g3['min_oof_rank1']}`",
                    "",
                    "## OOF Drop Attribution",
                    dataframe_to_markdown(drop_attr),
                    "",
                ]
            )
        )
    (output_root / "selection.md").write_text(
        "\n".join(
            [
                "# Rerank V8 Regularized Selection",
                "",
                f"Selected model: `{selected_name}`",
                "",
                "## Selected Metrics",
                dataframe_to_markdown(pd.DataFrame(selection_rows)),
                "",
                "## G3",
                dataframe_to_markdown(pd.DataFrame([g3])),
                "",
            ]
        )
    )
    report = [
        "# Unified Rerank V1",
        "",
        "## Summary",
        dataframe_to_markdown(pd.DataFrame([summary])),
        "",
        "## Model Selection",
        dataframe_to_markdown(pd.DataFrame(selection_rows)),
        "",
        "## G3",
        dataframe_to_markdown(pd.DataFrame([g3])),
        "",
        "## Feature Importances",
        dataframe_to_markdown(importances.head(20)),
        "",
        "## Metrics",
        dataframe_to_markdown(metrics_df),
        "",
    ]
    (output_root / "unified_rerank_report.md").write_text("\n".join(report))
    print(json.dumps(summary, indent=2, sort_keys=True))
    return summary


def main() -> None:
    """Parse CLI arguments and invoke :func:`run`."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", type=Path, default=Path("models/unified_rerank_v1"))
    parser.add_argument(
        "--unified-pool", type=Path, default=Path("models/unified_candidate_pool_v1/unified_candidates_per_control.csv")
    )
    parser.add_argument("--lib-id", type=int, default=209)
    parser.add_argument("--ion-mode", default="negative")
    parser.add_argument("--top-k", type=int, default=10)
    parser.add_argument(
        "--non-s0-pool",
        type=Path,
        default=Path("models/unified_candidate_pool_v1/unified_candidates_per_control_non_s0.csv"),
    )
    parser.add_argument("--v8-regularized", action="store_true", default=True)
    parser.add_argument("--model-selection-mode", choices=["auto", "gbm", "logistic"], default="auto")
    parser.add_argument("--training-control-scope", choices=["x_controls", "all"], default="x_controls")
    parser.add_argument("--v6-baseline-rank1", type=int, default=44)
    parser.add_argument("--wgcna-soft-prior", action="store_true")
    parser.add_argument("--class-match-score-feature", action="store_true")
    parser.add_argument("--drop-novelty-score", action="store_true")
    run(parser.parse_args())


if __name__ == "__main__":
    main()
