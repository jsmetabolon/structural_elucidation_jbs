#!/usr/bin/env python3
"""Shared helpers for protected 209 candidate fusion."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import pandas as pd

from denovo_improvement_utils import canonicalize_smiles
from x_domain_analog_utils import safe_float


PROTECTED_209_CONTROLS = ["X-17919", "X-18345", "X-21410", "X-21821", "X-24137", "X-24972", "X-25459"]
OPTIONAL_AUDIT_CONTROLS = ["X-25003"]
GOOD_TANIMOTO = 0.7


def parse_xid(value: Any) -> str:
    match = re.search(r"(X-\d+)", str(value or ""))
    return match.group(1) if match else ""


def canonical_key(value: Any) -> str:
    canonical = canonicalize_smiles(value)
    return canonical or str(value or "").strip()


def num(frame: pd.DataFrame, col: str, default: float = 0.0) -> pd.Series:
    if col not in frame.columns:
        return pd.Series(default, index=frame.index, dtype=float)
    return pd.to_numeric(frame[col], errors="coerce").fillna(default).astype(float)


def source_priority(label: str) -> float:
    text = str(label or "").lower()
    if "seed_source_acquisition" in text:
        return 1.00
    if "seed_source_local_edits" in text:
        return 0.92
    if "assembly_v4" in text:
        return 0.88
    if "formula_rescue" in text:
        return 0.86
    if "assembly_v3" in text:
        return 0.78
    return 0.60


def infer_source_label(path: Path) -> str:
    text = str(path)
    if "seed_source_acquisition_v1" in text:
        return "seed_source_acquisition_strict"
    if "209_seed_source_discovery_v1" in text:
        return "seed_source_local_edits"
    if "denovo_pipeline_v4" in text:
        return "assembly_v4"
    if "denovo_pipeline_v3" in text:
        return "assembly_v3"
    if "formula_rescue" in text:
        return "formula_rescue"
    if "ranking_scoring_improvement_v1" in text:
        return "ranking_scoring"
    return path.parent.name


def infer_profile(path: Path) -> str:
    parts = list(path.parts)
    if "selfies_profiles" in parts:
        idx = parts.index("selfies_profiles")
        if idx + 1 < len(parts):
            return parts[idx + 1]
    for part in parts:
        if part in {"strict_1da", "balanced"}:
            return part
    return "unknown"


def load_candidate_file(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    frame = pd.read_csv(path, low_memory=False)
    if frame.empty or "feature_key" not in frame.columns or "predicted_smiles" not in frame.columns:
        return pd.DataFrame()
    out = frame.copy()
    out["fusion_source_label"] = infer_source_label(path)
    out["fusion_source_profile"] = infer_profile(path)
    out["fusion_source_path"] = str(path)
    out["fusion_source_priority"] = source_priority(out["fusion_source_label"].iloc[0])
    out["control_id"] = out["feature_key"].map(parse_xid)
    out["_canonical_key"] = out["predicted_smiles"].map(canonical_key)
    out = out[out["_canonical_key"].astype(str).ne("")].copy()
    return out


def load_accepted_replay_scores(paths: list[Path]) -> pd.DataFrame:
    rows = []
    for path in paths:
        if not path.exists():
            continue
        frame = pd.read_csv(path, low_memory=False)
        if frame.empty or "feature_key" not in frame.columns or "predicted_smiles" not in frame.columns:
            continue
        rank = num(frame, "rank", 999999.0).clip(lower=1.0)
        tmp = pd.DataFrame(
            {
                "feature_key": frame["feature_key"].astype(str),
                "_canonical_key": frame["predicted_smiles"].map(canonical_key),
                "accepted_replay_score": 1.0,
                "accepted_replay_rank_score": (1.0 / rank).clip(0.0, 1.0),
                "accepted_replay_source": str(path),
            }
        )
        rows.append(tmp[tmp["_canonical_key"].astype(str).ne("")])
    if not rows:
        return pd.DataFrame(columns=["feature_key", "_canonical_key", "accepted_replay_score", "accepted_replay_source"])
    out = pd.concat(rows, ignore_index=True, sort=False)
    out = out.sort_values(["feature_key", "_canonical_key", "accepted_replay_score"], ascending=[True, True, False], kind="mergesort")
    return out.drop_duplicates(["feature_key", "_canonical_key"], keep="first")


def protected_summary(per_control: pd.DataFrame, policy: str, protected_controls: list[str] | None = None) -> dict[str, Any]:
    protected = set(protected_controls or PROTECTED_209_CONTROLS)
    scoped = per_control[per_control["policy"].astype(str).eq(str(policy))].copy() if not per_control.empty else pd.DataFrame()
    recovered = set(scoped[pd.to_numeric(scoped.get("top10_best_tanimoto", pd.Series(dtype=float)), errors="coerce").fillna(0.0).ge(GOOD_TANIMOTO)]["control_id"].astype(str))
    missing = sorted(protected - recovered)
    return {
        "policy": str(policy),
        "n_protected_controls": len(protected),
        "n_protected_preserved": len(protected) - len(missing),
        "protected_controls_preserved": len(missing) == 0,
        "missing_protected_controls": "|".join(missing),
        "protected_recovered_controls": "|".join(sorted(protected & recovered)),
    }


def choose_best_policy(decisions: pd.DataFrame) -> str:
    if decisions.empty:
        return "protected_replay_x_success"
    work = decisions.copy()
    work["_preserve"] = work.get("protected_controls_preserved", pd.Series(False, index=work.index)).astype(bool)
    work["_ge"] = pd.to_numeric(work.get("n_top10_ge_0p7", 0), errors="coerce").fillna(0)
    work["_mean"] = pd.to_numeric(work.get("mean_top10_tanimoto", 0.0), errors="coerce").fillna(0.0)
    work = work.sort_values(["_preserve", "_ge", "_mean", "policy"], ascending=[False, False, False, True], kind="mergesort")
    return str(work.iloc[0]["policy"]) if not work.empty else "protected_replay_x_success"
