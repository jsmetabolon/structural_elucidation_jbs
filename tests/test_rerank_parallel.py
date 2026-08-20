"""Parity tests for rerank labeling memoization and parallel LOOCV.

Both optimizations must be output-identical: memoized labeling reproduces the
row-wise loop exactly, and LOOCV scores are bit-equal at any worker count.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pandas as pd

from structural_elucidation.pipeline import rerank


def test_add_labels_matches_rowwise_reference():
    controls = [
        SimpleNamespace(control_id="X-1", smiles=["CCO", "CCN"]),
        SimpleNamespace(control_id="X-2", smiles=["c1ccccc1O"]),
    ]
    pool = pd.DataFrame(
        {
            "control_id": ["X-1", "X-1", "X-1", "X-2", "X-2", "X-3"],
            "predicted_smiles": ["CCO", "OCC", "CCCC", "c1ccccc1O", "CCO", None],
            "canonical_2d_smiles": ["", "", "", "", "", "CC"],
        }
    )

    out = rerank.add_labels(pool.copy(), controls)

    truth = {str(control.control_id): list(control.smiles) for control in controls}
    reference_tans: list[float] = []
    reference_labels: list[int] = []
    for _, row in pool.iterrows():
        canonical = rerank.canonicalize_smiles(row.get("predicted_smiles", "")) or str(
            row.get("canonical_2d_smiles", "") or ""
        )
        _exact, tan = rerank.best_tanimoto_for_smiles(canonical, truth.get(str(row.get("control_id", "")), []))
        reference_tans.append(float(tan))
        reference_labels.append(int(float(tan) >= 0.7))

    assert out["truth_tanimoto"].tolist() == reference_tans
    assert out["label_ge_0p7"].tolist() == reference_labels
    # OCC canonicalizes to the same structure as CCO, so the exact truth match labels it 1.
    assert out.loc[1, "label_ge_0p7"] == 1


def test_add_labels_empty_pool():
    out = rerank.add_labels(pd.DataFrame(columns=["control_id", "predicted_smiles"]), [])
    assert list(out["truth_tanimoto"]) == []
    assert list(out["label_ge_0p7"]) == []


def _synthetic_pool(n_controls: int = 6, rows_per_control: int = 30, seed: int = 7) -> pd.DataFrame:
    rng = np.random.RandomState(seed)
    frames = []
    for i in range(n_controls):
        frame = pd.DataFrame({col: rng.rand(rows_per_control) for col in rerank.FEATURE_COLUMNS})
        frame["control_id"] = f"X-{i}"
        if i == 0:
            labels = np.zeros(rows_per_control, dtype=int)  # heuristic-fallback control
        else:
            labels = (rng.rand(rows_per_control) > 0.8).astype(int)
        frame["label_ge_0p7"] = labels
        frames.append(frame)
    return pd.concat(frames, ignore_index=True)


def test_loocv_parallel_matches_serial():
    pool = _synthetic_pool()
    serial = rerank.loocv_scores(pool, n_jobs=1)
    parallel = rerank.loocv_scores(pool, n_jobs=2)
    pd.testing.assert_series_equal(serial, parallel)


def test_select_model_parallel_matches_serial():
    pool = _synthetic_pool(n_controls=5, rows_per_control=25, seed=11)
    serial = rerank.select_model(pool, "auto", n_jobs=1)
    parallel = rerank.select_model(pool, "auto", n_jobs=2)
    assert serial["name"] == parallel["name"]
    assert serial["full_rank1"] == parallel["full_rank1"]
    assert serial["oof_rank1"] == parallel["oof_rank1"]
    pd.testing.assert_frame_equal(serial["scores"], parallel["scores"])
