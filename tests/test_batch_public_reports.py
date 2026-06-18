from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from structural_elucidation.reports import batch_public


def write_predictions(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, index=False)


def fake_renderer(args: argparse.Namespace) -> dict[str, object]:
    args.out_dir.mkdir(parents=True, exist_ok=True)
    (args.out_dir / "index.html").write_text("index", encoding="utf-8")
    (args.out_dir / "public.html").write_text("public", encoding="utf-8")
    (args.out_dir / "guides.html").write_text("guides", encoding="utf-8")
    (args.out_dir / "manifest.json").write_text(json.dumps({"control_id": args.control_id}), encoding="utf-8")
    return {"html": str(args.out_dir / "index.html")}


def test_discover_and_index_prediction_csvs(tmp_path: Path) -> None:
    csv_path = tmp_path / "eq" / "rerank_lib209" / "predictions.csv"
    write_predictions(
        csv_path,
        [
            {"control_id": "X-1", "lib_id": 209, "rank": 1, "predicted_smiles": "CCO"},
            {"control_id": "KNOWN-1", "lib_id": 209, "rank": 1, "predicted_smiles": "CCC"},
        ],
    )

    discovered = batch_public.discover_prediction_csvs(tmp_path / "eq", [209], "predictions.csv")
    index = batch_public.load_prediction_index(discovered)
    jobs = batch_public.iter_report_jobs(index, out_root=tmp_path / "out")

    assert [item.path for item in discovered] == [csv_path]
    assert len(index) == 1
    assert len(jobs) == 1
    assert jobs[0].slug == "lib209_X-1"


def test_batch_run_writes_outputs(tmp_path: Path) -> None:
    csv_path = tmp_path / "eq" / "rerank_lib209" / "predictions.csv"
    write_predictions(
        csv_path,
        [
            {
                "control_id": "X-1",
                "lib_id": 209,
                "ion_mode": "negative",
                "rank": 1,
                "candidate_name": "ethanol",
                "predicted_smiles": "CCO",
                "predicted_formula": "C2H6O",
                "predicted_pathway": "Energy",
            }
        ],
    )
    out_root = tmp_path / "public"
    args = batch_public.parse_args(
        [
            "--equalization-root",
            str(tmp_path / "eq"),
            "--libs",
            "209",
            "--input-name",
            "predictions.csv",
            "--out-root",
            str(out_root),
            "--force",
        ]
    )

    summary = batch_public.run(args, renderer=fake_renderer)

    assert summary["n_generated"] == 1
    assert (out_root / "index.html").exists()
    assert (out_root / "manifest.csv").exists()
    assert (out_root / "manifest.json").exists()
    assert (out_root / "failures.csv").exists()
    assert (out_root / "compounds" / "lib209_X-1" / "index.html").exists()
    manifest = pd.read_csv(out_root / "manifest.csv")
    assert manifest.loc[0, "rank1_name"] == "ethanol"
    assert manifest.loc[0, "index_html"] == "compounds/lib209_X-1/index.html"


def test_failed_report_writes_failures(tmp_path: Path) -> None:
    csv_path = tmp_path / "predictions.csv"
    write_predictions(csv_path, [{"control_id": "X-2", "rank": 1, "predicted_smiles": "CCO"}])
    out_root = tmp_path / "public"
    args = batch_public.parse_args(
        [
            "--predictions",
            str(csv_path),
            "--out-root",
            str(out_root),
        ]
    )

    def failing_renderer(_: argparse.Namespace) -> dict[str, object]:
        raise RuntimeError("boom")

    summary = batch_public.run(args, renderer=failing_renderer)

    assert summary["n_failures"] == 1
    failures = pd.read_csv(out_root / "failures.csv")
    assert failures.loc[0, "control_id"] == "X-2"
    assert "boom" in failures.loc[0, "reason"]
