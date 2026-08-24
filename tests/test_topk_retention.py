"""Unit tests for the top-10 retention screen audit (tools/audit_topk_retention.py)."""
from __future__ import annotations

import importlib.util
import json
import math
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent


def load_tool(name: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "tools" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    return module


TK = load_tool("audit_topk_retention")

from structural_elucidation.reports import confidence as C  # noqa: E402

# A calibration whose arithmetic is checkable in your head: predicted = 2 * proxy + 100,
# and the residual is standardized by 100 index units.
CALIBRATION = {"209": {"slope": 2.0, "intercept": 100.0, "residual_rmse": 100.0}}

PROXIES = {"AAA": 400.0, "BBB": 650.0, "CCC": math.nan, "DDD": 300.0}


def fake_proxy(smiles: str) -> float:
    return PROXIES.get(smiles, math.nan)


def _rows() -> list[dict]:
    return [
        # z = (1500 - (2*400+100)) / 100 = 6.0 -> flagged; labeled incorrect.
        {"control_id": "X-1", "rank": "1", "rerank_policy": "unified_rerank_v1",
         "feature_key": "f_mz100.0_ri1500_X-1", "predicted_smiles": "AAA", "label_ge_0p7": "0"},
        # z = (1500 - 1400) / 100 = 1.0 -> kept; labeled correct.
        {"control_id": "X-1", "rank": "2", "rerank_policy": "unified_rerank_v1",
         "feature_key": "f_mz100.0_ri1500_X-1", "predicted_smiles": "BBB", "label_ge_0p7": "1"},
        # RI recorded as zero means "not measured": no axis.
        {"control_id": "X-2", "rank": "1", "rerank_policy": "unified_rerank_v1",
         "feature_key": "f_mz100.0_ri0_X-2", "predicted_smiles": "AAA", "label_ge_0p7": ""},
        # rank_frame's fallback key carries no index at all: no axis.
        {"control_id": "X-2", "rank": "2", "rerank_policy": "unified_rerank_v1",
         "feature_key": "unified_rerank_X-2", "predicted_smiles": "AAA", "label_ge_0p7": ""},
        # No SMILES: nothing to predict from.
        {"control_id": "X-3", "rank": "1", "rerank_policy": "unified_rerank_v1",
         "feature_key": "f_mz100.0_ri2000_X-3", "predicted_smiles": "", "label_ge_0p7": ""},
        # A proxy the descriptor model cannot compute: no axis.
        {"control_id": "X-3", "rank": "2", "rerank_policy": "unified_rerank_v1",
         "feature_key": "f_mz100.0_ri2000_X-3", "predicted_smiles": "CCC", "label_ge_0p7": ""},
        # Exactly at the bar: |z| = 3.0 flags.
        {"control_id": "X-4", "rank": "1", "rerank_policy": "unified_rerank_v1_non_s0",
         "feature_key": "f_mz100.0_ri1000_X-4", "predicted_smiles": "DDD", "label_ge_0p7": ""},
        # A row from a library that did not clear the feasibility audit: no axis.
        {"control_id": "X-5", "rank": "1", "rerank_policy": "unified_rerank_v1", "lib_id": "305",
         "feature_key": "f_mz100.0_ri1200_X-5", "predicted_smiles": "AAA", "label_ge_0p7": ""},
    ]


def test_score_rows_computes_the_confidence_models_z_exactly():
    scored = TK.score_rows(_rows(), "209", CALIBRATION, fake_proxy)
    by_key = {(row["control_id"], row["rank"]): row for row in scored}
    assert by_key[("X-1", 1)]["z"] == pytest.approx(6.0)
    assert by_key[("X-1", 2)]["z"] == pytest.approx(1.0)
    # Parity with the axis the report applies, not a reimplementation of it.
    direct = C.retention_residual_z(
        {"lib_id": "209", "observed_retention_index": 1500.0, "candidate_retention_proxy": 400.0},
        CALIBRATION,
    )
    assert by_key[("X-1", 1)]["z"] == pytest.approx(direct)


def test_rows_without_the_evidence_get_no_z_rather_than_a_zero():
    scored = TK.score_rows(_rows(), "209", CALIBRATION, fake_proxy)
    by_key = {(row["control_id"], row["rank"]): row for row in scored}
    for key in [("X-2", 1), ("X-2", 2), ("X-3", 1), ("X-3", 2), ("X-5", 1)]:
        assert by_key[key]["z"] is None
        assert by_key[key]["flagged"] is False


def test_the_flag_fires_at_the_bar_and_not_below_it():
    scored = TK.score_rows(_rows(), "209", CALIBRATION, fake_proxy)
    by_key = {(row["control_id"], row["rank"]): row for row in scored}
    assert by_key[("X-4", 1)]["z"] == pytest.approx(3.0)
    assert by_key[("X-4", 1)]["flagged"] is True
    assert by_key[("X-1", 2)]["flagged"] is False


def test_summary_counts_coverage_flags_and_false_kills():
    summary = TK.summarize_library(TK.score_rows(_rows(), "209", CALIBRATION, fake_proxy))
    assert summary["n_rows"] == 8
    assert summary["n_controls"] == 5
    assert summary["coverage"]["rows_scored"] == 3
    assert summary["coverage"]["rows_with_observed_ri"] == 6
    assert summary["coverage"]["controls_scored"] == 2
    flags = summary["flags"]
    assert flags["rows_flagged"] == 2
    assert flags["controls_with_flag"] == 2
    assert flags["share_controls_with_flag"] == pytest.approx(1.0)
    assert flags["rank1_flagged_controls"] == 2
    assert flags["flagged_rows_by_rank"] == {"1": 2}
    kills = summary["false_kills"]
    assert kills["labeled_rows_scored"] == 2
    assert kills["labeled_correct_rows_scored"] == 1
    assert kills["labeled_correct_rows_flagged"] == 0
    assert kills["max_abs_z_labeled_correct"] == pytest.approx(1.0)
    assert summary["abs_z_distribution"]["max"] == pytest.approx(6.0)
    # Worst offenders first in the examples.
    assert [row["z"] for row in summary["flagged_examples"]] == [6.0, 3.0]


def _evaluated(share: float | None, correct_flagged: int) -> dict:
    return {
        "status": "evaluated",
        "flags": {"share_controls_with_flag": share},
        "false_kills": {"labeled_correct_rows_flagged": correct_flagged},
    }


def test_gate_goes_only_when_it_fires_enough_and_kills_nothing_correct():
    go = TK.gate_verdict({"209": _evaluated(0.08, 0), "400": _evaluated(0.01, 0)})
    assert go["verdict"] == "GO"
    too_quiet = TK.gate_verdict({"209": _evaluated(0.01, 0)})
    assert too_quiet["verdict"] == "NO_GO"
    assert "under the" in too_quiet["reasoning"]
    false_kill = TK.gate_verdict({"209": _evaluated(0.50, 1)})
    assert false_kill["verdict"] == "NO_GO"
    assert "known-correct" in false_kill["reasoning"]
    nothing = TK.gate_verdict({"209": {"status": "missing_top10"}})
    assert nothing["verdict"] == "NOT_EVALUATED"


def test_parse_label_reads_the_shapes_the_csv_produces():
    assert TK.parse_label("1") is True
    assert TK.parse_label("1.0") is True
    assert TK.parse_label("0") is False
    assert TK.parse_label("0.0") is False
    assert TK.parse_label("True") is True
    assert TK.parse_label("False") is False
    assert TK.parse_label("") is None
    assert TK.parse_label("nan") is None
    assert TK.parse_label(None) is None


def test_percentile_is_nearest_rank():
    assert TK.percentile([1.0, 3.0, 6.0], 0.50) == 3.0
    assert TK.percentile([1.0, 3.0, 6.0], 0.90) == 6.0
    assert TK.percentile([5.0], 0.99) == 5.0
    assert TK.percentile([], 0.50) is None


def test_run_audit_reports_each_library_status_and_serializes(tmp_path):
    import csv as csv_module

    lib_dir = tmp_path / "rerank_lib209"
    lib_dir.mkdir()
    with (lib_dir / "unified_rerank_all_top10.csv").open("w", newline="") as handle:
        writer = csv_module.DictWriter(
            handle,
            fieldnames=["control_id", "rank", "rerank_policy", "feature_key", "predicted_smiles", "label_ge_0p7"],
        )
        writer.writeheader()
        for row in _rows():
            writer.writerow({key: row.get(key, "") for key in writer.fieldnames})
    calibration = {**CALIBRATION, "400": {"slope": 2.0, "intercept": 0.0, "residual_rmse": 100.0}}
    summary = TK.run_audit(
        rerank_root=tmp_path,
        libs=["209", "400", "305"],
        top10_name="unified_rerank_all_top10.csv",
        calibration=calibration,
        proxy_fn=fake_proxy,
    )
    assert summary["libraries"]["209"]["status"] == "evaluated"
    assert summary["libraries"]["400"]["status"] == "missing_top10"
    assert summary["libraries"]["305"]["status"] == "no_calibrated_axis"
    assert summary["gate"]["verdict"] == "GO"
    json.dumps(summary)  # every value must be serializable as committed audit data
    markdown = TK.render_markdown(summary)
    assert "GO" in markdown
    assert "| 209 |" in markdown


def test_the_tool_reads_the_committed_calibration_the_report_uses():
    """The audit and the report must screen with the same committed constants."""
    calibration = C.load_retention_calibration(ROOT / "docs/paper/data/retention/ri_calibration.json")
    assert set(calibration) == {"209", "400"}
    z = C.retention_residual_z(
        {"lib_id": "209", "observed_retention_index": 1000.0, "candidate_retention_proxy": 300.0},
        calibration,
    )
    # (1000 - (5.1967 * 300 - 590.7)) / 599.9, straight from the feasibility audit's fit.
    assert z == pytest.approx((1000.0 - (5.1967 * 300.0 - 590.7)) / 599.9)


def test_committed_topk_audit_is_consistent_with_its_own_gate():
    """Guard the committed audit artifact once pipeline outputs have produced it."""
    audit_path = ROOT / "docs/paper/data/retention/audit_topk_retention.json"
    if not audit_path.exists():
        pytest.skip("top-10 retention audit not yet generated (needs reranked pipeline outputs)")
    committed = json.loads(audit_path.read_text())
    if committed.get("status") == "not_evaluated":
        pytest.skip("committed audit records a not-evaluated run")
    assert committed["gate"] == TK.gate_verdict(committed["libraries"])
    assert committed["screen"]["implausible_abs_z"] == C.RETENTION_IMPLAUSIBLE_Z
