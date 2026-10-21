"""Offline driver contract: scored results survive optional rendering failures."""
from argparse import Namespace
import hashlib
import json
import logging
from pathlib import Path

import pandas as pd
import pytest

from structural_elucidation.adapters.untargeted_inputs import build_inputs_from_handoff
from structural_elucidation.pipeline import run_novel as driver
from structural_elucidation.reports import novel_scoring

CONTROL_ID = "X-900000001"
SMILES = "NC(Cc1ccc(O)cc1)C(=O)O"


def feature(feature_id=1):
    return {
        "feature_id": feature_id, "lib_id": 400, "ion_mode": "positive",
        "parent_mz": 182.081166, "neutral_mass": 181.07389, "adduct_type": "[M+H]+",
        "ri": 1500, "peaks": [[136.0757, 1.0], [119.0491, 0.5], [91.0542, 0.3]],
        "prior_match_level": "MS2", "prior_matched_name": "putative tyrosine",
        "compound_name": "putative tyrosine",
    }


def prediction(control_id=CONTROL_ID, **overrides):
    return {
        "control_id": control_id, "rank": 1, "predicted_smiles": SMILES,
        "rerank_score": 0.82, "pure_nce_score": 0.72, "nce_status": "scored",
        "forward_spectrum_score": 0.72, "predicted_spectrum_score_v2": 0.72,
        "candidate_name": "tyrosine", "mass_match_score": 1.0,
        "abs_mass_error_da": 0.0, "predicted_formula": "C9H11NO3", **overrides,
    }


def args_for(tmp_path, *, features=None, stages="rerank"):
    handoff = tmp_path / "handoff.json"
    handoff.write_text(json.dumps({"schema_version": 1, "source": "untargeted_workflow",
                                   "features": [feature(), feature(2)] if features is None else features}))
    return Namespace(handoff=handoff, out_root=tmp_path / "run", inputs_dir=None,
                     assets_root=tmp_path / "assets", stages=stages, libs="",
                     enable_denovo=False, denovo_extra="", top_k=10, jobs=1, dry_run=False)


def stub_stages(monkeypatch, rows, *, reports="success"):
    commands = []

    def run(self, module, args, **kwargs):
        commands.append((module, args))
        if module.endswith(".rerank"):
            root = Path(args[args.index("--output-root") + 1])
            root.mkdir(parents=True, exist_ok=True)
            pd.DataFrame(rows, columns=list(prediction())).to_csv(root / "unified_rerank_all_top10.csv", index=False)
        elif module.endswith("reports.batch"):
            if reports == "failure":
                raise RuntimeError("fixture report failed")
            root = Path(args[args.index("--out-dir") + 1])
            score_root = Path(args[args.index("--scoring-dir") + 1])
            root.mkdir(parents=True, exist_ok=True)
            for control in args[args.index("--control-ids") + 1].split(","):
                payload = json.loads((score_root / f"{control}.json").read_text())
                (root / control).mkdir(parents=True, exist_ok=True)
                (root / control / "index.html").write_text("<pre>" + json.dumps(payload["summary"]) + "</pre>")
            (root / "batch_summary.json").write_text(json.dumps({"failed": []}))
        return 0

    monkeypatch.setattr(driver.StageRunner, "run", run)
    return commands


def stub_scorer(monkeypatch, *, level="2*"):
    def score(candidates, metadata, **kwargs):
        summary = dict.fromkeys(driver.SCORE_COLUMNS)
        summary.update(se_msi_level=2, se_msi_label="Putative", se_msi_rationale="fixture",
                       se_confidence=0.8, se_evidence_score=0.0, se_confidence_tier="High",
                       se_reported_level=level, se_reported_identification="fixture identification",
                       se_candidate_identity_confidence=78.0, se_confidence_flags_json="[]")
        return {"schema_version": 1, "control_id": metadata["control_id"], "status": "scored",
                "summary": summary, "fields": {"fixture": True}, "errors": []}
    monkeypatch.setattr(novel_scoring, "score_novel_feature", score)
    return score


@pytest.mark.parametrize("report_mode", ["omitted", "success", "failure"])
def test_real_scoring_precedes_optional_reports(tmp_path, monkeypatch, report_mode):
    args = args_for(tmp_path, stages="rerank" if report_mode == "omitted" else "rerank,reports")
    commands = stub_stages(monkeypatch, [prediction()], reports=report_mode)
    # Supplied NCE scores must not cause model access during this fixture.
    monkeypatch.setattr("tools.rescore_pure_nce.score_key", lambda *_: pytest.fail("unexpected NCE model access"))
    metadata = driver.run(args)
    summary = pd.read_csv(metadata["summary"], dtype={"se_reported_level": "string"})
    assert summary["feature_id"].tolist() == [1]
    assert summary.loc[0, "se_rank1_smiles"] == SMILES
    assert metadata["scoring_errors"] == []
    canonical = json.loads((args.out_root / "scores" / f"{CONTROL_ID}.json").read_text())
    assert canonical["status"] == "scored"
    assert canonical["fields"].get("metabolon_level", "") == ""
    assert canonical["fields"].get("truth_match") is not True
    assert 0 <= summary.loc[0, "se_confidence"] <= 1
    assert 0 <= summary.loc[0, "se_evidence_score"] <= 1
    assert 0 <= summary.loc[0, "se_candidate_identity_confidence"] <= 100
    for column, value in canonical["summary"].items():
        if value is None:
            assert pd.isna(summary.loc[0, column])
        elif isinstance(value, float):
            assert summary.loc[0, column] == pytest.approx(value)
        else:
            assert summary.loc[0, column] == value
    if report_mode == "success":
        assert Path(summary.loc[0, "se_report_path"]).is_file()
        report_text = Path(summary.loc[0, "se_report_path"]).read_text()
        assert json.dumps(canonical["summary"]) in report_text
    else:
        assert pd.isna(summary.loc[0, "se_report_path"])
    assert bool(metadata["report_errors"]) == (report_mode == "failure")
    assert any(module.endswith("reports.batch") for module, _ in commands) == (report_mode != "omitted")


def test_summary_contains_only_current_handoff_predictions_and_preserves_level(tmp_path, monkeypatch):
    args = args_for(tmp_path)
    stub_stages(monkeypatch, [prediction(), prediction("X-900000002", predicted_smiles=""), prediction("X-900000099")])
    stub_scorer(monkeypatch)
    metadata = driver.run(args)
    summary = pd.read_csv(metadata["summary"], dtype={"se_reported_level": "string"})
    assert summary["feature_id"].tolist() == [1]
    assert summary.loc[0, "se_reported_level"] == "2*"
    assert summary.loc[0, "se_evidence_score"] == 0.0
    assert not (args.out_root / "scores" / "X-900000099.json").exists()
    assert metadata["n_feature_summaries"] == 1


@pytest.mark.parametrize("failure", ["no_output", "exception"])
def test_failed_rerank_never_reuses_stale_predictions_or_scores(tmp_path, monkeypatch, failure):
    args = args_for(tmp_path)
    root = args.out_root / "rerank_lib400"
    root.mkdir(parents=True)
    pd.DataFrame([prediction()]).to_csv(root / "unified_rerank_all_top10.csv", index=False)
    score_root = args.out_root / "scores"
    score_root.mkdir()
    (score_root / f"{CONTROL_ID}.json").write_text(json.dumps({"summary": {"se_confidence": 1.0}}))
    def stage(*_args, **_kwargs):
        if failure == "exception":
            raise RuntimeError("fixture failed")
        return 0
    monkeypatch.setattr(driver.StageRunner, "run", stage)
    monkeypatch.setattr(novel_scoring, "score_novel_feature", lambda *a, **k: pytest.fail("scoring stale predictions"))
    metadata = driver.run(args)
    assert pd.read_csv(metadata["summary"]).empty
    assert metadata["scopes"][0]["error"]
    assert metadata["scopes"][0]["rerank_top10"] is None


def test_scoring_error_keeps_prediction_but_no_invented_confidence(tmp_path, monkeypatch):
    args = args_for(tmp_path)
    stub_stages(monkeypatch, [prediction()])
    def fail(*args, **kwargs):
        raise ValueError("fixture evidence invalid")
    monkeypatch.setattr(novel_scoring, "score_novel_feature", fail)
    metadata = driver.run(args)
    summary = pd.read_csv(metadata["summary"])
    assert summary.loc[0, "se_rank1_smiles"] == SMILES
    assert summary[["se_confidence", "se_evidence_score", "se_candidate_identity_confidence"]].isna().all().all()
    assert "fixture evidence invalid" in metadata["scoring_errors"][0]["error"]


@pytest.mark.parametrize("features, rows", [([], []), ([feature()], [])])
def test_empty_handoff_or_predictions_produce_valid_empty_outputs(tmp_path, monkeypatch, features, rows):
    args = args_for(tmp_path, features=features)
    stub_stages(monkeypatch, rows)
    monkeypatch.setattr(novel_scoring, "score_novel_feature", lambda *a, **k: pytest.fail("scoring no predictions"))
    metadata = driver.run(args)
    summary = pd.read_csv(metadata["summary"])
    assert summary.empty
    assert set(driver.SCORE_COLUMNS) <= set(summary.columns)
    assert metadata["n_feature_summaries"] == 0
    assert json.loads((args.out_root / "run_novel_summary.json").read_text())["summary_schema_version"] == 2


def test_version_metadata_records_source_and_calibration_bytes(tmp_path, monkeypatch):
    calibration = tmp_path / "tools" / "identity_confidence_calibration.json"
    calibration.parent.mkdir()
    calibration.write_bytes(b'{"fixture":true}\n')
    monkeypatch.setattr(driver, "REPO_ROOT", tmp_path)
    def git(cmd, **kwargs):
        assert cmd[:3] == ["git", "-C", str(tmp_path)]
        return "a" * 40 + "\n" if cmd[3] == "rev-parse" else " M fixture.py\n"
    monkeypatch.setattr(driver.subprocess, "check_output", git)
    metadata = driver.version_metadata()
    assert metadata["companion_commit"] == "a" * 40
    assert metadata["companion_dirty"] is True
    assert metadata["calibration_checksum"] == hashlib.sha256(calibration.read_bytes()).hexdigest()
    assert metadata["calibration_path"] == str(calibration)


def test_dry_run_logs_every_stage_without_scoring_or_spawning(tmp_path, monkeypatch, caplog):
    args = args_for(tmp_path, stages=",".join(driver.ALL_STAGES))
    args.dry_run = True
    monkeypatch.setattr(driver, "score_predictions", lambda *a, **k: pytest.fail("dry run scored features"))
    monkeypatch.setattr(driver.subprocess, "run", lambda *a, **k: pytest.fail("dry run spawned a process"))
    monkeypatch.setattr(driver, "version_metadata", lambda: {"summary_schema_version": 2})
    with caplog.at_level(logging.INFO, logger=driver.LOGGER.name):
        metadata = driver.run(args)
    for stage in driver.ALL_STAGES:
        module = "structural_elucidation.reports.batch" if stage == "reports" else f"structural_elucidation.pipeline.{stage}"
        assert module in caplog.text
    assert metadata["n_feature_summaries"] == 0
    assert not (args.out_root / "scores").exists()


def test_denovo_success_without_new_output_does_not_use_stale_csv(tmp_path, monkeypatch):
    args = args_for(tmp_path, features=[feature()])
    manifest = build_inputs_from_handoff(args.handoff, tmp_path / "inputs")
    encoder = args.assets_root / "models" / "encoder" / "best_encoder.pt"
    encoder.parent.mkdir(parents=True)
    encoder.touch()
    stale = args.out_root / "denovo_lib400" / "novel_lib400_merged_predictions.csv"
    stale.parent.mkdir(parents=True)
    stale.write_text("stale fixture\n")
    monkeypatch.setattr(driver.StageRunner, "run", lambda *a, **k: 0)
    result = driver.run_denovo(driver.StageRunner(tmp_path), manifest.scopes[0], args.out_root, args.assets_root, [])
    assert result is None


@pytest.mark.parametrize("failure", ["serialization", "write"])
def test_one_score_persistence_failure_preserves_other_completed_scores(tmp_path, monkeypatch, failure):
    args = args_for(tmp_path)
    stub_stages(monkeypatch, [prediction(), prediction("X-900000002")])
    score = stub_scorer(monkeypatch)
    if failure == "serialization":
        def with_bad_payload(candidates, metadata, **kwargs):
            payload = score(candidates, metadata, **kwargs)
            if metadata["control_id"] == "X-900000002":
                payload["not_json"] = {object()}
            return payload
        monkeypatch.setattr(novel_scoring, "score_novel_feature", with_bad_payload)
    else:
        write = Path.write_text
        def fail_second(path, *args, **kwargs):
            if path.name == "X-900000002.json":
                raise OSError("fixture write failed")
            return write(path, *args, **kwargs)
        monkeypatch.setattr(Path, "write_text", fail_second)
    metadata = driver.run(args)
    summary = pd.read_csv(metadata["summary"]).set_index("feature_id")
    assert summary.index.tolist() == [1, 2]
    assert summary.loc[1, "se_confidence"] == 0.8
    assert summary.loc[1, "se_candidate_identity_confidence"] == 78.0
    assert (args.out_root / "scores" / f"{CONTROL_ID}.json").exists()
    assert metadata["scoring_errors"]


def test_failed_score_replacement_never_renders_stale_payload(tmp_path, monkeypatch):
    args = args_for(tmp_path, stages="rerank,reports")
    commands = stub_stages(monkeypatch, [prediction(), prediction("X-900000002")])
    stub_scorer(monkeypatch)
    stale = args.out_root / "scores" / "X-900000002.json"
    stale.parent.mkdir(parents=True)
    stale.write_text(json.dumps({"summary": {"se_confidence": 0.99, "se_reported_identification": "stale fixture"}}))
    write = Path.write_text
    def fail_replacement(path, *args, **kwargs):
        if path == stale:
            raise OSError("fixture replacement failed")
        return write(path, *args, **kwargs)
    monkeypatch.setattr(Path, "write_text", fail_replacement)
    metadata = driver.run(args)
    summary = pd.read_csv(metadata["summary"]).set_index("feature_id")
    assert summary.loc[1, "se_confidence"] == 0.8
    assert summary.loc[2, "se_confidence"] == 0.8
    assert Path(summary.loc[1, "se_report_path"]).is_file()
    assert pd.isna(summary.loc[2, "se_report_path"])
    assert not (args.out_root / "reports" / "X-900000002" / "index.html").exists()
    report_args = next(arguments for module, arguments in commands if module.endswith("reports.batch"))
    assert report_args[report_args.index("--control-ids") + 1] == CONTROL_ID
    assert metadata["scoring_errors"]


@pytest.mark.parametrize("through_run", [False, True])
def test_unreadable_calibration_is_recorded_without_losing_outputs(tmp_path, monkeypatch, through_run):
    calibration = tmp_path / "tools" / "identity_confidence_calibration.json"
    calibration.parent.mkdir()
    calibration.write_text("{}")
    monkeypatch.setattr(driver, "REPO_ROOT", tmp_path)
    monkeypatch.setattr(driver.subprocess, "check_output", lambda *a, **k: "a" * 40 + "\n")
    read = Path.read_bytes
    def denied(path):
        if path == calibration:
            raise PermissionError("fixture calibration unreadable")
        return read(path)
    monkeypatch.setattr(Path, "read_bytes", denied)
    if through_run:
        args = args_for(tmp_path)
        stub_stages(monkeypatch, [prediction()])
        stub_scorer(monkeypatch)
        metadata = driver.run(args)
        saved = json.loads((args.out_root / "run_novel_summary.json").read_text())
        assert saved["calibration_checksum"] is None
        summary = pd.read_csv(metadata["summary"])
        assert summary.loc[0, "se_rank1_smiles"] == SMILES
        assert summary.loc[0, "se_confidence"] == 0.8
        assert any("fixture calibration unreadable" in str(error) for error in metadata["scoring_errors"])
    else:
        metadata = driver.version_metadata()
    assert metadata["calibration_checksum"] is None
    assert "fixture calibration unreadable" in metadata["calibration_error"]
