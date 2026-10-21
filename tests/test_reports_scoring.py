"""Batch rendering forwards canonical headless confidence without recomputing it."""
import json
from pathlib import Path
from types import SimpleNamespace

from structural_elucidation.reports import batch


def test_batch_forwards_matching_precomputed_score_json(tmp_path, monkeypatch):
    control = "X-900000001"
    scores = tmp_path / "scores"
    scores.mkdir()
    scoring_path = scores / f"{control}.json"
    scoring_path.write_text(json.dumps({"schema_version": 1, "control_id": control}))
    out = tmp_path / "reports"
    predictions = tmp_path / "predictions.csv"
    args = batch.build_parser().parse_args(["--predictions", str(predictions), "--scoring-dir", str(scores)])
    def run(cmd, **kwargs):
        assert cmd[cmd.index("--scoring-json") + 1] == str(scoring_path)
        assert cmd[cmd.index("--control-id") + 1] == control
        compound_out = Path(cmd[cmd.index("--out-dir") + 1])
        compound_out.mkdir(parents=True)
        (compound_out / "index.html").write_text("fixture rendered")
        return SimpleNamespace(returncode=0, stdout="", stderr="")
    monkeypatch.setattr("subprocess.run", run)
    result = batch.build_one(control, "400", predictions, out, args)
    assert result.status == "generated"


def test_batch_missing_scoring_json_fails_without_rendering(tmp_path, monkeypatch):
    args = batch.build_parser().parse_args(["--predictions", "fixture.csv", "--scoring-dir", str(tmp_path)])
    monkeypatch.setattr("subprocess.run", lambda *a, **k: (_ for _ in ()).throw(AssertionError("unexpected subprocess")))
    result = batch.build_one("X-900000001", "400", tmp_path / "fixture.csv", tmp_path / "reports", args)
    assert result.status == "failed"
    assert "missing precomputed scoring JSON" in result.error_tail
