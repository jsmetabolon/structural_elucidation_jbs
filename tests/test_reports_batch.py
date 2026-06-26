"""Smoke tests for the batch report driver (no RDKit / heavy assets required).

Exercise the discovery, enumeration, and dry-run paths of
``structural_elucidation.reports.batch`` against tiny synthetic CSVs.
"""

from __future__ import annotations

import sys
from pathlib import Path

try:
    import pytest
except ImportError:  # pragma: no cover
    pytest = None  # type: ignore[assignment]

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


def _smoke(func):
    return pytest.mark.smoke(func) if pytest is not None else func


def _write_csv(path: Path, id_col: str, ids: list[str]) -> None:
    lines = [f"{id_col},predicted_smiles,rank"]
    for i, cid in enumerate(ids):
        lines.append(f"{cid},CCO,{(i % 3) + 1}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


@_smoke
def test_lib_from_path() -> None:
    from structural_elucidation.reports import batch

    assert batch.lib_from_path(Path("a/rerank_lib209/unified_rerank_all_top10.csv")) == "209"
    assert batch.lib_from_path(Path("a/b/c.csv")) == ""


@_smoke
def test_enumerate_dedup_and_id_fallback(tmp_path: Path) -> None:
    from structural_elucidation.reports import batch

    a = tmp_path / "rerank_lib209" / "unified_rerank_all_top10.csv"
    b = tmp_path / "rerank_lib400" / "unified_rerank_all_top10.csv"
    a.parent.mkdir(parents=True)
    b.parent.mkdir(parents=True)
    _write_csv(a, "control_id", ["X-1", "X-1", "X-2", "nan", ""])
    _write_csv(b, "x_compound_id", ["X-3", "X-2"])  # X-2 dup -> first (lib 209) wins

    rows = batch.enumerate_control_ids([a, b])
    ids = [r[0] for r in rows]
    libs = {r[0]: r[1] for r in rows}
    assert ids == ["X-1", "X-2", "X-3"]
    assert libs["X-2"] == "209"
    assert libs["X-3"] == "400"


@_smoke
def test_dry_run_counts(tmp_path: Path) -> None:
    from structural_elucidation.reports import batch

    csv = tmp_path / "rerank_lib305" / "unified_rerank_all_top10.csv"
    csv.parent.mkdir(parents=True)
    _write_csv(csv, "control_id", ["X-10", "X-11", "X-12", "X-10"])

    args = batch.build_parser().parse_args(
        ["--predictions", str(csv), "--out-dir", str(tmp_path / "reports"), "--dry-run"]
    )
    result = batch.run_batch(args)
    assert result["dry_run"] is True
    assert result["total_control_ids"] == 3
    assert result["selected"] == 3
    assert result["per_lib"].get("305") == 3


@_smoke
def test_limit_and_control_ids_filter(tmp_path: Path) -> None:
    from structural_elucidation.reports import batch

    csv = tmp_path / "rerank_lib402" / "unified_rerank_all_top10.csv"
    csv.parent.mkdir(parents=True)
    _write_csv(csv, "control_id", ["X-a", "X-b", "X-c", "X-d"])

    rows = batch.enumerate_control_ids([csv])
    assert batch.filter_control_ids(rows, None, 2) == rows[:2]
    only = batch.filter_control_ids(rows, {"X-b", "X-d"}, None)
    assert [r[0] for r in only] == ["X-b", "X-d"]


if __name__ == "__main__":
    test_lib_from_path()
    print("batch smoke: lib_from_path PASS")
