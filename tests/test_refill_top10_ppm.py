"""Unit tests for ppm-gated top-10 refill (no EC2 / RDKit)."""

from __future__ import annotations

import csv
import importlib.util
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def load_tool():
    path = ROOT / "tools" / "refill_top10_ppm.py"
    spec = importlib.util.spec_from_file_location("refill_top10_ppm", path)
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    sys.modules["refill_top10_ppm"] = module
    spec.loader.exec_module(module)
    return module


R = load_tool()


def _write(path: Path, rows: list[dict], fieldnames: list[str] | None = None) -> None:
    names = fieldnames or list(rows[0].keys())
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=names)
        writer.writeheader()
        writer.writerows(rows)


def test_ppm_conversion_matches_exec_summary():
    import pandas as pd

    error = pd.Series([0.002, None])
    mass = pd.Series([200.0, 200.0])
    ppm = R.abs_mass_error_ppm(error, mass)
    assert abs(float(ppm.iloc[0]) - 10.0) < 1e-9
    assert pd.isna(ppm.iloc[1])


def test_refill_drops_high_ppm_and_promotes_deeper_ranks(tmp_path: Path):
    lib = tmp_path / "rerank_lib209"
    lib.mkdir()
    rows = []
    # Ranks 1-2 fail ppm; 3-12 pass. After gate, new rank 1 is old rank 3.
    for rank in range(1, 13):
        ppm = 50.0 if rank <= 2 else 1.0
        mass = 200.0
        error_da = ppm * mass / 1e6
        rows.append(
            {
                "control_id": "X-1",
                "predicted_smiles": f"C{rank}",
                "unified_final_score": 1.0 - rank * 0.01,
                "abs_mass_error_da": error_da,
                "predicted_exact_mass": mass,
                "predicted_spectrum_score_v2": 0.5,
                "mass_match_score": 0.9,
                "formula_prior_score": 0.1,
                "min_source_rank": rank,
            }
        )
    _write(lib / "unified_candidates_scored.csv", rows)
    (lib / "non_s0_candidates_scored.csv").write_text("")
    _write(
        lib / "unified_rerank_all_top10.csv",
        [{**row, "rank": i} for i, row in enumerate(rows[:10], start=1)],
    )
    summary = R.refill_lib(
        lib, ppm_ceiling=10.0, top_k=10, out_name="unified_rerank_all_top10_ppm10.csv"
    )
    import pandas as pd

    out = pd.read_csv(lib / "unified_rerank_all_top10_ppm10.csv")
    assert list(out["rank"]) == list(range(1, 11))
    assert list(out["predicted_smiles"]) == [f"C{r}" for r in range(3, 13)]
    assert (out["abs_mass_error_ppm"] <= 10).all()
    assert int(summary["n_rows_from_rank_gt_top_k"]) == 2
    assert int(summary["n_rank1_smiles_changed"]) == 1
    assert int(summary["n_empty"]) == 0


def test_blank_mass_error_on_a_proton_hit_is_filled(tmp_path: Path):
    lib = tmp_path / "rerank_lib400"
    lib.mkdir()
    exact = 559.1748486
    hypothesis = 559.174114
    rows = [
        {
            "control_id": "X-24228",
            "predicted_smiles": "PROTON",
            "unified_final_score": 0.4,
            "abs_mass_error_da": "",
            "predicted_exact_mass": exact,
            "mass_hypothesis_neutral_mass": hypothesis,
            "mass_hypothesis_adduct": "[M+H]+",
            "predicted_spectrum_score_v2": 0.029,
            "mass_match_score": 0.5,
            "formula_prior_score": 0.1,
            "min_source_rank": 1,
        },
        {
            "control_id": "X-24228",
            "predicted_smiles": "SODIUM",
            "unified_final_score": 0.9,
            "abs_mass_error_da": "",
            "predicted_exact_mass": exact,
            "mass_hypothesis_neutral_mass": exact - 22.0,
            "mass_hypothesis_adduct": "[M+Na]+",
            "predicted_spectrum_score_v2": 0.9,
            "mass_match_score": 0.9,
            "formula_prior_score": 0.1,
            "min_source_rank": 2,
        },
    ]
    _write(lib / "unified_candidates_scored.csv", rows)
    summary = R.refill_lib(
        lib,
        ppm_ceiling=10.0,
        top_k=10,
        out_name="unified_rerank_all_top10_ppm10.csv",
        min_spectrum_score=0.0,
    )
    import pandas as pd

    out = pd.read_csv(lib / "unified_rerank_all_top10_ppm10.csv")
    assert list(out["predicted_smiles"]) == ["PROTON"]
    assert float(out.iloc[0]["abs_mass_error_ppm"]) < 10
    assert int(summary["n_passing_rows"]) == 1


def test_spectrum_floor_drops_a_mass_only_candidate(tmp_path: Path):
    lib = tmp_path / "rerank_lib209"
    lib.mkdir()
    mass = 200.0
    rows = [
        {
            "control_id": "X-1",
            "predicted_smiles": "MASS",
            "unified_final_score": 0.9,
            "abs_mass_error_da": 0.001,
            "predicted_exact_mass": mass,
            "predicted_spectrum_score_v2": 0.0,
            "mass_match_score": 0.99,
            "formula_prior_score": 0.1,
            "min_source_rank": 1,
        },
        {
            "control_id": "X-1",
            "predicted_smiles": "IONS",
            "unified_final_score": 0.2,
            "abs_mass_error_da": 0.001,
            "predicted_exact_mass": mass,
            "predicted_spectrum_score_v2": 0.4,
            "mass_match_score": 0.5,
            "formula_prior_score": 0.1,
            "min_source_rank": 2,
        },
    ]
    _write(lib / "unified_candidates_scored.csv", rows)
    summary = R.refill_lib(
        lib,
        ppm_ceiling=10.0,
        top_k=10,
        out_name="unified_rerank_all_top10_ppm10.csv",
        min_spectrum_score=0.0,
    )
    import pandas as pd

    out = pd.read_csv(lib / "unified_rerank_all_top10_ppm10.csv")
    assert list(out["predicted_smiles"]) == ["IONS"]
    assert int(summary["n_compounds_written"]) == 1


def test_missing_ppm_is_dropped(tmp_path: Path):
    lib = tmp_path / "rerank_lib305"
    lib.mkdir()
    rows = [
        {
            "control_id": "X-2",
            "predicted_smiles": "CC",
            "unified_final_score": 0.9,
            "abs_mass_error_da": "",
            "predicted_exact_mass": 100.0,
            "predicted_spectrum_score_v2": 0.1,
            "mass_match_score": 0.1,
            "formula_prior_score": 0.1,
            "min_source_rank": 1,
        },
        {
            "control_id": "X-2",
            "predicted_smiles": "CCC",
            "unified_final_score": 0.8,
            "abs_mass_error_da": 0.0005,
            "predicted_exact_mass": 100.0,
            "predicted_spectrum_score_v2": 0.1,
            "mass_match_score": 0.1,
            "formula_prior_score": 0.1,
            "min_source_rank": 2,
        },
    ]
    _write(lib / "unified_candidates_scored.csv", rows)
    summary = R.refill_lib(
        lib, ppm_ceiling=10.0, top_k=10, out_name="unified_rerank_all_top10_ppm10.csv"
    )
    import pandas as pd

    out = pd.read_csv(lib / "unified_rerank_all_top10_ppm10.csv")
    assert list(out["predicted_smiles"]) == ["CCC"]
    assert int(out.iloc[0]["rank"]) == 1
    assert int(summary["n_short_of_top_k"]) == 1


def test_cli_writes_all_libs(tmp_path: Path):
    root = tmp_path / "eq"
    for lib in ("209", "305"):
        d = root / f"rerank_lib{lib}"
        d.mkdir(parents=True)
        _write(
            d / "unified_candidates_scored.csv",
            [
                {
                    "control_id": f"X-{lib}",
                    "predicted_smiles": "CCO",
                    "unified_final_score": 0.5,
                    "abs_mass_error_da": 0.0001,
                    "predicted_exact_mass": 100.0,
                    "predicted_spectrum_score_v2": 0.2,
                    "mass_match_score": 0.9,
                    "formula_prior_score": 0.1,
                    "min_source_rank": 1,
                }
            ],
        )
    assert R.main(["--equalization-root", str(root), "--libs", "209,305"]) == 0
    assert (root / "rerank_lib209" / "unified_rerank_all_top10_ppm10.csv").exists()
    assert (root / "rerank_lib305" / "unified_rerank_all_top10_ppm10.csv").exists()
    ungated = root / "rerank_lib209" / "unified_rerank_all_top10.csv"
    assert not ungated.exists()
