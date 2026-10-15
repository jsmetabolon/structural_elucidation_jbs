"""Recovery selection keeps only unscored compounds that already have inputs."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "tools" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


sys.path.insert(0, str(ROOT / "tools"))
assign_recovery = _load("prepare_ppm10_recovery_inputs").assign_recovery
merge_new_controls = _load("merge_ppm10_recovery").merge_new_controls
replace_ids = _load("merge_ppm10_recovery").replace_ids
classify = _load("select_ppm10_remaining").classify
denovo = _load("prepare_denovo_ppm10_inputs")


def _catalog() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {"lib_id": 209, "x_compound_id": "X-1", "level": "5"},
            {"lib_id": 209, "x_compound_id": "X-2", "level": "5"},
            {"lib_id": 400, "x_compound_id": "X-3", "level": "3"},
            {"lib_id": 209, "x_compound_id": "X-4", "level": "Exclude"},
            {"lib_id": 209, "x_compound_id": "X-5", "level": "5"},
            {"lib_id": 209, "x_compound_id": "X-6", "level": "5"},
            {"lib_id": 209, "x_compound_id": "X-7", "level": "Artifact"},
        ]
    )


def test_assign_recovery_uses_job_files_then_recovered_ms2() -> None:
    assignments, counts = assign_recovery(
        _catalog(),
        scored={"X-6"},
        job_features={209: {"X-1", "X-2", "X-4", "X-6"}, 305: set(), 400: set(), 402: set()},
        job_hypotheses={209: {"X-1", "X-2", "X-4", "X-6"}, 305: set(), 400: set(), 402: set()},
        job_manifest={209: {"X-1"}, 305: set(), 400: set(), 402: set()},
        recovered_features={209: set(), 305: set(), 400: {"X-3"}, 402: set()},
        recovered_hypotheses={209: set(), 305: set(), 400: {"X-3"}, 402: set()},
    )
    by_id = {row["control_id"]: row for row in assignments}
    assert set(by_id) == {"X-1", "X-2", "X-3"}
    assert by_id["X-1"]["origin"] == "job_manifest"
    assert by_id["X-2"]["origin"] == "job_features"
    assert by_id["X-3"]["origin"] == "recovered_ms2"
    assert by_id["X-3"]["lib"] == "400"
    assert counts["no_inputs"] == 1
    assert counts["skipped_level"] == 2
    assert counts["selected"] == 3


def test_classify_splits_unpublished_by_reason() -> None:
    groups = classify(
        kept=["X-1", "X-2", "X-3", "X-4", "X-5"],
        excluded={"X-9"},
        published={"X-1"},
        main_scored={"X-1", "X-2"},
        recovery_selected={"X-3", "X-4"},
        recovery_scored={"X-3"},
    )
    assert groups["library10"] == ["X-2"]
    assert groups["searched10_over"] == ["X-3"]
    assert groups["nohit10"] == ["X-4"]
    assert groups["no_inputs"] == ["X-5"]
    assert groups["excluded"] == ["X-9"]


def test_denovo_mass_error_uses_closest_derived_target() -> None:
    candidates = pd.DataFrame(
        {
            "control_id": ["X-1", "X-1", "X-2"],
            "predicted_smiles": ["CCO", "CCN", "CCC"],
            "predicted_exact_mass": [200.0, 250.0, 100.0],
            "abs_mass_error_da": [50.0, 0.0, 1.0],
        }
    )
    targets = pd.DataFrame({"control_id": ["X-1", "X-1"], "target_neutral_mass": [200.001, 260.0]})
    out = denovo.rescore_to_targets(candidates, targets, lambda error, tol: max(0.0, 1.0 - error / tol))
    assert out["control_id"].tolist() == ["X-1", "X-1"]
    assert abs(out.loc[0, "abs_mass_error_da"] - 0.001) < 1e-9
    assert out.loc[0, "target_neutral_mass"] == 200.001
    assert abs(out.loc[1, "abs_mass_error_da"] - 10.0) < 1e-9
    assert out.loc[1, "stored_abs_mass_error_da"] == 0.0


def test_denovo_selection_keeps_assigned_library() -> None:
    selection = {"assignments": [{"control_id": "X-1", "lib": "400"}, {"control_id": "X-2", "lib": "209"}]}
    assert denovo.selected_by_lib({"X-1"}, selection) == {400: {"X-1"}}


def test_replace_ids_drops_mass_only_compounds() -> None:
    main = pd.DataFrame({"control_id": ["X-1", "X-2"], "rank": [1, 1]})
    updates = pd.DataFrame({"control_id": ["X-2"], "rank": [1]})
    out = replace_ids(main, updates, {"X-2", "X-3"})
    assert out["control_id"].tolist() == ["X-1", "X-2"]


def test_merge_appends_only_new_controls() -> None:
    main = pd.DataFrame({"control_id": ["X-1"], "rank": [1]})
    extra = pd.DataFrame({"control_id": ["X-1", "X-2"], "rank": [1, 1]})
    merged = merge_new_controls(main, extra)
    assert merged["control_id"].tolist() == ["X-1", "X-2"]


def test_catalog_assigns_an_id_the_ppm_files_do_not_have() -> None:
    ion = _load("prepare_ion_denovo")
    catalog = pd.DataFrame(
        [
            {"lib_id": 400, "x_compound_id": "X-9"},
            {"lib_id": 209, "x_compound_id": "X-9"},
            {"lib_id": 309, "x_compound_id": "X-8"},
        ]
    )
    assigned = ion.assign_libraries({"X-1", "X-9"}, [(209, {"X-1"})], catalog)
    assert assigned == {"X-1": 209, "X-9": 209}


def test_gap_groups_do_not_overlap() -> None:
    gap = _load("write_ppm10_gap")
    missing = {"X-1", "X-2", "X-3", "X-4", "X-5"}
    groups = gap.disjoint_missing(
        missing,
        {
            "mass_only": {"X-1", "X-2"},
            "no_inputs": {"X-2", "X-3"},
            "nohit10": {"X-4"},
            "searched10_over": {"X-5", "X-4"},
            "library10": {"X-5", "X-1"},
        },
    )
    assert groups["mass_only"] == ["X-1", "X-2"]
    assert groups["no_inputs"] == ["X-3"]
    assert groups["nohit_decoded"] == ["X-4"]
    assert groups["searched10_over"] == ["X-5"]
    assert groups["library10_still"] == []
    assert groups["unlisted"] == []


def test_stored_mass_gate_sends_only_unsupported_compounds_to_fallback(tmp_path: Path) -> None:
    fallback = _load("select_adduct_fallback")
    scored = tmp_path / "scored"
    mass = 200.0
    error = 5.0 * mass / 1e6
    rows = []
    for cid, spectrum in (("X-keep", 0.4), ("X-mass", 0.0)):
        rows.append(
            {
                "control_id": cid,
                "predicted_smiles": cid,
                "unified_final_score": 0.5,
                "abs_mass_error_da": error,
                "predicted_exact_mass": mass,
                "predicted_spectrum_score_v2": spectrum,
                "mass_match_score": 0.9,
                "formula_prior_score": 0.1,
                "min_source_rank": 1,
            }
        )
    for lib in (209, 305, 400, 402):
        lib_dir = scored / f"rerank_lib{lib}"
        lib_dir.mkdir(parents=True)
        frame = pd.DataFrame(rows if lib == 209 else [])
        if frame.empty:
            frame = pd.DataFrame(columns=list(rows[0].keys()))
        frame.to_csv(lib_dir / "unified_candidates_scored.csv", index=False)
    catalog = tmp_path / "catalog.csv"
    pd.DataFrame(
        [
            {"lib_id": 209, "x_compound_id": "X-keep", "level": "5"},
            {"lib_id": 209, "x_compound_id": "X-mass", "level": "5"},
            {"lib_id": 209, "x_compound_id": "X-absent", "level": "5"},
            {"lib_id": 209, "x_compound_id": "X-skip", "level": "Exclude"},
        ]
    ).to_csv(catalog, index=False)
    out = tmp_path / "summary"
    assert fallback.main(["--scored-root", str(scored), "--catalog", str(catalog), "--out-root", str(out)]) == 0
    settled = set((out / "settled_ids.txt").read_text().split())
    missed = (out / "fallback_ids.txt").read_text().split()
    assert settled == {"X-keep"}
    assert set(missed) == {"X-mass", "X-absent"}
    published = pd.read_csv(out / "rerank_lib209" / "unified_rerank_all_top10_ppm10.csv")
    extra = pd.DataFrame({"control_id": ["X-keep", "X-new"], "rank": [1, 1]})
    merged = merge_new_controls(published, extra)
    assert merged["control_id"].tolist() == ["X-keep", "X-new"]
