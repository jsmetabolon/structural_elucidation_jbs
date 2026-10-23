"""Unit tests for the leave-one-out vote tuning grid (T2)."""
from __future__ import annotations

import importlib.util
import sqlite3
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent


def load_tool(name: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "tools" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    # Register before exec so dataclass string-annotation resolution can find
    # the defining module in sys.modules.
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


T = load_tool("pathway_vote_tuning")


def test_mw_prior_is_gaussian_and_off_at_zero_sigma():
    assert T.mw_prior(500.0, 0.0) == 1.0
    assert T.mw_prior(0.0, 50.0) == 1.0
    assert T.mw_prior(50.0, 50.0) == pytest.approx(0.6065, abs=1e-4)
    assert T.mw_prior(-50.0, 50.0) == T.mw_prior(50.0, 50.0)
    assert T.mw_prior(100.0, 50.0) < T.mw_prior(50.0, 50.0)


def test_inverse_frequency_weighting_can_flip_a_majority_class_vote():
    # Truth for item 0 is Amino Acid; its closest neighbor is a Lipid from the
    # dominant class. Unweighted, Lipid wins; inverse-frequency flips it.
    supers = ["Amino Acid"] + ["Lipid"] * 9 + ["Amino Acid"]
    population = {
        "supers": supers,
        "subs": [""] * len(supers),
        "mws": [100.0] * len(supers),
    }
    cache = [[] for _ in supers]
    cache[0] = [(0.9, 1), (0.6, 10)]  # strong Lipid neighbor, weaker AA neighbor

    base = T.VoteConfig("morgan2_bits", k=10, floor=0.2, power=2.0, mw_sigma=0.0, class_weighting="none")
    balanced = T.VoteConfig(
        "morgan2_bits", k=10, floor=0.2, power=2.0, mw_sigma=0.0, class_weighting="inverse_frequency"
    )

    assert T.evaluate_config(base, population, cache)["super_accuracy"] == 0.0
    assert T.evaluate_config(balanced, population, cache)["super_accuracy"] == 1.0


def test_mw_prior_downweights_distant_mass_neighbors():
    supers = ["Amino Acid", "Lipid", "Amino Acid"]
    population = {
        "supers": supers,
        "subs": [""] * 3,
        "mws": [100.0, 600.0, 105.0],  # the Lipid neighbor is 500 Da away
    }
    cache = [[(0.8, 1), (0.7, 2)], [], []]
    no_prior = T.VoteConfig("morgan2_bits", 10, 0.2, 2.0, mw_sigma=0.0, class_weighting="none")
    with_prior = T.VoteConfig("morgan2_bits", 10, 0.2, 2.0, mw_sigma=50.0, class_weighting="none")
    assert T.evaluate_config(no_prior, population, cache)["super_accuracy"] == 0.0
    assert T.evaluate_config(with_prior, population, cache)["super_accuracy"] == 1.0


@pytest.mark.slow
def test_grid_runs_end_to_end_and_resumes(tmp_path: Path):
    pytest.importorskip("rdkit")
    db_path = tmp_path / "ref.sqlite"
    con = sqlite3.connect(db_path)
    con.execute("CREATE TABLE metabolites (smiles TEXT, metabolon_super TEXT, metabolon_sub TEXT)")
    con.executemany(
        "INSERT INTO metabolites VALUES (?, ?, ?)",
        [
            ("NCC(=O)O", "Amino Acid", "Glycine Metabolism"),
            ("CC(N)C(=O)O", "Amino Acid", "Alanine Metabolism"),
            ("OCC(N)C(=O)O", "Amino Acid", "Serine Metabolism"),
            ("CCCCCC(=O)O", "Lipid", "Medium Chain Fatty Acid"),
            ("CCCCCCCC(=O)O", "Lipid", "Medium Chain Fatty Acid"),
        ],
    )
    con.commit()
    con.close()

    out = tmp_path / "grid.csv"
    args = type(
        "Args",
        (),
        {
            "reference_db": db_path,
            "fingerprints": ["morgan2_bits"],
            "max_configs": 6,
            "out": out,
        },
    )()
    assert T.run_grid(args) == 0
    with out.open() as handle:
        first = handle.read().strip().splitlines()
    assert len(first) == 1 + 6  # header + six configs

    # Resume: already-done configs are skipped, so max_configs=6 adds new ones.
    args.max_configs = 0
    assert T.run_grid(args) == 0
    with out.open() as handle:
        done = T.load_done_keys(out)
    total_configs = len(T.grid_configs(("morgan2_bits",)))
    assert len(done) == total_configs
