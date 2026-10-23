#!/usr/bin/env python3
"""Leave-one-out tuning grid for the neighbor-consensus pathway vote (T2).

Searches vote hyperparameters against the reference DB's own labels via
leave-one-out — the tuning signal that cannot leak known_X (which stays a
one-shot final check through ``tools/pathway_variant_benchmark.py``).

Baseline (committed, ``docs/paper/data/pathway/audit_pathway_loo.json``):
Morgan r2 bits, k=10, floor 0.30, power 2 → 86.7% super / 71.1% sub, with the
weak [0.3,0.5) similarity bin at 68.6%.

Searched axes:
- fingerprint: ``morgan2_bits`` (baseline), ``morgan3_bits``, ``morgan2_counts``,
  ``rdkit_fp``
- k: 5 / 10 / 20
- similarity floor: 0.2 / 0.3 / 0.4
- weight power: 1 / 2 / 3
- molecular-weight prior sigma (Da): off / 50 / 100 — multiplies each vote by
  ``exp(-0.5 * (dMW / sigma)^2)``, restoring the production MW-prior idea in a
  documented form
- class weighting: ``none`` / ``inverse_frequency`` (down-weights supers that
  dominate the reference DB)

Similarities are computed once per fingerprint variant (top-64 neighbors at
floor >= 0.2 cached), so the parameter grid itself is cheap; expect a few
minutes per fingerprint variant on the ~5.2k-metabolite reference DB. Results
append to a checkpoint CSV and completed configs are skipped on resume.

Gated on the metabolite synonym SQLite (``./fetch_assets.sh``).

Usage::

    python tools/pathway_vote_tuning.py run \
        --out docs/paper/data/pathway/vote_tuning_grid.csv
    python tools/pathway_vote_tuning.py report \
        --grid docs/paper/data/pathway/vote_tuning_grid.csv
"""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import math
import sqlite3
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parent.parent
REFERENCE_DB_CANDIDATES = (
    REPO_ROOT / "Metabolite Synonym Database" / "metabolite_synonyms_merged_with_mtx_giga.sqlite",
    REPO_ROOT / "vendor" / "metabolite_synonyms.sqlite",
)
CACHE_MIN_FLOOR = 0.2
CACHE_TOP_N = 64

FINGERPRINTS = ("morgan2_bits", "morgan3_bits", "morgan2_counts", "rdkit_fp")
K_GRID = (5, 10, 20)
FLOOR_GRID = (0.2, 0.3, 0.4)
POWER_GRID = (1.0, 2.0, 3.0)
MW_SIGMA_GRID = (0.0, 50.0, 100.0)  # 0.0 = prior off
CLASS_WEIGHTING = ("none", "inverse_frequency")


@dataclass(frozen=True)
class VoteConfig:
    fingerprint: str
    k: int
    floor: float
    power: float
    mw_sigma: float
    class_weighting: str

    def key(self) -> str:
        return (
            f"{self.fingerprint}|k={self.k}|floor={self.floor}|power={self.power}"
            f"|mw={self.mw_sigma}|class={self.class_weighting}"
        )


def norm(value: Any) -> str:
    text = str("" if value is None else value).strip().lower()
    return "" if text in {"nan", "none", "null", "<na>"} else " ".join(text.split())


def find_reference_db() -> Path | None:
    for candidate in REFERENCE_DB_CANDIDATES:
        if candidate.exists():
            return candidate
    return None


def load_reference(db_path: Path) -> list[dict[str, str]]:
    con = sqlite3.connect(str(db_path))
    rows = con.execute(
        "SELECT smiles, metabolon_super, metabolon_sub FROM metabolites "
        "WHERE TRIM(COALESCE(smiles,'')) <> '' AND TRIM(COALESCE(metabolon_super,'')) <> ''"
    ).fetchall()
    con.close()
    return [
        {"smiles": s, "metabolon_super": (sup or "").strip(), "metabolon_sub": (sub or "").strip()}
        for s, sup, sub in rows
    ]


def make_fingerprinter(variant: str):
    from rdkit.Chem import rdFingerprintGenerator

    if variant == "morgan2_bits":
        gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
        return gen.GetFingerprint
    if variant == "morgan3_bits":
        gen = rdFingerprintGenerator.GetMorganGenerator(radius=3, fpSize=2048)
        return gen.GetFingerprint
    if variant == "morgan2_counts":
        gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
        return gen.GetCountFingerprint
    if variant == "rdkit_fp":
        gen = rdFingerprintGenerator.GetRDKitFPGenerator(fpSize=2048)
        return gen.GetFingerprint
    raise ValueError(f"unknown fingerprint variant: {variant}")


def prepare_population(reference: list[dict[str, str]], variant: str) -> dict[str, list[Any]]:
    """Parse molecules once; fingerprints, labels and exact weights aligned by index."""

    from rdkit import Chem
    from rdkit.Chem.Descriptors import ExactMolWt

    fingerprinter = make_fingerprinter(variant)
    fps: list[Any] = []
    supers: list[str] = []
    subs: list[str] = []
    weights: list[float] = []
    for row in reference:
        mol = Chem.MolFromSmiles(str(row["smiles"]).split("|")[0].strip())
        if mol is None:
            continue
        fps.append(fingerprinter(mol))
        supers.append(row["metabolon_super"])
        subs.append(row["metabolon_sub"])
        weights.append(ExactMolWt(mol))
    return {"fps": fps, "supers": supers, "subs": subs, "mws": weights}


def build_neighbor_cache(fps: list[Any]) -> list[list[tuple[float, int]]]:
    """Top-``CACHE_TOP_N`` neighbors >= ``CACHE_MIN_FLOOR`` for every metabolite.

    One full similarity pass per fingerprint variant; every (k, floor, power,
    prior) configuration re-votes from this cache, which keeps the grid cheap.
    """

    from rdkit.DataStructs import BulkTanimotoSimilarity

    cache: list[list[tuple[float, int]]] = []
    n = len(fps)
    for i in range(n):
        sims = BulkTanimotoSimilarity(fps[i], fps)
        neighbors = [(sims[j], j) for j in range(n) if j != i and sims[j] >= CACHE_MIN_FLOOR]
        neighbors.sort(reverse=True)
        cache.append(neighbors[:CACHE_TOP_N])
    return cache


def mw_prior(delta_mw: float, sigma: float) -> float:
    """Gaussian molecular-weight proximity prior; 1.0 when the prior is off."""

    if sigma <= 0:
        return 1.0
    return math.exp(-0.5 * (delta_mw / sigma) ** 2)


def evaluate_config(
    config: VoteConfig,
    population: dict[str, list[Any]],
    cache: list[list[tuple[float, int]]],
) -> dict[str, Any]:
    """Leave-one-out accuracy of one configuration using the neighbor cache."""

    supers = population["supers"]
    subs = population["subs"]
    mws = population["mws"]
    class_counts = Counter(supers)
    n = len(supers)

    evaluated = correct_super = correct_sub = sub_evaluable = 0
    weak_bin = [0, 0]  # correct, n for best similarity in [0.3, 0.5)
    for i in range(n):
        neighbors = [(s, j) for s, j in cache[i] if s >= config.floor][: config.k]
        if not neighbors:
            continue
        evaluated += 1
        votes: dict[str, float] = defaultdict(float)
        for sim, j in neighbors:
            weight = sim**config.power * mw_prior(mws[j] - mws[i], config.mw_sigma)
            if config.class_weighting == "inverse_frequency":
                weight /= class_counts[supers[j]]
            votes[supers[j]] += weight
        pred_super = max(votes, key=lambda name: votes[name])
        ok = norm(pred_super) == norm(supers[i])
        correct_super += int(ok)
        best_sim = neighbors[0][0]
        if 0.3 <= best_sim < 0.5:
            weak_bin[0] += int(ok)
            weak_bin[1] += 1
        if subs[i]:
            sub_votes: dict[str, float] = defaultdict(float)
            for sim, j in neighbors:
                if supers[j] == pred_super and subs[j]:
                    weight = sim**config.power * mw_prior(mws[j] - mws[i], config.mw_sigma)
                    if config.class_weighting == "inverse_frequency":
                        weight /= class_counts[supers[j]]
                    sub_votes[subs[j]] += weight
            pred_sub = max(sub_votes, key=lambda name: sub_votes[name]) if sub_votes else ""
            sub_evaluable += 1
            correct_sub += int(norm(pred_sub) == norm(subs[i]))

    return {
        **asdict(config),
        "n_evaluated": evaluated,
        "coverage": round(evaluated / n, 4) if n else 0.0,
        "super_accuracy": round(correct_super / evaluated, 4) if evaluated else 0.0,
        "sub_accuracy": round(correct_sub / sub_evaluable, 4) if sub_evaluable else 0.0,
        "weak_bin_accuracy": round(weak_bin[0] / weak_bin[1], 4) if weak_bin[1] else None,
        "weak_bin_n": weak_bin[1],
    }


def grid_configs(fingerprints: tuple[str, ...]) -> list[VoteConfig]:
    return [
        VoteConfig(fp, k, floor, power, mw, cw)
        for fp, k, floor, power, mw, cw in itertools.product(
            fingerprints, K_GRID, FLOOR_GRID, POWER_GRID, MW_SIGMA_GRID, CLASS_WEIGHTING
        )
    ]


def load_done_keys(out_path: Path) -> set[str]:
    if not out_path.exists():
        return set()
    with out_path.open(newline="") as handle:
        return {
            VoteConfig(
                row["fingerprint"],
                int(row["k"]),
                float(row["floor"]),
                float(row["power"]),
                float(row["mw_sigma"]),
                row["class_weighting"],
            ).key()
            for row in csv.DictReader(handle)
        }


def run_grid(args: argparse.Namespace) -> int:
    db_path = args.reference_db or find_reference_db()
    if db_path is None or not Path(db_path).exists():
        print(
            "[tuning] SKIPPED: metabolite synonym SQLite not present "
            f"(looked for {[str(p) for p in REFERENCE_DB_CANDIDATES]}). "
            "Fetch it with ./fetch_assets.sh and re-run."
        )
        return 0

    reference = load_reference(Path(db_path))
    fingerprints = tuple(args.fingerprints) if args.fingerprints else FINGERPRINTS
    configs = grid_configs(fingerprints)
    done = load_done_keys(args.out)
    todo = [c for c in configs if c.key() not in done]
    print(f"[tuning] {len(reference)} reference metabolites; {len(todo)}/{len(configs)} configs to run")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(asdict(configs[0]).keys()) + [
        "n_evaluated",
        "coverage",
        "super_accuracy",
        "sub_accuracy",
        "weak_bin_accuracy",
        "weak_bin_n",
    ]
    write_header = not args.out.exists()
    ran = 0
    for fingerprint in fingerprints:
        fp_todo = [c for c in todo if c.fingerprint == fingerprint]
        if not fp_todo:
            continue
        print(f"[tuning] fingerprinting + similarity pass: {fingerprint}")
        population = prepare_population(reference, fingerprint)
        cache = build_neighbor_cache(population["fps"])
        with args.out.open("a", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            if write_header:
                writer.writeheader()
                write_header = False
            for config in fp_todo:
                writer.writerow(evaluate_config(config, population, cache))
                handle.flush()
                ran += 1
                if args.max_configs and ran >= args.max_configs:
                    print(f"[tuning] stopping at --max-configs={args.max_configs}")
                    return 0
    print(f"[tuning] done; wrote {ran} new rows to {args.out}")
    return 0


def run_report(args: argparse.Namespace) -> int:
    if not args.grid.exists():
        print(f"[tuning] no grid results at {args.grid}; run the grid first")
        return 1
    with args.grid.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    for row in rows:
        row["super_accuracy"] = float(row["super_accuracy"])
        row["sub_accuracy"] = float(row["sub_accuracy"])
        row["coverage"] = float(row["coverage"])
    # Rank by super accuracy; require coverage within 2 points of the baseline
    # config so a high floor cannot win by abstaining on hard cases.
    baseline = next(
        (
            r
            for r in rows
            if r["fingerprint"] == "morgan2_bits"
            and r["k"] == "10"
            and float(r["floor"]) == 0.3
            and float(r["power"]) == 2.0
            and float(r["mw_sigma"]) == 0.0
            and r["class_weighting"] == "none"
        ),
        None,
    )
    floor_coverage = (baseline["coverage"] - 0.02) if baseline else 0.0
    eligible = [r for r in rows if r["coverage"] >= floor_coverage]
    top = sorted(eligible, key=lambda r: (-r["super_accuracy"], -r["sub_accuracy"]))[: args.top]
    print(json.dumps({"baseline": baseline, "top_configs": top}, indent=2))
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)

    pr = sub.add_parser("run", help="run the leave-one-out tuning grid (needs reference DB)")
    pr.add_argument("--reference-db", type=Path, default=None)
    pr.add_argument("--fingerprints", nargs="*", default=None, choices=list(FINGERPRINTS))
    pr.add_argument("--max-configs", type=int, default=0)
    pr.add_argument("--out", type=Path, default=REPO_ROOT / "docs/paper/data/pathway/vote_tuning_grid.csv")

    pp = sub.add_parser("report", help="rank completed grid configs")
    pp.add_argument("--grid", type=Path, default=REPO_ROOT / "docs/paper/data/pathway/vote_tuning_grid.csv")
    pp.add_argument("--top", type=int, default=10)

    args = parser.parse_args()
    if args.command == "run":
        return run_grid(args)
    if args.command == "report":
        return run_report(args)
    raise AssertionError(args.command)


if __name__ == "__main__":
    raise SystemExit(main())
