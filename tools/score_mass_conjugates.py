"""Build one conjugate per seed and keep products within 10 ppm that score above 0."""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from structural_elucidation.chemistry.conjugates import CONJUGATE_REACTIONS, apply_reaction
from structural_elucidation.nce.verifier_v2 import predicted_spectrum_score_v2


def _ppm(mass: float, target: float) -> float:
    if target <= 0:
        return math.inf
    return abs(mass - target) / target * 1e6


def products_for(seed: str, correction: str, target: float, max_products: int) -> list[dict[str, object]]:
    """Return conjugate products whose exact mass is within 10 ppm of the proton mass."""
    kept: list[dict[str, object]] = []
    seen: set[str] = set()
    for reaction in CONJUGATE_REACTIONS.get(correction, ()):
        for smiles in apply_reaction(seed, reaction, max_products=max_products):
            if smiles in seen:
                continue
            seen.add(smiles)
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                continue
            mass = float(rdMolDescriptors.CalcExactMolWt(mol))
            ppm = _ppm(mass, target)
            if ppm > 10:
                continue
            kept.append(
                {
                    "candidate_smiles": smiles,
                    "candidate_formula": rdMolDescriptors.CalcMolFormula(mol),
                    "predicted_exact_mass": mass,
                    "abs_mass_error_da": abs(mass - target),
                    "abs_mass_error_ppm": ppm,
                    "reaction": reaction,
                }
            )
    return kept


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--targets", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--max-products", type=int, default=12)
    args = parser.parse_args()
    reference_dir = os.environ.get("NCE_REFERENCE_DIR", "").strip() or None
    frame = pd.read_csv(args.targets, dtype=str, keep_default_na=False)
    seeds = frame[frame.action.eq("build_conjugate")].copy()
    rows: list[dict[str, object]] = []
    for record in seeds.itertuples(index=False):
        target = float(record.proton_neutral)
        built = products_for(record.seed_smiles, record.correction, target, args.max_products)
        best_score = 0.0
        best: dict[str, object] | None = None
        for product in built:
            score = float(
                predicted_spectrum_score_v2(
                    product["candidate_smiles"],
                    record.ion_json or "[]",
                    record.loss_json or "[]",
                    lib_id=int(float(record.lib_id)),
                    reference_dir=reference_dir,
                )
            )
            product["predicted_spectrum_score_v2"] = score
            if score > best_score:
                best_score = score
                best = product
        if best is None or best_score <= 0:
            rows.append(
                {
                    "x_compound_id": record.x_compound_id,
                    "passed": False,
                    "n_on_mass": len(built),
                    "best_score": best_score,
                }
            )
            continue
        rows.append(
            {
                "x_compound_id": record.x_compound_id,
                "passed": True,
                "lib_id": record.lib_id,
                "ion_mode": record.ion_mode,
                "proton_neutral": target,
                "correction": record.correction,
                "seed_smiles": record.seed_smiles,
                "n_on_mass": len(built),
                **best,
            }
        )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.out, index=False)
    passed = sum(1 for row in rows if row["passed"])
    print(json.dumps({"compounds": len(rows), "passed": passed, "out": str(args.out)}))


if __name__ == "__main__":
    main()
