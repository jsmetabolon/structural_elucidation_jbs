"""Score plausible-formula PubChem structures, then isomers of a score-0 seed."""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import sys
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from structural_elucidation.nce.verifier_v2 import predicted_spectrum_score_v2

TOOLS = Path(__file__).resolve().parent
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))
from spectrum_pilot import enumerate_constitutional_isomers

ELEMENT = re.compile(r"([A-Z][a-z]?)(\d*)")


def element_counts(formula: str) -> dict[str, int]:
    counts: dict[str, int] = {}
    for element, number in ELEMENT.findall(formula):
        counts[element] = counts.get(element, 0) + (int(number) if number else 1)
    return counts


def ppm(mass: float, target: float) -> float:
    if target <= 0:
        return math.inf
    return abs(mass - target) / target * 1e6


def score_smiles(smiles: str, record: pd.Series, reference_dir: str | None) -> float:
    return float(
        predicted_spectrum_score_v2(
            smiles,
            record.ion_json or "[]",
            record.loss_json or "[]",
            lib_id=int(float(record.library)),
            reference_dir=reference_dir,
        )
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cohort", type=Path, required=True)
    parser.add_argument("--formulas", type=Path, required=True)
    parser.add_argument("--structures", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--isomer-cap", type=int, default=50)
    args = parser.parse_args()
    reference_dir = os.environ.get("NCE_REFERENCE_DIR", "").strip() or None
    cohort = pd.read_csv(args.cohort, dtype=str, keep_default_na=False).set_index("x_compound_id")
    formulas = pd.read_csv(args.formulas, dtype=str, keep_default_na=False)
    structures = pd.read_csv(args.structures, dtype=str, keep_default_na=False)
    joined = formulas.merge(structures, on="formula", how="left")
    rows: list[dict[str, object]] = []
    from_formula = 0
    from_isomer = 0
    for cid, group in joined.groupby("x_compound_id"):
        record = cohort.loc[str(cid)]
        target = float(record.proton_neutral_mass)
        best = None
        best_score = 0.0
        closest = None
        seen: set[str] = set()
        for product in group.itertuples(index=False):
            smiles = str(getattr(product, "smiles", "") or "")
            if not smiles or smiles == "nan" or smiles in seen:
                continue
            seen.add(smiles)
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                continue
            formula = rdMolDescriptors.CalcMolFormula(mol)
            if element_counts(formula) != element_counts(str(product.formula)):
                continue
            mass = float(rdMolDescriptors.CalcExactMolWt(mol))
            error = ppm(mass, target)
            if error > 10:
                continue
            candidate = {
                "candidate_smiles": smiles,
                "candidate_formula": formula,
                "pubchem_cid": product.cid,
                "predicted_exact_mass": mass,
                "abs_mass_error_da": abs(mass - target),
                "abs_mass_error_ppm": error,
            }
            if closest is None or error < float(closest["abs_mass_error_ppm"]):
                closest = candidate
            score = score_smiles(smiles, record, reference_dir)
            if score > best_score:
                best_score = score
                best = {**candidate, "predicted_spectrum_score_v2": score, "structure_origin": "plausible_formula"}
        if best is None and closest is not None:
            for smiles in enumerate_constitutional_isomers(str(closest["candidate_smiles"]), limit=args.isomer_cap):
                score = score_smiles(smiles, record, reference_dir)
                if score > best_score:
                    best_score = score
                    mol = Chem.MolFromSmiles(smiles)
                    mass = float(rdMolDescriptors.CalcExactMolWt(mol)) if mol is not None else float(closest["predicted_exact_mass"])
                    best = {
                        **closest,
                        "candidate_smiles": smiles,
                        "pubchem_cid": "",
                        "predicted_exact_mass": mass,
                        "abs_mass_error_da": abs(mass - target),
                        "abs_mass_error_ppm": ppm(mass, target),
                        "predicted_spectrum_score_v2": score,
                        "structure_origin": "isomer",
                    }
        if best is None or best_score <= 0 or float(best["abs_mass_error_ppm"]) > 10:
            rows.append({"x_compound_id": cid, "passed": False, "best_score": best_score, "had_on_mass_seed": closest is not None})
            continue
        if best["structure_origin"] == "isomer":
            from_isomer += 1
        else:
            from_formula += 1
        rows.append({"x_compound_id": cid, "passed": True, "library": record.library, "ion_mode": record.ion_mode, "proton_neutral_mass": target, **best})
    scored = {row["x_compound_id"] for row in rows}
    for cid in cohort.index:
        if cid not in scored:
            rows.append({"x_compound_id": cid, "passed": False, "best_score": 0.0, "had_on_mass_seed": False})
    args.out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.out, index=False)
    passed = sum(1 for row in rows if row["passed"])
    print(json.dumps({"compounds": len(rows), "passed": passed, "from_formula": from_formula, "from_isomer": from_isomer}))


if __name__ == "__main__":
    main()
