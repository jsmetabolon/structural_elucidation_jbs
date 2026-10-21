"""List formulas closest to each compound's proton mass within 10 ppm."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

# REFERENCE_DIR is read when the formula-finder loader is imported.
os.environ.setdefault("REFERENCE_DIR", "/data/josh/spec2mol/shared/references/production")

import pandas as pd
from elucidation.formulaFinder import FormulaFinder


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cohort", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--ppm", type=float, default=10.0)
    parser.add_argument("--keep", type=int, default=8)
    parser.add_argument("--min-harmony", type=float, default=0.0)
    args = parser.parse_args()
    harmony = None
    if args.min_harmony > 0:
        from structural_elucidation.nce.evidence import formula_harmony_score

        harmony = formula_harmony_score
    cohort = pd.read_csv(args.cohort)
    finder = FormulaFinder.get_instance(209)
    rows: list[dict[str, object]] = []
    for record in cohort.itertuples(index=False):
        found = finder.fmlas_for_neutral_mass(float(record.proton_neutral_mass), ppm=args.ppm)
        if found.empty:
            continue
        found = found.assign(absppm=found["ppm"].abs()).sort_values("absppm")
        if harmony is not None:
            target = float(record.proton_neutral_mass)
            found = found[found["cand_fmla"].map(lambda formula: harmony(formula, target) >= args.min_harmony)]
        found = found.head(args.keep)
        for hit in found.itertuples(index=False):
            rows.append(
                {
                    "x_compound_id": record.x_compound_id,
                    "library": record.library,
                    "ion_mode": record.ion_mode,
                    "proton_neutral_mass": float(record.proton_neutral_mass),
                    "formula": hit.cand_fmla,
                    "formula_mass": float(hit.mass),
                    "ppm": float(hit.ppm),
                }
            )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    frame = pd.DataFrame(rows)
    frame.to_csv(args.out, index=False)
    compounds = frame["x_compound_id"].nunique() if not frame.empty else 0
    print(f"rows={len(frame)} compounds_with_formula={compounds} cohort={len(cohort)}")


if __name__ == "__main__":
    main()
