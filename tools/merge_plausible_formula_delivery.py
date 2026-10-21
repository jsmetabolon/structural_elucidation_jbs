"""Append plausible-formula hits onto a copy of the formula delivery."""

from __future__ import annotations

import argparse
import math
import shutil
from pathlib import Path

import pandas as pd

from merge_ppm10_recovery import merge_new_controls
from merge_inventory_delivery import refresh_rollups


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--scores", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    for name in (
        "exec_summary_top10.csv",
        "exec_summary_top10_truncated.csv",
        "exec_summary_by_compound.csv",
    ):
        destination = args.out / name
        if not destination.exists():
            shutil.copy2(args.source / name, destination)
    top10 = pd.read_csv(args.out / "exec_summary_top10.csv", low_memory=False, dtype=str)
    by_compound = pd.read_csv(args.out / "exec_summary_by_compound.csv", low_memory=False, dtype=str)
    scores = pd.read_csv(args.scores, dtype=str, keep_default_na=False)
    passed = scores[scores["passed"].eq("True")].copy()
    meta = by_compound.set_index("x_compound_id")
    rows: list[dict[str, object]] = []
    for hit in passed.itertuples(index=False):
        cid = str(hit.x_compound_id)
        base = meta.loc[cid].to_dict() if cid in meta.index else {"x_compound_id": cid}
        formula = str(hit.candidate_formula)
        origin = str(hit.structure_origin)
        note = f"PubChem structure for {formula}" if origin != "isomer" else f"isomer of PubChem {formula}"
        error = float(hit.abs_mass_error_da)
        score = float(hit.predicted_spectrum_score_v2)
        row = {column: base.get(column, "") for column in top10.columns}
        row.update(
            {
                "x_compound_id": cid,
                "control_id": cid,
                "prediction_source": "formula_library",
                "delivery_note": note,
                "rank": "1",
                "candidate_smiles": hit.candidate_smiles,
                "candidate_formula": formula,
                "evidence_confidence": f"{score:.6g}",
                "evidence_confidence_basis": "predicted_spectrum_score_v2",
                "mass_match_score": f"{math.exp(-abs(error) / 0.05):.6g}",
                "abs_mass_error_da": f"{error:.6g}",
                "abs_mass_error_ppm": f"{float(hit.abs_mass_error_ppm):.6g}",
                "predicted_spectrum_score_v2": f"{score:.6g}",
            }
        )
        rows.append(row)
    extra = pd.DataFrame(rows)
    rank1_before = top10.copy()
    rank1_before["_rank"] = pd.to_numeric(rank1_before["rank"], errors="coerce")
    rank1_before = rank1_before.sort_values("_rank", kind="mergesort").drop_duplicates("x_compound_id")
    before_smiles = rank1_before.set_index("x_compound_id")["candidate_smiles"].to_dict() if len(rank1_before) else {}
    before_source = rank1_before.set_index("x_compound_id")["prediction_source"].to_dict() if len(rank1_before) else {}
    before_rank = rank1_before.set_index("x_compound_id")["rank"].to_dict() if len(rank1_before) else {}
    merged = merge_new_controls(top10.assign(control_id=top10["x_compound_id"]), extra)
    merged = merged.drop(columns=["control_id"])
    merged.to_csv(args.out / "exec_summary_top10.csv", index=False)
    truncated_columns = list(pd.read_csv(args.out / "exec_summary_top10_truncated.csv", nrows=0).columns)
    merged.reindex(columns=truncated_columns).to_csv(args.out / "exec_summary_top10_truncated.csv", index=False)
    new_ids = set(passed["x_compound_id"]) - set(before_smiles)
    for cid in new_ids:
        note = extra.loc[extra["x_compound_id"].eq(cid), "delivery_note"].iloc[0]
        by_compound.loc[by_compound["x_compound_id"].eq(cid), "prediction_source"] = "formula_library"
        by_compound.loc[by_compound["x_compound_id"].eq(cid), "delivery_note"] = note
    by_compound = refresh_rollups(by_compound, merged)
    by_compound.to_csv(args.out / "exec_summary_by_compound.csv", index=False)
    after_ranked = merged.copy()
    after_ranked["_rank"] = pd.to_numeric(after_ranked["rank"], errors="coerce")
    after = after_ranked.sort_values("_rank", kind="mergesort").drop_duplicates("x_compound_id").set_index("x_compound_id")
    changed = 0
    for cid, smiles in before_smiles.items():
        if cid not in after.index:
            changed += 1
            continue
        if (
            str(after.loc[cid, "candidate_smiles"]) != str(smiles)
            or str(after.loc[cid, "prediction_source"]) != str(before_source.get(cid))
            or str(after.loc[cid, "rank"]) != str(before_rank.get(cid))
        ):
            changed += 1
    print(f"published={len(new_ids)} prior_changed={changed} rows={len(merged)}")


if __name__ == "__main__":
    main()
