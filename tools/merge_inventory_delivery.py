#!/usr/bin/env python3
"""Build a new delivery that lists every eligible compound once.

Stored-mass rows win when the current summary has the compound. Otherwise the
earlier adduct-derived summary is kept. Eligible compounds in neither summary
are listed with no candidate rows. Neither source directory is modified.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from prepare_ppm10_recovery_inputs import LIBS, catalog_targets

STORED_MASS = "stored_mass"
ADDUCT_FALLBACK = "adduct_fallback"
EARLIER = "earlier_adduct_derive"
NO_PREDICTION = "no_prediction"
NO_INPUTS = "no feature row or hypothesis"
NO_SUPPORTED = "no candidate within 10 ppm with a spectrum score above 0"


def assign_source(
    compound_id: str,
    *,
    current_ids: set[str],
    settled_ids: set[str],
    earlier_ids: set[str],
    no_input_ids: set[str],
) -> tuple[str, str]:
    """Return the winning source and a note for compounds with no rows."""

    if compound_id in current_ids:
        source = STORED_MASS if compound_id in settled_ids else ADDUCT_FALLBACK
        return source, ""
    if compound_id in earlier_ids:
        return EARLIER, ""
    if compound_id in no_input_ids:
        return NO_PREDICTION, NO_INPUTS
    return NO_PREDICTION, NO_SUPPORTED


def _read(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, dtype=str, keep_default_na=False, encoding="utf-8")
    frame["x_compound_id"] = frame["x_compound_id"].astype(str)
    return frame


def _take(frame: pd.DataFrame, ids: set[str]) -> pd.DataFrame:
    return frame[frame["x_compound_id"].isin(ids)].copy()


def _known_index(path: Path) -> pd.DataFrame:
    known = pd.read_csv(path, dtype=str, keep_default_na=False, encoding="utf-8")
    known["Unknown"] = known["Unknown"].astype(str)
    known["LIB_ID"] = pd.to_numeric(known["LIB_ID"], errors="coerce")
    known = known[known["LIB_ID"].isin(LIBS)].copy()
    order = {lib: index for index, lib in enumerate(LIBS)}
    known["ord"] = known["LIB_ID"].map(order)
    return known.sort_values(["ord", "Unknown"]).drop_duplicates("Unknown").set_index("Unknown")


def stub_rows(ids: list[str], columns: list[str], notes: dict[str, str], known: pd.DataFrame) -> pd.DataFrame:
    """One blank candidate rollup per compound, with the known-X level filled in."""

    columns = list(columns)
    for column in ("metabolon_level", "metabolon_sublevel", "known_x_tentative_information", "known_x_cleaned_external_info"):
        if column not in columns:
            columns.append(column)
    rows = []
    for compound_id in ids:
        row = {column: "" for column in columns}
        row["x_compound_id"] = compound_id
        row["n_candidates"] = "0"
        row["prediction_source"] = NO_PREDICTION
        row["delivery_note"] = notes[compound_id]
        if compound_id in known.index:
            record = known.loc[compound_id]
            row["metabolon_level"] = record.get("Level", "")
            row["metabolon_sublevel"] = record.get("SubLevel", "")
            row["known_x_tentative_information"] = record.get("Tentative information", "")
            row["known_x_cleaned_external_info"] = record.get("Cleaned External Info", "")
        rows.append(row)
    return pd.DataFrame(rows, columns=columns)


def fill_missing_candidates(top10: pd.DataFrame, by_compound: pd.DataFrame, gated: pd.DataFrame) -> pd.DataFrame:
    """Add gated candidates for compounds the executive summary exported with none."""

    present = set(top10["x_compound_id"]) if len(top10) else set()
    owners = by_compound.loc[
        ~by_compound["prediction_source"].eq(NO_PREDICTION) & ~by_compound["x_compound_id"].isin(present),
        ["x_compound_id", "prediction_source"],
    ]
    pool_for = {STORED_MASS: "current", ADDUCT_FALLBACK: "current", EARLIER: "earlier"}
    owners = owners.assign(pool=owners["prediction_source"].map(pool_for))
    extra = gated[gated["control_id"].astype(str).isin(set(owners["x_compound_id"]))].copy()
    if extra.empty:
        return top10
    extra["control_id"] = extra["control_id"].astype(str)
    if "pool" in extra.columns:
        allowed = set(zip(owners["x_compound_id"], owners["pool"]))
        extra = extra[extra.apply(lambda row: (row["control_id"], row["pool"]) in allowed, axis=1)].copy()
    if "lib" in extra.columns and not extra.empty:
        order = {lib: index for index, lib in enumerate(LIBS)}
        extra["ord"] = extra["lib"].map(order)
        first = extra.groupby("control_id")["ord"].transform("min")
        extra = extra[extra["ord"].eq(first)].copy()
    if extra.empty:
        return top10
    smiles = extra["predicted_smiles"] if "predicted_smiles" in extra.columns else extra.get("candidate_smiles", "")
    formula = extra["candidate_formula"] if "candidate_formula" in extra.columns else extra.get("predicted_formula", "")
    built = pd.DataFrame({
        "x_compound_id": extra["control_id"].astype(str),
        "rank": extra["rank"].astype(str),
        "candidate_smiles": smiles.astype(str),
        "candidate_formula": formula.astype(str) if hasattr(formula, "astype") else "",
        "mass_match_score": extra["mass_match_score"].astype(str) if "mass_match_score" in extra.columns else "",
        "abs_mass_error_da": extra["abs_mass_error_da"].astype(str) if "abs_mass_error_da" in extra.columns else "",
        "abs_mass_error_ppm": extra["abs_mass_error_ppm"].astype(str) if "abs_mass_error_ppm" in extra.columns else "",
        "evidence_confidence": extra["evidence_confidence"].astype(str) if "evidence_confidence" in extra.columns else "",
        "evidence_confidence_basis": extra["evidence_confidence_basis"].astype(str) if "evidence_confidence_basis" in extra.columns else "",
    })
    shared = [column for column in by_compound.columns if column not in built.columns and not column.startswith("rank1_") and column != "n_candidates"]
    built = built.merge(by_compound[["x_compound_id", *shared]], on="x_compound_id", how="left")
    out = pd.concat([top10, built.reindex(columns=top10.columns)], ignore_index=True)
    out["_rank"] = pd.to_numeric(out["rank"], errors="coerce")
    return out.sort_values(["x_compound_id", "_rank"], kind="mergesort").drop(columns="_rank").reset_index(drop=True)


def refresh_rollups(by_compound: pd.DataFrame, top10: pd.DataFrame) -> pd.DataFrame:
    """Fill a blank compound rollup from the rank-1 candidate row."""

    out = by_compound.copy()
    counts = top10.groupby("x_compound_id").size() if len(top10) else pd.Series(dtype=int)
    out["n_candidates"] = out["x_compound_id"].map(counts).fillna(0).astype(int).astype(str)
    ranked = top10.copy()
    ranked["_rank"] = pd.to_numeric(ranked["rank"], errors="coerce")
    rank1 = ranked.sort_values("_rank", kind="mergesort").drop_duplicates("x_compound_id").set_index("x_compound_id")
    for dest, src in (
        ("rank1_candidate_smiles", "candidate_smiles"),
        ("rank1_candidate_formula", "candidate_formula"),
        ("rank1_mass_match_score", "mass_match_score"),
        ("rank1_abs_mass_error_da", "abs_mass_error_da"),
        ("rank1_abs_mass_error_ppm", "abs_mass_error_ppm"),
        ("rank1_evidence_confidence", "evidence_confidence"),
        ("rank1_evidence_confidence_basis", "evidence_confidence_basis"),
    ):
        if dest not in out.columns or src not in rank1.columns:
            continue
        blank = out[dest].fillna("").eq("")
        out.loc[blank, dest] = out.loc[blank, "x_compound_id"].map(rank1[src]).fillna("")
    return out


def _with_source(frame: pd.DataFrame, source_by_id: dict[str, str], note_by_id: dict[str, str]) -> pd.DataFrame:
    out = frame.copy()
    out.insert(1, "prediction_source", out["x_compound_id"].map(source_by_id))
    out.insert(2, "delivery_note", out["x_compound_id"].map(note_by_id).fillna(""))
    return out


def build_delivery(
    *,
    eligible: list[str],
    current_by_compound: pd.DataFrame,
    earlier_by_compound: pd.DataFrame,
    current_top10: pd.DataFrame,
    earlier_top10: pd.DataFrame,
    settled_ids: set[str],
    no_input_ids: set[str],
    known: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, int]]:
    current_ids = set(current_by_compound["x_compound_id"])
    earlier_ids = set(earlier_by_compound["x_compound_id"])
    chosen: dict[str, str] = {}
    notes: dict[str, str] = {}
    for compound_id in eligible:
        source, note = assign_source(
            compound_id,
            current_ids=current_ids,
            settled_ids=settled_ids,
            earlier_ids=earlier_ids,
            no_input_ids=no_input_ids,
        )
        chosen[compound_id] = source
        notes[compound_id] = note

    def ids_for(source: str) -> set[str]:
        return {compound_id for compound_id, value in chosen.items() if value == source}

    stored_ids = ids_for(STORED_MASS)
    fallback_ids = ids_for(ADDUCT_FALLBACK)
    earlier_kept = ids_for(EARLIER)
    missing = [compound_id for compound_id in eligible if chosen[compound_id] == NO_PREDICTION]
    current_kept = stored_ids | fallback_ids

    columns = list(dict.fromkeys(["x_compound_id", *current_by_compound.columns, *earlier_by_compound.columns]))
    by_compound = pd.concat(
        [
            _take(current_by_compound, current_kept).reindex(columns=columns),
            _take(earlier_by_compound, earlier_kept).reindex(columns=columns),
        ],
        ignore_index=True,
    )
    by_compound = _with_source(by_compound, chosen, notes)
    stubs = stub_rows(missing, list(by_compound.columns), notes, known)
    by_compound = pd.concat([by_compound, stubs], ignore_index=True)
    by_compound = by_compound.sort_values("x_compound_id", kind="mergesort").reset_index(drop=True)

    top_columns = list(dict.fromkeys(["x_compound_id", *current_top10.columns, *earlier_top10.columns]))
    top10 = pd.concat(
        [
            _take(current_top10, current_kept).reindex(columns=top_columns),
            _take(earlier_top10, earlier_kept).reindex(columns=top_columns),
        ],
        ignore_index=True,
    )
    top10 = _with_source(top10, chosen, notes)
    top10["_rank"] = pd.to_numeric(top10.get("rank"), errors="coerce")
    top10 = top10.sort_values(["x_compound_id", "_rank"], kind="mergesort").drop(columns="_rank").reset_index(drop=True)

    counts = {
        "eligible": len(eligible),
        "compounds": int(by_compound["x_compound_id"].nunique()),
        STORED_MASS: len(stored_ids),
        ADDUCT_FALLBACK: len(fallback_ids),
        EARLIER: len(earlier_kept),
        NO_PREDICTION: len(missing),
        "candidate_rows": int(len(top10)),
        "candidate_compounds": int(top10["x_compound_id"].nunique()) if len(top10) else 0,
    }
    return by_compound, top10, counts


def _load_gated(roots: list[tuple[Path | None, str]]) -> pd.DataFrame:
    wanted = {
        "control_id", "rank", "predicted_smiles", "candidate_smiles", "candidate_formula",
        "predicted_formula", "mass_match_score", "abs_mass_error_da", "abs_mass_error_ppm",
        "evidence_confidence", "evidence_confidence_basis",
    }
    frames = []
    for root, root_pool in roots:
        if root is None:
            continue
        for lib in LIBS:
            path = root / f"rerank_lib{lib}" / "unified_rerank_all_top10_ppm10.csv"
            if not path.exists() or path.stat().st_size <= 12:
                continue
            frame = pd.read_csv(path, dtype=str, keep_default_na=False, encoding="utf-8", usecols=lambda column: column in wanted)
            if "control_id" in frame.columns:
                frame["lib"] = lib
                frame["pool"] = root_pool
                frames.append(frame)
    if not frames:
        return pd.DataFrame(columns=["control_id"])
    return pd.concat(frames, ignore_index=True)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--current-exec", type=Path, required=True)
    parser.add_argument("--earlier-exec", type=Path, required=True)
    parser.add_argument("--catalog", type=Path, required=True)
    parser.add_argument("--settled-ids", type=Path, required=True)
    parser.add_argument("--fallback-selection", type=Path, required=True)
    parser.add_argument("--known-x", type=Path, required=True)
    parser.add_argument("--current-gated", type=Path, default=None)
    parser.add_argument("--earlier-gated", type=Path, default=None)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args(argv)

    eligible, _skipped = catalog_targets(pd.read_csv(args.catalog, low_memory=False))
    settled = set(args.settled_ids.read_text(encoding="utf-8").split())
    selection = json.loads(args.fallback_selection.read_text(encoding="utf-8"))
    searched = {str(item["control_id"]) for item in selection["assignments"]}
    fallback_ids = set(Path(args.settled_ids).with_name("fallback_ids.txt").read_text(encoding="utf-8").split())
    no_input_ids = fallback_ids - searched
    known = _known_index(args.known_x)

    by_compound, top10, counts = build_delivery(
        eligible=eligible,
        current_by_compound=_read(args.current_exec / "exec_summary_by_compound.csv"),
        earlier_by_compound=_read(args.earlier_exec / "exec_summary_by_compound.csv"),
        current_top10=_read(args.current_exec / "exec_summary_top10.csv"),
        earlier_top10=_read(args.earlier_exec / "exec_summary_top10.csv"),
        settled_ids=settled,
        no_input_ids=no_input_ids,
        known=known,
    )
    truncated_name = "exec_summary_top10_truncated.csv"
    current_truncated = args.current_exec / truncated_name
    earlier_truncated = args.earlier_exec / truncated_name
    if current_truncated.exists() and earlier_truncated.exists():
        _by, truncated, _counts = build_delivery(
            eligible=eligible,
            current_by_compound=_read(args.current_exec / "exec_summary_by_compound.csv"),
            earlier_by_compound=_read(args.earlier_exec / "exec_summary_by_compound.csv"),
            current_top10=_read(current_truncated),
            earlier_top10=_read(earlier_truncated),
            settled_ids=settled,
            no_input_ids=no_input_ids,
            known=known,
        )
    else:
        truncated = top10

    gated = _load_gated([(args.current_gated, "current"), (args.earlier_gated, "earlier")])
    if not gated.empty:
        top10 = fill_missing_candidates(top10, by_compound, gated)
        truncated = fill_missing_candidates(truncated, by_compound, gated)
        by_compound = refresh_rollups(by_compound, top10)
        counts["candidate_rows"] = int(len(top10))
        counts["candidate_compounds"] = int(top10["x_compound_id"].nunique())

    args.out_dir.mkdir(parents=True, exist_ok=True)
    by_compound.to_csv(args.out_dir / "exec_summary_by_compound.csv", index=False)
    top10.to_csv(args.out_dir / "exec_summary_top10.csv", index=False)
    truncated.to_csv(args.out_dir / truncated_name, index=False)
    (args.out_dir / "delivery_counts.json").write_text(json.dumps(counts, indent=2) + "\n", encoding="utf-8")
    print("[delivery]", json.dumps(counts))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
