#!/usr/bin/env python3
"""Build candidate-generation inputs for compounds the 10 ppm run never scored.

Selects catalog compounds on libraries 209, 305, 400, and 402 that are absent
from the scored pools, drops Exclude and Artifact, and keeps a compound only
when a feature row and a neutral-mass hypothesis already exist. Job feature
files win. ``pipeline_equalization_v1_recovered_ms2`` fills compounds the job
files do not contain.

The written gap manifests are the control lists. Candidate generation should
run them with ``--ppm-tolerance 10`` and ``--derive-neutral-mass-from-adduct``.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

LIBS = (209, 305, 400, 402)
ION = {209: "negative", 305: "negative", 400: "positive", 402: "positive"}
SKIP_LEVELS = {"exclude", "artifact"}


def _ids(frame: pd.DataFrame) -> set[str]:
    if frame is None or frame.empty:
        return set()
    column = "control_id" if "control_id" in frame.columns else "x_compound_id"
    if column not in frame.columns:
        return set()
    return {str(value).strip() for value in frame[column].dropna() if str(value).strip()}


def catalog_targets(catalog: pd.DataFrame) -> tuple[list[str], set[str]]:
    """Return run-library catalog IDs (first library wins) and the Exclude/Artifact IDs."""

    work = catalog.copy()
    work["lib_id"] = pd.to_numeric(work["lib_id"], errors="coerce")
    work = work[work["lib_id"].isin(LIBS)]
    work["x_compound_id"] = work["x_compound_id"].astype(str)
    order = {lib: index for index, lib in enumerate(LIBS)}
    work["ord"] = work["lib_id"].map(order)
    work = work.sort_values(["ord", "x_compound_id"]).drop_duplicates("x_compound_id")
    level = work.get("level", pd.Series("", index=work.index)).astype(str).str.strip().str.lower()
    skipped = set(work.loc[level.isin(SKIP_LEVELS), "x_compound_id"])
    kept = [cid for cid in work["x_compound_id"] if cid not in skipped]
    return kept, skipped


def assign_inputs(
    ids: list[str],
    job_features: dict[int, set[str]],
    job_hypotheses: dict[int, set[str]],
    job_manifest: dict[int, set[str]],
    recovered_features: dict[int, set[str]],
    recovered_hypotheses: dict[int, set[str]],
) -> list[dict[str, str]]:
    """Place each ID on the earliest run library that has a feature row and a hypothesis.

    Job files are searched before the recovered-MS2 tables. IDs with neither are
    left out of the result.
    """

    assignments: list[dict[str, str]] = []
    for cid in ids:
        placed = False
        for lib in LIBS:
            if cid in job_features.get(lib, set()) and cid in job_hypotheses.get(lib, set()):
                origin = "job_manifest" if cid in job_manifest.get(lib, set()) else "job_features"
                assignments.append({"control_id": cid, "lib": str(lib), "ion": ION[lib], "origin": origin})
                placed = True
                break
        if placed:
            continue
        for lib in LIBS:
            if cid in recovered_features.get(lib, set()) and cid in recovered_hypotheses.get(lib, set()):
                assignments.append(
                    {"control_id": cid, "lib": str(lib), "ion": ION[lib], "origin": "recovered_ms2"}
                )
                break
    return assignments


def assign_recovery(
    catalog: pd.DataFrame,
    scored: set[str],
    job_features: dict[int, set[str]],
    job_hypotheses: dict[int, set[str]],
    job_manifest: dict[int, set[str]],
    recovered_features: dict[int, set[str]],
    recovered_hypotheses: dict[int, set[str]],
) -> tuple[list[dict[str, str]], dict[str, int]]:
    """Return one assignment per recoverable compound, plus skip counts."""

    kept, skipped = catalog_targets(catalog)
    skipped_level = len(skipped)
    never = [cid for cid in kept if cid not in scored]
    assignments = assign_inputs(
        never, job_features, job_hypotheses, job_manifest, recovered_features, recovered_hypotheses
    )
    counts = {
        "skipped_level": skipped_level,
        "never_scored": len(never),
        "job_manifest": sum(1 for row in assignments if row["origin"] == "job_manifest"),
        "job_features": sum(1 for row in assignments if row["origin"] == "job_features"),
        "recovered_ms2": sum(1 for row in assignments if row["origin"] == "recovered_ms2"),
        "no_inputs": len(never) - len(assignments),
        "selected": len(assignments),
    }
    return assignments, counts


def _read_ids(path: Path) -> set[str]:
    if not path.exists():
        return set()
    header = pd.read_csv(path, nrows=0).columns
    column = "control_id" if "control_id" in header else "x_compound_id"
    if column not in header:
        return set()
    return _ids(pd.read_csv(path, usecols=[column]))


def _filter_rows(path: Path, keep: set[str]) -> pd.DataFrame:
    frame = pd.read_csv(path, low_memory=False)
    column = "control_id" if "control_id" in frame.columns else None
    if column is None:
        return frame.iloc[0:0]
    return frame[frame[column].astype(str).isin(keep)].copy()


def write_inputs(
    assignments: list[dict[str, str]],
    feature_dir: Path,
    recovered_dir: Path,
    out_dir: Path,
) -> list[str]:
    """Write per-library feature, hypothesis, and manifest slices. Return libs used."""

    out_dir.mkdir(parents=True, exist_ok=True)
    used: list[str] = []
    for lib in LIBS:
        chosen = [row for row in assignments if int(row["lib"]) == lib]
        if not chosen:
            continue
        ion = ION[lib]
        keep = {row["control_id"] for row in chosen}
        job_mask = {row["control_id"] for row in chosen if row["origin"] != "recovered_ms2"}
        recovered_mask = keep - job_mask
        frames = []
        hyps = []
        if job_mask:
            frames.append(_filter_rows(feature_dir / f"lib{lib}_{ion}_feature_rows.csv", job_mask))
            hyps.append(_filter_rows(feature_dir / f"lib{lib}_{ion}_neutral_mass_hypotheses.csv", job_mask))
        if recovered_mask:
            frames.append(_filter_rows(recovered_dir / f"lib{lib}_{ion}_feature_rows.csv", recovered_mask))
            hyps.append(_filter_rows(recovered_dir / f"lib{lib}_{ion}_neutral_mass_hypotheses.csv", recovered_mask))
        features = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
        hypotheses = pd.concat(hyps, ignore_index=True) if hyps else pd.DataFrame()
        manifest_path = feature_dir / f"lib{lib}_{ion}_gap_manifest.csv"
        if manifest_path.exists():
            manifest = _filter_rows(manifest_path, job_mask)
        else:
            manifest = pd.DataFrame(columns=["control_id"])
        present = set(manifest["control_id"].astype(str)) if not manifest.empty else set()
        missing = [{"control_id": cid, "conjugate_hints": ""} for cid in sorted(keep - present)]
        if missing:
            manifest = pd.concat([manifest, pd.DataFrame(missing)], ignore_index=True)
        features.to_csv(out_dir / f"lib{lib}_{ion}_feature_rows.csv", index=False)
        hypotheses.to_csv(out_dir / f"lib{lib}_{ion}_neutral_mass_hypotheses.csv", index=False)
        manifest.to_csv(out_dir / f"lib{lib}_{ion}_gap_manifest.csv", index=False)
        used.append(str(lib))
    return used


def scored_ids(root: Path, libs: tuple[int, ...] = LIBS) -> set[str]:
    found: set[str] = set()
    for lib in libs:
        for name in ("unified_candidates_scored.csv", "non_s0_candidates_scored.csv"):
            found |= _read_ids(root / f"rerank_lib{lib}" / name)
    return found


def load_input_ids(feature_dir: Path, recovered_dir: Path) -> tuple[dict[int, set[str]], ...]:
    """Return job features, hypotheses, manifest, then recovered features and hypotheses."""

    def per_lib(directory: Path, suffix: str) -> dict[int, set[str]]:
        return {lib: _read_ids(directory / f"lib{lib}_{ION[lib]}_{suffix}.csv") for lib in LIBS}

    return (
        per_lib(feature_dir, "feature_rows"),
        per_lib(feature_dir, "neutral_mass_hypotheses"),
        per_lib(feature_dir, "gap_manifest"),
        per_lib(recovered_dir, "feature_rows"),
        per_lib(recovered_dir, "neutral_mass_hypotheses"),
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog", type=Path, required=True)
    parser.add_argument("--scored-root", type=Path, required=True)
    parser.add_argument("--feature-dir", type=Path, required=True)
    parser.add_argument("--recovered-dir", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args(argv)

    catalog = pd.read_csv(args.catalog, low_memory=False)
    inputs = load_input_ids(args.feature_dir, args.recovered_dir)
    assignments, counts = assign_recovery(catalog, scored_ids(args.scored_root), *inputs)
    used = write_inputs(assignments, args.feature_dir, args.recovered_dir, args.out_dir)
    payload = {"counts": counts, "run_libs": ",".join(used), "assignments": assignments}
    args.out_dir.mkdir(parents=True, exist_ok=True)
    (args.out_dir / "selection.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(
        "[recover] "
        f"selected={counts['selected']} manifest={counts['job_manifest']} "
        f"features_only={counts['job_features']} recovered={counts['recovered_ms2']} "
        f"no_inputs={counts['no_inputs']} skipped_level={counts['skipped_level']} "
        f"libs={payload['run_libs'] or 'none'}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
