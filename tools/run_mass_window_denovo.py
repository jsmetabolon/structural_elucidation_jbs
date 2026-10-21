"""Decode unpublished compounds in mass bands and keep structures within 10 ppm.

The encoder's spectrum axis starts at m/z 50 with 0.5 Da bins, which is the
checkpoint width. Each batch raises only the upper edge to that batch's
corrected masses, so peaks stay in the bins the encoder was trained on and the
call is not the 50–950 sweep.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
from pathlib import Path

import pandas as pd

ION = {209: "negative", 305: "negative", 400: "positive", 402: "positive"}
PROTON = 1.007276466583
PAIR = {"negative": "neg_209_305", "positive": "pos_400_402"}


def assign_libraries(targets: pd.DataFrame, bundle: Path) -> pd.DataFrame:
    """Fill a library for compounds that have an encoder row but no local hypothesis."""
    out = targets.copy()
    missing = set(out.loc[out.lib_id.eq(""), "x_compound_id"])
    if not missing:
        return out
    found: dict[str, int] = {}
    for lib, ion in ION.items():
        table = bundle / f"lib_{lib}_{ion}" / "feature_table_prepared.csv"
        if not table.exists():
            continue
        frame = pd.read_csv(table, usecols=["compound_name"], low_memory=False)
        names = frame["compound_name"].astype(str)
        for cid in list(missing):
            if names.str.contains(cid, regex=False).any():
                found[cid] = lib
                missing.discard(cid)
    out.loc[out.x_compound_id.isin(found), "lib_id"] = out.loc[out.x_compound_id.isin(found), "x_compound_id"].map(found).astype(str)
    return out


def write_selection(frame: pd.DataFrame, path: Path) -> None:
    assignments = [
        {"control_id": row.x_compound_id, "lib": int(row.lib_id)}
        for row in frame.itertuples(index=False)
        if str(row.lib_id).strip()
    ]
    path.write_text(json.dumps({"assignments": assignments}), encoding="utf-8")


def band_feature_file(features: pd.DataFrame, width: float) -> list[pd.DataFrame]:
    """Split a feature table into groups whose corrected masses span at most ``width`` Da."""
    meta = features.groupby("feature_key", sort=False).agg(parent_mz=("parent_mz", "first"), target=("target_neutral_mass", "first"))
    meta = meta.sort_values("target")
    bands: list[list[str]] = []
    current: list[str] = []
    start = None
    for key, row in meta.iterrows():
        mass = float(row.target)
        if start is None or mass - start <= width:
            current.append(str(key))
            start = mass if start is None else start
        else:
            bands.append(current)
            current = [str(key)]
            start = mass
    if current:
        bands.append(current)
    return [features[features.feature_key.astype(str).isin(keys)].copy() for keys in bands]


def decode_band(
    *,
    features: Path,
    lib: int,
    ion: str,
    maxmass: float,
    tolerance_da: float,
    out_dir: Path,
    workspace: Path,
    python: str,
    encoder: Path,
    samples: int,
) -> None:
    pair = PAIR[ion]
    out_dir.mkdir(parents=True, exist_ok=True)
    command = [
        python,
        "generate_selfies_denovo.py",
        "--model_ckpt",
        f"models/selfies_decoder_{pair}/selfies_decoder_{pair}.pt",
        "--vocab_json",
        f"models/selfies_decoder_{pair}/selfies_vocab_{pair}.json",
        "--feature_table",
        str(features.resolve()),
        "--encoder_ckpt",
        str(encoder),
        "--encoder_kind",
        "merged",
        "--lib_id",
        str(lib),
        "--control_lib_id",
        str(lib),
        "--ion_mode",
        ion,
        "--output_dir",
        str(out_dir.resolve()),
        "--artifact_prefix",
        f"mass_{lib}",
        "--num_samples",
        str(samples),
        "--device",
        "cpu",
        "--minmass",
        "50",
        "--maxmass",
        f"{maxmass:.4f}",
        "--parent_mz_min",
        "50",
        "--parent_mz_max",
        "950",
        "--mass_tolerance_da",
        f"{tolerance_da:.6f}",
        "--disable_adduct_hypotheses",
        "--disable_ri_conditioning",
        "--skip_control_eval",
    ]
    env = os.environ.copy()
    env.pop("PYTHONPATH", None)
    subprocess.run(command, cwd=workspace, env=env, check=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--targets", type=Path, required=True)
    parser.add_argument("--bundle-dir", type=Path, required=True)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--encoder", type=Path, required=True)
    parser.add_argument("--python", required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--band-width", type=float, default=80.0)
    parser.add_argument("--samples", type=int, default=60)
    args = parser.parse_args()
    raw = pd.read_csv(args.targets, dtype=str, keep_default_na=False)
    wanted = raw[raw.action.isin(["decode", "build_conjugate"])].copy()
    # Conjugate rows that already passed are removed by the caller; failures stay marked build_conjugate.
    assigned = assign_libraries(wanted, args.bundle_dir)
    assigned = assigned[assigned.lib_id.astype(str).str.strip().ne("")]
    args.work.mkdir(parents=True, exist_ok=True)
    selection = args.work / "selection.json"
    ids = args.work / "ids.txt"
    write_selection(assigned, selection)
    ids.write_text("\n".join(assigned.x_compound_id) + "\n", encoding="utf-8")
    slice_cmd = [
        args.python,
        str(args.repo / "tools" / "prepare_denovo_ppm10_inputs.py"),
        "slice",
        "--bundle-dir",
        str(args.bundle_dir),
        "--ids",
        str(ids),
        "--selection",
        str(selection),
        "--out-dir",
        str(args.work / "features"),
    ]
    env = os.environ.copy()
    env["PYTHONPATH"] = str(args.repo)
    subprocess.run(slice_cmd, cwd=args.repo, env=env, check=True)
    mass_by_id = {}
    for row in assigned.itertuples(index=False):
        if str(row.proton_neutral).strip():
            mass_by_id[row.x_compound_id] = float(row.proton_neutral)
    manifest = []
    for feature_path in sorted((args.work / "features").glob("lib*_features.csv")):
        frame = pd.read_csv(feature_path, low_memory=False)
        name = feature_path.name
        lib = int(name.split("_")[0].replace("lib", ""))
        ion = "negative" if "negative" in name else "positive"
        cids = frame["compound_name"].astype(str).str.extract(r"(X-\d+)", expand=False)
        parent = pd.to_numeric(frame["parent_mz"], errors="coerce")
        proton = cids.map(mass_by_id)
        fallback = parent + (PROTON if ion == "negative" else -PROTON)
        frame["target_neutral_mass"] = proton.fillna(fallback)
        for index, band in enumerate(band_feature_file(frame, args.band_width)):
            target = pd.to_numeric(band["target_neutral_mass"], errors="coerce")
            parent_max = float(pd.to_numeric(band["parent_mz"], errors="coerce").max())
            ceiling = max(float(target.max()), parent_max)
            maxmass = min(950.0, ceiling + 30.0)
            tolerance = max(ceiling * 10.0 / 1e6, 0.002)
            band_path = args.work / "bands" / f"lib{lib}_{ion}_band{index:02d}.csv"
            band_path.parent.mkdir(parents=True, exist_ok=True)
            band.drop(columns=["target_neutral_mass"]).to_csv(band_path, index=False)
            out_dir = args.work / "decode" / f"lib_{lib}_{ion}" / f"band{index:02d}" / "selfies_denovo"
            n_compounds = int(band["compound_name"].astype(str).str.extract(r"(X-\d+)", expand=False).nunique())
            print(json.dumps({"lib": lib, "band": index, "compounds": n_compounds, "maxmass": round(maxmass, 3), "tolerance_da": round(tolerance, 5)}), flush=True)
            decode_band(
                features=band_path,
                lib=lib,
                ion=ion,
                maxmass=maxmass,
                tolerance_da=tolerance,
                out_dir=out_dir,
                workspace=args.workspace,
                python=args.python,
                encoder=args.encoder,
                samples=args.samples,
            )
            manifest.append({"lib": lib, "ion": ion, "band": index, "out": str(out_dir / "candidate_metrics_all.csv")})
    (args.work / "decode_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps({"bands": len(manifest)}))


if __name__ == "__main__":
    main()
