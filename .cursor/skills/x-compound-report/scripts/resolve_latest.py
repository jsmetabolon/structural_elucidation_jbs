#!/usr/bin/env python3
"""Inventory structural-elucidation result sources and pick the latest of each kind.

Three kinds of source can feed an X-compound report:

- ``predictions``: per-library rerank CSVs, ``<root>/rerank_lib{209,305,400,402}/<name>.csv``.
- ``scores``: canonical per-compound scoring JSONs (``schema_version`` 1, written by
  ``structural_elucidation.pipeline.run_novel`` / ``reports.novel_scoring``).
- ``reports``: an already-rendered report bundle (``<root>/X-*/manifest.json``).

Nothing is modified. The output is a JSON inventory with the chosen source per
kind, every candidate it was chosen from, and warnings for stale or conflicting
results (calibration checksum drift, same compound with a different rank-1
structure in two prediction roots, the June 3,897-compound bundle).

Usage:
    python resolve_latest.py --search output models /data/deliveries \
        [--compounds ids.txt] [--json-out latest.json]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

REPO_ROOT = Path(__file__).resolve().parents[4]
CALIBRATION = REPO_ROOT / "tools" / "identity_confidence_calibration.json"
LIBS = ("209", "305", "400", "402")
# The 10 ppm refill is the delivered gate since the October naming-v2 release;
# unrefilled CSVs are the June bundle lineage.
PREDICTION_NAMES = (
    "unified_rerank_all_top10_ppm10.csv",
    "unified_rerank_all_top10_with_pathways.csv",
    "unified_rerank_all_top10.csv",
)
# Hashes of the June 26 report-inputs handoff (branch handoff/se-report-inputs).
JUNE_BUNDLE_SIZES = {"209": 42609879, "305": 3260580, "400": 14690763, "402": 11935632}
ID_RE = re.compile(r"^X-\d+$")
SKIP_DIRS = {".git", ".venv", "node_modules", "__pycache__", "assets", "vendor"}


def sha256(path: Path, limit: int | None = None) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        remaining = limit
        while True:
            chunk = handle.read(1 << 20 if remaining is None else min(1 << 20, remaining))
            if not chunk:
                break
            digest.update(chunk)
            if remaining is not None:
                remaining -= len(chunk)
                if remaining <= 0:
                    break
    return digest.hexdigest()


def iso(ts: float) -> str:
    return datetime.fromtimestamp(ts, tz=timezone.utc).isoformat(timespec="seconds")


def git_head() -> str | None:
    try:
        return subprocess.run(["git", "-C", str(REPO_ROOT), "rev-parse", "HEAD"], capture_output=True,
                              text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def is_ancestor(commit: str, head: str) -> bool | None:
    try:
        proc = subprocess.run(["git", "-C", str(REPO_ROOT), "merge-base", "--is-ancestor", commit, head],
                              capture_output=True)
    except OSError:
        return None
    return {0: True, 1: False}.get(proc.returncode)


def walk(roots: Iterable[Path], max_depth: int = 6) -> Iterable[Path]:
    """Directories under each root, skipping heavy or irrelevant trees."""
    for root in roots:
        if not root.is_dir():
            continue
        stack = [(root, 0)]
        while stack:
            directory, depth = stack.pop()
            yield directory
            if depth >= max_depth:
                continue
            try:
                children = sorted(p for p in directory.iterdir() if p.is_dir())
            except OSError:
                continue
            for child in children:
                if child.name in SKIP_DIRS or ID_RE.match(child.name):
                    continue
                stack.append((child, depth + 1))


def read_ids(path: Path) -> set[str]:
    import csv

    with path.open(newline="", encoding="utf-8", errors="replace") as handle:
        reader = csv.reader(handle)
        header = next(reader, [])
        column = next((i for i, name in enumerate(header) if name in ("control_id", "x_compound_id")), None)
        if column is None:
            return set()
        return {row[column] for row in reader if len(row) > column and row[column]}


def rank1_smiles(path: Path, wanted: set[str]) -> dict[str, str]:
    import csv

    out: dict[str, tuple[float, str]] = {}
    with path.open(newline="", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle)
        id_col = "control_id" if "control_id" in (reader.fieldnames or []) else "x_compound_id"
        for row in reader:
            cid = row.get(id_col, "")
            if wanted and cid not in wanted:
                continue
            try:
                rank = float(row.get("rank") or "inf")
            except ValueError:
                rank = float("inf")
            smiles = row.get("predicted_smiles") or row.get("smiles") or ""
            if cid not in out or rank < out[cid][0]:
                out[cid] = (rank, smiles)
    return {cid: smiles for cid, (_, smiles) in out.items()}


def find_predictions(roots: list[Path]) -> list[dict[str, Any]]:
    found = []
    for directory in walk(roots):
        libs = {lib: directory / f"rerank_lib{lib}" for lib in LIBS}
        if not any(path.is_dir() for path in libs.values()):
            continue
        for name in PREDICTION_NAMES:
            files = {lib: path / name for lib, path in libs.items() if (path / name).is_file()}
            if not files:
                continue
            ids: set[str] = set()
            for path in files.values():
                ids |= read_ids(path)
            june = all(JUNE_BUNDLE_SIZES.get(lib) == path.stat().st_size for lib, path in files.items())
            found.append({
                "kind": "predictions", "root": str(directory), "predictions_name": name,
                "libs": sorted(files), "n_compounds": len(ids),
                "newest_mtime": iso(max(p.stat().st_mtime for p in files.values())),
                "ppm10_gated": name.endswith("_ppm10.csv"),
                "june_2026_bundle": june,
                "files": {lib: {"path": str(p), "bytes": p.stat().st_size, "sha256_first_mb": sha256(p, 1 << 20)}
                          for lib, p in files.items()},
                "_ids": ids,
            })
    return found


def scoring_payload_summary(path: Path) -> dict[str, Any] | None:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    if not isinstance(payload, dict) or payload.get("schema_version") != 1 or "control_id" not in payload:
        return None
    return {"control_id": payload["control_id"], "status": payload.get("status"),
            "calibration_checksum": payload.get("calibration_checksum")}


def find_scores(roots: list[Path], current_calibration: str | None, head: str | None) -> list[dict[str, Any]]:
    found = []
    for directory in walk(roots):
        files = sorted(directory.glob("X-*.json"))
        if not files:
            continue
        payloads = [s for s in (scoring_payload_summary(p) for p in files) if s]
        if not payloads:
            continue
        checksums = sorted({p["calibration_checksum"] for p in payloads if p["calibration_checksum"]})
        statuses: dict[str, int] = {}
        for p in payloads:
            statuses[str(p["status"])] = statuses.get(str(p["status"]), 0) + 1
        ids = {p["control_id"] for p in payloads}
        synthetic = sum(1 for cid in ids if int(cid.split("-")[1]) >= 900000000)
        entry: dict[str, Any] = {
            "kind": "scores", "root": str(directory), "n_payloads": len(payloads), "statuses": statuses,
            "n_synthetic_ids": synthetic, "calibration_checksums": checksums,
            "calibration_matches_repo": bool(checksums) and checksums == [current_calibration],
            "newest_mtime": iso(max(p.stat().st_mtime for p in files)), "_ids": ids,
        }
        summary = directory.parent / "run_novel_summary.json"
        if summary.is_file():
            try:
                meta = json.loads(summary.read_text(encoding="utf-8"))
            except ValueError:
                meta = {}
            commit = meta.get("companion_commit")
            entry["run_novel_summary"] = {
                "path": str(summary), "commit": commit, "dirty": meta.get("companion_dirty"),
                "calibration_sha256": meta.get("calibration_checksum"),
                "commit_in_current_history": is_ancestor(commit, head) if commit and head else None,
            }
        found.append(entry)
    return found


def find_reports(roots: list[Path]) -> list[dict[str, Any]]:
    found = []
    for directory in walk(roots):
        manifests = list(directory.glob("X-*/manifest.json"))
        if not manifests:
            continue
        ids, versions = set(), {}
        newest = 0.0
        for path in manifests:
            ids.add(path.parent.name)
            newest = max(newest, path.stat().st_mtime)
            try:
                version = json.loads(path.read_text(encoding="utf-8")).get("naming_policy_version") or "v1"
            except ValueError:
                version = "unreadable"
            versions[version] = versions.get(version, 0) + 1
        found.append({"kind": "reports", "root": str(directory), "n_compounds": len(ids),
                      "naming_policy_versions": versions, "newest_mtime": iso(newest), "_ids": ids})
    return found


def choose(candidates: list[dict[str, Any]], key) -> dict[str, Any] | None:
    return max(candidates, key=key) if candidates else None


def resolve(search: list[Path], compounds: set[str] | None = None) -> dict[str, Any]:
    head = git_head()
    calibration = sha256(CALIBRATION) if CALIBRATION.is_file() else None
    predictions = find_predictions(search)
    scores = find_scores(search, calibration, head)
    reports = find_reports(search)
    warnings: list[str] = []

    chosen_predictions = choose(predictions, lambda c: (c["ppm10_gated"], not c["june_2026_bundle"], c["newest_mtime"]))
    chosen_scores = choose(scores, lambda c: (c["calibration_matches_repo"], c["n_payloads"] - c["n_synthetic_ids"] > 0,
                                              c["newest_mtime"]))
    chosen_reports = choose(reports, lambda c: (c["naming_policy_versions"].get("v2", 0) > 0, c["newest_mtime"]))

    for entry in predictions:
        if entry["june_2026_bundle"]:
            warnings.append(f"{entry['root']}: June 26 2026 report-inputs bundle (pre 10 ppm gate, pre naming v2, "
                            "includes Exclude/Artifact compounds); use only if no newer delivery is available")
    if chosen_predictions and not chosen_predictions["ppm10_gated"]:
        warnings.append("No 10 ppm-refilled predictions (unified_rerank_all_top10_ppm10.csv) found; the October "
                        "delivery used them. Ask for the merged delivery root before sending reports to a client")
    for entry in scores:
        if not entry["calibration_matches_repo"]:
            warnings.append(f"{entry['root']}: scoring JSON calibration checksum {entry['calibration_checksums']} "
                            f"does not match tools/identity_confidence_calibration.json ({calibration}); stale")
        if entry["n_synthetic_ids"] == entry["n_payloads"]:
            warnings.append(f"{entry['root']}: every scoring JSON uses a synthetic untargeted id (X-9xxxxxxxx); "
                            "they cannot back a registered X-compound list")
        rn = entry.get("run_novel_summary") or {}
        if rn.get("commit_in_current_history") is False:
            warnings.append(f"{entry['root']}: produced by commit {rn['commit']}, "
                            "which is not in this checkout's history")
        if rn.get("dirty"):
            warnings.append(f"{entry['root']}: produced from a dirty checkout")

    conflicts: list[dict[str, Any]] = []
    if compounds and len(predictions) > 1:
        per_root = {}
        for entry in predictions:
            smiles: dict[str, str] = {}
            for info in entry["files"].values():
                for cid, smi in rank1_smiles(Path(info["path"]), compounds).items():
                    smiles.setdefault(cid, smi)
            per_root[f"{entry['root']}::{entry['predictions_name']}"] = smiles
        for cid in sorted(compounds):
            values = {root: s[cid] for root, s in per_root.items() if cid in s}
            if len(set(values.values())) > 1:
                conflicts.append({"control_id": cid, "rank1_by_source": values})
        if conflicts:
            warnings.append(f"{len(conflicts)} requested compounds have a different rank-1 structure in different "
                            "prediction roots; the chosen root wins, confirm it is the intended delivery")

    def coverage(entry: dict[str, Any] | None) -> dict[str, Any] | None:
        if entry is None or not compounds:
            return None
        have = entry["_ids"] & compounds
        return {"requested": len(compounds), "covered": len(have), "missing": sorted(compounds - have)}

    def public(entry: dict[str, Any] | None) -> dict[str, Any] | None:
        return None if entry is None else {k: v for k, v in entry.items() if not k.startswith("_")}

    return {
        "generated_at": iso(datetime.now(tz=timezone.utc).timestamp()),
        "repo_head": head, "repo_calibration_sha256": calibration,
        "search_roots": [str(p) for p in search],
        "chosen": {
            "predictions": public(chosen_predictions), "predictions_coverage": coverage(chosen_predictions),
            "scores": public(chosen_scores), "scores_coverage": coverage(chosen_scores),
            "reports": public(chosen_reports), "reports_coverage": coverage(chosen_reports),
        },
        "candidates": [public(e) for e in (*predictions, *scores, *reports)],
        "rank1_conflicts": conflicts,
        "warnings": warnings,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--search", nargs="+", type=Path, default=[REPO_ROOT / "output", REPO_ROOT / "models"],
                        help="Directories to scan (default: <repo>/output <repo>/models)")
    parser.add_argument("--compounds", type=Path, default=None, help="Normalized id list (one X-id per line)")
    parser.add_argument("--json-out", type=Path, default=None)
    args = parser.parse_args(argv)
    compounds = None
    if args.compounds:
        compounds = {line.strip() for line in args.compounds.read_text().splitlines() if ID_RE.match(line.strip())}
    result = resolve([p.expanduser().resolve() for p in args.search], compounds)
    text = json.dumps(result, indent=2)
    if args.json_out:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(text + "\n", encoding="utf-8")
    print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
