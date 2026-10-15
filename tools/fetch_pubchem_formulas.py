"""Fetch up to 25 PubChem structures for each formula, with an on-disk cache."""

from __future__ import annotations

import argparse
import csv
import io
import json
import time
import urllib.error
import urllib.request
from pathlib import Path

USER_AGENT = "structural-elucidation-formula-search"
PROPERTY_URL = (
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastformula/{formula}"
    "/property/CanonicalSMILES/CSV?MaxRecords={limit}"
)


def cache_path(cache: Path, formula: str) -> Path:
    return cache / f"{formula}.csv"


def fetch_formula(formula: str, limit: int) -> list[tuple[str, str]]:
    url = PROPERTY_URL.format(formula=formula, limit=limit)
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    try:
        with urllib.request.urlopen(request, timeout=60) as response:
            payload = response.read().decode("utf-8")
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            return []
        raise
    rows = list(csv.DictReader(io.StringIO(payload)))
    kept: list[tuple[str, str]] = []
    for row in rows:
        smiles = row.get("CanonicalSMILES") or row.get("ConnectivitySMILES") or row.get("SMILES") or ""
        cid = row.get("CID", "")
        if smiles:
            kept.append((cid, smiles))
        if len(kept) >= limit:
            break
    return kept


def read_cache(path: Path) -> list[tuple[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return [(row["cid"], row["smiles"]) for row in csv.DictReader(handle) if row.get("smiles")]


def write_cache(path: Path, rows: list[tuple[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["cid", "smiles"])
        writer.writeheader()
        for cid, smiles in rows:
            writer.writerow({"cid": cid, "smiles": smiles})


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--formulas", type=Path, required=True)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--limit", type=int, default=25)
    parser.add_argument("--pause", type=float, default=0.2)
    args = parser.parse_args()
    args.cache.mkdir(parents=True, exist_ok=True)
    frame_formulas = []
    with args.formulas.open(newline="", encoding="utf-8") as handle:
        seen: set[str] = set()
        for row in csv.DictReader(handle):
            formula = row["formula"]
            if formula not in seen:
                seen.add(formula)
                frame_formulas.append(formula)
    failures = 0
    for index, formula in enumerate(frame_formulas, start=1):
        path = cache_path(args.cache, formula)
        if path.exists():
            continue
        for attempt in range(5):
            try:
                rows = fetch_formula(formula, args.limit)
                write_cache(path, rows)
                break
            except Exception as exc:
                if attempt == 4:
                    failures += 1
                    path.with_suffix(".error").write_text(str(exc), encoding="utf-8")
                    break
                if "429" in str(exc):
                    time.sleep(20 * (attempt + 1))
                else:
                    time.sleep(2 ** attempt)
        time.sleep(args.pause)
        if index % 100 == 0:
            print(json.dumps({"fetched": index, "of": len(frame_formulas), "failures": failures}), flush=True)
    with args.out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["formula", "cid", "smiles"])
        writer.writeheader()
        structures = 0
        formulas_with_hits = 0
        for formula in frame_formulas:
            path = cache_path(args.cache, formula)
            if not path.exists():
                continue
            rows = read_cache(path)
            if rows:
                formulas_with_hits += 1
            for cid, smiles in rows:
                writer.writerow({"formula": formula, "cid": cid, "smiles": smiles})
                structures += 1
    print(json.dumps({
        "formulas": len(frame_formulas),
        "formulas_with_structures": formulas_with_hits,
        "structures": structures,
        "failures": failures,
    }))


if __name__ == "__main__":
    main()
