#!/usr/bin/env python3
"""Batch-generate IUPAC names with round-trip verification.

Uses STOUT when installed; otherwise falls back to PubChem PUG REST lookups.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from structural_elucidation.reports.tentative_naming_ext import connectivity_block  # noqa: E402

try:
    from STOUT import translate_forward, translate_reverse

    HAS_STOUT = True
except ImportError:  # pragma: no cover
    HAS_STOUT = False


def _load_cache(path: Path) -> dict[str, dict]:
    if not path.is_file():
        return {}
    out: dict[str, dict] = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            try:
                row = json.loads(line)
            except json.JSONDecodeError:
                continue
            smiles = str(row.get("smiles") or "").strip()
            if smiles:
                out[smiles] = row
    return out


def _write_cache(path: Path, rows: dict[str, dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for smiles in sorted(rows):
            payload = dict(rows[smiles])
            payload["smiles"] = smiles
            handle.write(json.dumps(payload, ensure_ascii=False) + "\n")


def _pubchem_request(url: str) -> dict | None:
    request = urllib.request.Request(url, headers={"User-Agent": "structural-elucidation/1.0"})
    try:
        with urllib.request.urlopen(request, timeout=30) as response:
            payload = json.loads(response.read().decode("utf-8"))
    except (urllib.error.URLError, TimeoutError, json.JSONDecodeError, ValueError):
        return None
    return payload if isinstance(payload, dict) else None


def _pubchem_properties(smiles: str) -> dict[str, str]:
    encoded = urllib.parse.quote(smiles, safe="")
    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/"
        f"{encoded}/property/IUPACName,InChIKey,CanonicalSMILES/JSON"
    )
    payload = _pubchem_request(url)
    if not payload:
        return {}
    props = payload.get("PropertyTable", {}).get("Properties", [])
    if not props:
        return {}
    row = props[0]
    return {
        "iupac_name": str(row.get("IUPACName") or "").strip(),
        "inchikey": str(row.get("InChIKey") or "").strip(),
        "canonical_smiles": str(row.get("CanonicalSMILES") or "").strip(),
    }


def _pubchem_round_trip_ok(smiles: str, iupac_name: str, inchikey: str) -> bool:
    if not iupac_name or not inchikey:
        return False
    encoded = urllib.parse.quote(iupac_name, safe="")
    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
        f"{encoded}/property/InChIKey/JSON"
    )
    payload = _pubchem_request(url)
    if not payload:
        return False
    props = payload.get("PropertyTable", {}).get("Properties", [])
    if not props:
        return False
    back_key = str(props[0].get("InChIKey") or "").strip()
    if not back_key:
        return False
    return connectivity_block(back_key) == connectivity_block(inchikey)


def _generate_stout(smiles: str) -> dict:
    iupac_name = str(translate_forward(smiles) or "").strip()
    ok = False
    if iupac_name:
        try:
            back_smiles = translate_reverse(iupac_name)
            from structural_elucidation.reports.tentative_naming_ext import inchikey_from_smiles

            ok = connectivity_block(inchikey_from_smiles(back_smiles)) == connectivity_block(
                inchikey_from_smiles(smiles)
            )
        except Exception:
            ok = False
    return {
        "smiles": smiles,
        "iupac_name": iupac_name,
        "round_trip_ok": ok,
        "source_version": "stout-pypi",
    }


def _generate_pubchem(smiles: str) -> dict:
    props = _pubchem_properties(smiles)
    iupac_name = props.get("iupac_name", "")
    inchikey = props.get("inchikey", "")
    ok = _pubchem_round_trip_ok(smiles, iupac_name, inchikey) if iupac_name else False
    time.sleep(0.25)
    return {
        "smiles": smiles,
        "iupac_name": iupac_name,
        "round_trip_ok": ok,
        "source_version": "pubchem-pug-rest",
        "inchikey": inchikey,
    }


def _generate(smiles: str) -> dict:
    if HAS_STOUT:
        return _generate_stout(smiles)
    return _generate_pubchem(smiles)


def collect_residual_smiles(reports_root: Path) -> list[str]:
    smiles: list[str] = []
    seen: set[str] = set()
    for report_dir in sorted(path for path in reports_root.iterdir() if path.is_dir()):
        slug = report_dir.name.lower()
        confidence_path = report_dir / "data" / f"{slug}_confidence_inputs.json"
        if not confidence_path.is_file():
            continue
        try:
            fields = json.loads(confidence_path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            continue
        if fields.get("name_resolved") is not False:
            continue
        display_name = str(fields.get("display_name") or "")
        if not display_name.startswith("De novo structure"):
            continue
        candidate = str(fields.get("smiles") or "").strip()
        if not candidate or candidate in seen:
            continue
        seen.add(candidate)
        smiles.append(candidate)
    return smiles


def build_cache(reports_root: Path, output_path: Path, *, limit: int | None = None) -> dict[str, int]:
    rows = _load_cache(output_path)
    smiles_list = collect_residual_smiles(reports_root)
    if limit is not None:
        smiles_list = smiles_list[:limit]
    generated = 0
    passed = 0
    for smiles in smiles_list:
        if smiles in rows:
            continue
        rows[smiles] = _generate(smiles)
        generated += 1
        if rows[smiles].get("round_trip_ok"):
            passed += 1
        if generated % 25 == 0:
            _write_cache(output_path, rows)
    _write_cache(output_path, rows)
    total_ok = sum(1 for row in rows.values() if row.get("round_trip_ok"))
    return {
        "backend": "stout" if HAS_STOUT else "pubchem",
        "unique_smiles": len(smiles_list),
        "cache_rows": len(rows),
        "generated_this_run": generated,
        "round_trip_ok_total": total_ok,
    }


def smoke_test() -> None:
    samples = {
        "aspirin": "CC(=O)Oc1ccccc1C(=O)O",
        "caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "glucose": "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O",
        "quercetin": "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12",
        "cholesterol": "CC(C)CCCC(C)C1CCC2C1CCC3C2CCC4C3(C)CCC(O)C4(C)C",
    }
    passed = 0
    for label, smiles in samples.items():
        row = _generate(smiles)
        ok = bool(row.get("round_trip_ok"))
        passed += int(ok)
        print(f"[smoke] {label}: ok={ok} name={row.get('iupac_name', '')[:80]}")
    print(f"[smoke] round_trip_pass_rate={passed}/{len(samples)}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reports-root", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--smoke-test", action="store_true")
    args = parser.parse_args()
    if args.smoke_test:
        smoke_test()
        return
    if not args.reports_root or not args.output:
        raise SystemExit("--reports-root and --output are required unless --smoke-test is set")
    stats = build_cache(args.reports_root.resolve(), args.output.resolve(), limit=args.limit)
    print(f"[build_iupac_cache] wrote {args.output}")
    for key, value in stats.items():
        print(f"  {key}: {value}")


if __name__ == "__main__":
    main()
