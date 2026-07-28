"""Second-stage tentative naming layers for residual de novo structures."""
from __future__ import annotations

import csv
import json
import sqlite3
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

try:
    from structural_elucidation.reports.tentative_naming import (
        HAS_RDKIT,
        TentativeName,
        _clean,
        _format_label,
        _num,
        _sentence_case_label,
    )
except ImportError:  # pragma: no cover - direct module import in tests
    from tentative_naming import (  # type: ignore[no-redef]
        HAS_RDKIT,
        TentativeName,
        _clean,
        _format_label,
        _num,
        _sentence_case_label,
    )

try:
    from rdkit import Chem

    _RDKitChem = Chem
except ImportError:  # pragma: no cover
    _RDKitChem = None

ANALOG_THRESHOLD = 0.50
ANALOG_HIGH_THRESHOLD = 0.85
IUPAC_TITLE_MAX_LEN = 90
_GENERIC_CLASS_NAMES = {
    "organic compounds",
    "organoheterocyclic compounds",
    "organic oxygen compounds",
    "organic nitrogen compounds",
    "organosulfur compounds",
    "lipids and lipid-like molecules",
}


@dataclass(frozen=True)
class ExtendedNamingContext:
    """Optional resources for layers 6-9."""

    external_db_path: Path | None = None
    classyfire_db_path: Path | None = None
    iupac_cache_path: Path | None = None
    chemont_model_dir: Path | None = None
    lib_id: int | None = None


def default_extended_context(
    *,
    reports_root: Path | None = None,
    lib_id: int | None = 400,
) -> ExtendedNamingContext | None:
    """Discover second-stage naming resources from standard workspace locations."""

    project_root = Path(__file__).resolve().parents[2]
    workspace_root = project_root.parents[1]
    reports_root = reports_root or (workspace_root / "reports 4")
    external_db = project_root / "reference_data" / "external_name_index.sqlite"
    classyfire_db = project_root / "reference_data" / "classyfire_hmdb.sqlite"
    iupac_cache = reports_root / "iupac_cache.jsonl"
    chemont_model_dir = workspace_root / "references" / "production"
    if not any(path.is_file() for path in (external_db, classyfire_db, iupac_cache)):
        return None
    return ExtendedNamingContext(
        external_db_path=external_db if external_db.is_file() else None,
        classyfire_db_path=classyfire_db if classyfire_db.is_file() else None,
        iupac_cache_path=iupac_cache if iupac_cache.is_file() else None,
        chemont_model_dir=chemont_model_dir if chemont_model_dir.is_dir() else None,
        lib_id=lib_id,
    )


def inchikey_from_smiles(smiles: str) -> str:
    if not HAS_RDKIT or not _RDKitChem:
        return ""
    mol = _RDKitChem.MolFromSmiles(smiles)
    if mol is None:
        return ""
    return _RDKitChem.MolToInchiKey(mol)


def connectivity_block(inchikey: str) -> str:
    return _clean(inchikey).split("-")[0]


def _label_with_formula(base: str, formula: str, *, tentative: bool = True) -> str:
    text = _format_label(base)
    formula_text = _clean(formula)
    if formula_text and formula_text not in text:
        text = f"{text} ({formula_text})"
    return _format_label(text)


def truncate_for_title(label: str, *, max_len: int = IUPAC_TITLE_MAX_LEN) -> str:
    text = _format_label(label)
    if len(text) <= max_len:
        return text
    return text[: max_len - 1].rstrip() + "…"


def load_iupac_cache(path: Path | str | None) -> dict[str, dict[str, Any]]:
    cache_path = Path(path) if path else None
    if cache_path is None or not cache_path.is_file():
        return {}
    out: dict[str, dict[str, Any]] = {}
    with cache_path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            try:
                row = json.loads(line)
            except json.JSONDecodeError:
                continue
            smiles = _clean(row.get("smiles"))
            if smiles:
                out[smiles] = row
    return out


def _lookup_external_db(
    connection: sqlite3.Connection,
    inchikey: str,
) -> tuple[str, str] | None:
    full = _clean(inchikey)
    skeleton = connectivity_block(full)
    if full:
        row = connection.execute(
            "SELECT canonical_name FROM names "
            "WHERE inchikey=? AND match_kind='full' "
            "ORDER BY priority LIMIT 1",
            (full,),
        ).fetchone()
        if row:
            return str(row[0]), "external_db_exact"
    if skeleton:
        row = connection.execute(
            "SELECT canonical_name FROM names "
            "WHERE inchikey=? AND match_kind='skeleton' "
            "ORDER BY priority LIMIT 1",
            (skeleton,),
        ).fetchone()
        if row:
            return str(row[0]), "external_db_skeleton"
    return None


def resolve_layer6_external_db(
    inchikey: str,
    formula: str,
    index_conn: sqlite3.Connection | None,
) -> TentativeName | None:
    if index_conn is None:
        return None
    hit = _lookup_external_db(index_conn, inchikey)
    if not hit:
        return None
    name, kind = hit
    if kind == "external_db_exact":
        label = _label_with_formula(_sentence_case_label(name), formula, tentative=False)
        confidence = "high"
    else:
        label = _label_with_formula(
            f"{_sentence_case_label(name)} (stereo unspecified)",
            formula,
            tentative=True,
        )
        confidence = "medium"
    return TentativeName(
        kind=kind,
        label=label,
        parent_name=name,
        parent_inchikey=_clean(inchikey) or None,
        confidence=confidence,
        source="external_db",
    )


def _read_neighbors_csv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    try:
        with path.open(newline="", encoding="utf-8") as handle:
            return list(csv.DictReader(handle))
    except (OSError, csv.Error):
        return []


def resolve_layer7_analog(neighbors_csv_path: Path | str | None) -> TentativeName | None:
    path = Path(neighbors_csv_path) if neighbors_csv_path else None
    if path is None:
        return None
    rows = _read_neighbors_csv(path)
    if not rows:
        return None

    def sort_key(row: dict[str, str]) -> float:
        return _num(row.get("similarity_to_target")) or -1.0

    best = max(rows, key=sort_key)
    tanimoto = _num(best.get("similarity_to_target"))
    neighbor_name = _clean(best.get("canonical_name"))
    if tanimoto is None or tanimoto < ANALOG_THRESHOLD or not neighbor_name:
        return None
    neighbor_ik = _clean(best.get("inchikeys")).split("|")[0].split(",")[0].strip()
    if tanimoto >= ANALOG_HIGH_THRESHOLD:
        base = f"≈ {neighbor_name} (Tanimoto {tanimoto:.2f})"
        confidence = "medium"
    else:
        base = f"Putative analog of {neighbor_name} (Tanimoto {tanimoto:.2f})"
        confidence = "low"
    return TentativeName(
        kind="analog",
        label=_label_with_formula(base, "", tentative=True),
        parent_name=neighbor_name,
        parent_inchikey=neighbor_ik or None,
        confidence=confidence,
        source="similarity_neighbors",
        analog_neighbor_name=neighbor_name,
        analog_neighbor_inchikey=neighbor_ik or None,
        analog_neighbor_tanimoto=tanimoto,
    )


def _lookup_classyfire(
    connection: sqlite3.Connection,
    inchikey: str,
) -> tuple[str, str, str] | None:
    full = _clean(inchikey)
    skeleton = connectivity_block(full)
    for key in (full, skeleton):
        if not key:
            continue
        row = connection.execute(
            "SELECT superclass, class_name, subclass FROM classyfire "
            "WHERE inchikey=? OR inchikey_skeleton=? LIMIT 1",
            (key, key),
        ).fetchone()
        if row:
            return str(row[0] or ""), str(row[1] or ""), str(row[2] or "")
    return None


def _class_label(superclass: str, class_name: str, subclass: str) -> str | None:
    candidates = [subclass, class_name, superclass]
    for candidate in candidates:
        text = _clean(candidate)
        if text and text.lower() not in _GENERIC_CLASS_NAMES:
            return text
    for candidate in candidates:
        text = _clean(candidate)
        if text:
            return text
    return None


def _predict_chemont_ms2(
    *,
    bundle_dir: Path,
    lib_id: int | None,
    chemont_model_dir: Path | None,
) -> tuple[str, str, str] | None:
  # MS2 fallback is optional and only used when raw spectrum data is available.
    confidence_path = bundle_dir / "data" / f"{bundle_dir.name.lower()}_confidence_inputs.json"
    if not confidence_path.is_file():
        return None
    try:
        fields = json.loads(confidence_path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    if not fields.get("ms2_spectrum_present"):
        return None
    if lib_id is None or chemont_model_dir is None:
        return None
    spectrum_path = bundle_dir / "data" / f"{bundle_dir.name.lower()}_ms2_spectrum.json"
    if not spectrum_path.is_file():
        return None
    try:
        spectrum_payload = json.loads(spectrum_path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    peaks = spectrum_payload.get("peaks") or spectrum_payload.get("spectrum")
    if not isinstance(peaks, list) or not peaks:
        return None
    parent_mass = _num(fields.get("observed_neutral_mass") or fields.get("observed_mz"))
    if parent_mass is None:
        return None
    try:
        import sys

        vendor_root = Path(__file__).resolve().parents[2] / "vendor" / "nce"
        if str(vendor_root) not in sys.path:
            sys.path.insert(0, str(vendor_root))
        from elucidation.chemontPredictor import ChemontPredictor  # type: ignore
    except Exception:  # pragma: no cover - optional heavy dependency chain
        return None
    try:
        predictor = ChemontPredictor.get_instance(int(lib_id))
        spectrum = {
            float(peak.get("mz") if isinstance(peak, dict) else peak[0]): float(
                peak.get("intensity") if isinstance(peak, dict) else peak[1]
            )
            for peak in peaks
            if peak
        }
        frame = predictor.predict_for_spectrum(parent_mass, spectrum)
    except Exception:  # pragma: no cover
        return None
    if frame is None or frame.empty:
        return None
    top = frame.iloc[0]
    name = _clean(top.get("chemont_name"))
    if not name:
        return None
    return "", name, name


def resolve_layer8_class(
    inchikey: str,
    formula: str,
    classyfire_conn: sqlite3.Connection | None,
    *,
    bundle_dir: Path | None = None,
    ctx: ExtendedNamingContext | None = None,
) -> TentativeName | None:
    superclass = class_name = subclass = ""
    source = "classyfire_hmdb"
    if classyfire_conn is not None:
        hit = _lookup_classyfire(classyfire_conn, inchikey)
        if hit:
            superclass, class_name, subclass = hit
    if not any((superclass, class_name, subclass)) and bundle_dir is not None and ctx is not None:
        ms2_hit = _predict_chemont_ms2(
            bundle_dir=bundle_dir,
            lib_id=ctx.lib_id,
            chemont_model_dir=ctx.chemont_model_dir,
        )
        if ms2_hit:
            superclass, class_name, subclass = ms2_hit
            source = "chemont_ms2"
    label_class = _class_label(superclass, class_name, subclass)
    if not label_class:
        return None
    label = _format_label(f"{label_class} — chemical class (tentative)")
    if formula:
        label = _format_label(f"{label} ({formula})")
    confidence = "medium" if source == "classyfire_hmdb" else "low"
    return TentativeName(
        kind="chemical_class",
        label=label,
        confidence=confidence,
        source=source,
        class_superclass=superclass or None,
        class_class=class_name or None,
        class_subclass=subclass or None,
    )


def resolve_layer9_iupac(
    smiles: str,
    formula: str,
    iupac_cache: Mapping[str, Mapping[str, Any]] | None,
) -> TentativeName | None:
    if not iupac_cache:
        return None
    row = iupac_cache.get(_clean(smiles))
    if not row:
        return None
    if not bool(row.get("round_trip_ok")):
        return None
    iupac_name = _clean(row.get("iupac_name"))
    if not iupac_name:
        return None
    base = f"{iupac_name} — IUPAC (generated)"
    label = _label_with_formula(base, formula, tentative=False)
    return TentativeName(
        kind="iupac_generated",
        label=label,
        confidence="medium",
        source="stout",
        iupac_name=iupac_name,
        iupac_round_trip_ok=True,
    )


def apply_extended_layers(
    smiles: str,
    formula: str,
    inchikey: str,
    bundle_dir: Path,
    ctx: ExtendedNamingContext | None,
) -> TentativeName | None:
    """Evaluate layers 6-9 in priority order."""

    if ctx is None:
        return None
    smiles_text = _clean(smiles)
    formula_text = _clean(formula)
    inchikey_text = _clean(inchikey) or inchikey_from_smiles(smiles_text)

    external_conn = None
    classyfire_conn = None
    try:
        if ctx.external_db_path and ctx.external_db_path.is_file():
            external_conn = sqlite3.connect(str(ctx.external_db_path))
            hit = resolve_layer6_external_db(inchikey_text, formula_text, external_conn)
            if hit:
                return hit

        neighbors_csv = bundle_dir / "data" / "similarity_neighbors.csv"
        hit = resolve_layer7_analog(neighbors_csv)
        if hit:
            return hit

        if ctx.classyfire_db_path and ctx.classyfire_db_path.is_file():
            classyfire_conn = sqlite3.connect(str(ctx.classyfire_db_path))
            hit = resolve_layer8_class(
                inchikey_text,
                formula_text,
                classyfire_conn,
                bundle_dir=bundle_dir,
                ctx=ctx,
            )
            if hit:
                return hit

        cache = load_iupac_cache(ctx.iupac_cache_path)
        return resolve_layer9_iupac(smiles_text, formula_text, cache)
    finally:
        if external_conn is not None:
            external_conn.close()
        if classyfire_conn is not None:
            classyfire_conn.close()
