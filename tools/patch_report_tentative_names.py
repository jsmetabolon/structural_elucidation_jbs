#!/usr/bin/env python3
"""Reconcile display-only structural descriptors on existing report bundles.

Dry-run is the default and writes ``tentative_naming_report.csv``.  ``--apply``
updates confidence/evidence sidecars and report HTML without changing any
confidence, MSI, tier, or identity values.  The pass is idempotent: a second
run should change zero files when names are already current.
"""
from __future__ import annotations

import argparse
import csv
import html
import json
import re
import sys
from pathlib import Path
from typing import Any

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from structural_elucidation.reports import tentative_naming as TN  # noqa: E402
from structural_elucidation.reports.tentative_naming_ext import (  # noqa: E402
    ExtendedNamingContext,
    truncate_for_title,
)

csv.field_size_limit(10**8)

TITLE_RE = re.compile(r"<title>(.*?)</title>", re.I | re.S)
H1_RE = re.compile(r"<h1>(.*?)</h1>", re.I | re.S)
CAPTION_RE = re.compile(r"<p><strong>(.*?)</strong></p>", re.I | re.S)
REPORT_TITLE_SUFFIX_RE = re.compile(r"(\s+Structure Resolution \((?:Public|Internal)\))$", re.I)
BADGE_RE = re.compile(
    r"<!--TENTATIVE_DESCRIPTOR_BADGE-->.*?<!--/TENTATIVE_DESCRIPTOR_BADGE-->",
    re.S,
)
BADGE_ANCHOR = '<span class="badge ok">top-1 candidate</span>'
BADGE_TITLE = "Descriptor specificity only; not identification confidence"
ALT_RE = re.compile(r"""alt\s*=\s*['"]([^'"]+?)\s+structure['"]""", re.I)
_TENTATIVE_SUFFIX_RE = re.compile(r"\s*—\s*tentative\s*$", re.I)

SUMMARY_COLUMNS = [
    "control_id",
    "descriptor_kind",
    "old_title",
    "new_title",
    "parent_name",
    "parent_formula",
    "conjugate_count",
    "confidence",
]
PROTECTED_CONFIDENCE_KEYS = {
    "overall",
    "identity_confidence",
    "identity_status",
    "msi_level",
    "tier",
}


def _read_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}
    return value if isinstance(value, dict) else {}


def _html_alt_fallback(report_dir: Path) -> str:
    for name in ("public.html", "index.html"):
        path = report_dir / name
        if not path.is_file():
            continue
        try:
            text = path.read_text(encoding="utf-8")
        except OSError:
            continue
        match = ALT_RE.search(text)
        if match:
            return html.unescape(match.group(1).strip())
    return ""


def _read_rank1(path: Path) -> dict[str, str]:
    if not path.exists():
        return {}
    try:
        with path.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
    except (OSError, csv.Error):
        return {}

    def rank(row: dict[str, str]) -> float:
        try:
            return float(row.get("rank") or 1e12)
        except (TypeError, ValueError):
            return 1e12

    return min(rows, key=rank) if rows else {}


def _write_json_if_changed(path: Path, payload: dict[str, Any], *, apply: bool) -> bool:
    rendered = json.dumps(payload, indent=2)
    try:
        current = path.read_text(encoding="utf-8")
    except OSError:
        current = ""
    if current == rendered:
        return False
    if apply:
        path.write_text(rendered, encoding="utf-8")
    return True


def patch_html(
    text: str,
    control_id: str,
    new_title: str,
    descriptor: TN.TentativeName | None,
) -> str:
    """Replace report names and idempotently add/remove the descriptor badge."""

    display_title = truncate_for_title(new_title) if descriptor and descriptor.kind == "iupac_generated" else new_title
    escaped = html.escape(display_title)
    escaped_full = html.escape(new_title)
    report_title = (
        display_title
        if descriptor or (new_title and new_title != TN.UNNAMED_STRUCTURE_LABEL and new_title != TN.LEGACY_UNNAMED_PREFIX)
        else control_id
    )
    escaped_report_title = html.escape(report_title)

    def title_replacement(match: re.Match[str]) -> str:
        content = html.unescape(match.group(1)).strip()
        suffix_match = REPORT_TITLE_SUFFIX_RE.search(content)
        suffix = suffix_match.group(1) if suffix_match else ""
        return f"<title>{escaped_report_title}{html.escape(suffix)}</title>"

    patched = TITLE_RE.sub(title_replacement, text, count=1)
    patched = H1_RE.sub(lambda _match: f"<h1>{escaped_report_title}</h1>", patched, count=1)

    if descriptor:
        confidence = html.escape(descriptor.confidence.title())
        badge = (
            "<!--TENTATIVE_DESCRIPTOR_BADGE-->"
            f'<span class="badge" title="{html.escape(BADGE_TITLE)}">'
            f"Tentative descriptor \N{MIDDLE DOT} {confidence}</span>"
            "<!--/TENTATIVE_DESCRIPTOR_BADGE-->"
        )
        if BADGE_RE.search(patched):
            patched = BADGE_RE.sub(lambda _match: badge, patched)
        elif BADGE_ANCHOR in patched:
            patched = patched.replace(BADGE_ANCHOR, f"{BADGE_ANCHOR}\n          {badge}", 1)
    else:
        patched = BADGE_RE.sub("", patched)

    return CAPTION_RE.sub(
        lambda _match: (
            f"<p><strong title='{escaped_full}'>{escaped}</strong></p>"
            if descriptor and descriptor.kind == "iupac_generated" and display_title != new_title
            else f"<p><strong>{escaped}</strong></p>"
        ),
        patched,
        count=1,
    )


def _is_report_html(path: Path, text: str) -> bool:
    return path.name != "guides.html" and (
        "Structure-resolution report" in text or "Structure Resolution (" in text
    )


def _is_library_named_fields(fields: dict[str, Any]) -> bool:
    """Whether a report name was resolved out of the named-metabolite library."""

    source = str(fields.get("source") or "").strip().lower()
    return fields.get("name_resolved") is True or source == "exact"


def _is_known_x_named(report_dir: Path, fields: dict[str, Any]) -> bool:
    """Whether the displayed name came from the compound's own known_X row."""

    source = str(fields.get("descriptor_source") or "").strip().lower()
    return source == "known_x_row" or bool(TN.read_known_x_name(report_dir))


def _in_isomer_scope(report_dir: Path, fields: dict[str, Any]) -> bool:
    """Reports the ``Isomer of`` policy decides the name for, either way.

    A report already carrying the prefix stays in scope so that a name later
    found to be the compound's own known_X identity loses it again.
    """

    display_name = str(fields.get("display_name") or "").strip()
    return _is_library_named_fields(fields) or display_name.lower().startswith(
        TN.ISOMER_PREFIX.lower()
    )


def _patch_titles(text: str, report_title: str, caption: str) -> str:
    """Set the page title, heading, and caption without touching anything else."""

    escaped_title = html.escape(report_title)

    def title_replacement(match: re.Match[str]) -> str:
        content = html.unescape(match.group(1)).strip()
        suffix_match = REPORT_TITLE_SUFFIX_RE.search(content)
        suffix = suffix_match.group(1) if suffix_match else ""
        return f"<title>{escaped_title}{html.escape(suffix)}</title>"

    patched = TITLE_RE.sub(title_replacement, text, count=1)
    patched = H1_RE.sub(lambda _match: f"<h1>{escaped_title}</h1>", patched, count=1)
    return CAPTION_RE.sub(
        lambda _match: f"<p><strong>{html.escape(caption)}</strong></p>",
        patched,
        count=1,
    )


def _patch_named_report(
    report_dir: Path,
    fields: dict[str, Any],
    *,
    apply: bool,
) -> dict[str, Any]:
    """Apply the ``Isomer of`` policy to a report that already carries a name.

    The name is rebuilt from ``resolved_name`` the way the report build does it,
    so a compound whose name belongs to its own known_X row is asserted plainly
    while a name borrowed from the metabolite library is qualified.
    """

    data_dir = report_dir / "data"
    slug = report_dir.name.lower()
    confidence_path = data_dir / f"{slug}_confidence_inputs.json"
    evidence_path = data_dir / f"{slug}_evidence_summary.json"
    manifest_path = report_dir / "manifest.json"
    evidence = _read_json(evidence_path)
    target = evidence.get("target") if isinstance(evidence.get("target"), dict) else {}
    old_title = str(
        fields.get("display_name") or target.get("name") or report_dir.name
    ).strip()
    resolved_name = TN.bare_library_name(
        str(fields.get("resolved_name") or old_title).strip()
    )
    # ``sheet_feature_name`` carries the guards that keep the qualifier off
    # synthesized descriptors, so the report and the summary sheet agree.
    label = (
        resolved_name
        if _is_known_x_named(report_dir, fields)
        else TN.sheet_feature_name(
            resolved_name,
            name_resolved=bool(fields.get("name_resolved")),
            source=str(fields.get("source") or ""),
            descriptor_kind=str(fields.get("descriptor_kind") or ""),
        )
    )
    level = TN.msi_int(fields.get("msi_level"))
    new_title = TN.display_name_for_level(
        label,
        msi_level=fields.get("msi_level"),
        formula=str(target.get("formula") or "").strip(),
        observed_mz=fields.get("observed_mz"),
        name_is_reference=False,
    )
    # Matches the report renderer: levels 4-5 withhold the structure, and a
    # compound with no descriptor has no name to head the page with either.
    withholding = level is not None and level >= 4
    report_title = (
        label
        if not withholding and fields.get("descriptor_kind")
        else report_dir.name
    )
    summary = {
        "control_id": report_dir.name,
        "descriptor_kind": "library_isomer",
        "old_title": old_title,
        "new_title": new_title,
        "parent_name": resolved_name,
        "parent_formula": "",
        "conjugate_count": 0,
        "confidence": "",
    }

    updated_fields = dict(fields)
    updated_fields["display_name"] = new_title
    changed = _write_json_if_changed(confidence_path, updated_fields, apply=apply)
    if evidence_path.exists():
        updated_evidence = dict(evidence)
        updated_target = dict(target)
        updated_target["name"] = new_title
        updated_evidence["target"] = updated_target
        changed = _write_json_if_changed(evidence_path, updated_evidence, apply=apply) or changed
    if manifest_path.exists():
        manifest = _read_json(manifest_path)
        manifest_target = manifest.get("target") if isinstance(manifest.get("target"), dict) else {}
        updated_manifest = dict(manifest)
        updated_manifest_target = dict(manifest_target)
        updated_manifest_target["name"] = new_title
        updated_manifest["target"] = updated_manifest_target
        changed = _write_json_if_changed(manifest_path, updated_manifest, apply=apply) or changed

    html_changed = False
    for path in sorted(report_dir.glob("*.html")):
        try:
            text = path.read_text(encoding="utf-8")
        except OSError:
            continue
        if not _is_report_html(path, text):
            continue
        patched = _patch_titles(text, report_title, new_title)
        if patched != text:
            html_changed = True
            if apply:
                path.write_text(patched, encoding="utf-8")
    changed = changed or html_changed
    return {
        "status": "library_isomer",
        "changed": changed,
        "title_changed": old_title != new_title or html_changed,
        "summary": summary,
    }


def patch_report(
    report_dir: Path,
    db_names: dict[str, Any],
    *,
    apply: bool,
    extended_ctx: ExtendedNamingContext | None = None,
) -> dict[str, Any]:
    """Patch one report bundle and return a structured outcome."""

    data_dir = report_dir / "data"
    slug = report_dir.name.lower()
    confidence_path = data_dir / f"{slug}_confidence_inputs.json"
    evidence_path = data_dir / f"{slug}_evidence_summary.json"
    manifest_path = report_dir / "manifest.json"
    if not confidence_path.exists():
        return {"status": "no_fields", "changed": False, "summary": None}

    fields = _read_json(confidence_path)
    if _in_isomer_scope(report_dir, fields):
        return _patch_named_report(report_dir, fields, apply=apply)
    # Scope is deliberately strict: missing/null is not an unnamed report.
    if fields.get("name_resolved") is not False:
        return {"status": "named_or_unscoped", "changed": False, "summary": None}

    display_name = str(fields.get("display_name") or "").strip()
    stripped_display = _TENTATIVE_SUFFIX_RE.sub("", display_name).strip()

    evidence = _read_json(evidence_path)
    target = evidence.get("target") if isinstance(evidence.get("target"), dict) else {}
    rank1 = _read_rank1(data_dir / f"{slug}_top_candidates.csv")
    if not rank1.get("lipid_class") and stripped_display in TN.GENERIC_CLASS_LABELS:
        inferred = TN.infer_lipid_class_from_generic_label(stripped_display)
        if inferred:
            rank1 = dict(rank1)
            rank1["lipid_class"] = inferred
    smiles = str(target.get("smiles") or fields.get("smiles") or rank1.get("predicted_smiles") or "").strip()
    formula = str(target.get("formula") or rank1.get("predicted_formula") or "").strip()
    old_title = str(
        fields.get("display_name")
        or target.get("name")
        or (f"{TN.UNNAMED_STRUCTURE_LABEL} ({formula})" if formula else report_dir.name)
    ).strip()
    inchikey = str(target.get("inchikey") or fields.get("inchikey") or "").strip()
    msi_level = str(fields.get("msi_level") or "").strip()
    truth_match = bool(fields.get("truth_match"))
    html_alt_name = _html_alt_fallback(report_dir)

    tentative = TN.resolve(
        smiles,
        formula,
        rank1,
        db_names,
        name_resolved=False,
        bundle_dir=report_dir,
        extended_ctx=extended_ctx,
        inchikey=inchikey,
        control_id=report_dir.name,
        msi_level=msi_level,
        truth_match=truth_match,
        html_alt_name=html_alt_name,
    )
    raw_current = display_name or old_title
    current_name = TN.normalize_descriptor_label(raw_current)
    if tentative:
        new_title = TN.normalize_descriptor_label(tentative.label)
    elif current_name != raw_current:
        new_title = current_name
    else:
        new_title = TN.UNNAMED_STRUCTURE_LABEL

    current_rank = TN.specificity_rank(str(fields.get("descriptor_kind") or ""), current_name)
    proposed_rank = TN.specificity_rank(tentative.kind if tentative else "", new_title)

    summary = {
        "control_id": report_dir.name,
        "descriptor_kind": tentative.kind if tentative else "",
        "old_title": old_title,
        "new_title": new_title,
        "parent_name": (tentative.parent_name or "") if tentative else "",
        "parent_formula": (tentative.parent_formula or "") if tentative else "",
        "conjugate_count": tentative.conjugate_count if tentative else 0,
        "confidence": tentative.confidence if tentative else "",
    }

    suffix_only_cleanup = current_name == new_title and raw_current != new_title
    if proposed_rank < current_rank and not suffix_only_cleanup:
        return {
            "status": "preserved",
            "changed": False,
            "title_changed": False,
            "summary": summary,
        }
    if current_name == new_title and raw_current == new_title:
        return {
            "status": "no_change",
            "changed": False,
            "title_changed": False,
            "summary": summary,
        }

    protected_before = {key: fields.get(key) for key in PROTECTED_CONFIDENCE_KEYS}
    updated_fields = dict(fields)
    updated_fields["smiles"] = smiles
    updated_fields["display_name"] = new_title
    updated_fields.update(TN.descriptor_fields(tentative))
    protected_after = {key: updated_fields.get(key) for key in PROTECTED_CONFIDENCE_KEYS}
    if protected_before != protected_after:  # pragma: no cover - defensive invariant
        raise RuntimeError(f"protected confidence values changed for {report_dir.name}")

    changed = _write_json_if_changed(confidence_path, updated_fields, apply=apply)
    if evidence_path.exists():
        updated_evidence = dict(evidence)
        updated_target = dict(target)
        updated_target["name"] = new_title
        if tentative:
            updated_target["tentative_descriptor"] = tentative.as_dict()
        else:
            updated_target.pop("tentative_descriptor", None)
        updated_evidence["target"] = updated_target
        changed = _write_json_if_changed(evidence_path, updated_evidence, apply=apply) or changed

    if manifest_path.exists():
        manifest = _read_json(manifest_path)
        manifest_target = manifest.get("target") if isinstance(manifest.get("target"), dict) else {}
        updated_manifest = dict(manifest)
        updated_manifest_target = dict(manifest_target)
        updated_manifest_target["name"] = new_title
        if tentative:
            updated_manifest_target["tentative_descriptor"] = tentative.as_dict()
        else:
            updated_manifest_target.pop("tentative_descriptor", None)
        updated_manifest["target"] = updated_manifest_target
        changed = _write_json_if_changed(manifest_path, updated_manifest, apply=apply) or changed

    html_changed = False
    for path in sorted(report_dir.glob("*.html")):
        try:
            text = path.read_text(encoding="utf-8")
        except OSError:
            continue
        if not _is_report_html(path, text):
            continue
        patched = patch_html(text, report_dir.name, new_title, tentative)
        if patched != text:
            html_changed = True
            if apply:
                path.write_text(patched, encoding="utf-8")
    changed = changed or html_changed

    title_changed = old_title != new_title or html_changed
    return {
        "status": "descriptor" if tentative else "fallback",
        "changed": changed,
        "title_changed": title_changed,
        "summary": summary,
    }


def _discover_default(path_name: str, candidates: list[Path]) -> Path | None:
    return next((path for path in candidates if path.is_file()), None)


def _discover_db(reports_root: Path) -> Path | None:
    workspace_root = PROJECT_ROOT.parents[1]
    candidates = [
        reports_root / "metabolite_synonyms.sqlite",
        PROJECT_ROOT / "vendor" / "metabolite_synonyms.sqlite",
        workspace_root / "spec2mol" / "Metabolite Synonym Database" / "metabolite_synonyms.sqlite",
        workspace_root
        / "spec2mol_mass_expansion_fullrange"
        / "Metabolite Synonym Database"
        / "metabolite_synonyms.sqlite",
    ]
    return next((path for path in candidates if path.is_file()), None)


def _build_extended_context(args: argparse.Namespace, reports_root: Path) -> ExtendedNamingContext | None:
    external_db = args.external_db.resolve() if args.external_db else _discover_default(
        "external_db",
        [
            PROJECT_ROOT / "reference_data" / "external_name_index.sqlite",
            reports_root / "reference_data" / "external_name_index.sqlite",
        ],
    )
    classyfire_db = args.classyfire_db.resolve() if args.classyfire_db else _discover_default(
        "classyfire_db",
        [
            PROJECT_ROOT / "reference_data" / "classyfire_hmdb.sqlite",
            reports_root / "reference_data" / "classyfire_hmdb.sqlite",
        ],
    )
    iupac_cache = args.iupac_cache.resolve() if args.iupac_cache else _discover_default(
        "iupac_cache",
        [reports_root / "iupac_cache.jsonl"],
    )
    chemont_model_dir = (
        args.chemont_model_dir.resolve()
        if args.chemont_model_dir
        else _discover_default(
            "chemont_model_dir",
            [PROJECT_ROOT.parents[1] / "references" / "production"],
        )
    )
    if not any((external_db, classyfire_db, iupac_cache)):
        return None
    return ExtendedNamingContext(
        external_db_path=external_db,
        classyfire_db_path=classyfire_db,
        iupac_cache_path=iupac_cache,
        chemont_model_dir=chemont_model_dir,
        lib_id=args.chemont_model_lib,
    )


def _write_summary(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reports_root_positional", nargs="?", type=Path)
    parser.add_argument("--reports-root", type=Path)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--apply", action="store_true", help="Write sidecar and HTML changes")
    mode.add_argument("--dry-run", action="store_true", help="Preview only (default)")
    parser.add_argument(
        "--named-library-only",
        action="store_true",
        help="Patch only reports the Isomer-of naming policy decides the name for",
    )
    parser.add_argument("--db", type=Path, help="Metabolite synonym SQLite database")
    parser.add_argument("--external-db", type=Path, help="Merged external name index SQLite")
    parser.add_argument("--classyfire-db", type=Path, help="ClassyFire HMDB SQLite index")
    parser.add_argument("--iupac-cache", type=Path, help="STOUT IUPAC cache JSONL")
    parser.add_argument("--chemont-model-dir", type=Path, help="ChemOnt model directory")
    parser.add_argument("--chemont-model-lib", type=int, default=400, help="ChemOnt library id")
    parser.add_argument("--summary", type=Path, help="Summary CSV output path")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    reports_root = args.reports_root or args.reports_root_positional or Path("_audit_pull/reports")
    reports_root = reports_root.resolve()
    if not reports_root.is_dir():
        raise SystemExit(f"missing reports root: {reports_root}")
    db_path = args.db.resolve() if args.db else _discover_db(reports_root)
    db_names = {} if args.named_library_only else TN.load_db_names(db_path)
    extended_ctx = (
        None
        if args.named_library_only
        else _build_extended_context(args, reports_root)
    )
    summary_path = args.summary or (reports_root / "tentative_naming_report.csv")

    counts: dict[str, int] = {}
    summaries: list[dict[str, Any]] = []
    changed = 0
    title_changes = 0
    for report_dir in sorted(path for path in reports_root.iterdir() if path.is_dir()):
        if args.named_library_only:
            slug = report_dir.name.lower()
            fields = _read_json(
                report_dir / "data" / f"{slug}_confidence_inputs.json"
            )
            if not _in_isomer_scope(report_dir, fields):
                continue
        outcome = patch_report(
            report_dir,
            db_names,
            apply=bool(args.apply),
            extended_ctx=extended_ctx,
        )
        status = str(outcome["status"])
        counts[status] = counts.get(status, 0) + 1
        changed += int(bool(outcome["changed"]))
        title_changes += int(bool(outcome.get("title_changed")))
        if outcome["summary"]:
            summaries.append(outcome["summary"])
    _write_summary(summary_path, summaries)
    mode = "apply" if args.apply else "dry-run"
    print(
        f"[tentative-naming] mode={mode} reports={sum(counts.values())} "
        f"changed={changed} title_changes={title_changes} db_names={len(db_names)} "
        f"counts={counts} summary={summary_path}"
    )


if __name__ == "__main__":
    main()
