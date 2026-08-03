#!/usr/bin/env python3
"""Inject the hybrid identification confidence into already-generated reports.

For each report, reads ``<slug>_confidence_inputs.json`` (from
tools/enrich_report_evidence.py), computes the confidence with the shared
``structural_elucidation/reports/confidence.py`` model, and injects into each
HTML variant:
  - a confidence badge in the header,
  - an 'Overall identification confidence' metric card,
  - a scored-parameters table replacing the old text 'Parameters...' body,
  - the real observed parent m/z (fixing the ``nan`` display).

Compliance: public variants may show the aggregate MS2 score that contributes
to confidence, but omit ion-level and per-fragment detail. Idempotent via HTML
comment markers.

Usage:
    python3 tools/patch_report_confidence.py <reports_root> [--dry-run]
"""
from __future__ import annotations

import argparse
import csv
import html
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "structural_elucidation" / "reports"))
import confidence as C  # noqa: E402

BADGE_ANCHOR = '<span class="badge ok">top-1 candidate</span>'
CARD_ANCHOR = '<div class="metric-card"><div class="metric-label">Top candidate structural similarity</div>'

BADGE_RE = re.compile(r"<!--CONF_BADGE-->.*?<!--/CONF_BADGE-->", re.S)
CARD_RE = re.compile(r"<!--CONF_CARD-->.*?<!--/CONF_CARD-->", re.S)
PARAMS_MARK_RE = re.compile(r"<!--CONF_PARAMS-->.*?<!--/CONF_PARAMS-->", re.S)
PARAMS_UL_RE = re.compile(
    r"(Parameters Contributing to (?:Identification|ID) Confidence.*?</h2>)\s*<ul>.*?</ul>",
    re.S | re.I,
)
ID_PARAMETER_CARDS_RE = re.compile(
    r"<!--ID_PARAMETER_CARDS-->.*?<!--/ID_PARAMETER_CARDS-->", re.S
)
OLD_PARAMETER_CARD_RE = re.compile(
    r'<div class="metric-card"><div class="metric-label">'
    r"(?:Top candidate structural similarity|Candidate Chemical Pathway|"
    r"Candidate Chemical SubPathway|Library neighbors)"
    r"</div>.*?</div></div>",
    re.S | re.I,
)
STRUCTURE_NOTE_RE = re.compile(
    r'<p class="note">(?:Exact stereochemistry is not overclaimed; '
    r'this is the top-ranked candidate structure for this feature\.|'
    r'This is the top ranked candidate structure for this feature\. '
    r'Note that exact stereochemistry is not known from available evidence\.)</p>'
)
CONSENSUS_HEADING_RE = re.compile(
    r"<h3>(?:Final Consensus Atom-Position Confidence|"
    r"Atom-position confidence mapped to top candidate)</h3>",
    re.I,
)
MASS_ROW_RE = re.compile(
    r"(<tr><td>[0-9.]+</td><td>)[^<]*(</td><td>(?:Negative|Positive)[^<]*</td><td>)[^<]*(</td></tr>)"
)
PARENT_MZ_CARD_RE = re.compile(
    r'(Observed parent m/z</div><div class="metric-value">)([^<]*)(</div>)'
)
CORR_NOTE_RE = re.compile(
    r'(Top Correlating Metabolites.*?</h2>\s*<p class="note">).*?(</p>)', re.S | re.I
)
CORR_TBODY_RE = re.compile(
    r"(<th>Correlation \(R\)</th><th>)Pathway(</th></tr></thead>\s*<tbody>).*?(</tbody>)", re.S
)
MS_FIT_RE = re.compile(
    r"(<h2>(?:MS Fit Scores|Feature Match Scores)\b.*?</h2>)(.*?)(</section>)", re.S
)
FRAGMENT_SCORES_RE = re.compile(
    r"(<h2>(?:Fragment Scores|Substructure Scores)\b.*?</h2>.*?<tbody>)(.*?)(</tbody>)",
    re.S,
)
FRAGMENT_SUMMARY_MARK_RE = re.compile(r"<!--FRAGMENT_SUMMARY-->.*?<!--/FRAGMENT_SUMMARY-->", re.S)
FRAGMENT_DIAGNOSTICS_MARK_RE = re.compile(r"<!--FRAGMENT_DIAGNOSTICS-->.*?<!--/FRAGMENT_DIAGNOSTICS-->", re.S)


def _fmt(v, digits=4):
    try:
        return f"{float(v):.{digits}f}"
    except (TypeError, ValueError):
        return None


def patch_html(
    text: str,
    conf: C.Confidence,
    fields: dict,
    control_id: str,
    parameter_cards: str,
    fragment_candidate_rows: list[dict] | None = None,
) -> str:
    public = "Dashboard \u2022 Internal" not in text
    identity_badge = C.identity_badge_html(conf)
    identity_badge_block = ("\n          " + identity_badge) if identity_badge else ""
    badge = (
        f"<!--CONF_BADGE-->{C.confidence_badge_html(conf)}\n          "
        f"{C.msi_badge_html(conf)}"
        f"{identity_badge_block}<!--/CONF_BADGE-->"
    )
    card = (
        f"<!--CONF_CARD-->{C.confidence_metric_card_html(conf)}\n        "
        f"{C.msi_metric_card_html(conf)}\n        "
        f"{C.identity_metric_card_html(conf)}<!--/CONF_CARD-->"
    )
    params = f"<!--CONF_PARAMS-->{C.scored_parameters_section_html(conf, public=public)}<!--/CONF_PARAMS-->"

    # 1. header badge (idempotent)
    if BADGE_RE.search(text):
        text = BADGE_RE.sub(lambda _m: badge, text)
    elif BADGE_ANCHOR in text:
        text = text.replace(BADGE_ANCHOR, BADGE_ANCHOR + "\n          " + badge, 1)

    # 2. confidence metric card (idempotent)
    if CARD_RE.search(text):
        text = CARD_RE.sub(lambda _m: card, text)
    elif CARD_ANCHOR in text:
        text = text.replace(CARD_ANCHOR, card + "\n        " + CARD_ANCHOR, 1)

    # 3. relocate high-level evidence cards and update the section name.
    text = ID_PARAMETER_CARDS_RE.sub("", text)
    text = OLD_PARAMETER_CARD_RE.sub("", text)
    text = re.sub(
        r"Parameters Contributing To? (?:Identification|ID) Confidence",
        "Parameters Contributing to ID Confidence",
        text,
        flags=re.I,
    )
    text = text.replace(
        "parameters-contributing-to-identification-confidence",
        "parameters-contributing-to-id-confidence",
    )

    # 4. scored-parameters section (idempotent: marker first, else old <ul>)
    if PARAMS_MARK_RE.search(text):
        text = PARAMS_MARK_RE.sub(lambda _m: params, text)
    else:
        text = PARAMS_UL_RE.sub(lambda m: m.group(1) + "\n      " + params, text)
    if "<!--CONF_PARAMS-->" in text:
        text = re.sub(r"\n[ \t]*\n([ \t]*<!--CONF_PARAMS-->)", r"\n\1", text)
        text = text.replace("<!--CONF_PARAMS-->", parameter_cards + "\n      <!--CONF_PARAMS-->", 1)

    # 5. requested structure and consensus explanatory copy.
    structure_note = (
        '<p class="note">This is the top ranked candidate structure for this feature. '
        "Note that exact stereochemistry is not known from available evidence.</p>"
    )
    text = STRUCTURE_NOTE_RE.sub(structure_note, text)
    text = re.sub(
        r"\s*<!--CONSENSUS_NOTE-->.*?<!--/CONSENSUS_NOTE-->",
        "",
        text,
        flags=re.S,
    )
    text = CONSENSUS_HEADING_RE.sub(
        "<h3>Atom-position confidence mapped to top candidate</h3>",
        text,
        count=1,
    )

    # 6. observed parent m/z (MS1, allowed in public) - fix nan display
    observed = _fmt(fields.get("observed_mz"))
    err = _fmt(fields.get("abs_mass_error_da"))
    obs_disp = observed if observed is not None else "not recorded"
    err_disp = err if err is not None else "not recorded"
    text = MASS_ROW_RE.sub(lambda m: m.group(1) + obs_disp + m.group(2) + err_disp + m.group(3), text)
    text = PARENT_MZ_CARD_RE.sub(lambda m: m.group(1) + obs_disp + m.group(3), text)

    # 7. co-abundance correlations (gap 1): replace proxy with WGCNA placeholders
    corr = fields.get("correlations") or []
    if corr:
        note = ("Illustrative co-abundance neighbors from reference-cohort WGCNA modules "
                "(placeholder; client-specific co-abundance correlations pending).")
        text = CORR_TBODY_RE.sub(
            lambda m: m.group(1) + "Co-abundance module" + m.group(2)
            + C.correlations_table_html(corr) + m.group(3),
            text,
        )
    else:
        note = ("Co-abundance correlations pending (client-specific dataset); rows below are "
                "structural-proxy neighbors, not measured co-abundance.")
    text = CORR_NOTE_RE.sub(lambda m: m.group(1) + note + m.group(2), text)

    fragment_rows = fragment_candidate_rows or []
    diagnostics = C.diagnostic_substructure_rows_html(fragment_rows[0] if fragment_rows else {}, {})
    if FRAGMENT_DIAGNOSTICS_MARK_RE.search(text):
        text = FRAGMENT_DIAGNOSTICS_MARK_RE.sub(lambda _m: diagnostics, text)
    else:
        text = FRAGMENT_SCORES_RE.sub(lambda m: m.group(1) + diagnostics + m.group(3), text)

    summary = C.fragment_summary_section_html(fields, fragment_rows, public=public)
    if FRAGMENT_SUMMARY_MARK_RE.search(text):
        text = FRAGMENT_SUMMARY_MARK_RE.sub(lambda _m: summary, text)
    else:
        text = MS_FIT_RE.sub(
            lambda m: (
                m.group(1)
                + '\n      <p class="note">Metabolon’s elucidation pipeline uses mass fragment (MS2) matching as a key parameter for structural elucidation. The tables below summarize observed fragment counts and how well each candidate explains the spectrum.</p>\n      '
                + summary
                + "\n    "
                + m.group(3)
            ),
            text,
        )

    # 8. display name (gap 2): replace the structure-card caption showing the control id
    name = fields.get("display_name")
    if name and control_id:
        text = text.replace(
            f"<p><strong>{control_id}</strong></p>",
            f"<p><strong>{html.escape(name)}</strong></p>",
        )
    return text


def _read_csv_rows(path: Path) -> list[dict]:
    if not path.exists():
        return []
    try:
        with path.open(newline="", encoding="utf-8-sig") as handle:
            return list(csv.DictReader(handle))
    except (OSError, csv.Error):
        return []


def _report_context(data: Path, slug: str, fields: dict, conf: C.Confidence) -> str:
    """Build the shared parameter-card block from persisted report sidecars."""

    evidence_path = data / f"{slug}_evidence_summary.json"
    if evidence_path.exists():
        try:
            evidence = json.loads(evidence_path.read_text())
            target = evidence.get("target") or {}
            if target.get("exact_mass") is not None:
                fields["candidate_exact_mass"] = target["exact_mass"]
        except (ValueError, OSError):
            pass

    pathway_rows = _read_csv_rows(data / f"{slug}_pathway_prediction.csv")
    pathway_row = pathway_rows[0] if pathway_rows else {}
    pathway_info = {
        "super": pathway_row.get("predicted_pathway", ""),
        "sub": pathway_row.get("predicted_subpathway", ""),
        "source": pathway_row.get("pathway_source", ""),
        "tier": pathway_row.get("pathway_confidence_tier", ""),
        "similarity": pathway_row.get("pathway_best_reference_similarity", 0.0),
        "status": pathway_row.get("pathway_assignment_status", ""),
        "note": (
            "From reference DB record"
            if pathway_row.get("pathway_source") == "exact"
            else "Predicted from reference-library neighbors"
        ),
    }
    try:
        top_similarity = float(pathway_info["similarity"] or fields.get("structure_similarity") or 0.0)
    except (TypeError, ValueError):
        top_similarity = 0.0
    candidate_scores = _read_csv_rows(data / f"{slug}_top_candidates.csv")[:10]
    neighbor_count = len(_read_csv_rows(data / "similarity_neighbors.csv"))
    return C.identification_parameter_cards_html(
        conf,
        candidate_scores=candidate_scores,
        pathway_info=pathway_info,
        top_similarity=top_similarity,
        neighbor_count=neighbor_count,
    ), candidate_scores[:3]


def patch_report(report_dir: Path, apply: bool, identity_overrides: dict | None = None) -> str:
    data = report_dir / "data"
    if not data.is_dir():
        return "no_data"
    slug = report_dir.name.lower()
    fields = C.read_fields_from_sidecars(data, slug)
    if not fields:
        return "no_fields"
    # Inject curated known-identity override (curator allowlist) so the library-
    # identity metric can render >= 90% without rerunning enrichment. Never
    # touches the calibrated ``overall`` P(correct) or its tier.
    if identity_overrides:
        override = identity_overrides.get(report_dir.name)
        if override:
            fields["known_identity_override"] = override
    # Load persisted target mass before confidence computation so the mass card
    # can compare the observed neutral mass with the predicted candidate mass.
    evidence_path = data / f"{slug}_evidence_summary.json"
    if evidence_path.exists():
        try:
            target = (json.loads(evidence_path.read_text()).get("target") or {})
            if target.get("exact_mass") is not None:
                fields["candidate_exact_mass"] = target["exact_mass"]
        except (ValueError, OSError):
            pass
    conf = C.compute_confidence(fields)
    parameter_cards, fragment_candidate_rows = _report_context(data, slug, fields, conf)
    # Persist identity fields into the sidecar so rebuild_master_index can surface them.
    if apply:
        enriched = data / f"{slug}_confidence_inputs.json"
        if enriched.exists():
            try:
                sidecar = json.loads(enriched.read_text())
                sidecar["identity_status"] = conf.identity_status
                sidecar["identity_confidence"] = conf.identity_confidence
                sidecar["identity_basis"] = conf.identity_basis
                sidecar["candidate_exact_mass"] = conf.candidate_exact_mass
                if fields.get("known_identity_override"):
                    sidecar["known_identity_override"] = fields["known_identity_override"]
                enriched.write_text(json.dumps(sidecar, indent=2))
            except (ValueError, OSError):
                pass
    changed = 0
    for hf in sorted(report_dir.glob("*.html")):
        text = hf.read_text()
        if hf.name == "guides.html":
            new = re.sub(
                r"Parameters Contributing to Identification Confidence",
                "Parameters Contributing to ID Confidence",
                text,
                flags=re.I,
            ).replace(
                "parameters-contributing-to-identification-confidence",
                "parameters-contributing-to-id-confidence",
            )
            if new != text:
                if apply:
                    hf.write_text(new)
                changed += 1
            continue
        if BADGE_ANCHOR not in text and "<!--CONF_BADGE-->" not in text:
            continue
        new = patch_html(text, conf, fields, report_dir.name, parameter_cards, fragment_candidate_rows)
        if new != text:
            if apply:
                hf.write_text(new)
            changed += 1
    return "patched" if changed else "unchanged"


DEFAULT_IDENTITY_OVERRIDES = (
    Path(__file__).resolve().parent.parent / "data" / "known_identity_overrides.csv"
)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("reports_root", type=Path)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument(
        "--identity-overrides",
        type=Path,
        default=DEFAULT_IDENTITY_OVERRIDES,
        help=(
            "CSV allowlist of curated known-identity control_ids "
            "(columns: control_id, identity_confidence, reason). "
            "Defaults to data/known_identity_overrides.csv; pass an empty path or "
            "a missing file to disable."
        ),
    )
    args = ap.parse_args()
    identity_overrides = C.load_identity_overrides(args.identity_overrides)
    if identity_overrides:
        print(f"[patch-confidence] loaded {len(identity_overrides)} known-identity overrides "
              f"from {args.identity_overrides}", file=sys.stderr)
    counts: dict[str, int] = {}
    for d in sorted(p for p in args.reports_root.iterdir() if p.is_dir()):
        st = patch_report(d, apply=not args.dry_run, identity_overrides=identity_overrides)
        counts[st] = counts.get(st, 0) + 1
    print(f"[patch-confidence] {'(dry-run) ' if args.dry_run else ''}{counts}")


if __name__ == "__main__":
    main()
