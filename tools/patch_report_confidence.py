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

Compliance: public variants omit MS2 itemization and ion-level detail; only the
internal variant breaks out the MS2 axis. The single overall confidence (which
uses MS2 where available) is shown in both. Idempotent via HTML comment markers.

Usage:
    python3 tools/patch_report_confidence.py <reports_root> [--dry-run]
"""
from __future__ import annotations

import argparse
import html
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
    r"(Parameters Contributing to Identification Confidence.*?</h2>)\s*<ul>.*?</ul>", re.S | re.I
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


def _fmt(v, digits=4):
    try:
        return f"{float(v):.{digits}f}"
    except (TypeError, ValueError):
        return None


def patch_html(text: str, conf: C.Confidence, fields: dict, control_id: str) -> str:
    public = "Dashboard \u2022 Internal" not in text
    badge = f"<!--CONF_BADGE-->{C.confidence_badge_html(conf)}<!--/CONF_BADGE-->"
    card = f"<!--CONF_CARD-->{C.confidence_metric_card_html(conf)}<!--/CONF_CARD-->"
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

    # 3. scored-parameters section (idempotent: marker first, else old <ul>)
    if PARAMS_MARK_RE.search(text):
        text = PARAMS_MARK_RE.sub(lambda _m: params, text)
    else:
        text = PARAMS_UL_RE.sub(lambda m: m.group(1) + "\n      " + params, text)

    # 4. observed parent m/z (MS1, allowed in public) - fix nan display
    observed = _fmt(fields.get("observed_mz"))
    err = _fmt(fields.get("abs_mass_error_da"))
    obs_disp = observed if observed is not None else "not recorded"
    err_disp = err if err is not None else "not recorded"
    text = MASS_ROW_RE.sub(lambda m: m.group(1) + obs_disp + m.group(2) + err_disp + m.group(3), text)
    text = PARENT_MZ_CARD_RE.sub(lambda m: m.group(1) + obs_disp + m.group(3), text)

    # 5. co-abundance correlations (gap 1): replace proxy with WGCNA placeholders
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

    # 5b. MS Fit / Feature Match section: pipeline MS2 evidence (internal) or note (public)
    text = MS_FIT_RE.sub(
        lambda m: m.group(1) + "\n      " + C.ms2_evidence_section_html(fields, public=public) + "\n    " + m.group(3),
        text,
    )

    # 6. display name (gap 2): replace the structure-card caption showing the control id
    name = fields.get("display_name")
    if name and control_id:
        text = text.replace(
            f"<p><strong>{control_id}</strong></p>",
            f"<p><strong>{html.escape(name)}</strong></p>",
        )
    return text


def patch_report(report_dir: Path, apply: bool) -> str:
    data = report_dir / "data"
    if not data.is_dir():
        return "no_data"
    slug = report_dir.name.lower()
    fields = C.read_fields_from_sidecars(data, slug)
    if not fields:
        return "no_fields"
    conf = C.compute_confidence(fields)
    changed = 0
    for hf in sorted(report_dir.glob("*.html")):
        text = hf.read_text()
        if BADGE_ANCHOR not in text and "<!--CONF_BADGE-->" not in text:
            continue  # not a compound report page (e.g. guides.html)
        new = patch_html(text, conf, fields, report_dir.name)
        if new != text:
            if apply:
                hf.write_text(new)
            changed += 1
    return "patched" if changed else "unchanged"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("reports_root", type=Path)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    counts: dict[str, int] = {}
    for d in sorted(p for p in args.reports_root.iterdir() if p.is_dir()):
        st = patch_report(d, apply=not args.dry_run)
        counts[st] = counts.get(st, 0) + 1
    print(f"[patch-confidence] {'(dry-run) ' if args.dry_run else ''}{counts}")


if __name__ == "__main__":
    main()
