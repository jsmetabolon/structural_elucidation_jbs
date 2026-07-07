#!/usr/bin/env python3
"""Patch the top-10 candidate widget in already-generated reports.

Surfaces the pipeline's ``evidence_confidence`` + ``mass_match_score`` (the real
identification signals, read from each report's ``confidence_inputs.json``
``candidate_evidence`` array) and relabels the raw ``candidate_report_score`` as a
relative tie-breaker. Four in-place edits per HTML variant:

  1. the ``<select>`` option labels (``... score X`` -> ``... conf E``),
  2. the server-rendered initial summary card,
  3. the embedded ``candidate-consensus-data`` JSON (add ev + mass_match per rank),
  4. the ``renderCandidateSummary`` JS builder,

plus a one-line "reading the scores" note. Idempotent (sentinel-guarded). Pure
standard library.

Usage:
    python3 tools/patch_candidate_scores.py <reports_root>
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

SENTINEL = "Rank score (relative)"

JSON_RE = re.compile(
    r'(<script type="application/json" id="candidate-consensus-data">)(.*?)(</script>)', re.S
)
SUMMARY_DIV_RE = re.compile(r'(<div id="candidate-consensus-summary">)(.*?)(</div>\s*<div class="consensus-box")', re.S)
JS_RE = re.compile(r"      function renderCandidateSummary\(rank\) \{.*?\n      \}", re.S)
OPTION_RE = re.compile(r'<option value="(\d+)"([^>]*)>(.*?)</option>', re.S)
SCORE_SUFFIX_RE = re.compile(r"\s*(?:&middot;|\u00b7)\s*score\s*[-0-9.eE]+\s*$")
NOTE_ANCHOR = "atom circles show maximum-common-substructure retention across the top-10 candidate pool.</p>"
NOTE_HTML = (
    '\n      <p class="note"><strong>Reading the scores:</strong> "Evidence confidence" '
    "(library + mass + spectral corroboration) and \"Mass match\" are the identification "
    'signals. The "Rank score (relative)" is a small, uncalibrated value the reranker uses '
    "only to order otherwise-equivalent candidates (e.g. mass-consistent isomers), so a low "
    "absolute value does not indicate a poor match.</p>"
)

NEW_JS = """      function renderCandidateSummary(rank) {
        const summary = consensusSummaries[String(rank)] || {};
        const bits = [];
        if (summary.rank) {
          bits.push('<p><strong>Rank ' + summary.rank + '</strong> \u00b7 '
            + (summary.predicted_formula || '') + '</p>');
        }
        const fmt = (v, nd) => (v !== undefined && v !== ''
          && Number.isFinite(Number(v))) ? Number(v).toFixed(nd) : null;
        const seg = [];
        const ev = fmt(summary.evidence_confidence, 3);
        const mm = fmt(summary.mass_match_score, 3);
        const rs = fmt(summary.candidate_report_score, 4);
        const me = fmt(summary.abs_mass_error_da, 6);
        if (ev !== null) seg.push('Evidence confidence: <strong>' + ev + '</strong>');
        if (mm !== null) seg.push('Mass match: <strong>' + mm + '</strong>');
        if (rs !== null) seg.push('Rank score (relative): <strong>' + rs + '</strong>');
        if (me !== null) seg.push('Mass error (Da): <strong>' + me + '</strong>');
        if (seg.length) bits.push('<p>' + seg.join(' \u00b7 ') + '</p>');
        if (summary.predicted_smiles) bits.push('<p><code>' + summary.predicted_smiles + '</code></p>');
        candidateConsensusSummary.innerHTML = bits.join('') || '<p class=\"note\">No candidate summary available.</p>';
      }"""


def _fmt(v, nd):
    try:
        return f"{float(v):.{nd}f}"
    except (TypeError, ValueError):
        return None


def _score_line(summary: dict) -> str:
    parts = []
    ev = _fmt(summary.get("evidence_confidence"), 3)
    mm = _fmt(summary.get("mass_match_score"), 3)
    rs = _fmt(summary.get("candidate_report_score"), 4)
    me = _fmt(summary.get("abs_mass_error_da"), 6)
    if ev is not None:
        parts.append(f"Evidence confidence: <strong>{ev}</strong>")
    if mm is not None:
        parts.append(f"Mass match: <strong>{mm}</strong>")
    if rs is not None:
        parts.append(f"Rank score (relative): <strong>{rs}</strong>")
    if me is not None:
        parts.append(f"Mass error (Da): <strong>{me}</strong>")
    return ("<p>" + " \u00b7 ".join(parts) + "</p>") if parts else ""


def _summary_card(summary: dict) -> str:
    from html import escape
    return (
        "<div class='candidate-summary'>"
        f"<p><strong>Rank {escape(str(summary.get('rank', '')))}</strong> \u00b7 "
        f"{escape(str(summary.get('predicted_formula', '')))}</p>"
        f"{_score_line(summary)}"
        f"<p><code>{escape(str(summary.get('predicted_smiles', '')))}</code></p>"
        "</div>"
    )


def _evidence_by_rank(data_dir: Path, slug: str) -> dict[int, dict]:
    p = data_dir / f"{slug}_confidence_inputs.json"
    if not p.exists():
        return {}
    try:
        rows = json.loads(p.read_text()).get("candidate_evidence") or []
    except (ValueError, OSError):
        return {}
    out = {}
    for r in rows:
        try:
            out[int(r["rank"])] = r
        except (KeyError, TypeError, ValueError):
            continue
    return out


def patch_html(text: str, ev_by_rank: dict[int, dict]) -> str | None:
    if SENTINEL in text:
        return None  # already patched (idempotent)

    m = JSON_RE.search(text)
    if not m:
        return None
    try:
        payload = json.loads(m.group(2))
    except ValueError:
        return None
    summaries = payload.get("summaries") or {}
    if not summaries:
        return None
    for rank_s, summ in summaries.items():
        try:
            ev = ev_by_rank.get(int(rank_s), {})
        except (TypeError, ValueError):
            ev = {}
        summ["evidence_confidence"] = ev.get("evidence_confidence", "")
        summ["mass_match_score"] = ev.get("mass_match_score", "")
    new_json = json.dumps(payload).replace("</", "<\\/")
    text = text[: m.start()] + m.group(1) + new_json + m.group(3) + text[m.end():]

    # server-rendered initial card -> rebuild from the default (selected) rank
    sel = re.search(r'<option value="(\d+)"[^>]*\bselected\b', text)
    default_rank = sel.group(1) if sel else "1"
    default_summary = summaries.get(default_rank) or summaries.get("1") or {}
    if default_summary:
        text = SUMMARY_DIV_RE.sub(
            lambda mm: mm.group(1) + _summary_card(default_summary) + mm.group(3),
            text,
            count=1,
        )

    # relabel <select> options: "... score X" -> "... conf E"
    def _opt(mm: re.Match) -> str:
        rank_s, attrs, label = mm.group(1), mm.group(2), mm.group(3)
        base = SCORE_SUFFIX_RE.sub("", label)
        try:
            ev = ev_by_rank.get(int(rank_s), {})
        except (TypeError, ValueError):
            ev = {}
        conf = _fmt(ev.get("evidence_confidence"), 3)
        suffix = f" \u00b7 conf {conf}" if conf is not None else ""
        return f'<option value="{rank_s}"{attrs}>{base}{suffix}</option>'

    text = OPTION_RE.sub(_opt, text)

    # JS builder
    text = JS_RE.sub(lambda _m: NEW_JS, text, count=1)

    # reading-the-scores note (once)
    if NOTE_ANCHOR in text and "Reading the scores" not in text:
        text = text.replace(NOTE_ANCHOR, NOTE_ANCHOR + NOTE_HTML, 1)
    return text


def patch_report(report_dir: Path) -> str:
    data = report_dir / "data"
    if not data.is_dir():
        return "no_data"
    slug = report_dir.name.lower()
    ev_by_rank = _evidence_by_rank(data, slug)
    if not ev_by_rank:
        return "no_evidence"
    changed = False
    seen = False
    for name in ("index.html", "internal.html", "public.html"):
        f = report_dir / name
        if not f.exists():
            continue
        seen = True
        text = f.read_text()
        new = patch_html(text, ev_by_rank)
        if new is not None and new != text:
            f.write_text(new)
            changed = True
    if not seen:
        return "no_html"
    return "patched" if changed else "unchanged"


def main() -> None:
    root = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("_audit_pull/reports")
    counts: dict[str, int] = {}
    for d in sorted(p for p in root.iterdir() if p.is_dir() and p.name.upper().startswith("X-")):
        r = patch_report(d)
        counts[r] = counts.get(r, 0) + 1
    print(f"[patch-candidate-scores] {counts}")


if __name__ == "__main__":
    main()
