#!/usr/bin/env python3
"""Patch the top-10 candidate widget in already-generated reports (v2).

Features a per-candidate ``candidate_confidence`` that meshes with the 1..10
ranking (monotone non-increasing; rank 1 equals the report overall confidence),
read from each report's ``confidence_inputs.json`` ``candidate_evidence`` array
(written by ``tools/enrich_report_evidence.py``). The class-level
``evidence_confidence`` is relabeled "Formula/class evidence" and shown as
supporting context alongside "Mass match"; the tiny, uncalibrated raw rerank
score is no longer featured. Tied same-formula isomers are flagged.

This supersedes the v1 patch (which featured the saturated, isomer-shared
``evidence_confidence`` and labeled the raw score "Rank score (relative)"). It
transforms both original and v1-patched reports into v2 and is idempotent via an
HTML-comment sentinel embedded in the summary card.

Edits per HTML variant:
  1. the ``<select>`` option labels (``... score X`` / ``... conf X`` -> ``... conf C``),
  2. the server-rendered initial summary card (rebuilt from the default rank),
  3. the embedded ``candidate-consensus-data`` JSON (candidate_confidence + tie +
     evidence + mass_match per rank),
  4. the ``renderCandidateSummary`` JS builder,
  5. the "Reading the scores" note (replaced or injected).

Pure standard library.

Usage:
    python3 tools/patch_candidate_scores.py <reports_root>
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "structural_elucidation" / "reports"))
import confidence as C  # noqa: E402

# Idempotency marker embedded in the rebuilt summary card (always present after v2).
V2_SENTINEL = "<!--cc-v2-->"

JSON_RE = re.compile(
    r'(<script type="application/json" id="candidate-consensus-data">)(.*?)(</script>)', re.S
)
SUMMARY_DIV_RE = re.compile(r'(<div id="candidate-consensus-summary">)(.*?)(</div>\s*<div class="consensus-box")', re.S)
JS_RE = re.compile(r"      function renderCandidateSummary\(rank\) \{.*?\n      \}", re.S)
OPTION_RE = re.compile(r'<option value="(\d+)"([^>]*)>(.*?)</option>', re.S)
# Strip an existing " · score 0.03" or " · conf 0.99 (tied)" suffix from an option.
SUFFIX_RE = re.compile(r"\s*(?:&middot;|\u00b7)\s*(?:score|conf)\s*[-0-9.eE]+(?:\s*\(tied\))?\s*$")
NOTE_ANCHOR = "atom circles show maximum-common-substructure retention across the top-10 candidate pool.</p>"
NOTE_BODY = C.CANDIDATE_SCORES_NOTE
NOTE_HTML = "\n      " + NOTE_BODY
NOTE_RE = re.compile(r'<p class="note"><strong>Reading the scores:</strong>.*?</p>', re.S)

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
        const cc = fmt(summary.candidate_confidence, 3);
        const mm = fmt(summary.mass_match_score, 3);
        const ev = fmt(summary.evidence_confidence, 3);
        const me = fmt(summary.abs_mass_error_da, 6);
        if (cc !== null) seg.push('Candidate confidence: <strong>' + cc + '</strong>');
        if (mm !== null) seg.push('Mass match: <strong>' + mm + '</strong>');
        if (ev !== null) seg.push('Formula/class evidence: <strong>' + ev + '</strong>');
        if (me !== null) seg.push('Mass error (Da): <strong>' + me + '</strong>');
        if (seg.length) bits.push('<p>' + seg.join(' \u00b7 ') + '</p>');
        if (summary.candidate_confidence_tied) {
          bits.push('<p class="note">Tied with other same-formula candidates: the '
            + 'formula/class is supported but the specific isomer is not resolved.</p>');
        }
        if (summary.predicted_smiles) bits.push('<p><code>' + summary.predicted_smiles + '</code></p>');
        candidateConsensusSummary.innerHTML = bits.join('') || '<p class=\"note\">No candidate summary available.</p>';
      }"""


def _fmt(v, nd):
    try:
        return f"{float(v):.{nd}f}"
    except (TypeError, ValueError):
        return None


def _summary_card(summary: dict) -> str:
    from html import escape
    return (
        "<div class='candidate-summary'>" + V2_SENTINEL
        + f"<p><strong>Rank {escape(str(summary.get('rank', '')))}</strong> \u00b7 "
        f"{escape(str(summary.get('predicted_formula', '')))}</p>"
        f"{C.candidate_score_line_html(summary)}"
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
    if V2_SENTINEL in text:
        return None  # already v2 (idempotent)

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
        summ["candidate_confidence"] = ev.get("candidate_confidence", "")
        summ["candidate_confidence_tied"] = bool(ev.get("candidate_confidence_tied"))
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

    # relabel <select> options: "... score X" / "... conf X" -> "... conf C"
    def _opt(mm: re.Match) -> str:
        rank_s, attrs, label = mm.group(1), mm.group(2), mm.group(3)
        base = SUFFIX_RE.sub("", label)
        try:
            ev = ev_by_rank.get(int(rank_s), {})
        except (TypeError, ValueError):
            ev = {}
        conf = _fmt(ev.get("candidate_confidence"), 3)
        suffix = f" \u00b7 conf {conf}" if conf is not None else ""
        if ev.get("candidate_confidence_tied"):
            suffix += " (tied)"
        return f'<option value="{rank_s}"{attrs}>{base}{suffix}</option>'

    text = OPTION_RE.sub(_opt, text)

    # JS builder
    text = JS_RE.sub(lambda _m: NEW_JS, text, count=1)

    # "N of 10" note when the reranker generated fewer than 10 candidates
    n_shown = len(summaries)
    count_note = (
        ""
        if n_shown >= 10
        else (
            f'<p class="note">Showing <strong>{n_shown}</strong> of up to 10 candidate structures; '
            "the reranker generated fewer than 10 distinct candidates for this feature.</p>"
        )
    )

    # reading-the-scores note: replace an existing note, else inject after the anchor
    if NOTE_RE.search(text):
        text = NOTE_RE.sub(lambda _m: NOTE_BODY + count_note, text, count=1)
    elif NOTE_ANCHOR in text:
        text = text.replace(NOTE_ANCHOR, NOTE_ANCHOR + NOTE_HTML + count_note, 1)
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
    # Patch every HTML variant in the report dir that carries the candidate widget
    # (index/internal/public and the canonical X_<id>.html); files without the
    # widget - e.g. guides.html - are no-ops because patch_html finds no JSON block.
    for f in sorted(report_dir.glob("*.html")):
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
