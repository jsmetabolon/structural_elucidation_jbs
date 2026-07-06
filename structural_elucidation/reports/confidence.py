"""Hybrid identification-confidence model for X-compound reports.

Pure standard library so the renderer and the delivered-set patcher share one
implementation and produce identical output.

Design (see docs/report_presentation_gap_findings.md for the evidence behind it):

- Confidence is built from up to four evidence axes, each a 0..1 sub-score with
  an availability flag: structure (nearest-reference Tanimoto / exact catalog
  match), mass (candidate vs target neutral-mass fit), MS2 (fragment/neutral-loss
  spectral fit), and pathway (neighbor pathway-vote consensus).
- ``overall`` is a *calibrated* estimate of P(top structure correct), fitted with
  1-D Platt scaling on the 92 truth-verified controls in the delivered set
  (tools/calibrate_confidence.py). It is deliberately honest: with MS2 absent and
  the observed m/z missing for ~85% of features, the available evidence is weakly
  discriminative, so most predicted compounds land in a narrow mid range.
- ``tier`` is a corroboration-aware rubric, NOT a rename of the score. "Exact
  reference match" is NOT treated as high confidence: on truth-verified controls
  exact matches are correct only ~38%, so an exact catalog hit alone caps at
  Moderate; High requires independent spectral corroboration.
- ``completeness`` and ``flags`` make missing evidence explicit.
"""
from __future__ import annotations

import csv
import html
import json
import math
from dataclasses import dataclass, field
from pathlib import Path

# --- calibration constants (from tools/calibrate_confidence.py) --------------
WEIGHTS = {"structure": 0.40, "mass": 0.30, "ms2": 0.20, "pathway": 0.10}
# Platt scaling of the completeness-weighted raw evidence score -> P(correct).
# Fit by tools/calibrate_confidence.py on the 92 controls that have a known true
# structure (the complete structure-truth set available). The evidence is weakly
# discriminative (MS2 absent; observed m/z missing for ~85%), so probabilities
# compress to a narrow mid range - this is intentionally honest.
PLATT_A = 1.887084
PLATT_B = -1.520154
CALIBRATION_N = 92          # truth-verified controls used to fit/validate
CALIBRATION_BASE_RATE = 0.478
# Tier cutoffs on the calibrated probability (validated to be monotonic in
# observed correctness on the truth set). "High" is intentionally set above the
# range the current evidence reaches: with MS2 present for <=44% and no authentic
# standard, essentially nothing is confidently confirmable yet - honest by design.
TIER_HIGH = 0.65
TIER_MODERATE = 0.48
TIER_LOW = 0.38

AXIS_LABELS = {
    "structure": "Structure similarity",
    "mass": "Mass / formula match",
    "ms2": "MS2 fragment fit",
    "pathway": "Pathway consensus",
}

# Corroboration thresholds
MASS_STRONG = 0.90
MASS_WEAK = 0.50
STRUCT_STRONG = 0.85
STRUCT_MODERATE = 0.70
STRUCT_WEAK = 0.45
MS2_PRESENT = 0.30


@dataclass
class Axis:
    key: str
    label: str
    score: float | None
    weight: float
    available: bool
    detail: str = ""

    def as_dict(self) -> dict:
        return {
            "key": self.key,
            "label": self.label,
            "score": self.score,
            "weight": self.weight,
            "available": self.available,
            "detail": self.detail,
        }


@dataclass
class Confidence:
    overall: float
    tier: str
    tier_note: str
    completeness: float
    descriptor: str
    axes: list[Axis] = field(default_factory=list)
    flags: list[str] = field(default_factory=list)
    source: str = ""
    # display extras
    observed_mz: float | None = None
    observed_neutral_mass: float | None = None
    ms2_n_pairs: int = 0

    def as_dict(self) -> dict:
        return {
            "overall": self.overall,
            "tier": self.tier,
            "tier_note": self.tier_note,
            "completeness": self.completeness,
            "descriptor": self.descriptor,
            "source": self.source,
            "observed_mz": self.observed_mz,
            "observed_neutral_mass": self.observed_neutral_mass,
            "ms2_n_pairs": self.ms2_n_pairs,
            "axes": [a.as_dict() for a in self.axes],
            "flags": self.flags,
        }


def _clamp01(x: float) -> float:
    return 0.0 if x < 0 else 1.0 if x > 1 else x


def _num(v) -> float | None:
    if v is None:
        return None
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(f) else f


def _sigmoid(z: float) -> float:
    if z < -60:
        return 0.0
    if z > 60:
        return 1.0
    return 1.0 / (1.0 + math.exp(-z))


def build_axes(fields: dict) -> list[Axis]:
    """Construct the four evidence axes from a normalized ``fields`` dict."""
    source = str(fields.get("source", "")).strip()

    # Structure: nearest-reference Tanimoto; exact catalog hit pins to 1.0. This
    # measures catalog proximity/plausibility, NOT proof of correct identity.
    struct = 1.0 if source == "exact" else _clamp01(_num(fields.get("structure_similarity")) or 0.0)
    if source == "exact":
        s_detail = "Top structure exactly matches a cataloged reference metabolite"
    elif struct >= STRUCT_STRONG:
        s_detail = "Top structure is very close to a cataloged metabolite"
    elif struct >= STRUCT_WEAK:
        s_detail = "Top structure is moderately similar to the nearest cataloged metabolite"
    else:
        s_detail = "Top structure is distant from any cataloged metabolite (de novo)"

    # Mass: candidate neutral mass vs the target neutral mass hypothesis.
    mass = _num(fields.get("mass_match_score"))
    abs_err = _num(fields.get("abs_mass_error_da"))
    mass_available = mass is not None
    if not mass_available:
        m_detail = "No mass-match score available"
    elif mass >= MASS_STRONG:
        m_detail = "Candidate formula/mass matches the target neutral mass"
        if abs_err is not None:
            m_detail += f" (|error| {abs_err:.3f} Da)"
    elif mass >= MASS_WEAK:
        m_detail = "Partial mass agreement"
    else:
        m_detail = "Candidate mass does NOT match the target neutral mass"
        if abs_err is not None:
            m_detail += f" (|error| {abs_err:.3f} Da)"

    # MS2: pipeline spectral fit of the candidate against the observed MS2
    # spectrum. Availability = an observed MS2 spectrum exists for the feature;
    # the score is the pipeline's forward/predicted spectral match (0 when a
    # spectrum exists but the candidate does not explain it - an informative low).
    ms2_present = bool(fields.get("ms2_spectrum_present"))
    ms2_score = _num(fields.get("ms2_score"))
    ms2 = _clamp01(ms2_score if ms2_score is not None else 0.0) if ms2_present else None
    ms2_available = ms2_present
    if not ms2_present:
        m2_detail = "No observed MS2 spectrum for this feature"
    elif (ms2 or 0.0) >= MS2_PRESENT:
        m2_detail = f"Observed MS2 spectrum explained by candidate (fit {ms2:.2f})"
    else:
        m2_detail = f"Observed MS2 spectrum poorly explained by candidate (fit {ms2 or 0.0:.2f})"

    # Pathway: neighbor pathway-vote consensus strength (falls back to the
    # persisted pathway confidence tier when vote share is not in the sidecar).
    vote = _num(fields.get("pathway_vote_share"))
    if vote is not None:
        pathway = _clamp01(vote)
        p_detail = f"Neighbor pathway consensus {pathway:.0%}"
    else:
        tier_map = {"high": 0.90, "medium": 0.60, "low": 0.30}
        pathway = tier_map.get(str(fields.get("pathway_tier", "")).strip().lower())
        p_detail = (
            f"Pathway confidence tier: {fields.get('pathway_tier')}"
            if pathway is not None
            else "No pathway consensus"
        )
    pathway_available = pathway is not None

    return [
        Axis("structure", AXIS_LABELS["structure"], struct, WEIGHTS["structure"], True, s_detail),
        Axis("mass", AXIS_LABELS["mass"], mass, WEIGHTS["mass"], mass_available, m_detail),
        Axis("ms2", AXIS_LABELS["ms2"], ms2, WEIGHTS["ms2"], ms2_available, m2_detail),
        Axis("pathway", AXIS_LABELS["pathway"], pathway, WEIGHTS["pathway"], pathway_available, p_detail),
    ]


def raw_score(axes: list[Axis]) -> tuple[float, float]:
    """Return (completeness-weighted raw score, completeness fraction)."""
    num = den = avail = 0.0
    total = sum(WEIGHTS.values())
    for ax in axes:
        if ax.available and ax.score is not None:
            num += ax.weight * ax.score
            den += ax.weight
            avail += ax.weight
    raw = num / den if den else 0.0
    return raw, (avail / total if total else 0.0)


def _tier(prob: float, axes: list[Axis], source: str) -> tuple[str, str]:
    """Tier from the calibrated probability, with a mass-mismatch override.

    Tiering off the calibrated probability keeps tiers monotonic in observed
    correctness (unlike a hand rubric). A top candidate whose mass disagrees with
    the observed feature is always demoted, regardless of probability.
    """
    by = {a.key: a for a in axes}
    mass = by["mass"].score if by["mass"].available else None
    if mass is not None and mass < MASS_WEAK:
        return ("Tentative", "Top candidate mass does not match the observed feature - treat with caution")
    if prob >= TIER_HIGH:
        return ("High", "Strong, corroborated evidence across structure, mass, and MS2")
    if prob >= TIER_MODERATE:
        return ("Moderate", "Plausible identification with partial corroboration; not standard-confirmed")
    if prob >= TIER_LOW:
        return ("Low", "Weak-to-moderate support; corroborating evidence is limited")
    return ("Tentative", "Little corroborating evidence beyond structure proximity")


def _flags(axes: list[Axis], fields: dict) -> list[str]:
    by = {a.key: a for a in axes}
    flags: list[str] = []
    if not fields.get("observed_mz_available", False):
        flags.append("No observed m/z recorded for this feature; mass fit uses the hypothesized neutral mass.")
    if not by["ms2"].available:
        flags.append("No matched MS2 fragments; independent spectral confirmation is unavailable.")
    if by["mass"].available and (by["mass"].score or 0.0) < MASS_WEAK:
        flags.append("Top candidate mass disagrees with the target neutral mass.")
    if not fields.get("correlations_real", True):
        flags.append("Co-abundance correlations are proxy (nearest structural neighbors), not measured.")
    margin = _num(fields.get("candidate_margin"))
    if margin is not None and margin < 0.02:
        flags.append("Top candidates are near-ties (isomeric ambiguity).")
    if str(fields.get("source", "")) == "predicted" and (by["structure"].score or 0.0) < STRUCT_WEAK:
        flags.append("Proposed structure is a de novo prediction with no close cataloged analog.")
    return flags


def compute_confidence(fields: dict) -> Confidence:
    """Compute the hybrid identification confidence from normalized ``fields``."""
    source = str(fields.get("source", "")).strip()
    axes = build_axes(fields)
    raw, completeness = raw_score(axes)
    overall = _sigmoid(PLATT_A * raw + PLATT_B)
    tier, tier_note = _tier(overall, axes, source)
    if source == "exact":
        descriptor = "Cataloged reference metabolite (exact structural match)"
    elif source == "predicted":
        descriptor = "De novo / predicted structure"
    else:
        descriptor = "Structure proposal"
    return Confidence(
        overall=round(overall, 4),
        tier=tier,
        tier_note=tier_note,
        completeness=round(completeness, 4),
        descriptor=descriptor,
        axes=axes,
        flags=_flags(axes, fields),
        source=source,
        observed_mz=_num(fields.get("observed_mz")),
        observed_neutral_mass=_num(fields.get("observed_neutral_mass")),
        ms2_n_pairs=int(_num(fields.get("ms2_n_pairs")) or 0),
    )


# --- shared input normalization (renderer + patcher + audit + calibration) ---

def fields_from_parts(
    css_row: dict | None,
    css_next: dict | None,
    pw_row: dict | None,
    *,
    observed_mz_available: bool,
    ms2_spectrum_present: bool = False,
    ms2_score: float | None = None,
    library_neighbors_ge_050: int = 0,
    correlations_real: bool = True,
) -> dict:
    """Normalize a candidate-summary row + pathway row into confidence ``fields``.

    ``css_row``/``css_next`` are rank-1/rank-2 rows from
    ``<slug>_candidate_score_summary.csv`` (or the equivalent compact dicts the
    renderer already builds); ``pw_row`` is the ``<slug>_pathway_prediction.csv``
    row. ``ms2_score``/``ms2_spectrum_present`` come from the pipeline spectral
    scores + observed-MS2 lookup; ``observed_mz_available`` from feature rows.
    This single normalizer keeps the renderer and the patcher identical.
    """
    css_row = css_row or {}
    pw_row = pw_row or {}
    margin = None
    r1 = _num((css_row or {}).get("candidate_report_score"))
    r2 = _num((css_next or {}).get("candidate_report_score")) if css_next else None
    if r1 is not None and r2 is not None:
        margin = abs(r1 - r2)
    return {
        "source": str(pw_row.get("pathway_source", "")).strip(),
        "structure_similarity": _num(pw_row.get("pathway_best_reference_similarity")),
        "mass_match_score": _num(css_row.get("mass_match_score")),
        "abs_mass_error_da": _num(css_row.get("abs_mass_error_da")),
        "ms2_spectrum_present": bool(ms2_spectrum_present),
        "ms2_score": _num(ms2_score),
        "pathway_vote_share": _num(pw_row.get("pathway_vote_share")),
        "pathway_tier": str(pw_row.get("pathway_confidence_tier", "")).strip(),
        "observed_mz_available": bool(observed_mz_available),
        "candidate_margin": margin,
        "library_neighbors_ge_050": library_neighbors_ge_050,
        "correlations_real": correlations_real,
    }


# --- shared HTML builders (renderer + patcher; variant-aware) ----------------
# Compliance: observed MS2 fragment ions/peaks are internal-only. Public reports
# never itemize MS2 or show ion-level detail; the MS2 fit still contributes to the
# single overall confidence value shown in both variants.

_TIER_CLASS = {"High": "badge ok", "Moderate": "badge ok", "Low": "badge", "Tentative": "badge"}


def _pct(x: float | None) -> str:
    return f"{round(100 * (x or 0.0))}%"


def confidence_badge_html(conf: Confidence) -> str:
    cls = _TIER_CLASS.get(conf.tier, "badge")
    return f'<span class="{cls}">Confidence {_pct(conf.overall)} \u00b7 {html.escape(conf.tier)}</span>'


def confidence_metric_card_html(conf: Confidence) -> str:
    return (
        '<div class="metric-card"><div class="metric-label">Overall identification confidence</div>'
        f'<div class="metric-value">{_pct(conf.overall)} \u00b7 {html.escape(conf.tier)}</div>'
        f'<div class="metric-note">{html.escape(conf.tier_note)}</div></div>'
    )


def scored_parameters_section_html(conf: Confidence, *, public: bool) -> str:
    """Body HTML for the 'Parameters Contributing to Identification Confidence' section."""
    rows = []
    for ax in conf.axes:
        if ax.key == "ms2" and public:  # public: no MS2 itemization at all
            continue
        score = "n/a" if (ax.score is None or not ax.available) else f"{ax.score:.2f}"
        avail = "yes" if ax.available else "no"
        rows.append(
            f"<tr><td>{html.escape(ax.label)}</td><td>{score}</td>"
            f"<td>{ax.weight:.2f}</td><td>{avail}</td><td>{html.escape(ax.detail)}</td></tr>"
        )
    table = (
        '<div class="table-wrap"><table><thead><tr>'
        "<th>Evidence axis</th><th>Sub-score (0-1)</th><th>Weight</th>"
        "<th>Available</th><th>Interpretation</th></tr></thead>"
        f"<tbody>{''.join(rows)}</tbody></table></div>"
    )
    intro = (
        '<p class="note">The overall identification confidence is a calibrated estimate of the '
        "probability that the top candidate structure is correct, combining structural similarity, "
        "mass/formula match, MS2 spectral fit, and pathway consensus (calibrated on "
        f"{CALIBRATION_N} structure-verified controls; base rate {CALIBRATION_BASE_RATE:.0%}). "
        f"Proposed structure: <strong>{html.escape(conf.descriptor)}</strong>. "
        f"Assessment: <strong>{html.escape(conf.tier)}</strong> - {html.escape(conf.tier_note)}.</p>"
    )
    if public:
        tail = (
            '<p class="note">Additional MS2 spectral evidence is incorporated into the overall '
            "confidence where available; per-ion detail is provided in the internal report.</p>"
        )
    else:
        tail = (
            f'<p class="note">Evidence completeness: {_pct(conf.completeness)} of the weighted '
            "evidence was available for this feature."
        )
        if conf.ms2_n_pairs:
            tail += f" Observed MS2 spectrum: {int(conf.ms2_n_pairs)} fragment peaks."
        tail += "</p>"
    flags = [f for f in conf.flags if not (public and "MS2" in f)]
    flag_html = ""
    if flags:
        items = "".join(f"<li>{html.escape(f)}</li>" for f in flags)
        flag_html = f'<p class="note"><strong>Caveats:</strong></p><ul>{items}</ul>'
    return intro + table + tail + flag_html


def _read_csv_rows(path: Path) -> list[dict]:
    if not path.exists():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def _rank1_row(rows: list[dict]) -> dict:
    best, best_rank = {}, None
    for row in rows:
        r = _num(row.get("rank"))
        r = 1e9 if r is None else r
        if best_rank is None or r < best_rank:
            best_rank, best = r, row
    return best


def read_fields_from_sidecars(data_dir: Path, slug: str) -> dict | None:
    """Build confidence ``fields`` from a report's on-disk sidecars.

    Prefers the enriched ``<slug>_confidence_inputs.json`` (written by
    tools/enrich_report_evidence.py with observed m/z + MS2 folded in); otherwise
    derives what it can from the base sidecars, taking the MS2 score from the
    pipeline spectral columns in ``<slug>_top_candidates.csv``.
    """
    enriched = data_dir / f"{slug}_confidence_inputs.json"
    if enriched.exists():
        try:
            return json.loads(enriched.read_text())
        except (ValueError, OSError):
            pass

    css = _read_csv_rows(data_dir / f"{slug}_candidate_score_summary.csv")
    pw = _read_csv_rows(data_dir / f"{slug}_pathway_prediction.csv")
    if not css and not pw:
        return None
    css_row = css[0] if css else {}
    css_next = css[1] if len(css) > 1 else None
    pw_row = pw[0] if pw else {}

    observed = False
    ev = data_dir / f"{slug}_evidence_summary.csv"
    if ev.exists():
        evd = {r.get("field"): r.get("value") for r in _read_csv_rows(ev)}
        observed = _num(evd.get("parent_mz")) is not None

    ms2_score = None
    ms2_present = False
    tc = _read_csv_rows(data_dir / f"{slug}_top_candidates.csv")
    if tc:
        row = _rank1_row(tc)
        cands = [
            _num(row.get("predicted_spectrum_score_v2")),
            _num(row.get("forward_spectrum_score")),
            _num(row.get("ms2_match_score")),
        ]
        cands = [c for c in cands if c is not None]
        if cands:
            ms2_score = max(cands)
            ms2_present = ms2_score > 0

    neigh = 0
    for row in _read_csv_rows(data_dir / "similarity_neighbors.csv"):
        s = _num(row.get("similarity_to_target"))
        if s is not None and s >= 0.50:
            neigh += 1

    correlations_real = True
    status = data_dir / f"{slug}_top_correlating_metabolites_status.json"
    if status.exists():
        try:
            blocker = json.loads(status.read_text()).get("blocker")
            correlations_real = not blocker
        except (ValueError, OSError):
            correlations_real = True

    return fields_from_parts(
        css_row,
        css_next,
        pw_row,
        observed_mz_available=observed,
        ms2_spectrum_present=ms2_present,
        ms2_score=ms2_score,
        library_neighbors_ge_050=neigh,
        correlations_real=correlations_real,
    )
