"""Hybrid identification-confidence model for X-compound reports.

Pure standard library so the renderer and the delivered-set patcher share one
implementation and produce identical output.

Design (see docs/report_presentation_gap_findings.md for the evidence behind it):

- Confidence is built from up to four evidence axes, each a 0..1 sub-score with
  an availability flag: structure (nearest-reference Tanimoto / exact catalog
  match), mass (candidate vs target neutral-mass fit), MS2 (fragment/neutral-loss
  spectral fit), and pathway (neighbor pathway-vote consensus).
- ``overall`` is an exploratory Platt-scaled evidence score fitted on the 92
  truth-verified controls (tools/calibrate_confidence.py). Publication reanalysis
  shows weak leave-one-out calibration (see tools/audit_structure_predictions.py),
  so it must not be interpreted as an externally validated absolute probability.
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
from typing import Any

# --- calibration constants (from tools/calibrate_confidence.py) --------------
WEIGHTS = {"structure": 0.40, "mass": 0.30, "ms2": 0.20, "pathway": 0.10}
# Platt scaling of the completeness-weighted raw evidence score -> P(correct).
# Fit in-sample by tools/calibrate_confidence.py on the 92 controls with known
# structures.  The delivered evidence is weakly discriminative (MS2 coverage
# 43.6%); leave-one-out publication validation shows poor absolute calibration.
PLATT_A = 1.887084
PLATT_B = -1.520154
CALIBRATION_N = 92          # truth-verified controls used to fit/validate
CALIBRATION_BASE_RATE = 0.478
# Operational tier cutoffs on the in-sample Platt score. They were monotonic on
# the small truth set but are not externally validated; publication-facing text
# therefore labels the aggregate score exploratory.
TIER_HIGH = 0.65
TIER_MODERATE = 0.48
TIER_LOW = 0.38
# Policy floor for authentic-standard-confirmed reports (MSI Level 1). This is a
# display value, not a recalibrated probability: a run authentic standard is the
# gold standard, so it reads 100%. Curator-known-identity promotions are clamped
# strictly below this (see IDENTITY_PROMOTED_MAX) to preserve the ordering
# authentic-standard > known-identity-pending.
CONFIRMED_OVERALL = 1.0

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

# --- Library identity ("authentic-standard pending") metric ------------------
# Separate, clearly-labeled from the calibrated ``overall`` P(correct). Reflects
# a library / database identity assessment for compounds whose top structure is
# a named reference-DB metabolite with corroborating accurate mass (and, when
# available, observed MS2). It is a curator-driven rubric value, not a
# validated probability. Confirmation-by-authentic-standard remains pending.
IDENTITY_STATUS_KNOWN = "known_pending"
IDENTITY_STATUS_PUTATIVE = "putative"
IDENTITY_STATUS_UNKNOWN = "unknown"

# Auto-rule: promote every named exact-match library metabolite with strong
# mass corroboration. OFF by default so behavior stays curator-driven; when
# ON the auto-rule uses IDENTITY_CONF_AUTO_MS2 / _NO_MS2 below.
IDENTITY_AUTO_RULE = False

# Rubric confidence values for the library-identity metric (policy, not P(correct)).
IDENTITY_CONF_OVERRIDE_DEFAULT = 0.98  # allowlist override without a per-row value
IDENTITY_CONF_AUTO_MS2 = 0.95          # auto-promoted, observed MS2 corroborates
IDENTITY_CONF_AUTO_NO_MS2 = 0.92       # auto-promoted, no observed MS2

# Upper bound on the overall headline value when a curator known-identity
# promotion applies. 100% is reserved for authentic-standard-confirmed reports
# (see CONFIRMED_OVERALL); a curator override can never reach the same ceiling.
IDENTITY_PROMOTED_MAX = 0.99


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
    raw_evidence: float = 0.0
    axes: list[Axis] = field(default_factory=list)
    flags: list[str] = field(default_factory=list)
    source: str = ""
    # display extras
    observed_mz: float | None = None
    observed_neutral_mass: float | None = None
    candidate_exact_mass: float | None = None
    ms2_n_pairs: int = 0
    metabolon_level: str = ""
    mass_alt: dict | None = None
    msi_level: int = 5
    msi_label: str = ""
    msi_rationale: str = ""
    msi_sublevel: str = ""          # derived fine sublevel (only 3A is defensible)
    metabolon_sublevel: str = ""    # curated known_X SubLevel (authoritative)
    # Library identity ("authentic-standard pending") metric - separate from ``overall``.
    identity_status: str = IDENTITY_STATUS_UNKNOWN
    identity_confidence: float | None = None
    identity_basis: str = ""

    def as_dict(self) -> dict:
        return {
            "overall": self.overall,
            "tier": self.tier,
            "tier_note": self.tier_note,
            "completeness": self.completeness,
            "descriptor": self.descriptor,
            "raw_evidence": self.raw_evidence,
            "source": self.source,
            "observed_mz": self.observed_mz,
            "observed_neutral_mass": self.observed_neutral_mass,
            "candidate_exact_mass": self.candidate_exact_mass,
            "ms2_n_pairs": self.ms2_n_pairs,
            "metabolon_level": self.metabolon_level,
            "mass_alt": self.mass_alt,
            "msi_level": self.msi_level,
            "msi_label": self.msi_label,
            "msi_rationale": self.msi_rationale,
            "msi_sublevel": self.msi_sublevel,
            "metabolon_sublevel": self.metabolon_sublevel,
            "identity_status": self.identity_status,
            "identity_confidence": self.identity_confidence,
            "identity_basis": self.identity_basis,
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
    observed_neutral_mass = _num(fields.get("observed_neutral_mass"))
    candidate_exact_mass = _num(fields.get("candidate_exact_mass"))
    mass_available = mass is not None
    if not mass_available:
        m_detail = "No mass-match score available"
    elif mass >= MASS_STRONG:
        if (
            observed_neutral_mass is not None
            and candidate_exact_mass is not None
        ):
            m_detail = (
                f"Strong: observed neutral mass {observed_neutral_mass:.6f} Da "
                f"vs predicted {candidate_exact_mass:.6f} Da"
            )
            if abs_err is not None:
                m_detail += f" (absolute error {abs_err:.6f} Da)"
        else:
            m_detail = "Strong: candidate formula/mass matches the target neutral mass"
            if abs_err is not None:
                m_detail += f" (absolute error {abs_err:.6f} Da)"
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
    """Tier from the exploratory Platt score, with a mass-mismatch override.

    The tiers were monotonic on the small in-sample truth set but are not externally
    calibrated. A top candidate whose mass disagrees with the observed feature is
    always demoted, regardless of score.
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


MSI_LABELS = {
    1: "Confirmed (authentic standard)",
    2: "Probable structure",
    3: "Tentative candidate structure",
    4: "Molecular formula",
    5: "Exact mass of interest",
}


# Curated Metabolon SubLevel tokens (authoritative; surfaced as-is).
CURATED_SUBLEVELS = {"1", "2A", "2B", "2C", "3A", "3B", "3C", "3D", "4", "4?", "5", "KNOWN", "Released"}

# Derived fine-sublevel labels (only the defensible subset; see
# docs/msi_sublevel_findings.md). Fine 2C/3B/3C/3D and 2A/2B are NOT derived
# because current evidence neither reproduces the curated letters (9.5% exact)
# nor orders them by correctness; they require observed RT, a reference MS2
# spectral library, and diagnostic-fragment matching (upstream handoff).
MSI_SUBLEVEL_LABELS = {
    "3A": "Positional-isomer candidates (multiple same-formula structures)",
}


def curated_sublevel(fields: dict) -> str:
    """Authoritative Metabolon SubLevel from known_X, or '' if absent/non-metabolite."""
    curated = str(fields.get("metabolon_sublevel", "") or "").strip()
    return curated if curated in CURATED_SUBLEVELS else ""


def assign_msi_sublevel(fields: dict, conf: "Confidence") -> str:
    """Derived fine sublevel for our evidence-based MSI level.

    Only ``3A`` (positional-isomer candidates) is assigned, since same-formula
    multiplicity is the one fine distinction our evidence supports defensibly. The
    curated Metabolon SubLevel is surfaced separately as an authoritative
    cross-reference; it is not overwritten by this derived value.
    """
    n_isomers = _num(fields.get("n_isomers")) or 0
    if conf.msi_level == 3 and n_isomers >= 2:
        return "3A"
    return ""


def assign_msi_level(fields: dict, conf: "Confidence") -> tuple[int, str, str]:
    """Assign a Schymanski/MSI identification-confidence level (1-5) from evidence.

    L1 confirmed (authentic standard) / L2 probable structure (strong structure +
    observed MS2 + mass) / L3 tentative candidate structure (hypothesis) / L4
    molecular formula only (structure not supported or mass-rejected) / L5 exact
    mass only. Structures are only claimed at L1-L3; L4-L5 report formula/mass.
    """
    by = {a.key: a for a in conf.axes}
    struct = by["structure"].score or 0.0
    mass = by["mass"].score if by["mass"].available else None
    ms2_ok = by["ms2"].available and (by["ms2"].score or 0.0) >= MS2_PRESENT
    mass_ok = mass is not None and mass >= MASS_STRONG
    mass_mismatch = mass is not None and mass < MASS_WEAK
    observed = bool(fields.get("observed_mz_available"))

    if conf.tier == "Confirmed":
        return 1, MSI_LABELS[1], "Metabolon Level-1 standard match and the top candidate matches it"
    if struct >= STRUCT_MODERATE and ms2_ok and mass_ok:
        return 2, MSI_LABELS[2], "Strong structure similarity with observed-MS2 and mass corroboration"
    if (not mass_mismatch) and struct >= STRUCT_WEAK:
        return 3, MSI_LABELS[3], "Structure candidate consistent with mass but not spectrally confirmed - a hypothesis"
    if observed or mass_ok:
        note = (
            "Observed mass is known but the top candidate mass does not fit - structure rejected"
            if mass_mismatch
            else "Observed mass/formula is supported but no structure is adequately supported"
        )
        return 4, MSI_LABELS[4], note
    return 5, MSI_LABELS[5], "Only an exact mass is available; no confident formula or structure"


def assign_identity(fields: dict, conf: "Confidence") -> tuple[str, float | None, str]:
    """Assign the library-identity ("authentic-standard pending") status.

    This is DISTINCT from ``conf.overall`` (the calibrated P(top structure
    correct)). It reflects a library / database identity claim:

      - ``known_pending`` = top structure is a named reference-DB metabolite
        with corroborating accurate mass (and, if available, observed MS2);
        confirmation by authentic standard is pending. Reached either by an
        explicit curator override (``fields["known_identity_override"]``) or by
        the auto-rule when ``IDENTITY_AUTO_RULE`` is enabled.
      - ``putative`` = library-plausible but does not meet the corroboration bar
        (weak mass, mass mismatch, or unresolved name).
      - ``unknown`` = de novo / no library correspondence.

    Confidence values are a curator rubric, NOT a validated probability.
    """
    source = str(fields.get("source", "")).strip()
    name_resolved = bool(fields.get("name_resolved"))
    by = {a.key: a for a in conf.axes}
    mass = by["mass"].score if by["mass"].available else None
    mass_ok = mass is not None and mass >= MASS_STRONG
    mass_mismatch = mass is not None and mass < MASS_WEAK
    ms2_ok = by["ms2"].available and (by["ms2"].score or 0.0) >= MS2_PRESENT

    override = fields.get("known_identity_override")
    if override:
        # Accept truthy scalar (bool) or a dict {confidence, reason}
        if isinstance(override, dict):
            cval = _num(override.get("confidence")) or _num(override.get("identity_confidence"))
            reason = str(override.get("reason") or "").strip()
        else:
            cval = None
            reason = ""
        identity_confidence = cval if cval is not None else IDENTITY_CONF_OVERRIDE_DEFAULT
        basis_bits = ["Curated known identity (authentic standard pending)"]
        if reason:
            basis_bits.append(reason)
        if name_resolved:
            basis_bits.append("named reference-DB metabolite")
        if mass_ok:
            basis_bits.append("accurate-mass corroborated")
        return IDENTITY_STATUS_KNOWN, round(float(identity_confidence), 4), "; ".join(basis_bits)

    auto_named_exact = (
        IDENTITY_AUTO_RULE
        and source == "exact"
        and name_resolved
        and mass_ok
        and not mass_mismatch
    )
    if auto_named_exact:
        cval = IDENTITY_CONF_AUTO_MS2 if ms2_ok else IDENTITY_CONF_AUTO_NO_MS2
        basis = (
            "Named reference-library metabolite with accurate-mass corroboration"
            + (" and observed-MS2 fit" if ms2_ok else "")
            + "; authentic-standard confirmation pending"
        )
        return IDENTITY_STATUS_KNOWN, round(float(cval), 4), basis

    if source == "predicted":
        return IDENTITY_STATUS_UNKNOWN, None, "De novo / predicted structure; no library correspondence"

    # source == "exact" or other library-adjacent but not meeting the bar
    if source == "exact":
        if mass_mismatch:
            basis = "Library-match structure but mass disagrees with the observed feature"
        elif not mass_ok:
            basis = "Library-match structure with partial mass corroboration; identity not established"
        elif not name_resolved:
            basis = "Library-match structure but no canonical name resolved from the reference DB"
        else:
            basis = "Library-match structure; identity not curator-confirmed"
        return IDENTITY_STATUS_PUTATIVE, None, basis
    return IDENTITY_STATUS_UNKNOWN, None, "No library correspondence for this feature"


def compute_confidence(fields: dict) -> Confidence:
    """Compute the hybrid identification confidence from normalized ``fields``."""
    source = str(fields.get("source", "")).strip()
    axes = build_axes(fields)
    raw, completeness = raw_score(axes)
    overall = _sigmoid(PLATT_A * raw + PLATT_B)
    tier, tier_note = _tier(overall, axes, source)
    # Authentic-standard override (gap 6): a Metabolon Level-1 control whose top
    # candidate matches the confirmed structure is genuinely confirmable. Display
    # policy: an authentic-standard-confirmed identification is the gold standard,
    # so the headline overall reads 100% (CONFIRMED_OVERALL). The underlying
    # exploratory calibrated evidence score is still shown in the scored-parameters
    # section text; this is a policy display value, not a recalibrated probability.
    if str(fields.get("metabolon_level", "")).strip() == "1" and bool(fields.get("truth_match")):
        tier = "Confirmed"
        tier_note = "Authentic-standard confirmed (Metabolon Level 1); top candidate matches the confirmed structure"
        overall = max(overall, CONFIRMED_OVERALL)
    if source == "exact":
        descriptor = "Cataloged reference metabolite (exact structural match)"
    elif source == "predicted":
        descriptor = "De novo / predicted structure"
    else:
        descriptor = "Structure proposal"
    result = Confidence(
        overall=round(overall, 4),
        tier=tier,
        tier_note=tier_note,
        completeness=round(completeness, 4),
        descriptor=descriptor,
        raw_evidence=round(raw, 4),
        axes=axes,
        flags=_flags(axes, fields),
        source=source,
        observed_mz=_num(fields.get("observed_mz")),
        observed_neutral_mass=_num(fields.get("observed_neutral_mass")),
        candidate_exact_mass=_num(fields.get("candidate_exact_mass")),
        ms2_n_pairs=int(_num(fields.get("ms2_n_pairs")) or 0),
        metabolon_level=str(fields.get("metabolon_level", "") or ""),
        mass_alt=fields.get("mass_alt") if isinstance(fields.get("mass_alt"), dict) else None,
        metabolon_sublevel=curated_sublevel(fields),
    )
    result.msi_level, result.msi_label, result.msi_rationale = assign_msi_level(fields, result)
    result.msi_sublevel = assign_msi_sublevel(fields, result)
    result.identity_status, result.identity_confidence, result.identity_basis = assign_identity(fields, result)
    # Curator-known library identity: promote the headline overall confidence to the
    # identity rubric value (>= 90%) so reports reflect "we know what it is,
    # authentic standard pending" without upgrading to MSI Level 1 / Confirmed.
    if (
        result.identity_status == IDENTITY_STATUS_KNOWN
        and result.identity_confidence is not None
        and result.identity_confidence >= 0.90
        and result.tier != "Confirmed"
    ):
        calibrated = result.overall
        # Clamp below the authentic-standard-confirmed floor so a curator override
        # can never match the Level 1 headline (100%).
        promoted = min(float(result.identity_confidence), IDENTITY_PROMOTED_MAX)
        result.overall = round(promoted, 4)
        result.tier = "High"
        result.tier_note = (
            "This is a High Confidence probable structure based on corroborating "
            "evidence derived from Metabolon\u2019s validated structural elucidation workflow."
        )
        if calibrated < result.overall - 0.05:
            result.flags.append(
                f"Overall confidence reflects curator-known library identity "
                f"(underlying calibrated evidence score was {round(100 * calibrated)}%)."
            )
    return result


def load_identity_overrides(path: Path | str | None) -> dict[str, dict]:
    """Load the curated known-identity allowlist CSV.

    Columns (case-insensitive): ``control_id`` (required),
    ``identity_confidence`` (optional; defaults to ``IDENTITY_CONF_OVERRIDE_DEFAULT``),
    ``reason`` (optional; free-text rationale). Returns a mapping
    ``{control_id: {"confidence": float, "reason": str}}``. Missing / unreadable
    files return an empty dict so callers can pass in an optional path.
    """
    if path is None:
        return {}
    p = Path(path)
    if not p.exists():
        return {}
    out: dict[str, dict] = {}
    try:
        with p.open(newline="") as handle:
            # Skip full-line comments starting with '#' (they must not appear inside quoted fields).
            filtered = (line for line in handle if not line.lstrip().startswith("#"))
            reader = csv.DictReader(filtered)
            for row in reader:
                cid = str(row.get("control_id") or row.get("Control_ID") or "").strip()
                if not cid or cid.startswith("#"):
                    continue
                cval = _num(row.get("identity_confidence") or row.get("confidence"))
                reason = str(row.get("reason") or row.get("note") or "").strip()
                out[cid] = {
                    "confidence": cval if cval is not None else IDENTITY_CONF_OVERRIDE_DEFAULT,
                    "reason": reason,
                }
    except (OSError, ValueError):
        return {}
    return out


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
# may show the aggregate MS2 fit that contributes to confidence, but never raw
# ion-level or per-fragment detail.

_TIER_CLASS = {
    "Confirmed": "badge ok", "High": "badge ok", "Moderate": "badge ok",
    "Low": "badge", "Tentative": "badge",
}


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


def msi_badge_html(conf: Confidence) -> str:
    cls = "badge ok" if conf.msi_level <= 2 else "badge"
    level_str = conf.msi_sublevel or str(conf.msi_level)
    return f'<span class="{cls}">MSI Level {level_str} \u00b7 {html.escape(conf.msi_label)}</span>'


def msi_metric_card_html(conf: Confidence) -> str:
    level_str = conf.msi_sublevel or str(conf.msi_level)
    note = conf.msi_label
    if conf.msi_sublevel:
        note += " - " + MSI_SUBLEVEL_LABELS.get(conf.msi_sublevel, "")
    return (
        '<div class="metric-card"><div class="metric-label">MSI identification level</div>'
        f'<div class="metric-value">Level {level_str}</div>'
        f'<div class="metric-note">{html.escape(note)}</div></div>'
    )


# --- Library-identity ("authentic-standard pending") HTML builders -----------
# These render a SEPARATE metric from ``overall`` (which stays untouched). The
# metric expresses a curator-driven library / database identity claim and is
# clearly labeled as such - it is NOT a validated P(structure correct).

_IDENTITY_LABEL = {
    IDENTITY_STATUS_KNOWN: "Known - authentic standard pending",
    IDENTITY_STATUS_PUTATIVE: "Putative (library match)",
    IDENTITY_STATUS_UNKNOWN: "Not established",
}


def _identity_pct(conf: Confidence) -> str:
    if conf.identity_confidence is None:
        return "not established"
    return f"{round(100 * conf.identity_confidence)}%"


def identity_badge_html(conf: Confidence) -> str:
    """Header badge for the library-identity metric.

    Only shown when the status is ``known_pending`` (>=90% rubric value).
    """
    if conf.identity_status != IDENTITY_STATUS_KNOWN:
        return ""
    return (
        f'<span class="badge ok">Identity: Known \u00b7 standard pending \u00b7 '
        f"{_identity_pct(conf)}</span>"
    )


def identity_metric_card_html(conf: Confidence) -> str:
    """Metric card for the library-identity metric.

    Emits a card for all statuses so the reader always sees an explicit answer;
    only the ``known_pending`` case reports a numeric confidence.
    """
    if conf.identity_status == IDENTITY_STATUS_KNOWN:
        return (
            '<div class="metric-card"><div class="metric-label">Library identity status</div>'
            '<div class="metric-value">High Confidence ID, Authentic Standard Confirmation '
            "Pending to Promote to Level 1 Library</div></div>"
        )
    label = _IDENTITY_LABEL.get(conf.identity_status, "Not established")
    return (
        '<div class="metric-card"><div class="metric-label">Library identity status</div>'
        f'<div class="metric-value">{html.escape(label)}</div>'
        f'<div class="metric-note">{html.escape(conf.identity_basis or "")}</div></div>'
    )


def candidate_pathway_agreement_note(
    candidates: list[dict[str, Any]],
    expected: Any,
    field: str,
    label: str,
) -> str:
    """Summarize agreement among up to ten ranked candidate pathway calls."""

    expected_text = str(expected or "").strip()
    if not expected_text:
        return ""
    top = candidates[:10]
    if not top:
        return ""
    matches = sum(
        str(row.get(field) or "").strip().casefold() == expected_text.casefold()
        for row in top
    )
    share = matches / len(top)
    level = "High" if share >= 0.80 else ("Moderate" if share >= 0.60 else "Low")
    return (
        f"{level} Confidence: {matches}/{len(top)} top candidate structures "
        f"map to this {label}"
    )


def identification_parameter_cards_html(
    confidence: Confidence,
    *,
    candidate_scores: list[dict[str, Any]],
    pathway_info: dict[str, Any],
    top_similarity: float,
    neighbor_count: int,
) -> str:
    """Render report-level evidence summaries from shared persisted inputs."""

    axes = {axis.key: axis for axis in confidence.axes}
    mass_axis = axes.get("mass")
    fragment_axis = axes.get("ms2")
    if mass_axis is not None and mass_axis.available and mass_axis.score is not None:
        mass_strength = (
            "Strong"
            if mass_axis.score >= MASS_STRONG
            else ("Moderate" if mass_axis.score >= MASS_WEAK else "Weak")
        )
        mass_score = f"{mass_strength} \N{MIDDLE DOT} {mass_axis.score:.2f}"
    else:
        mass_score = "Not available"
    mass_note = mass_axis.detail if mass_axis is not None else "Mass evidence was not available"
    fragment_score = (
        f"{fragment_axis.score:.2f}"
        if fragment_axis is not None and fragment_axis.available and fragment_axis.score is not None
        else "Not available"
    )
    fragment_note = (
        fragment_axis.detail
        if fragment_axis is not None
        else "Aggregate fragment-match evidence was not available"
    )
    pathway_super = str(pathway_info.get("super") or "").strip()
    pathway_sub = str(pathway_info.get("sub") or "").strip()
    super_note = candidate_pathway_agreement_note(
        candidate_scores, pathway_super, "predicted_pathway", "Pathway"
    )
    sub_note = candidate_pathway_agreement_note(
        candidate_scores, pathway_sub, "predicted_subpathway", "SubPathway"
    )
    if top_similarity >= 0.999:
        similarity_note = (
            "The top candidate exactly matches a Metabolon library standard, "
            "enabling a higher confidence structural candidate"
        )
    elif top_similarity > 0:
        similarity_note = (
            f"The top candidate’s nearest Metabolon library-standard similarity is "
            f"{top_similarity:.2f}; closer matches support higher confidence"
        )
    else:
        similarity_note = "No Metabolon library-standard structure match was available"
    neighbor_word = "standard" if neighbor_count == 1 else "standards"
    neighbor_note = (
        f"This metabolite is similar to {neighbor_count} Metabolon library {neighbor_word}, "
        "enabling a higher confidence structural candidate"
        if neighbor_count
        else "No Metabolon library standards were above the similarity floor"
    )

    def card(label: str, value: Any, note: str) -> str:
        return (
            '<div class="metric-card">'
            f'<div class="metric-label">{html.escape(label)}</div>'
            f'<div class="metric-value">{html.escape(str(value))}</div>'
            f'<div class="metric-note">{html.escape(note)}</div></div>'
        )

    cards = [
        card("Mass Match", mass_score, mass_note),
        card("Library Neighbors", neighbor_count, neighbor_note),
        card("Top Candidate Structural Similarity", f"{top_similarity:.2f}", similarity_note),
        card(
            "Candidate Chemical Pathway",
            pathway_super or "Not available",
            super_note or "Candidate-level pathway agreement was not available",
        ),
        card(
            "Candidate Chemical SubPathway",
            pathway_sub or "Not available",
            sub_note or "Candidate-level SubPathway agreement was not available",
        ),
        card("Fragment Match Score", fragment_score, fragment_note),
    ]
    return (
        '<!--ID_PARAMETER_CARDS--><div class="metrics">'
        + "".join(cards)
        + "</div><!--/ID_PARAMETER_CARDS-->"
    )


def identity_note_html(conf: Confidence, *, public: bool) -> str:
    """Explanatory paragraph placed alongside the scored-parameters section.

    Publicly, the note must not mention observed MS2 (mirrors the public/internal
    split used elsewhere in this module). The value shown is explicitly a curator
    rubric, not the calibrated structure-prediction probability.
    """
    if conf.identity_status == IDENTITY_STATUS_KNOWN:
        pieces = [
            "<strong>Library identity:</strong> ",
            html.escape(_IDENTITY_LABEL[IDENTITY_STATUS_KNOWN]),
            f" ({_identity_pct(conf)}). ",
            "This reflects a named reference-library match with corroborating "
            "accurate mass",
        ]
        if not public and conf.ms2_n_pairs:
            pieces.append(" and observed-MS2 fit")
        pieces.append(
            "; confirmation by authentic standard is pending. "
            "This is a library/database identity assessment (curator rubric), "
            "distinct from the calibrated structure-prediction confidence shown above."
        )
        body = "".join(pieces)
    elif conf.identity_status == IDENTITY_STATUS_PUTATIVE:
        body = (
            "<strong>Library identity:</strong> "
            + html.escape(_IDENTITY_LABEL[IDENTITY_STATUS_PUTATIVE])
            + ". "
            + html.escape(conf.identity_basis or "")
        )
    else:
        body = (
            "<strong>Library identity:</strong> "
            + html.escape(_IDENTITY_LABEL[IDENTITY_STATUS_UNKNOWN])
            + ". "
            + html.escape(conf.identity_basis or "")
        )
    return f'<p class="note">{body}</p>'


def scored_parameters_section_html(conf: Confidence, *, public: bool) -> str:
    """Body HTML for the 'Parameters Contributing to ID Confidence' section."""
    rows = []
    for ax in conf.axes:
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
    level_str = conf.msi_sublevel or str(conf.msi_level)
    sublevel_note = ""
    if conf.msi_sublevel:
        sublevel_note = f" Sub-level {conf.msi_sublevel}: {MSI_SUBLEVEL_LABELS.get(conf.msi_sublevel, '')}."
    msi_line = (
        f'<p class="note"><strong>MSI identification level: Level {level_str} - '
        f"{html.escape(conf.msi_label)}.</strong> {html.escape(conf.msi_rationale)}.{html.escape(sublevel_note)} "
        "Levels 1-3 present a structure (Level 3 as a hypothesis); Levels 4-5 report formula/mass "
        "without a specific structure.</p>"
    )
    intro = msi_line + (
        '<p class="note">The overall identification confidence is an exploratory evidence score '
        "combining structural similarity, mass/formula match, MS2 spectral fit, and pathway "
        f"consensus (fit in-sample on {CALIBRATION_N} structure-verified controls; "
        f"base rate {CALIBRATION_BASE_RATE:.0%}; not externally calibrated). "
        f"Proposed structure: <strong>{html.escape(conf.descriptor)}</strong>. "
        f"Assessment: <strong>{html.escape(conf.tier)}</strong> - "
        f"{html.escape(conf.tier_note.rstrip('.'))}.</p>"
    ) + identity_note_html(conf, public=public)
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
    extra = ""
    if conf.metabolon_level or conf.metabolon_sublevel:
        sub = f" (SubLevel {html.escape(conf.metabolon_sublevel)})" if conf.metabolon_sublevel else ""
        extra += (
            f'<p class="note"><strong>Metabolon identification Level (curated):</strong> '
            f"{html.escape(conf.metabolon_level or 'n/a')}{sub}"
            + (" - authentic-standard confirmed" if conf.metabolon_level == "1" else "")
            + ". This is Metabolon's authoritative curated assignment.</p>"
        )
    if conf.mass_alt:
        alt = conf.mass_alt
        rank = alt.get("rank")
        rank_txt = f"rank {int(rank)}" if isinstance(rank, (int, float)) else "a lower-ranked candidate"
        formula = html.escape(str(alt.get("formula", "")).strip() or "n/a")
        ms = alt.get("mass_match_score")
        err = alt.get("abs_mass_error_da")
        fit_bits = []
        if isinstance(ms, (int, float)):
            fit_bits.append(f"mass match {ms:.3f}")
        if isinstance(err, (int, float)):
            fit_bits.append(f"\u0394 {err:.4g} Da")
        fit_txt = (" (" + ", ".join(fit_bits) + ")") if fit_bits else ""
        smiles = str(alt.get("smiles", "")).strip()
        smiles_txt = f" SMILES <code>{html.escape(smiles)}</code>." if smiles else ""
        extra += (
            '<p class="note"><strong>Mass-consistent alternative (recoverable structure):</strong> '
            "the top candidate's mass does not fit the observed feature, but "
            f"{rank_txt} \u2014 formula <strong>{formula}</strong>{fit_txt} \u2014 fits the observed mass."
            f"{smiles_txt} A mass-aware reranker would promote this candidate; treat it as the "
            "stronger mass-based structure hypothesis for this feature.</p>"
        )
    flags = [f for f in conf.flags if not (public and "MS2" in f)]
    flag_html = ""
    if flags:
        items = "".join(f"<li>{html.escape(f)}</li>" for f in flags)
        flag_html = f'<p class="note"><strong>Caveats:</strong></p><ul>{items}</ul>'
    return intro + table + tail + extra + flag_html


def ms2_evidence_section_html(fields: dict, *, public: bool) -> str:
    """Body HTML for the MS Fit / Feature Match section.

    Public: no per-fragment/ion detail (compliance). Internal: pipeline spectral
    scores + observed peak counts (aggregate; still not raw ion m/z, which is not
    decodably available)."""
    d = fields.get("ms2_detail") or {}
    present = bool(d.get("present"))
    if public:
        return (
            '<p class="note">MS2 spectral fit is incorporated into the overall identification '
            "confidence. Per-fragment spectral detail is provided in the internal report.</p>"
        )
    if not present:
        return '<p class="note">No observed MS2 spectrum was available for this feature.</p>'

    def s(v):
        return f"{float(v):.3f}" if isinstance(v, (int, float)) else "n/a"

    return (
        '<p class="note">Observed MS2 spectral evidence for this feature (internal only). '
        "Scores are the pipeline\'s spectral match for the top candidate; peak counts describe "
        "the observed spectrum.</p>"
        '<div class="table-wrap"><table><thead><tr><th>Metric</th><th>Value</th></tr></thead><tbody>'
        f"<tr><td>Forward spectral score</td><td>{s(d.get('forward'))}</td></tr>"
        f"<tr><td>Predicted-spectrum score (v2)</td><td>{s(d.get('predicted_v2'))}</td></tr>"
        f"<tr><td>MS2 match score</td><td>{s(d.get('match'))}</td></tr>"
        f"<tr><td>Observed fragment peaks</td><td>{int(d.get('n_pairs') or 0)}</td></tr>"
        f"<tr><td>Observed fragment-ion bins</td><td>{int(d.get('n_ion_bins') or 0)}</td></tr>"
        f"<tr><td>Observed neutral-loss bins</td><td>{int(d.get('n_loss_bins') or 0)}</td></tr>"
        "</tbody></table></div>"
    )


def _present_field(value: Any) -> bool:
    if value is None:
        return False
    text = str(value).strip()
    return text.lower() not in {"", "nan", "none"}


def _score_display(value: Any) -> str:
    num = _num(value)
    if num is not None:
        return f"{num:.3f}"
    return html.escape(str(value).strip())


def _candidate_ms2_fit(row: dict[str, Any]) -> float | None:
    scores = [
        _num(row.get(key))
        for key in ("forward_spectrum_score", "predicted_spectrum_score_v2", "ms2_match_score")
    ]
    scores = [score for score in scores if score is not None]
    return max(scores) if scores else None


def _candidate_label(row: dict[str, Any]) -> str:
    formula = str(row.get("predicted_formula") or "").strip() or "n/a"
    smiles = str(row.get("predicted_smiles") or "").strip()
    if not smiles:
        return formula
    short = smiles if len(smiles) <= 40 else f"{smiles[:40]}…"
    return f"{formula} ({short})"


def diagnostic_substructure_rows_html(rank1: dict[str, Any], evidence: dict[str, Any] | None = None) -> str:
    """Render diagnostic substructure / biotransformation rows for the Fragment Scores table."""

    evidence = evidence or {}
    rows: list[tuple[str, Any]] = []

    lactoyl = rank1.get("lactoyl_hypothesis_score")
    if not _present_field(lactoyl):
        lactoyl = evidence.get("lactoyl_hypothesis_score")
    if _present_field(lactoyl):
        rows.append(("Lactoyl", lactoyl))

    adduct_label = str(rank1.get("adduct_label") or "").strip()
    adduct_score = rank1.get("adduct_score") or rank1.get("nce__score_adduct")
    if _present_field(adduct_score):
        label = f"Adduct ({adduct_label})" if adduct_label else "Adduct"
        rows.append((label, adduct_score))

    transform = str(rank1.get("transform_name") or rank1.get("conjugate_name") or "").strip()
    bio_score = rank1.get("biotransform_v2_score")
    if transform and _present_field(bio_score):
        rows.append((f"Biotransformation ({transform})", bio_score))
    elif transform:
        rows.append((f"Biotransformation ({transform})", "supported"))
    elif _present_field(bio_score):
        rows.append(("Biotransformation", bio_score))

    recursive = rank1.get("recursive__recursive_substructure_weighted_score")
    if _present_field(recursive):
        rows.append(("Recursive substructure", recursive))

    lib_overlap = rank1.get("library_fragment_overlap_score")
    if _present_field(lib_overlap):
        rows.append(("Library fragment overlap", lib_overlap))

    neutral_loss = rank1.get("spectrum_neutral_loss_score")
    if _present_field(neutral_loss):
        rows.append(("Neutral-loss fit", neutral_loss))

    if not rows:
        body = (
            '<tr><td colspan="2">No diagnostic substructure hypotheses were scored for this feature.</td></tr>'
        )
    else:
        body = "\n".join(
            f"<tr><td>{html.escape(label)}</td><td>{_score_display(score)}</td>"
            for label, score in rows
        )
    return f"<!--FRAGMENT_DIAGNOSTICS-->{body}<!--/FRAGMENT_DIAGNOSTICS-->"


def fragment_summary_section_html(
    fields: dict[str, Any],
    candidate_rows: list[dict[str, Any]],
    *,
    public: bool,
    top_n: int = 3,
) -> str:
    """Non-raw MS2 fragment summary for the MS Fit / Feature Match section."""

    detail = fields.get("ms2_detail") or {}
    present = bool(detail.get("present"))
    parts: list[str] = []

    if present:
        profile_rows = [
            ("Observed fragment peaks", int(detail.get("n_pairs") or 0)),
            ("Fragment-ion bins", int(detail.get("n_ion_bins") or 0)),
            ("Neutral-loss bins", int(detail.get("n_loss_bins") or 0)),
        ]
        profile_body = "".join(
            f"<tr><td>{html.escape(label)}</td><td>{value}</td></tr>" for label, value in profile_rows
        )
        parts.append(
            '<p class="note">Summary of the observed MS2 spectrum (counts only; no ion m/z).</p>'
            '<div class="table-wrap"><table><thead><tr><th>Metric</th><th>Count</th></tr></thead>'
            f"<tbody>{profile_body}</tbody></table></div>"
        )
    else:
        parts.append('<p class="note">No observed MS2 spectrum was available for this feature.</p>')

    fit_by_rank = {
        int(_num(item.get("rank")) or 0): item
        for item in (fields.get("fragment_candidate_fits") or [])
        if _num(item.get("rank")) is not None
    }
    candidate_body: list[str] = []
    for row in candidate_rows[:top_n]:
        rank_num = _num(row.get("rank"))
        rank_key = int(rank_num) if rank_num is not None else 0
        fit = fit_by_rank.get(rank_key, {})
        forward = _num(row.get("forward_spectrum_score"))
        predicted_v2 = _num(row.get("predicted_spectrum_score_v2"))
        best_fit = _candidate_ms2_fit(row)
        neutral = _num(row.get("spectrum_neutral_loss_score"))
        matched_ions = fit.get("matched_ion_bins")
        matched_losses = fit.get("matched_loss_bins")
        matched_ions_txt = (
            str(int(matched_ions))
            if isinstance(matched_ions, (int, float)) and not math.isnan(float(matched_ions))
            else "n/a"
        )
        matched_losses_txt = (
            str(int(matched_losses))
            if isinstance(matched_losses, (int, float)) and not math.isnan(float(matched_losses))
            else "n/a"
        )
        candidate_body.append(
            "<tr>"
            f"<td>{html.escape(str(row.get('rank', '')))}</td>"
            f"<td>{html.escape(_candidate_label(row))}</td>"
            f"<td>{_score_display(forward) if forward is not None else 'n/a'}</td>"
            f"<td>{_score_display(predicted_v2) if predicted_v2 is not None else 'n/a'}</td>"
            f"<td>{_score_display(best_fit) if best_fit is not None else 'n/a'}</td>"
            f"<td>{_score_display(neutral) if neutral is not None else 'n/a'}</td>"
            f"<td>{html.escape(matched_ions_txt)}</td>"
            f"<td>{html.escape(matched_losses_txt)}</td>"
            "</tr>"
        )

    if candidate_body:
        parts.append(
            '<p class="note">Per-candidate spectral fit scores (0–1 fractions; higher means more of '
            "the observed spectrum is explained by the candidate). Matched bin counts are shown when "
            "available.</p>"
            '<div class="table-wrap"><table><thead><tr>'
            "<th>Rank</th><th>Candidate</th><th>Forward fit</th><th>Predicted fit (v2)</th>"
            "<th>Best fit</th><th>Neutral-loss fit</th><th>Matched ion bins</th>"
            "<th>Matched loss bins</th></tr></thead><tbody>"
            + "".join(candidate_body)
            + "</tbody></table></div>"
        )
    elif present:
        parts.append('<p class="note">Per-candidate spectral fit scores were not available for this report.</p>')

    if public:
        parts.append(
            '<p class="note">MS2 spectral fit is incorporated into the overall identification '
            "confidence. Per-ion m/z detail is provided in the internal report.</p>"
        )
    else:
        parts.append(ms2_evidence_section_html(fields, public=False))

    body = "\n      ".join(parts)
    return f"<!--FRAGMENT_SUMMARY-->{body}<!--/FRAGMENT_SUMMARY-->"


def correlations_table_html(correlations: list[dict]) -> str:
    """Table body (rows) for the 'Top Correlating Metabolites' section from WGCNA."""
    rows = []
    for c in correlations[:6]:
        score = c.get("score")
        score_txt = f"{float(score):.2f}" if isinstance(score, (int, float)) else ""
        rows.append(
            f"<tr><td>{html.escape(str(c.get('id', '')))}</td>"
            f"<td>{html.escape(str(c.get('name', '')))}</td>"
            f"<td>{score_txt}</td>"
            f"<td>{html.escape(str(c.get('source', '')))}</td></tr>"
        )
    return "".join(rows)


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


def mass_consistent_alt(css: list[dict]) -> dict | None:
    """If the rank-1 candidate mass mismatches but a lower rank fits, return that alt.

    Shared by ``tools/enrich_report_evidence.py`` (post-hoc) and
    ``compound.build_report`` (from source) so both surface the same callout.
    """
    if not css:
        return None
    r1 = _num(css[0].get("mass_match_score"))
    if r1 is None or r1 >= MASS_WEAK:
        return None
    r1_err = _num(css[0].get("abs_mass_error_da"))
    r1_formula = str(css[0].get("predicted_formula", "") or "").strip()
    for row in css[1:]:
        ms = _num(row.get("mass_match_score"))
        if ms is not None and ms >= MASS_STRONG:
            return {
                "rank": _num(row.get("rank")),
                "formula": str(row.get("predicted_formula", "") or "").strip(),
                "smiles": str(row.get("predicted_smiles", "") or "").strip(),
                "mass_match_score": ms,
                "abs_mass_error_da": _num(row.get("abs_mass_error_da")),
                "rank1_formula": r1_formula,
                "rank1_abs_mass_error_da": r1_err,
            }
    return None


# --- per-candidate display confidence (renderer + enrich + patcher + audit) --
# The top-10 candidate widget must show a per-candidate confidence that meshes
# with the 1..10 ranking. The pipeline's ``evidence_confidence`` cannot do this:
# it saturates near 1.0 and is shared across mass-consistent isomers (so every
# candidate looks identical), and it is a reranker *input feature*, not the
# ranking objective, so it is frequently higher for a lower-ranked candidate.
#
# ``candidate_confidences`` derives, per candidate, an interpretable 0..1 value
# that is monotone non-increasing with rank by construction: the rank-1 candidate
# is anchored to the report's overall identification confidence (the headline the
# reader already sees), and each lower-ranked candidate receives that ceiling
# scaled by the reranker ordering-score ratio to rank 1, with a running-minimum so
# nothing can invert. Genuinely tied candidates (equal ordering score / confidence
# within ``CANDIDATE_TIE_EPS``) are flagged so the report presents them as an
# unresolved isomer group rather than fake-distinct numbers.
CANDIDATE_TIE_EPS = 0.005


def candidate_confidences(rows: list[dict], feature_overall: float | None) -> list[dict]:
    """Monotone, rank-consistent per-candidate confidence for the top-N widget.

    Args:
        rows: rank-ordered candidate dicts, each with ``rank`` and
            ``candidate_report_score`` (= ``unified_final_score``, the reranker
            ordering objective from which ``rank`` is derived).
        feature_overall: the rank-1 ceiling in [0,1] - the report's evidence
            match strength (``Confidence.raw_evidence``: completeness-weighted
            structure/mass/MS2/pathway evidence). Raw evidence is used instead of
            the Platt ``overall`` because the latter is compressed (caps ~0.59 for
            non-authentic-standard matches) and deflates strong best matches; raw
            evidence spans a fuller range and still tracks truth (AUC ~0.79).

    Returns:
        A list aligned to ``rows``; each item is
        ``{"rank": int, "candidate_confidence": float, "tie_group": int,
        "tied": bool}`` with ``candidate_confidence`` non-increasing across rank.
    """
    ceiling = _clamp01(_num(feature_overall) or 0.0)
    n = len(rows)
    if n == 0:
        return []

    def _rank_of(row: dict, position: int) -> int:
        value = _num(row.get("rank"))
        return int(value) if value is not None else position + 1

    order = sorted(range(n), key=lambda i: _rank_of(rows[i], i))
    scores: list[float] = []
    for i in order:
        s = _num(rows[i].get("candidate_report_score"))
        scores.append(max(0.0, s) if s is not None else 0.0)
    top = scores[0] if scores else 0.0

    running = 1.0
    confidences: list[float] = []
    for s in scores:
        ratio = 1.0 if top <= 0 else _clamp01(s / top)
        running = min(running, ratio)  # enforce monotonic non-increase with rank
        confidences.append(round(_clamp01(ceiling * running), 4))

    groups: list[int] = []
    group = 0
    prev: float | None = None
    for cc in confidences:
        if prev is None or (prev - cc) > CANDIDATE_TIE_EPS:
            group += 1
        groups.append(group)
        prev = cc
    sizes: dict[int, int] = {}
    for g in groups:
        sizes[g] = sizes.get(g, 0) + 1

    result: list[dict] = [None] * n  # type: ignore[list-item]
    for pos, i in enumerate(order):
        result[i] = {
            "rank": _rank_of(rows[i], i),
            "candidate_confidence": confidences[pos],
            "tie_group": groups[pos],
            "tied": sizes[groups[pos]] > 1,
        }
    return result


CANDIDATE_TIE_NOTE = (
    "<p class='note'>Tied with other same-formula candidates: the formula/class is "
    "supported but the specific isomer is not resolved.</p>"
)

CANDIDATE_SCORES_NOTE = (
    '<p class="note"><strong>Reading the scores:</strong> "Candidate confidence" is a 0-1 '
    "evidence-based match strength for this specific structure (completeness-weighted structure, "
    "mass, MS2, and pathway support); it decreases down the ranking, and the rank-1 value is the "
    'feature\'s overall evidence match strength. "Mass match" and "Formula/class evidence" '
    "corroborate the molecular formula/class but are shared across same-formula isomers, so they "
    "do not distinguish individual candidates. When candidates are marked tied, the formula/class "
    'is supported but the exact isomer is not resolved (check the "Mass-consistent alternative" '
    "note for any better-fitting structure). The separate header shows the calibrated overall "
    "identification confidence, which is intentionally conservative.</p>"
)


def candidate_score_line_html(row: dict) -> str:
    """Shared server-rendered score line for one candidate (renderer + patcher).

    Leads with ``candidate_confidence`` (monotone with rank), then ``Mass match``
    and the class-level ``Formula/class evidence`` (the pipeline
    ``evidence_confidence``), then ``Mass error``. Flags tied isomers. Both the
    from-source renderer and the delivered-set patcher call this so their cards
    are byte-identical.
    """

    def fnum(value, nd):
        try:
            return f"{float(value):.{nd}f}"
        except (TypeError, ValueError):
            return None

    parts = []
    cc = fnum(row.get("candidate_confidence"), 3)
    mm = fnum(row.get("mass_match_score"), 3)
    ev = fnum(row.get("evidence_confidence"), 3)
    me = fnum(row.get("abs_mass_error_da"), 6)
    if cc is not None:
        parts.append(f"Candidate confidence: <strong>{cc}</strong>")
    if mm is not None:
        parts.append(f"Mass match: <strong>{mm}</strong>")
    if ev is not None:
        parts.append(f"Formula/class evidence: <strong>{ev}</strong>")
    if me is not None:
        parts.append(f"Mass error (Da): <strong>{me}</strong>")
    line = ("<p>" + " &middot; ".join(parts) + "</p>") if parts else ""
    if row.get("candidate_confidence_tied"):
        line += CANDIDATE_TIE_NOTE
    return line


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
