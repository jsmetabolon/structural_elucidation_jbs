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
# Bar an exact catalog hit must clear to claim a probable structure (MSI L2).
# A catalog hit pins both the structure axis and the pathway axis to 1.0, so half
# the evidence weight rests on one observation - and that observation is correct
# only 38.1% of the time on the truth set, against 56.0% for de novo predictions
# (docs/paper/data/structure/audit_structure_summary.json, source_subgroups). The
# tier rubric already refuses to call an exact match High for this reason; the
# MSI rule needs the same guard, so an exact-source claim of "probable structure"
# has to rest on spectral evidence the catalog did not supply. The value is
# provisional: tools/audit_thresholds.py sweeps it once axis-level evidence for
# the cohort is available.
MS2_PROBABLE_EXACT = 0.50

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

# --- Assurance: human / laboratory verdicts ----------------------------------
# ``evidence_score`` is model-only and says nothing about human review. Human and
# laboratory findings are a separate, discrete, auditable state so that a report
# never presents a policy display value as if it were a fitted probability, and
# so a *refuting* result has somewhere to live. Before this existed, a standard
# run that came back negative could not be recorded at all: the only human inputs
# were promotions.
ASSURANCE_STANDARD_CONFIRMED = "standard_confirmed"
ASSURANCE_STANDARD_REFUTED = "standard_refuted"
ASSURANCE_CURATOR_KNOWN = "curator_known"
ASSURANCE_CURATOR_REJECTED = "curator_rejected"
ASSURANCE_NONE = "none"

ASSURANCE_VERDICTS = (
    ASSURANCE_STANDARD_CONFIRMED,
    ASSURANCE_STANDARD_REFUTED,
    ASSURANCE_CURATOR_KNOWN,
    ASSURANCE_CURATOR_REJECTED,
    ASSURANCE_NONE,
)

# Verdicts that argue against the proposed structure. These cap the MSI level and
# suppress the name; no evidence score can override them.
REFUTING_ASSURANCE = frozenset({ASSURANCE_STANDARD_REFUTED, ASSURANCE_CURATOR_REJECTED})
# Verdicts that establish identity by means outside the model.
SUPPORTING_ASSURANCE = frozenset({ASSURANCE_STANDARD_CONFIRMED, ASSURANCE_CURATOR_KNOWN})

ASSURANCE_LABELS = {
    ASSURANCE_STANDARD_CONFIRMED: "Confirmed against an authentic standard",
    ASSURANCE_STANDARD_REFUTED: "Ruled out by an authentic-standard run",
    ASSURANCE_CURATOR_KNOWN: "Curator-established identity (authentic standard pending)",
    ASSURANCE_CURATOR_REJECTED: "Rejected on curator review",
    ASSURANCE_NONE: "No human review recorded",
}


def assurance_of(fields: dict) -> str:
    """Return the recorded human/laboratory verdict for a report.

    Reads the explicit ``assurance`` field when the curation ledger supplied one,
    and otherwise derives it from the legacy inputs so existing bundles keep
    behaving as before: a Metabolon Level-1 control whose top candidate matches
    the confirmed structure is ``standard_confirmed``, and an entry on the
    curator allowlist is ``curator_known``.
    """
    explicit = str(fields.get("assurance", "") or "").strip().lower()
    if explicit in ASSURANCE_VERDICTS:
        return explicit
    if str(fields.get("metabolon_level", "")).strip() == "1" and bool(fields.get("truth_match")):
        return ASSURANCE_STANDARD_CONFIRMED
    if fields.get("known_identity_override"):
        return ASSURANCE_CURATOR_KNOWN
    return ASSURANCE_NONE


# --- how many candidate structures a report may display ----------------------
# A ten-deep list was shown to every reader. On the truth-verified cohort that
# depth surfaces the correct structure for 73.9% of compounds against 64.1% at
# three, and buys those 9.8 points by displaying 70.6 additional incorrect
# structures for each compound it helps (docs/paper/data/thresholds/). Ten
# near-tied isomers also imply a resolution the evidence does not support, so the
# public list is short and the internal one stays deep for triage.
PUBLIC_CANDIDATE_N = 3
INTERNAL_CANDIDATE_N = 10
# Public reports at Level 1-2 assert one structure; showing alternates alongside
# an asserted identification contradicts the claim.
ASSERTED_CANDIDATE_N = 1
# Relative, not absolute. The rank score is an ordering value whose rank-1 median
# is about 0.035 across the delivered set, so any fixed cutoff - 0.1 included -
# would discard most rank-1 candidates. A candidate an order of magnitude below
# rank 1 is not a live alternative and only pads the list.
CANDIDATE_MARGIN_RATIO = 0.10


def candidates_to_display(
    candidates: list,
    *,
    msi_level: int,
    public: bool,
    limit: int | None = None,
    score_key: str = "candidate_report_score",
) -> list:
    """Return the candidate structures this report variant may show.

    Levels 4 and 5 return nothing: those reports withhold the structure, and a
    list of candidate structures is a structural claim whatever it is captioned.
    """
    if msi_level >= 4:
        return []
    if limit is None:
        limit = PUBLIC_CANDIDATE_N if public else INTERNAL_CANDIDATE_N
    if public and msi_level <= 2:
        limit = min(limit, ASSERTED_CANDIDATE_N)
    kept = list(candidates)[: max(0, int(limit))]
    if len(kept) <= 1:
        return kept
    lead = _num(kept[0].get(score_key) if isinstance(kept[0], dict) else None)
    if lead is None or lead <= 0:
        return kept
    floor = CANDIDATE_MARGIN_RATIO * lead
    out = [kept[0]]
    for row in kept[1:]:
        score = _num(row.get(score_key)) if isinstance(row, dict) else None
        if score is None or score >= floor:
            out.append(row)
    return out


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
    """Identification confidence for one report.

    ``evidence_score`` and ``overall`` are deliberately two numbers.
    ``evidence_score`` is the model's own output and nothing is allowed to
    change it: it is what the calibration was fitted on and what the audits
    measure, so overwriting it would make every downstream statistic a
    statement about policy rather than about evidence. ``overall`` is the
    headline the report displays, which folds in ``assurance`` - the human or
    laboratory verdict. When the two differ the report says so.
    """

    overall: float
    tier: str
    tier_note: str
    completeness: float
    descriptor: str
    axes: list[Axis] = field(default_factory=list)
    flags: list[str] = field(default_factory=list)
    source: str = ""
    # Model-only evidence score, never promoted, pinned, or otherwise overwritten.
    evidence_score: float = 0.0
    # Human / laboratory verdict; see ASSURANCE_VERDICTS.
    assurance: str = ASSURANCE_NONE
    assurance_label: str = ASSURANCE_LABELS[ASSURANCE_NONE]
    assurance_basis: str = ""
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
    # Library identity ("authentic-standard pending") metric - separate from ``overall``.
    identity_status: str = IDENTITY_STATUS_UNKNOWN
    identity_confidence: float | None = None
    identity_basis: str = ""

    def as_dict(self) -> dict:
        return {
            "overall": self.overall,
            "evidence_score": self.evidence_score,
            "assurance": self.assurance,
            "assurance_label": self.assurance_label,
            "assurance_basis": self.assurance_basis,
            "tier": self.tier,
            "tier_note": self.tier_note,
            "completeness": self.completeness,
            "descriptor": self.descriptor,
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


@dataclass(frozen=True)
class MsiThresholds:
    """Boundary parameters for :func:`assign_msi_level`.

    Exposed as a value object so tools/audit_thresholds.py can sweep the
    boundaries against the truth set without re-implementing the rule. Drift
    between a swept rule and the delivered rule would make the sweep useless.
    """

    struct_moderate: float = STRUCT_MODERATE
    struct_weak: float = STRUCT_WEAK
    mass_strong: float = MASS_STRONG
    mass_weak: float = MASS_WEAK
    ms2_present: float = MS2_PRESENT
    ms2_probable_exact: float = MS2_PROBABLE_EXACT


DEFAULT_MSI_THRESHOLDS = MsiThresholds()


def assign_msi_level(
    fields: dict,
    conf: "Confidence",
    thresholds: MsiThresholds = DEFAULT_MSI_THRESHOLDS,
) -> tuple[int, str, str]:
    """Assign a Schymanski/MSI identification-confidence level (1-5) from evidence.

    L1 confirmed (authentic standard) / L2 probable structure (strong structure +
    observed MS2 + mass) / L3 tentative candidate structure (hypothesis) / L4
    molecular formula only (structure not supported or mass-rejected) / L5 exact
    mass only. Structures are only claimed at L1-L3; L4-L5 report formula/mass.

    A refuting human verdict (a standard run that came back negative, or a
    curator rejection) caps the level at 4 no matter what the evidence axes say:
    a failed authentic-standard run is stronger evidence against the structure
    than any model score is for it.

    An exact catalog hit faces a higher spectral bar for L2 than a de novo
    prediction does. Provenance pins the structure axis to 1.0 for those
    candidates while predicting correctness worse than the de novo ones do, so
    "probable structure" has to rest on evidence the catalog did not supply.
    """
    by = {a.key: a for a in conf.axes}
    struct = by["structure"].score or 0.0
    mass = by["mass"].score if by["mass"].available else None
    ms2_score = (by["ms2"].score or 0.0) if by["ms2"].available else None
    mass_ok = mass is not None and mass >= thresholds.mass_strong
    mass_mismatch = mass is not None and mass < thresholds.mass_weak
    observed = bool(fields.get("observed_mz_available"))
    from_catalog = str(fields.get("source", "")).strip() == "exact"
    ms2_probable_bar = (
        thresholds.ms2_probable_exact if from_catalog else thresholds.ms2_present
    )

    refuted = assurance_of(fields) in REFUTING_ASSURANCE
    if refuted:
        note = (
            "An authentic-standard run ruled this structure out"
            if assurance_of(fields) == ASSURANCE_STANDARD_REFUTED
            else "A curator rejected this structure on review"
        )
        return 4, MSI_LABELS[4], f"{note} - reporting formula/mass only"
    if conf.tier == "Confirmed":
        return 1, MSI_LABELS[1], "Metabolon Level-1 standard match and the top candidate matches it"
    if struct >= thresholds.struct_moderate and ms2_score is not None and ms2_score >= ms2_probable_bar and mass_ok:
        rationale = "Strong structure similarity with observed-MS2 and mass corroboration"
        if from_catalog:
            rationale += " that is independent of the catalog match"
        return 2, MSI_LABELS[2], rationale
    if (not mass_mismatch) and struct >= thresholds.struct_weak:
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

    assurance = assurance_of(fields)
    if assurance in REFUTING_ASSURANCE:
        return (
            IDENTITY_STATUS_UNKNOWN,
            None,
            "Structure ruled out on review; no identity is claimed",
        )

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
    """Compute the hybrid identification confidence from normalized ``fields``.

    Order matters. The model score is computed first and then frozen. Human and
    laboratory verdicts are folded into the headline afterwards, and only the
    ones that come from an authentic-standard run are allowed to move the MSI
    level: a curator saying "we know what this is" is grounds for a confident
    headline but not for claiming Level 1, which means a standard was run.
    """
    source = str(fields.get("source", "")).strip()
    axes = build_axes(fields)
    raw, completeness = raw_score(axes)
    # Model-only, and frozen from here on. Everything below adjusts the headline
    # ``overall`` instead, so the audits keep measuring evidence rather than policy.
    evidence = _sigmoid(PLATT_A * raw + PLATT_B)
    assurance = assurance_of(fields)
    tier, tier_note = _tier(evidence, axes, source)
    if source == "exact":
        descriptor = "Cataloged reference metabolite (exact structural match)"
    elif source == "predicted":
        descriptor = "De novo / predicted structure"
    else:
        descriptor = "Structure proposal"
    result = Confidence(
        overall=round(evidence, 4),
        tier=tier,
        tier_note=tier_note,
        completeness=round(completeness, 4),
        descriptor=descriptor,
        axes=axes,
        flags=_flags(axes, fields),
        source=source,
        evidence_score=round(evidence, 4),
        assurance=assurance,
        assurance_label=ASSURANCE_LABELS.get(assurance, ASSURANCE_LABELS[ASSURANCE_NONE]),
        assurance_basis=str(fields.get("assurance_reason", "") or "").strip(),
        observed_mz=_num(fields.get("observed_mz")),
        observed_neutral_mass=_num(fields.get("observed_neutral_mass")),
        candidate_exact_mass=_num(fields.get("candidate_exact_mass")),
        ms2_n_pairs=int(_num(fields.get("ms2_n_pairs")) or 0),
        metabolon_level=str(fields.get("metabolon_level", "") or ""),
        mass_alt=fields.get("mass_alt") if isinstance(fields.get("mass_alt"), dict) else None,
    )
    _apply_standard_run(result)
    result.msi_level, result.msi_label, result.msi_rationale = assign_msi_level(fields, result)
    result.identity_status, result.identity_confidence, result.identity_basis = assign_identity(fields, result)
    _apply_curator_identity(result)
    return result


def _apply_standard_run(result: Confidence) -> None:
    """Fold an authentic-standard outcome into the headline, before MSI is set.

    Both directions run through here. A confirming standard is the gold standard
    and reads 100%; a refuting one means the proposed structure is wrong, which
    no evidence score should be allowed to soften.
    """
    if result.assurance == ASSURANCE_STANDARD_CONFIRMED:
        result.tier = "Confirmed"
        result.tier_note = (
            "Authentic-standard confirmed (Metabolon Level 1); top candidate matches the confirmed structure"
        )
        # Display policy, not a recalibrated probability. ``evidence_score`` keeps
        # the model's own value and the scored-parameters section reports it.
        result.overall = max(result.overall, CONFIRMED_OVERALL)
        if result.evidence_score < CONFIRMED_OVERALL - 0.05:
            result.flags.append(
                f"Headline confidence reflects authentic-standard confirmation "
                f"(the model's own evidence score was {round(100 * result.evidence_score)}%)."
            )
    elif result.assurance in REFUTING_ASSURANCE:
        result.tier = "Refuted"
        result.tier_note = (
            "An authentic-standard run ruled this structure out"
            if result.assurance == ASSURANCE_STANDARD_REFUTED
            else "A curator rejected this structure on review"
        )
        if result.assurance_basis:
            result.tier_note += f" ({result.assurance_basis})"
        result.overall = 0.0
        result.flags.append(
            "This structure has been ruled out by review; the report is limited to formula and mass."
        )


def _apply_curator_identity(result: Confidence) -> None:
    """Promote the headline for a curator-established identity, after MSI is set.

    Runs after :func:`assign_msi_level` on purpose: a curator knowing the
    compound justifies a confident headline but not MSI Level 1, which asserts
    that an authentic standard was run.
    """
    if result.tier == "Confirmed" or result.assurance in REFUTING_ASSURANCE:
        return
    if (
        result.identity_status != IDENTITY_STATUS_KNOWN
        or result.identity_confidence is None
        or result.identity_confidence < 0.90
    ):
        return
    # Clamp below the authentic-standard-confirmed floor so a curator override
    # can never match the Level 1 headline (100%).
    result.overall = round(min(float(result.identity_confidence), IDENTITY_PROMOTED_MAX), 4)
    result.tier = "High"
    result.tier_note = (
        "This is a High Confidence probable structure based on corroborating "
        "evidence derived from Metabolon\u2019s validated structural elucidation workflow."
    )
    if result.evidence_score < result.overall - 0.05:
        result.flags.append(
            f"Overall confidence reflects curator-known library identity "
            f"(underlying calibrated evidence score was {round(100 * result.evidence_score)}%)."
        )


@dataclass(frozen=True)
class CurationEntry:
    """One reviewed verdict from the curation ledger."""

    control_id: str
    verdict: str
    identity_confidence: float | None = None
    evidence_type: str = ""
    reviewer: str = ""
    date: str = ""
    standard_run_id: str = ""
    reason: str = ""

    def note(self) -> str:
        """Human-readable provenance for the report card."""
        bits = [bit for bit in (self.reason, self.evidence_type) if bit]
        if self.standard_run_id:
            bits.append(f"standard run {self.standard_run_id}")
        if self.reviewer:
            bits.append(f"reviewed by {self.reviewer}")
        if self.date:
            bits.append(self.date)
        return "; ".join(bits)


def load_curation_ledger(path: Path | str | None) -> dict[str, CurationEntry]:
    """Load the reviewed verdicts keyed by control id.

    Columns (case-insensitive): ``control_id`` (required), ``verdict`` (one of
    :data:`ASSURANCE_VERDICTS`), and the optional ``identity_confidence``,
    ``evidence_type``, ``reviewer``, ``date``, ``standard_run_id`` and ``reason``.

    Files predating the ledger carried an allowlist with no ``verdict`` column;
    every row in one of those means the same thing, so they are read as
    ``curator_known``.

    A missing file is not an error - callers pass an optional path - but an
    unrecognized verdict is, because a mistyped verdict would silently discard a
    refutation, which is the one input that must never be lost.
    """
    if path is None:
        return {}
    source = Path(path)
    if not source.exists():
        return {}
    out: dict[str, CurationEntry] = {}
    try:
        with source.open(newline="") as handle:
            # Skip full-line comments (they must not appear inside quoted fields).
            filtered = (line for line in handle if not line.lstrip().startswith("#"))
            rows = list(csv.DictReader(filtered))
    except OSError:
        return {}
    for row in rows:
        lowered = {str(key or "").strip().lower(): value for key, value in row.items()}
        control_id = str(lowered.get("control_id") or "").strip()
        if not control_id or control_id.startswith("#"):
            continue
        verdict = str(lowered.get("verdict") or "").strip().lower() or ASSURANCE_CURATOR_KNOWN
        if verdict not in ASSURANCE_VERDICTS:
            raise ValueError(
                f"{source}: unrecognized verdict {verdict!r} for {control_id}; "
                f"expected one of {', '.join(ASSURANCE_VERDICTS)}"
            )
        confidence = _num(lowered.get("identity_confidence") or lowered.get("confidence"))
        if verdict == ASSURANCE_CURATOR_KNOWN and confidence is None:
            confidence = IDENTITY_CONF_OVERRIDE_DEFAULT
        out[control_id] = CurationEntry(
            control_id=control_id,
            verdict=verdict,
            identity_confidence=confidence,
            evidence_type=str(lowered.get("evidence_type") or "").strip(),
            reviewer=str(lowered.get("reviewer") or "").strip(),
            date=str(lowered.get("date") or "").strip(),
            standard_run_id=str(lowered.get("standard_run_id") or "").strip(),
            reason=str(lowered.get("reason") or lowered.get("note") or "").strip(),
        )
    return out


def apply_curation_entry(fields: dict, entry: CurationEntry | None) -> None:
    """Fold a reviewed verdict into the confidence ``fields`` in place."""
    if entry is None:
        return
    fields["assurance"] = entry.verdict
    fields["assurance_reason"] = entry.note()
    if entry.verdict == ASSURANCE_CURATOR_KNOWN:
        fields["known_identity_override"] = {
            "confidence": entry.identity_confidence
            if entry.identity_confidence is not None
            else IDENTITY_CONF_OVERRIDE_DEFAULT,
            "reason": entry.reason,
        }


def load_identity_overrides(path: Path | str | None) -> dict[str, dict]:
    """Known-identity allowlist view of the ledger, for callers that predate it.

    Returns ``{control_id: {"confidence": float, "reason": str}}`` for the
    ``curator_known`` entries only.
    """
    try:
        ledger = load_curation_ledger(path)
    except ValueError:
        return {}
    return {
        control_id: {
            "confidence": entry.identity_confidence
            if entry.identity_confidence is not None
            else IDENTITY_CONF_OVERRIDE_DEFAULT,
            "reason": entry.reason,
        }
        for control_id, entry in ledger.items()
        if entry.verdict == ASSURANCE_CURATOR_KNOWN
    }


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
    "Low": "badge", "Tentative": "badge", "Refuted": "badge fail",
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


def withheld_structures_note(conf: Confidence) -> str:
    """Explain why a Level 4-5 report shows no candidate structures."""
    if conf.assurance in REFUTING_ASSURANCE:
        return (
            "The proposed structure was ruled out on review, so no candidate structure "
            "is shown. The molecular formula and observed mass below stand on their own."
        )
    if conf.msi_level >= 5:
        return (
            "Only an exact mass is established for this feature (MSI Level 5). No formula "
            "or structure is supported well enough to show."
        )
    return (
        "This feature is reported at the level of its molecular formula (MSI Level 4). "
        "No candidate structure is supported well enough to show; "
        f"{conf.msi_rationale[:1].lower()}{conf.msi_rationale[1:]}."
    )


def assurance_note_html(conf: Confidence) -> str:
    """State the human/laboratory verdict, and the model score whenever it differs.

    The headline and the model's own score are allowed to disagree - that is the
    point of recording a verdict - but a report that showed only the headline
    would be presenting a policy decision as a measurement.
    """
    if conf.assurance == ASSURANCE_NONE:
        return ""
    label = html.escape(ASSURANCE_LABELS.get(conf.assurance, ""))
    basis = f" {html.escape(conf.assurance_basis)}." if conf.assurance_basis else ""
    note = f'<p class="note"><strong>Human review:</strong> {label}.{basis}'
    if abs(conf.overall - conf.evidence_score) >= 0.005:
        note += (
            f" The headline confidence of {_pct(conf.overall)} reflects that review; "
            f"the model's own evidence score for this feature is {_pct(conf.evidence_score)}."
        )
    return note + "</p>"


def assurance_badge_html(conf: Confidence) -> str:
    """Header badge for a refuting verdict, which a reader must not miss."""
    if conf.assurance not in REFUTING_ASSURANCE:
        return ""
    return f'<span class="badge fail">{html.escape(ASSURANCE_LABELS[conf.assurance])}</span>'


def msi_badge_html(conf: Confidence) -> str:
    cls = "badge ok" if conf.msi_level <= 2 else "badge"
    return f'<span class="{cls}">MSI Level {conf.msi_level} \u00b7 {html.escape(conf.msi_label)}</span>'


def msi_metric_card_html(conf: Confidence) -> str:
    return (
        '<div class="metric-card"><div class="metric-label">MSI identification level</div>'
        f'<div class="metric-value">Level {conf.msi_level}</div>'
        f'<div class="metric-note">{html.escape(conf.msi_label)}</div></div>'
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
    msi_line = (
        f'<p class="note"><strong>MSI identification level: Level {conf.msi_level} - '
        f"{html.escape(conf.msi_label)}.</strong> {html.escape(conf.msi_rationale)}. "
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
    ) + assurance_note_html(conf) + identity_note_html(conf, public=public)
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
    if conf.metabolon_level:
        extra += (
            f'<p class="note"><strong>Metabolon identification Level:</strong> '
            f"{html.escape(conf.metabolon_level)}"
            + (" (authentic-standard confirmed)" if conf.metabolon_level == "1" else "")
            + ".</p>"
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
        # A public report that withholds the structure cannot name an alternative
        # one either; the formula and rank still point a curator at the fix.
        if public and conf.msi_level >= 4:
            smiles = ""
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
    structures_withheld: bool = False,
) -> str:
    """Non-raw MS2 fragment summary for the MS Fit / Feature Match section.

    ``structures_withheld`` only changes what the empty table says. A Level 4-5
    report has no candidates to score, which is different from having candidates
    whose scores went missing, and saying the wrong one invites a bug report.
    """

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
    elif structures_withheld:
        parts.append(
            '<p class="note">No candidate structure is reported for this feature, so there '
            "are no per-candidate spectral fits to show.</p>"
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
