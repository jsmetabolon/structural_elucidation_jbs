"""Reconstruct the neighbor vote behind a pathway assignment, for the report.

The vote that produced a compound's pathway is not recorded anywhere. The
enrichment stage keeps six scalars per compound and discards which reference
metabolites voted, how similar they were, and how the tally came out
(``PATHWAY_COLUMNS`` in ``structural_elucidation/pipeline/enrichment/pathway.py``).
The production builder that cast those votes is a legacy import that is not part
of this repository, so there is nothing to replay either.

What this module does instead is re-run the mechanism the manuscript documents -
Morgan radius-2 2048-bit fingerprints, k=10 neighbors, Tanimoto floor 0.30,
similarity-squared vote weights - over the same reference table, and report the
result as a reconstruction. That is honest and it is nearly free: the report
already fingerprints the target and scores it against every reference metabolite
to build its neighbor table, so all this adds is a sort and a tally.

It will not always agree with the delivered assignment. The production builder
also applied molecular-weight priors, SMARTS rules and amino-acid overrides that
are not documented here. A disagreement is worth showing rather than hiding: it
marks a compound whose assignment rests on something the published mechanism
does not account for.

Standard library only, so the renderer and the delivered-set patcher can share
one implementation. The caller supplies similarities; nothing here needs RDKit.
"""
from __future__ import annotations

import html
import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Sequence

# The documented mechanism. These match tools/audit_pathway_predictions.py, whose
# leave-one-out validation is the published measurement of how well this vote
# performs; a reconstruction that used different parameters could not borrow that
# measurement to describe itself.
VOTE_SIMILARITY_FLOOR = 0.30
VOTE_K = 10
VOTE_WEIGHT_POWER = 2.0

# Bin edges for the measured accuracy lookup, matching SIM_BINS in the audit so a
# compound lands in the same bin the audit scored.
SIMILARITY_BINS: tuple[tuple[float, float], ...] = ((0.0, 0.3), (0.3, 0.5), (0.5, 0.7), (0.7, 1.0001))

_PAPER_DATA = Path(__file__).resolve().parents[2] / "docs" / "paper" / "data"
PATHWAY_AUDIT_PATH = _PAPER_DATA / "pathway" / "audit_pathway_summary.json"
# The ceiling lives in the harmonization analysis rather than the pathway audit,
# because it is a property of the two label vocabularies and not of the predictor.
TAXONOMY_CEILING_PATH = _PAPER_DATA / "taxonomy_harmonization.json"

# How a report came to treat a compound as a catalogue lookup. All three
# currently collapse into source="exact" with tier "high" and similarity 1.0,
# which overstates the weaker two: a name match says a compound with this name is
# catalogued, not that this structure is that compound.
BASIS_INCHIKEY = "inchikey_skeleton"
BASIS_CANONICAL_NAME = "canonical_name"
BASIS_KNOWN_X_LABEL = "known_x_label"
BASIS_OPERATOR_SUPPLIED = "operator_supplied"
BASIS_NEIGHBOR_VOTE = "neighbor_vote"
BASIS_NONE = "none"

BASIS_LABELS = {
    BASIS_INCHIKEY: "Structure match to a catalogued compound",
    BASIS_CANONICAL_NAME: "Name match to a catalogued compound",
    BASIS_KNOWN_X_LABEL: "Curator-supplied library label from known_X",
    BASIS_OPERATOR_SUPPLIED: "Chemical id supplied when the report was built",
    BASIS_NEIGHBOR_VOTE: "Consensus of the most similar catalogued compounds",
    BASIS_NONE: "No catalogue match",
}

BASIS_DETAIL = {
    BASIS_INCHIKEY: (
        "The proposed structure's InChIKey skeleton matches a catalogued compound, so that "
        "compound's pathway was copied across. This is the strongest of the lookup routes: "
        "it is a structural match, not a name match."
    ),
    BASIS_CANONICAL_NAME: (
        "A catalogued compound shares this compound's name, and its pathway was copied across. "
        "A name match is weaker than a structure match: it says a compound with this name is "
        "catalogued, not that the proposed structure is that compound."
    ),
    BASIS_KNOWN_X_LABEL: (
        "A curator recorded a library label for this feature in known_X, and the catalogued "
        "compound behind that label supplied the pathway. The pathway is only as good as the "
        "label; no structure was compared."
    ),
    BASIS_OPERATOR_SUPPLIED: (
        "The chemical identifier was supplied directly when this report was built, and its "
        "catalogued pathway was copied across."
    ),
    BASIS_NEIGHBOR_VOTE: (
        "No catalogued compound matched, so the pathway was inferred from the catalogued "
        "compounds most similar to the proposed structure."
    ),
    BASIS_NONE: (
        "Nothing in the catalogue was close enough to support a pathway assignment."
    ),
}

# Only a structure match supports the claim that this compound is the catalogued
# one. The report should not let the name-based routes read as if they did.
BASIS_IS_STRUCTURAL = {
    BASIS_INCHIKEY: True,
    BASIS_CANONICAL_NAME: False,
    BASIS_KNOWN_X_LABEL: False,
    BASIS_OPERATOR_SUPPLIED: False,
    BASIS_NEIGHBOR_VOTE: False,
    BASIS_NONE: False,
}


@dataclass(frozen=True)
class Voter:
    """One catalogued metabolite that carried weight in the vote."""

    chemical_id: str
    name: str
    smiles: str
    similarity: float
    super_pathway: str
    sub_pathway: str
    weight: float

    def as_dict(self) -> dict[str, Any]:
        return {
            "chemical_id": self.chemical_id,
            "canonical_name": self.name,
            "smiles": self.smiles,
            "similarity_to_target": round(self.similarity, 6),
            "metabolon_super": self.super_pathway,
            "metabolon_sub": self.sub_pathway,
            "vote_weight": round(self.weight, 6),
        }


@dataclass
class PathwayVote:
    """The reconstructed vote and everything needed to explain it."""

    voters: list[Voter] = field(default_factory=list)
    super_tally: list[tuple[str, float]] = field(default_factory=list)
    sub_tally: list[tuple[str, float]] = field(default_factory=list)
    winning_super: str = ""
    winning_sub: str = ""
    super_share: float = 0.0
    sub_share: float = 0.0
    total_weight: float = 0.0
    n_above_floor: int = 0
    floor: float = VOTE_SIMILARITY_FLOOR
    k: int = VOTE_K
    power: float = VOTE_WEIGHT_POWER

    @property
    def available(self) -> bool:
        return bool(self.voters)

    @property
    def best_similarity(self) -> float | None:
        return self.voters[0].similarity if self.voters else None

    def as_dict(self) -> dict[str, Any]:
        return {
            "winning_super": self.winning_super,
            "winning_sub": self.winning_sub,
            "super_share": round(self.super_share, 4),
            "sub_share": round(self.sub_share, 4),
            "total_weight": round(self.total_weight, 6),
            "n_voters": len(self.voters),
            "n_above_floor": self.n_above_floor,
            "floor": self.floor,
            "k": self.k,
            "power": self.power,
            "super_tally": [[name, round(weight, 6)] for name, weight in self.super_tally],
            "sub_tally": [[name, round(weight, 6)] for name, weight in self.sub_tally],
            "voters": [voter.as_dict() for voter in self.voters],
        }


def _text(value: Any) -> str:
    text = str("" if value is None else value).strip()
    return "" if text.lower() in {"nan", "none", "null", "<na>"} else text


def _number(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(number) else number


def _ranked(tally: dict[str, float]) -> list[tuple[str, float]]:
    """Sort a tally by weight, keeping insertion order as the tie-break.

    Insertion order here is descending similarity, because the caller iterates
    voters in that order, so a tie resolves toward the class whose closest
    neighbor is closest. That matches how the audit's ``max()`` over a
    ``defaultdict`` resolves ties, and unlike ``max()`` it is explicit.
    """
    order = {name: index for index, name in enumerate(tally)}
    return sorted(tally.items(), key=lambda item: (-item[1], order[item[0]]))


def pathway_vote(
    reference_rows: Iterable[dict[str, Any]],
    *,
    floor: float = VOTE_SIMILARITY_FLOOR,
    k: int = VOTE_K,
    power: float = VOTE_WEIGHT_POWER,
) -> PathwayVote:
    """Run the documented neighbor-consensus vote over scored reference rows.

    ``reference_rows`` are catalogued metabolites already scored against the
    target, each carrying ``similarity_to_target``, ``metabolon_super`` and
    ``metabolon_sub``. Rows without a super-pathway label are dropped, matching
    the reference population the audit votes over.

    The subpathway is voted only among the neighbors that backed the winning
    super-pathway. Voting it globally would let a subpathway win under a
    super-pathway the compound was not assigned to.
    """
    scored: list[tuple[float, dict[str, Any]]] = []
    for row in reference_rows:
        similarity = _number(row.get("similarity_to_target"))
        if similarity is None or similarity < floor:
            continue
        if not _text(row.get("metabolon_super")):
            continue
        scored.append((similarity, row))
    scored.sort(key=lambda item: item[0], reverse=True)
    n_above_floor = len(scored)
    kept = scored[: max(0, int(k))]
    if not kept:
        return PathwayVote(n_above_floor=n_above_floor, floor=floor, k=k, power=power)

    voters: list[Voter] = []
    super_tally: dict[str, float] = {}
    for similarity, row in kept:
        weight = similarity**power
        super_pathway = _text(row.get("metabolon_super"))
        voters.append(
            Voter(
                chemical_id=_text(row.get("chemical_id")),
                name=_text(row.get("canonical_name")),
                smiles=_text(row.get("smiles")),
                similarity=similarity,
                super_pathway=super_pathway,
                sub_pathway=_text(row.get("metabolon_sub")),
                weight=weight,
            )
        )
        super_tally[super_pathway] = super_tally.get(super_pathway, 0.0) + weight

    ranked_super = _ranked(super_tally)
    winning_super, winning_weight = ranked_super[0]
    total_weight = sum(super_tally.values())

    sub_tally: dict[str, float] = {}
    for voter in voters:
        if voter.super_pathway == winning_super and voter.sub_pathway:
            sub_tally[voter.sub_pathway] = sub_tally.get(voter.sub_pathway, 0.0) + voter.weight
    ranked_sub = _ranked(sub_tally)
    winning_sub, winning_sub_weight = ranked_sub[0] if ranked_sub else ("", 0.0)
    sub_total = sum(sub_tally.values())

    return PathwayVote(
        voters=voters,
        super_tally=ranked_super,
        sub_tally=ranked_sub,
        winning_super=winning_super,
        winning_sub=winning_sub,
        super_share=(winning_weight / total_weight) if total_weight else 0.0,
        sub_share=(winning_sub_weight / sub_total) if sub_total else 0.0,
        total_weight=total_weight,
        n_above_floor=n_above_floor,
        floor=floor,
        k=k,
        power=power,
    )


def similarity_bin(value: Any) -> str:
    """Label the audit similarity bin a value falls in, or ``""``."""
    number = _number(value)
    if number is None:
        return ""
    for low, high in SIMILARITY_BINS:
        if low <= number < high:
            return f"[{low:.1f},{min(high, 1.0):.1f})"
    return ""


_AUDIT_CACHE: dict[str, dict] = {}


def _load_json(source: Path) -> dict:
    key = str(source)
    if key not in _AUDIT_CACHE:
        try:
            _AUDIT_CACHE[key] = json.loads(source.read_text())
        except (OSError, ValueError):
            _AUDIT_CACHE[key] = {}
    return _AUDIT_CACHE[key]


def load_pathway_audit(path: Path | str | None = None) -> dict:
    """Load the committed cohort audit, or an empty dict if it is unreadable.

    The report quotes measured agreement rates rather than asserting its own, so
    the numbers a reader sees stay tied to the audit that produced them. A
    missing file drops the reliability paragraph rather than failing the render.
    """
    return _load_json(Path(path) if path is not None else PATHWAY_AUDIT_PATH)


def measured_reliability(
    *,
    best_similarity: Any = None,
    source: str = "",
    audit: dict | None = None,
) -> dict[str, Any] | None:
    """Measured super/sub agreement for a compound's own evidence level.

    Prefers the similarity bin, which is specific to this compound, and falls
    back to the coarser provenance split when no similarity is recorded. Returns
    None when neither is available, so the caller can omit the claim instead of
    inventing one.
    """
    audit = load_pathway_audit() if audit is None else audit
    accuracy = audit.get("accuracy") or {}

    bin_label = similarity_bin(best_similarity)
    block = (accuracy.get("by_similarity_bin") or {}).get(bin_label) if bin_label else None
    basis = f"catalogued compounds whose nearest reference similarity fell in {bin_label}"
    if block is None:
        group = "exact_reference_lookup" if str(source).strip() == "exact" else "neighbor_consensus"
        block = (accuracy.get("by_assignment_source") or {}).get(group)
        basis = (
            "compounds assigned by direct catalogue lookup"
            if group == "exact_reference_lookup"
            else "compounds assigned by neighbor consensus"
        )
    if not block:
        return None

    super_block = block.get("super") or {}
    sub_block = block.get("sub") or {}
    if super_block.get("accuracy") is None:
        return None
    return {
        "basis": basis,
        "n": super_block.get("n"),
        "super_accuracy": super_block.get("accuracy"),
        "super_ci95": super_block.get("ci95"),
        "sub_accuracy": sub_block.get("accuracy"),
        "sub_ci95": sub_block.get("ci95"),
    }


def annotation_ceiling(path: Path | str | None = None) -> dict[str, Any] | None:
    """How often the two label sources agree on compounds of identical structure.

    Every agreement rate in the report is measured against ``known_X`` labels,
    and on compounds whose structure both sources agree on they still disagree
    about the pathway. That bounds any score quoted against ``known_X``, so
    quoting one without the other reads as prediction error when part of it is
    the two vocabularies disagreeing.
    """
    harmonization = _load_json(Path(path) if path is not None else TAXONOMY_CEILING_PATH)
    super_block = ((harmonization.get("super_pathway") or {}).get("raw") or {}).get("reference_ceiling") or {}
    sub_block = ((harmonization.get("subpathway") or {}).get("raw") or {}).get("reference_ceiling") or {}
    if super_block.get("accuracy") is None:
        return None
    return {
        "n": super_block.get("n"),
        "super_agreement": super_block.get("accuracy"),
        "sub_agreement": sub_block.get("accuracy"),
    }


def vote_agrees_with_delivered(vote: PathwayVote, delivered_super: Any, delivered_sub: Any) -> dict[str, bool]:
    """Compare the reconstruction with the assignment the report is showing."""

    def same(left: Any, right: Any) -> bool:
        return " ".join(_text(left).lower().split()) == " ".join(_text(right).lower().split())

    return {
        "super": same(vote.winning_super, delivered_super),
        "sub": same(vote.winning_sub, delivered_sub),
    }


def vote_rows(vote: PathwayVote) -> list[dict[str, Any]]:
    """Ledger rows for the ``<slug>_pathway_vote.csv`` sidecar."""
    rows = []
    for rank, voter in enumerate(vote.voters, start=1):
        row = voter.as_dict()
        row["rank"] = rank
        row["backed_winning_super"] = int(voter.super_pathway == vote.winning_super)
        rows.append(row)
    return rows


PUBLIC_VOTER_LIMIT = 5


def _pct(value: Any) -> str:
    number = _number(value)
    return "n/a" if number is None else f"{round(100 * number)}%"


def _assignment_cards_html(pathway_info: dict[str, Any], vote: PathwayVote) -> str:
    basis = str(pathway_info.get("basis") or BASIS_NONE)
    tier = _text(pathway_info.get("tier")) or "not recorded"
    similarity = _number(pathway_info.get("similarity"))
    similarity_text = "n/a" if similarity is None else f"{similarity:.2f}"
    cards = [
        (
            "Assigned pathway",
            _text(pathway_info.get("super")) or "Not assigned",
            BASIS_LABELS.get(basis, BASIS_LABELS[BASIS_NONE]),
        ),
        (
            "Assigned subpathway",
            _text(pathway_info.get("sub")) or "Not assigned",
            "The finer of the two labels, and the less reliable",
        ),
        (
            "Nearest catalogued compound",
            similarity_text,
            "Similarity to the closest catalogued match, as recorded with the assignment",
        ),
        (
            "Assignment confidence tier",
            tier.title(),
            f"{vote.n_above_floor} catalogued compounds cleared the {vote.floor:.2f} similarity floor"
            if vote.n_above_floor
            else "No catalogued compound cleared the similarity floor",
        ),
    ]
    return '<div class="metrics">' + "".join(
        '<div class="metric-card">'
        f'<div class="metric-label">{html.escape(label)}</div>'
        f'<div class="metric-value">{html.escape(str(value))}</div>'
        f'<div class="metric-note">{html.escape(note)}</div></div>'
        for label, value, note in cards
    ) + "</div>"


def _voter_table_html(vote: PathwayVote, *, limit: int | None) -> str:
    voters = vote.voters if limit is None else vote.voters[:limit]
    if not voters:
        return ""
    body = "".join(
        "<tr>"
        f"<td>{index}</td>"
        f"<td>{html.escape(voter.name or voter.chemical_id or 'unnamed')}</td>"
        f"<td>{voter.similarity:.3f}</td>"
        f"<td>{html.escape(voter.super_pathway)}</td>"
        f"<td>{html.escape(voter.sub_pathway or '-')}</td>"
        f"<td>{voter.weight:.3f}</td>"
        "</tr>"
        for index, voter in enumerate(voters, start=1)
    )
    omitted = len(vote.voters) - len(voters)
    tail = (
        f'<p class="note">{omitted} further voter(s) are listed in the internal report '
        f"and in the data sidecar.</p>"
        if omitted > 0
        else ""
    )
    return (
        '<div class="table-wrap"><table><thead><tr>'
        "<th>#</th><th>Catalogued compound</th><th>Similarity</th>"
        "<th>Pathway</th><th>Subpathway</th><th>Vote weight</th>"
        "</tr></thead><tbody>"
        + body
        + "</tbody></table></div>"
        + tail
    )


def _tally_table_html(vote: PathwayVote) -> str:
    if not vote.super_tally:
        return ""
    total = vote.total_weight or 1.0
    super_body = "".join(
        "<tr>"
        f"<td>{html.escape(name)}</td>"
        f"<td>{weight:.3f}</td>"
        f"<td>{100 * weight / total:.0f}%</td>"
        f"<td>{'winner' if name == vote.winning_super else ''}</td>"
        "</tr>"
        for name, weight in vote.super_tally
    )
    sub_total = sum(weight for _name, weight in vote.sub_tally) or 1.0
    sub_body = "".join(
        "<tr>"
        f"<td>{html.escape(name)}</td>"
        f"<td>{weight:.3f}</td>"
        f"<td>{100 * weight / sub_total:.0f}%</td>"
        f"<td>{'winner' if name == vote.winning_sub else ''}</td>"
        "</tr>"
        for name, weight in vote.sub_tally
    )
    out = (
        '<div class="table-wrap"><table><thead><tr>'
        "<th>Pathway</th><th>Weight</th><th>Share</th><th></th>"
        "</tr></thead><tbody>" + super_body + "</tbody></table></div>"
    )
    if sub_body:
        out += (
            f'<p class="note">Subpathway is decided only among the voters that backed '
            f"{html.escape(vote.winning_super)}. Voting it across all neighbors could return a "
            "subpathway belonging to a pathway this compound was not assigned to.</p>"
            '<div class="table-wrap"><table><thead><tr>'
            "<th>Subpathway</th><th>Weight</th><th>Share</th><th></th>"
            "</tr></thead><tbody>" + sub_body + "</tbody></table></div>"
        )
    return out


def _reliability_html(pathway_info: dict[str, Any]) -> str:
    measured = measured_reliability(
        best_similarity=pathway_info.get("similarity"),
        source=str(pathway_info.get("source") or ""),
    )
    if not measured:
        return ""
    ceiling = annotation_ceiling()
    text = (
        f'<p class="note"><strong>How often this kind of assignment has been right.</strong> '
        f"Across {measured['n']} {html.escape(measured['basis'])}, the pathway matched an "
        f"independent annotation {_pct(measured['super_accuracy'])} of the time"
    )
    if measured.get("sub_accuracy") is not None:
        text += f" and the subpathway {_pct(measured['sub_accuracy'])} of the time"
    text += "."
    if ceiling:
        text += (
            f" Both figures are bounded: on {ceiling['n']} compounds where the two label sources "
            f"describe the same structure, they still agree only {_pct(ceiling['super_agreement'])} "
            f"of the time on pathway and {_pct(ceiling['sub_agreement'])} on subpathway, so part of "
            "the gap is the two vocabularies disagreeing rather than a wrong call."
        )
    return text + "</p>"


def _agreement_html(vote: PathwayVote, pathway_info: dict[str, Any]) -> str:
    """State whether the reconstruction landed where the delivered assignment did."""
    if not vote.available:
        return ""
    agreement = vote_agrees_with_delivered(vote, pathway_info.get("super"), pathway_info.get("sub"))
    if agreement["super"]:
        return (
            '<p class="note">Re-running the published mechanism over the catalogue reproduces the '
            "assigned pathway.</p>"
        )
    return (
        '<p class="note"><strong>Re-running the published mechanism does not reproduce the assigned '
        f"pathway.</strong> The vote shown above lands on "
        f"<strong>{html.escape(vote.winning_super or 'no pathway')}</strong> while the report shows "
        f"<strong>{html.escape(_text(pathway_info.get('super')) or 'none')}</strong>. The production "
        "assignment also applied molecular-weight priors, substructure rules and amino-acid "
        "overrides that the published mechanism does not include, so a disagreement means the "
        "assignment rests on one of those rather than on the neighbors alone.</p>"
    )


def pathway_section_html(
    pathway_info: dict[str, Any],
    vote: PathwayVote,
    *,
    public: bool,
    msi_level: int = 3,
    candidate_agreement: str = "",
) -> str:
    """Body HTML for the pathway and subpathway assignment section.

    Shown at every MSI level, because a pathway is assigned to every feature
    including the ones whose structure is withheld. That is exactly why Levels
    4-5 need the warning: the pathway was inferred from the rank-1 structure, so
    a report that rejects that structure has nothing left holding its pathway up.
    """
    basis = str(pathway_info.get("basis") or BASIS_NONE)
    parts: list[str] = [_assignment_cards_html(pathway_info, vote)]

    if msi_level >= 4:
        parts.append(
            '<p class="note"><strong>This assignment rests on a structure this report does not '
            "stand behind.</strong> The pathway was taken from the top-ranked structure, and that "
            "structure did not survive the evidence checks, so treat the pathway as unsupported "
            "rather than as context that outlived the structure.</p>"
        )

    # The neighbor-vote blurb is subsumed by the vote paragraph below, which says
    # the same thing with this compound's actual numbers.
    if not (basis == BASIS_NEIGHBOR_VOTE and vote.available):
        parts.append(f'<p class="note">{html.escape(BASIS_DETAIL.get(basis, BASIS_DETAIL[BASIS_NONE]))}</p>')

    if basis in {BASIS_INCHIKEY, BASIS_CANONICAL_NAME, BASIS_KNOWN_X_LABEL, BASIS_OPERATOR_SUPPLIED}:
        parts.append(
            '<p class="note">A lookup copies the catalogued label exactly, so the label itself is '
            "not in doubt. What is in doubt is whether the proposed structure is really that "
            "compound: across the structure-verified controls, a top-ranked structure that matched "
            "the catalogue was the correct structure 38% of the time. The pathway is only as sound "
            "as that match.</p>"
        )
    elif vote.available:
        closest = vote.voters[0]
        parts.append(
            f'<p class="note">No catalogued compound matched this structure outright, so the '
            f"pathway comes from the company its structure keeps. The {len(vote.voters)} "
            f"catalogued compounds below are the closest matches that cleared a Tanimoto "
            f"similarity of {vote.floor:.2f}, the nearest being "
            f"<strong>{html.escape(closest.name or closest.chemical_id or 'an unnamed compound')}</strong> "
            f"at {closest.similarity:.2f}. Each votes for its own pathway with a weight of its "
            f"similarity squared, so a close neighbor counts for far more than a distant one, and "
            f"the pathway holding the most weight wins. Here that was "
            f"<strong>{html.escape(vote.winning_super or 'no pathway')}</strong> with "
            f"{_pct(vote.super_share)} of the weight.</p>"
        )
        parts.append(_voter_table_html(vote, limit=PUBLIC_VOTER_LIMIT if public else None))
        parts.append(_tally_table_html(vote))
        parts.append(_agreement_html(vote, pathway_info))
    else:
        parts.append(
            '<p class="note">No catalogued compound was similar enough to vote, so there is no '
            "neighbor evidence to show. Any pathway displayed came from the upstream enrichment "
            "and cannot be explained from this report's own catalogue comparison.</p>"
        )

    if candidate_agreement:
        parts.append(
            f'<p class="note"><strong>Across the candidate list.</strong> {html.escape(candidate_agreement)}. '
            "A pathway that holds across the alternative structures is more robust than one that "
            "depends on the top-ranked structure being right.</p>"
        )

    parts.append(_reliability_html(pathway_info))
    return "\n      ".join(part for part in parts if part)


def reference_rows_from_frame(records: Sequence[dict[str, Any]]) -> list[dict[str, Any]]:
    """Normalize scored reference records into the shape :func:`pathway_vote` wants."""
    return [
        {
            "chemical_id": _text(row.get("chemical_id")),
            "canonical_name": _text(row.get("canonical_name")),
            "smiles": _text(row.get("smiles")),
            "similarity_to_target": _number(row.get("similarity_to_target")),
            "metabolon_super": _text(row.get("metabolon_super")),
            "metabolon_sub": _text(row.get("metabolon_sub")),
        }
        for row in records
    ]
