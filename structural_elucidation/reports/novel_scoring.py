"""Headless scoring for the novel-feature driver and its report renderer.

The payload is the durable scientific result: rendering consumes its fields and
never recomputes a different confidence from presentation-specific artifacts.
"""
from __future__ import annotations

import hashlib
import html
import json
import math
from dataclasses import replace
from pathlib import Path
from typing import Any

import pandas as pd

from structural_elucidation.reports import confidence as confmod
from structural_elucidation.reports import pathway_explain as pathmod
from structural_elucidation.utils.domain import ri_from_feature_key
from structural_elucidation.nce.evidence import ri_proxy_from_smiles

SUMMARY_COLUMNS = (
    "se_msi_level", "se_msi_label", "se_msi_rationale", "se_confidence",
    "se_confidence_tier", "se_evidence_score", "se_reported_level",
    "se_reported_identification", "se_candidate_identity_confidence",
    "se_confidence_flags_json",
)


def number(value: Any) -> float | None:
    try:
        result = float(value)
        return result if math.isfinite(result) else None
    except (TypeError, ValueError):
        return None


def text(value: Any) -> str:
    result = str(value if value is not None else "").strip()
    return "" if result.lower() in {"nan", "none", "null", "<na>"} else result


def json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if hasattr(value, "item"):
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def prepare_confidence_fields(
    rank1: dict, candidate_scores: list[dict], pathway_info: dict, *,
    pathway_vote_share: float | None, neighbors: list[dict], n_peaks: int,
    target_smiles: str, target_exact_mass: float | None,
    known_x_row: dict | None = None, naming_decision: dict | None = None,
    curation_entry: Any = None, correlations_real: bool = False,
) -> dict:
    """The renderer's evidence normalization, shared with headless scoring."""
    ms2 = [number(rank1.get(key)) for key in (
        "predicted_spectrum_score_v2", "forward_spectrum_score", "ms2_match_score",
    )]
    ms2 = [value for value in ms2 if value is not None]
    ms2_score = max(ms2) if ms2 else None
    observed_mz = number(rank1.get("parent_mz"))
    fields = confmod.fields_from_parts(
        candidate_scores[0] if candidate_scores else {},
        candidate_scores[1] if len(candidate_scores) > 1 else None,
        {"pathway_source": pathway_info.get("source", ""),
         "pathway_best_reference_similarity": pathway_info.get("similarity"),
         "pathway_vote_share": pathway_vote_share,
         "pathway_confidence_tier": pathway_info.get("tier", "")},
        observed_mz_available=observed_mz is not None and observed_mz > 0,
        ms2_spectrum_present=n_peaks > 0 and ms2_score is not None,
        ms2_score=ms2_score,
        library_neighbors_ge_050=sum(
            1 for row in neighbors
            if (number(row.get("similarity_to_target")) or 0) >= confmod.NEIGHBOR_SIMILARITY_FLOOR
        ),
        correlations_real=correlations_real,
        lib_id=rank1.get("lib_id"),
        observed_retention_index=number(ri_from_feature_key(rank1.get("feature_key"))),
        candidate_retention_proxy=number(ri_proxy_from_smiles(target_smiles)) if target_smiles else None,
    )
    neutral_mass = number(rank1.get("target_neutral_mass"))
    if neutral_mass is None:
        neutral_mass = number(rank1.get("neutral_mass"))
    level = text((known_x_row or {}).get("Level"))
    fields.update({
        "observed_mz": observed_mz, "observed_neutral_mass": neutral_mass,
        "candidate_exact_mass": target_exact_mass,
        "metabolon_level": "" if level.lower() in {"artifact", "exclude", "4?"} else level,
        "mass_alt": confmod.mass_consistent_alt(candidate_scores),
        "ms2_detail": {
            "present": n_peaks > 0, "n_pairs": n_peaks, "n_ion_bins": 0, "n_loss_bins": 0,
            "forward": number(rank1.get("forward_spectrum_score")),
            "predicted_v2": number(rank1.get("predicted_spectrum_score_v2")),
            "match": number(rank1.get("ms2_match_score")),
        },
        "ms2_n_pairs": n_peaks,
    })
    if naming_decision and naming_decision.get("exception_basis") == "standard_confirmed":
        fields["truth_match"] = True
    if naming_decision and naming_decision.get("library_collision") and curation_entry:
        if curation_entry.verdict == "curator_known":
            fields["curation_note"] = curation_entry.note()
            curation_entry = None
    confmod.apply_curation_entry(fields, curation_entry)
    return fields


def confidence_from_fields(fields: dict):
    """Respect the existing naming-v2 unresolved level as well as scored MSI."""
    confidence = confmod.compute_confidence(fields)
    if fields.get("unresolved_msi_level"):
        level = int(fields["unresolved_msi_level"])
        confidence = replace(confidence, msi_level=level, msi_label=confmod.MSI_LABELS[level],
                             msi_rationale=fields["unresolved_msi_rationale"], tier="Tentative")
    return confidence


def summarize_scores(control_id: str, fields: dict, candidates: list[dict], *, flags=()) -> dict:
    """Apply the existing confidence and publication rules without changing thresholds."""
    from tools import reported_identification as publication

    confidence = confidence_from_fields(fields)
    fields.update({
        "msi_level": confidence.msi_level, "msi_label": confidence.msi_label,
        "msi_rationale": confidence.msi_rationale, "evidence_score": confidence.evidence_score,
        "overall": confidence.overall, "tier": confidence.tier,
        "identity_status": confidence.identity_status,
        "identity_confidence": confidence.identity_confidence,
        "identity_basis": confidence.identity_basis,
    })
    # The publication workbook vocabulary differs from the ranked CSV schema.
    candidates = [dict(row, smiles=text(row.get("smiles") or row.get("candidate_smiles") or row.get("predicted_smiles")),
                       formula=text(row.get("formula") or row.get("candidate_formula") or row.get("predicted_formula")))
                  for row in candidates]
    call = publication.reported_identification(control_id, candidates, fields)
    candidate_confidence = publication.candidate_identity_confidence(candidates, fields)
    fields["reported_identification"] = call
    fields["reported_level"] = publication.reported_level(call)
    fields["candidate_identity_confidence"] = candidate_confidence
    all_flags = list(dict.fromkeys([*confidence.flags, *flags]))
    fields["scoring_flags"] = all_flags
    return {
        "se_msi_level": confidence.msi_level, "se_msi_label": confidence.msi_label,
        "se_msi_rationale": confidence.msi_rationale, "se_confidence": confidence.overall,
        "se_confidence_tier": confidence.tier, "se_evidence_score": confidence.evidence_score,
        "se_reported_level": publication.reported_level(call),
        "se_reported_identification": call,
        "se_candidate_identity_confidence": candidate_confidence,
        "se_confidence_flags_json": json.dumps(all_flags),
    }


def publication_cards_html(fields: dict) -> str:
    """Publication output accompanies, rather than replaces, the MSI rubric."""
    call = text(fields.get("reported_identification"))
    if not call:
        return ""
    score = number(fields.get("candidate_identity_confidence"))
    label = f"{score:.1f}%" if score is not None else "Unavailable"
    return (
        '<div class="metric-card"><div class="metric-label">Publication identification</div>'
        f'<div class="metric-value">{html.escape(call)}</div></div>'
        '<div class="metric-card"><div class="metric-label">Candidate identity confidence</div>'
        f'<div class="metric-value">{label}</div>'
        '<div class="metric-note">Publication tier and candidate-family margin</div></div>'
    )


def read_scoring_payload(path: Path | str, control_id: str) -> dict:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    if payload.get("schema_version") != 1 or payload.get("control_id") != control_id:
        raise ValueError("Scoring payload schema or control_id does not match the report")
    if payload.get("status") != "scored" or not payload.get("fields"):
        raise ValueError("Scoring payload does not contain completed confidence fields")
    return payload


def _reference_evidence(db_path: Path | None, smiles: str, rank1: dict, excluded_ids: set[str]) -> tuple[dict, list[dict], float | None]:
    from structural_elucidation.reports import compound
    from structural_elucidation.utils.denovo import inchikey_from_smiles

    detail: dict = {}
    scored: list[dict] = []
    if db_path and db_path.exists():
        identity = compound.resolve_db_identity_by_inchikey(db_path, inchikey_from_smiles(smiles))
        if identity:
            rows = compound.query_df(db_path, "SELECT * FROM metabolites WHERE chemical_id = ?", (identity[0],))
            detail = rows.iloc[0].to_dict() if not rows.empty else {}
        references = compound.query_df(db_path, "SELECT * FROM metabolites WHERE smiles IS NOT NULL AND trim(smiles) <> ''")
        _, target_fp = compound.mol_fp(smiles)
        for row in references.to_dict(orient="records"):
            if text(row.get("chemical_id")) in excluded_ids:
                continue
            reference_smiles = text(row.get("smiles")).split("|")[0].strip()
            _, fp = compound.mol_fp(reference_smiles)
            if target_fp is None or fp is None:
                continue
            similarity = float(compound.TanimotoSimilarity(target_fp, fp))
            if similarity >= min(0.28, pathmod.VOTE_SIMILARITY_FLOOR):
                scored.append({**row, "smiles": reference_smiles, "similarity_to_target": similarity})
    pathway = compound.resolve_pathway_display(detail, pd.Series(rank1), chemical_id_basis=pathmod.BASIS_INCHIKEY if detail else pathmod.BASIS_NONE)
    vote = pathmod.pathway_vote(pathmod.reference_rows_from_frame(scored))
    neighbors = sorted((row for row in scored if row["similarity_to_target"] >= 0.28), key=lambda row: -row["similarity_to_target"])[:10]
    return pathway, neighbors, vote.super_share if vote.available else None


def score_novel_feature(
    candidates: list[dict] | pd.DataFrame, metadata: dict, *,
    db_path: Path | None = None, peaks: pd.DataFrame | None = None,
    curation_ledger: Path | None = None,
) -> dict:
    """Score one feature; the caller persists this result before attempting HTML.

    Handoff names and prior matches are provenance only. Authentic-standard
    assurance can enter solely via an explicit synthetic-ID curation entry.
    """
    from tools import reported_identification as publication
    from tools.rescore_pure_nce import score_key
    from structural_elucidation.adapters.nce_bins import observed_bins_for_spectrum
    from structural_elucidation.reports import compound, naming_v2, tentative_naming as naming
    from structural_elucidation.utils.denovo import exact_mass_from_smiles, formula_from_smiles, inchikey_from_smiles

    control_id = text(metadata.get("control_id"))
    payload = {"schema_version": 1, "control_id": control_id, "status": "no_predictions",
               "summary": dict.fromkeys(SUMMARY_COLUMNS), "fields": {}, "candidates": [],
               "errors": [], "calibration_checksum": None}
    frame = candidates.copy() if isinstance(candidates, pd.DataFrame) else pd.DataFrame(candidates)
    if frame.empty:
        return payload
    if not control_id:
        raise ValueError("A synthetic control_id is required for novel scoring")
    frame["control_id"] = control_id
    frame = compound.predictions_to_candidate_pool(frame, control_id)
    if frame.empty:
        return payload
    for key in ("parent_mz", "lib_id", "ion_mode", "feature_key", "neutral_mass"):
        if metadata.get(key) is not None:
            frame[key] = metadata[key]
    if number(metadata.get("neutral_mass")) is not None:
        frame["target_neutral_mass"] = float(metadata["neutral_mass"])
    frame = frame.sort_values("rank", kind="mergesort")
    raw_peaks = metadata.get("peaks") or []
    if peaks is not None and not raw_peaks:
        mz_column = next((key for key in ("mz", "mass") if key in peaks), None)
        intensity_column = next((key for key in ("intensity", "rel_intensity") if key in peaks), None)
        if mz_column:
            raw_peaks = list(zip(peaks[mz_column], peaks[intensity_column] if intensity_column else [1.0] * len(peaks)))
    if peaks is None:
        peaks = pd.DataFrame(raw_peaks, columns=["mz", "intensity"])
        peaks["parent_mz"] = number(metadata.get("parent_mz"))
    rows = json_safe(frame.to_dict(orient="records"))
    flags: list[str] = []
    try:
        calibration = publication.active_calibration()
    except (OSError, ValueError, TypeError) as exc:
        calibration = None
        flags.append(f"Identity calibration unavailable: {exc}")
    if publication.CALIBRATION_PATH.exists():
        try:
            payload["calibration_checksum"] = hashlib.sha256(publication.CALIBRATION_PATH.read_bytes()).hexdigest()
        except OSError as exc:
            flags.append(f"Identity calibration checksum unavailable: {exc}")
    if calibration is None:
        flags.append("Identity calibration unavailable; using existing confidence bands")
    if not db_path or not Path(db_path).exists():
        flags.append("Reference database unavailable; using ranked pathway evidence")
    ion_json, loss_json = "[]", "[]"
    if raw_peaks:
        bins = observed_bins_for_spectrum(raw_peaks, float(metadata["parent_mz"]), int(metadata["lib_id"]))
        ion_json, loss_json = json.dumps(bins[0]), json.dumps(bins[1])
    for row in rows:
        status = text(row.get("nce_status"))
        score = number(row.get("pure_nce_score"))
        if score is None and not status:
            score, status = score_key((text(row.get("predicted_smiles")), ion_json, loss_json, int(metadata.get("lib_id") or row.get("lib_id") or 0)))
        if status in {"scored", "zero"} and number(score) is not None and 0 <= float(score) <= 1:
            row["pure_nce_score"] = float(score)
        elif not status and score is not None and 0 <= score <= 1:
            row["pure_nce_score"] = score
            status = "zero" if score == 0 else "scored"
        else:
            row["pure_nce_score"] = None
            flags.append(f"Pure NCE unavailable ({status or 'missing'}); evidence-score fallback where needed")
            if status == "error":
                payload["errors"].append({"stage": "pure_nce", "rank": row.get("rank"), "error": "matcher failed or unavailable"})
        row["nce_status"] = status
    payload["raw_candidates"] = rows
    known = {"Unknown": control_id, "Level": "",
             "Formula": text(metadata.get("formula")),
             "Tentative information": text(metadata.get("known_x_tentative_information")),
             "Cleaned External Info": text(metadata.get("known_x_cleaned_external_info"))}
    entry = confmod.load_curation_ledger(curation_ledger).get(control_id) if curation_ledger else None
    if db_path and Path(db_path).exists():
        index = naming_v2.CollisionIndex(db_path, supplied_names=[text(row.get("library_name")) for row in rows])
    else:
        # The absence is explicitly flagged; a missing database never creates a
        # positive match or grants laboratory assurance.
        index = naming_v2.CollisionIndex.__new__(naming_v2.CollisionIndex)
        index.by_connectivity, index.by_name, index.by_chemical_id = {}, {}, {}
        flags.append("Naming library-collision check unavailable without reference database")
    decision = naming_v2.decide(rows, index, known, assurance=entry.verdict if entry else "none")
    payload["naming_decision"] = decision
    payload["known_x_row"] = known
    unresolved = decision["identity_candidate"] is None
    if not unresolved:
        allowed = {row["model_rank"] for row in decision["candidate_annotations"] if row["identity_eligible"]}
        selected = int(decision["identity_candidate_rank"])
        rows = [dict(row, model_rank=row["rank"]) for row in rows if int(row["rank"]) in allowed]
        rows.sort(key=lambda row: (int(row["rank"]) != selected, int(row["rank"])))
        for rank, row in enumerate(rows, 1):
            row["rank"] = rank
    frame = pd.DataFrame(rows)
    rank1 = rows[0]
    smiles = text(rank1.get("predicted_smiles"))
    if not unresolved and number(exact_mass_from_smiles(smiles)) is None:
        raise ValueError("Rank-1 candidate has no valid molecular structure")
    score_input = compound.score_smiles_spectrum_alignment(frame.head(20), peaks, tolerance_da=1.0)
    score_frame = compound.annotate_adduct_hypotheses(compound.compact_candidate_scores(score_input), number(rank1.get("parent_mz")) or 0.0, compound._ion_mode_from_row(rank1))
    scores = json_safe(score_frame.to_dict(orient="records"))
    excluded_ids = {text(match.get("chemical_id")) for item in decision["excluded_library_matches"] for match in item.get("matches", [])}
    pathway, neighbors, vote_share = _reference_evidence(Path(db_path) if db_path else None, smiles, rank1, excluded_ids)
    fields = prepare_confidence_fields(rank1, scores, pathway, pathway_vote_share=vote_share,
        neighbors=neighbors, n_peaks=len(peaks), target_smiles=smiles,
        target_exact_mass=exact_mass_from_smiles(smiles), curation_entry=entry, naming_decision=decision)
    # No prior putative assignment may become a laboratory confirmation.
    fields["metabolon_level"] = ""
    fields.setdefault("assurance", "none")
    fields["smiles"] = smiles
    fields["formula"] = formula_from_smiles(smiles)
    fields["feature_name"] = ""
    identity = compound.resolve_db_identity_by_inchikey(Path(db_path), inchikey_from_smiles(smiles)) if db_path else None
    name_resolved = bool(identity)
    resolved_name = text(rank1.get("candidate_name") or rank1.get("base_name")) or (identity[1] if identity else control_id)
    tentative = naming.resolve(smiles, fields["formula"], rank1,
        naming.SqliteNameLookup(Path(db_path)) if db_path else {}, name_resolved=name_resolved,
        control_id=control_id, inchikey=inchikey_from_smiles(smiles))
    if tentative:
        resolved_name = tentative.label
        fields.update(naming.descriptor_fields(tentative))
    elif not name_resolved:
        resolved_name = naming.UNNAMED_STRUCTURE_LABEL
    fields["resolved_name"] = resolved_name
    fields["name_resolved"] = name_resolved
    fields.update({key: decision[key] for key in (
        "naming_policy_version", "name_source", "identity_candidate_rank", "library_collision",
        "identity_eligible", "identity_exclusion_reason", "excluded_library_matches",
        "raw_rank1_candidate", "identity_candidate", "exception_basis", "review_verdict")})
    if unresolved:
        hypothesis = decision["curated_hypothesis"]
        formula = known["Formula"] or text(rank1.get("target_formula"))
        level = 3 if hypothesis else 4 if formula else 5
        rationale = ("Compound-specific tentative X annotation; molecular structure remains unassigned" if hypothesis
                     else "No identity-eligible structure; reporting the annotated formula" if formula
                     else "No identity-eligible structure or molecular formula")
        unresolved_observation = {key: rank1.get(key) for key in ("parent_mz", "lib_id", "feature_key", "neutral_mass", "target_neutral_mass")}
        fields = prepare_confidence_fields(unresolved_observation, [], {"source": "predicted", "tier": "low"},
            pathway_vote_share=None, neighbors=[], n_peaks=len(peaks), target_smiles="", target_exact_mass=None)
        fields.update({"assurance": "none", "metabolon_level": "", "smiles": "", "formula": formula,
                       "resolved_name": hypothesis, "name_resolved": False,
                       "unresolved_msi_level": level, "unresolved_msi_rationale": rationale,
                       **{key: decision[key] for key in ("exception_basis", "identity_exclusion_reason", "library_collision")}})
        resolved_name = hypothesis
        if decision["excluded_library_matches"]:
            matches = decision["excluded_library_matches"][0].get("matches") or []
            fields["excluded_library_reference"] = matches[0]["name"] if matches else ""
    confidence = confidence_from_fields(fields)
    display_name = naming.display_name_for_level(resolved_name, msi_level=confidence.msi_level,
        formula=fields["formula"], observed_mz=number(rank1.get("parent_mz")), name_is_reference=name_resolved)
    if not unresolved and confidence.msi_level <= 3:
        bare = naming.bare_library_name(resolved_name)
        display_name = bare if confidence.msi_level == 1 else f"Probable {bare}" if confidence.msi_level == 2 else f"Possible {bare}"
    fields["display_name"] = display_name
    fields["report_display_name"] = display_name
    # Only explicit feature annotations, never prior_matched_name/compound_name.
    for key in ("known_x_cleaned_external_info", "known_x_tentative_information"):
        fields[key] = text(metadata.get(key))
    payload.update(status="scored", fields=fields, candidates=rows,
                   candidate_scores=scores, pathway_info=pathway)
    payload["summary"] = summarize_scores(control_id, fields, [] if unresolved else rows, flags=flags)
    return json_safe(payload)
