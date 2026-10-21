from __future__ import annotations

import pandas as pd
import pytest

from tools.build_validated_x_purchase_list import (
    apply_vendor_validation,
    candidate_consistency,
    formula_from_smiles,
    inchikey_from_smiles,
    load_detection,
    peptide_consistency,
    rank_candidates,
    validate_vendor_identity,
)

HYDROXY_ISOLEUCINE = "CC(O)C(C)C(N)C(=O)O"
UNRELATED_ESTER = "CC(=O)OCCNCCO"
HYDROXY_TRYPTOPHAN = "NC(Cc1c[nH]c2ccc(O)cc12)C(=O)O"


def test_detection_eligibility_uses_blood_or_urine_not_feces(tmp_path) -> None:
    workbook = tmp_path / "detection.xlsx"
    pd.DataFrame(
        [
            {
                "X Compound ID": "X-URINE",
                "Name": "amino acid derivative",
                "Top Match Formula": "C6H13NO3",
                "MSI Level": 3,
                "MSI Sub-level": "3A",
                "Human Plasma": 0.0,
                "Human Serum": 0.0,
                "Human Urine": 0.25,
                "Human Feces": 0.90,
            },
            {
                "X Compound ID": "X-FECES",
                "Name": "amino acid derivative",
                "Top Match Formula": "C6H13NO3",
                "MSI Level": 3,
                "MSI Sub-level": "3A",
                "Human Plasma": 0.0,
                "Human Serum": 0.0,
                "Human Urine": 0.0,
                "Human Feces": 0.90,
            },
        ]
    ).to_excel(workbook, sheet_name="Latest_Unique_X", index=False)
    detection = load_detection(workbook, "Latest_Unique_X", 10.0, 10.0).set_index("x_compound_id")
    assert bool(detection.loc["X-URINE", "frequency_eligible"])
    assert not bool(detection.loc["X-FECES", "frequency_eligible"])
    assert "human_feces" not in detection.columns


def detection_frame(*x_ids: str, msi_sublevel: object = "3A") -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "x_compound_id": x_id,
                "x_annotation_name": (
                    "5-hydroxynorleucine OR 5-hydroxyisoleucine"
                    if x_id == "X-11787"
                    else "amino acid derivative"
                ),
                "x_annotation_formula": "C6H13NO3",
                "msi_level": 3,
                "msi_sublevel": msi_sublevel,
                "human_plasma": 80.0,
                "human_serum": 0.0,
                "human_urine": 40.0,
                "max_blood_detection": 80.0,
                "urine_detection": 40.0,
                "frequency_eligible": True,
                "msi_sublevel_1": str(msi_sublevel) == "1",
            }
            for x_id in x_ids
        ]
    )


def truth_frame(*x_ids: str) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "x_compound_id": x_id,
                "resolved_known_name": "",
                "known_formula": "",
                "x_annotation_smiles": "",
            }
            for x_id in x_ids
        ]
    )


def candidate_frame(rows: list[tuple[str, int, int, str, float]]) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "x_compound_id": x_id,
                "rank": rank,
                "recommended_lib_id": str(lib),
                "recommended_predicted_smiles": smiles,
                "selection_evidence_score": score,
            }
            for x_id, lib, rank, smiles, score in rows
        ]
    )


def test_x11787_prefers_lib402_amino_acid_but_stays_unconfirmed() -> None:
    """The lib209 ester (CID 170749) must lose to the lib402 hydroxy amino acid.

    The annotation itself names two alternatives, so the winner is only
    ``ambiguous`` and never reaches the purchase list.
    """
    candidates = candidate_frame(
        [
            ("X-11787", 209, 1, UNRELATED_ESTER, 0.50),
            ("X-11787", 402, 2, HYDROXY_ISOLEUCINE, 0.25),
        ]
    )
    audit = rank_candidates(
        candidates,
        detection_frame("X-11787"),
        truth_frame("X-11787"),
        metabolon_keys=set(),
        min_evidence_score=0.05,
    )
    best = audit.sort_values("selection_score", ascending=False).iloc[0]
    assert best["recommended_lib_id"] == "402"
    assert best["annotation_structure_consistency"] == "ambiguous"
    ester = audit[audit["recommended_lib_id"].eq("209")].iloc[0]
    assert ester["annotation_structure_consistency"] == "mismatch"
    assert not audit["selected_for_vendor_validation"].any()


def test_x12100_hydroxytryptophan_is_not_accepted_as_oxindolylalanine() -> None:
    result = candidate_consistency(
        "beta-3-oxindolylalanine",
        "C11H12N2O3",
        "",
        HYDROXY_TRYPTOPHAN,
    )
    assert result["formula_agreement"] is True
    assert result["annotation_structure_consistency"] == "mismatch"


@pytest.mark.parametrize(
    "x_id,annotation",
    [
        ("X-11795", "amino acid derivative"),
        ("X-17438", "amino acid derivative"),
        ("X-18059", "amino acid derivative"),
        ("X-18921", "amino acid derivative"),
        ("X-23173", "amino acid derivative"),
        ("X-24400", "amino acid derivative"),
        ("X-25519", "amino acid derivative"),
    ],
)
def test_audited_formula_isomer_mismatches_fail(x_id: str, annotation: str) -> None:
    result = candidate_consistency(annotation, "C6H13NO3", "", UNRELATED_ESTER)
    assert result["formula_agreement"] is True
    assert result["annotation_structure_consistency"] == "mismatch", x_id


def test_formula_equality_is_not_sufficient() -> None:
    result = candidate_consistency("5-hydroxyisoleucine", "C6H13NO3", "", UNRELATED_ESTER)
    assert result == {
        "annotation_structure_consistency": "mismatch",
        "consistency_reason": "missing_annotation_class:amino_acid",
        "formula_agreement": True,
        "annotation_tanimoto": 0.0,
        "class_agreement": False,
    }


def test_x11795_heterogeneous_name_is_not_treated_as_confirmed_identity() -> None:
    result = candidate_consistency(
        "pantonine, 3-methyoxyvaline OR 2,4-Dihydroxy-3,3-dimethylbutanamide",
        "C6H13NO3",
        "",
        "CC(CO)NC(C)C(=O)O",
    )
    assert result["annotation_structure_consistency"] == "ambiguous"
    assert result["consistency_reason"] == "heterogeneous_annotation_alternatives"


def test_sulfate_annotation_rejects_carbon_bonded_sulfonate_isomer() -> None:
    sulfonate_smiles = "O=C(NCCS(=O)(=O)O)C(O)c1ccc(O)cc1"
    sulfonate = candidate_consistency(
        "N-acetyldopamine sulfate",
        formula_from_smiles(sulfonate_smiles),
        "",
        sulfonate_smiles,
    )
    assert sulfonate["annotation_structure_consistency"] == "mismatch"
    sulfate_smiles = "O=S(=O)(O)Oc1ccccc1"
    sulfate_ester = candidate_consistency(
        "phenol sulfate",
        formula_from_smiles(sulfate_smiles),
        "",
        sulfate_smiles,
    )
    assert sulfate_ester["annotation_structure_consistency"] == "plausible_analog"


DECANOYL_GLUCURONIDE = "CCCCCCCCCC(=O)OC1OC(C(=O)O)C(O)C(O)C1O"
PRO_GLY_GLY = "O=C(O)CNC(=O)CNC(=O)C1CCCN1"


def test_declared_peptide_residues_must_be_present() -> None:
    result = candidate_consistency(
        "Tripeptide with Pro, Asn",
        formula_from_smiles(PRO_GLY_GLY),
        "",
        PRO_GLY_GLY,
    )
    assert result["formula_agreement"] is True
    assert result["annotation_structure_consistency"] == "mismatch"
    assert result["consistency_reason"] == "missing_declared_residues:asn"


def test_declared_peptide_length_must_match() -> None:
    ok, reason = peptide_consistency("Dipeptide with Pro, Gly", PRO_GLY_GLY)
    assert ok is False
    assert reason == "peptide_length_mismatch:declared=2;observed=3"


def test_glucoside_annotation_rejects_uronic_acid_structure() -> None:
    glucoside = candidate_consistency(
        "glucoside of C16H28O8",
        formula_from_smiles(DECANOYL_GLUCURONIDE),
        "",
        DECANOYL_GLUCURONIDE,
    )
    assert glucoside["annotation_structure_consistency"] == "mismatch"
    assert "uronic_acid_but_annotation_says_glycoside" in glucoside["consistency_reason"]
    glucuronide = candidate_consistency(
        "glucuronide of C10H20O2",
        formula_from_smiles(DECANOYL_GLUCURONIDE),
        "",
        DECANOYL_GLUCURONIDE,
    )
    assert glucuronide["annotation_structure_consistency"] == "plausible_analog"


def test_glucuronide_annotation_requires_uronic_acid() -> None:
    acetylated_disaccharide = "CC(=O)OCC1OC(OC2OCC(O)C(O)C2O)C(O)C(O)C1O"
    result = candidate_consistency(
        "glucuronide of C8H16O6",
        formula_from_smiles(acetylated_disaccharide),
        "",
        acetylated_disaccharide,
    )
    assert result["annotation_structure_consistency"] == "mismatch"
    assert "missing_uronic_acid" in result["consistency_reason"]


def test_glycine_conjugate_annotation_requires_glycine() -> None:
    without_glycine = candidate_consistency(
        "glycine conjugate of C10H12O2",
        formula_from_smiles("CCCCc1ccc(OCC(=O)O)cc1"),
        "",
        "CCCCc1ccc(OCC(=O)O)cc1",
    )
    assert without_glycine["annotation_structure_consistency"] == "mismatch"
    assert "missing_glycine_conjugate" in without_glycine["consistency_reason"]


def test_unannotated_carnitine_motif_is_a_conflict() -> None:
    carnitine_ester = "C[N+](C)(C)CC(OC(=O)c1ccccc1)CC(=O)[O-]"
    result = candidate_consistency("small peptide", "", "", carnitine_ester)
    assert result["annotation_structure_consistency"] == "mismatch"
    assert "unexpected_quaternary_ammonium" in result["consistency_reason"]


def test_or_annotations_are_ambiguous_regardless_of_class_match() -> None:
    result = candidate_consistency(
        "5-hydroxynorleucine OR 5-hydroxyisoleucine",
        "C6H13NO3",
        "",
        HYDROXY_ISOLEUCINE,
    )
    assert result["annotation_structure_consistency"] == "ambiguous"
    assert result["consistency_reason"] == "heterogeneous_annotation_alternatives"


def test_msi_1_and_existing_library_candidates_are_excluded() -> None:
    candidate = candidate_frame([("X-11787", 402, 1, HYDROXY_ISOLEUCINE, 0.4)])
    key14 = inchikey_from_smiles(HYDROXY_ISOLEUCINE).split("-")[0]
    audit = rank_candidates(
        candidate,
        detection_frame("X-11787", msi_sublevel=1),
        truth_frame("X-11787"),
        metabolon_keys={key14},
        min_evidence_score=0.05,
    )
    reason = audit.iloc[0]["pre_vendor_exclusion_reason"]
    assert "msi_sublevel_1" in reason
    assert "already_in_metabolon_library" in reason
    assert not audit.iloc[0]["selected_for_vendor_validation"]


def test_duplicate_inchikeys_across_x_ids_are_flagged() -> None:
    candidates = candidate_frame(
        [
            ("X-A", 402, 1, HYDROXY_ISOLEUCINE, 0.4),
            ("X-B", 402, 1, HYDROXY_ISOLEUCINE, 0.4),
        ]
    )
    detection = detection_frame("X-A", "X-B")
    detection["x_annotation_name"] = "5-hydroxyisoleucine"
    audit = rank_candidates(candidates, detection, truth_frame("X-A", "X-B"), set(), 0.05)
    assert audit["prediction_reused"].all()
    assert set(audit["prediction_reused_by_x_ids"]) == {"X-A|X-B"}


def test_pubchem_requires_full_inchikey_formula_name_and_vendor() -> None:
    full_key = inchikey_from_smiles(HYDROXY_ISOLEUCINE)
    candidate = {
        "recommended_predicted_inchikey": full_key,
        "recommended_predicted_formula": "C6H13NO3",
    }
    valid_record = {
        "pubchem_cid": "123",
        "pubchem_name": "5-hydroxyisoleucine",
        "pubchem_inchikey": full_key,
        "pubchem_molecular_formula": "C6H13NO3",
        "vendor_sources": "Example Vendor",
        "vendor_count": 1,
    }
    assert validate_vendor_identity(candidate, valid_record)[0]
    for missing in ("pubchem_name", "pubchem_inchikey", "pubchem_molecular_formula"):
        invalid = dict(valid_record)
        invalid[missing] = ""
        assert not validate_vendor_identity(candidate, invalid)[0]
    mismatch = dict(valid_record, pubchem_inchikey="AAAAAAAAAAAAAA-BBBBBBBBBB-C")
    assert validate_vendor_identity(candidate, mismatch) == (False, "pubchem_inchikey_mismatch")


def test_cached_connectivity_only_record_cannot_be_recommended() -> None:
    candidate = candidate_frame([("X-CACHE", 402, 1, HYDROXY_ISOLEUCINE, 0.4)])
    audit = rank_candidates(
        candidate,
        detection_frame("X-CACHE"),
        truth_frame("X-CACHE"),
        set(),
        0.05,
    )
    key14 = inchikey_from_smiles(HYDROXY_ISOLEUCINE).split("-")[0]
    old_cache = {
        key14: {
            "pubchem_cid": "123",
            "pubchem_compound_name": "5-hydroxyisoleucine",
            "vendor_sources": "Example Vendor",
            "vendor_count": 1,
        }
    }
    validated = apply_vendor_validation(audit, old_cache, "cached", 1.0, 0.0)
    assert validated.iloc[0]["pubchem_vendor_listed"]
    assert not validated.iloc[0]["purchase_recommended"]
    assert validated.iloc[0]["exclusion_reason"] == "cached_pubchem_inchikey_missing"
