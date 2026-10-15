"""Tests for display-only tentative structural descriptors."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from structural_elucidation.reports import tentative_naming as TN
from tools import patch_report_tentative_names as patcher


def test_peptide_ordered_branch() -> None:
    result = TN.resolve(
        "CC(=O)NCC(=O)NCC(=O)O",
        "C8H15N3O4",
        {
            "is_peptide_candidate": True,
            "candidate_ordered_sequence": "Ala-Gly-Ser",
            "peptide_connectivity_confidence": 0.91,
            "peptide_branch_or_cycle_penalty": 0.0,
        },
        {},
        name_resolved=False,
    )
    assert result is not None
    assert result.kind == "peptide"
    assert result.label == "Tripeptide (Ala–Gly–Ser)"
    assert result.confidence == "high"


def test_peptide_composition_wins_when_ordered_parse_is_incomplete() -> None:
    result = TN.resolve(
        "",
        "C8H14N4O5",
        {
            "is_peptide_candidate": True,
            "candidate_ordered_sequence": "Gly-Gly",
            "normalized_peptide_composition": "Asn-Gly-Gly",
            "sequence_matches_declared_composition": False,
            "peptide_connectivity_confidence": 1.0,
        },
        {},
        name_resolved=False,
    )
    assert result is not None
    assert result.label == "Tripeptide (Asn/Gly/Gly composition)"


def test_long_peptide_uses_compact_composition() -> None:
    result = TN.resolve(
        "",
        "C48H77N21O16",
        {
            "is_peptide_candidate": True,
            "normalized_peptide_composition": "Ala-Ala-Arg-Arg-Arg-Gln-Gly-Gly-Pyr-Pyr-Pyr",
            "peptide_connectivity_confidence": 0.9,
        },
        {},
        name_resolved=False,
    )
    assert result is not None
    assert result.label == "11-residue peptide (Ala/Arg/Gln/Gly/Pyr composition)"


@pytest.mark.skipif(not TN.HAS_RDKIT, reason="RDKit is not installed")
def test_lipid_branch() -> None:
    smiles = "CCCCCCCCCCCCCCCC(=O)OCC(COP(=O)(O)OCC[N+](C)(C)C)OC(=O)CCCCCCC/C=C\\CCCCCCCC"
    result = TN.resolve(smiles, "C42H82NO8P", {}, {}, name_resolved=False)
    assert result is not None
    assert result.kind == "lipid"
    assert result.label == "Phosphatidylcholine (34:1)"
    assert result.confidence == "high"


CAFFEIC_GLUCURONIDE = "O=C(O)/C=C/c1ccc(OC2OC(C(=O)O)C(O)C(O)C2O)cc1O"
TYROSOL_SULFATE = "O=S(=O)(O)OCCc1ccc(O)cc1"
O_CRESOL_DISULFATE = "Cc1ccccc1OS(=O)(=O)OS(=O)(=O)O"
MEVALONATE_DIPHOSPHATE = "CC(O)(CCOP(=O)(O)OP(=O)(O)O)CC(=O)O"
AMBIGUOUS_PHOSPHATE = "N=c1nc(O)c2c([nH]1)NC1OC3COP(=O)(O)OC3C(O)(O)C1N2"
BRIDGED_PHOSPHATE = "CCCCCC/C=C\\C/C=C\\C/C=C\\CCCCCCCC(=O)OC(CO)COP(=O)(O)OCCCCCC"


@pytest.mark.skipif(not TN.HAS_RDKIT, reason="RDKit is not installed")
def test_conjugate_resolved_branch() -> None:
    parent = TN._reverse_products(CAFFEIC_GLUCURONIDE, "glucuronide")[0]
    parent_inchikey = TN._mol_inchikey(parent)
    db_names = {parent_inchikey.split("-")[0]: {"name": "caffeic acid"}}
    result = TN.resolve(
        CAFFEIC_GLUCURONIDE,
        "C15H16O10",
        {},
        db_names,
        name_resolved=False,
    )
    assert result is not None
    assert result.kind == "conjugate_resolved"
    assert result.label == "Caffeic acid glucuronide"
    assert result.parent_inchikey == parent_inchikey
    assert result.parent_formula == "C9H8O4"
    assert result.conjugate_count == 1
    assert result.confidence == "high"


def test_resolved_conjugate_avoids_repeated_class() -> None:
    assert (
        TN._resolved_conjugate_base(("glucuronide",), "equol glucuronide")
        == "equol diglucuronide"
    )


@pytest.mark.skipif(not TN.HAS_RDKIT, reason="RDKit is not installed")
def test_conjugate_unresolved_branch() -> None:
    result = TN.resolve(
        CAFFEIC_GLUCURONIDE,
        "C15H16O10",
        {},
        {},
        name_resolved=False,
    )
    assert result is not None
    assert result.kind == "conjugate_unresolved"
    assert result.label == "Glucuronide conjugate of C9H8O4"
    assert result.parent_formula == "C9H8O4"
    assert result.conjugate_count == 1
    assert result.confidence == "medium"


@pytest.mark.skipif(not TN.HAS_RDKIT, reason="RDKit is not installed")
def test_sulfate_resolves_equal_carbon_parent() -> None:
    parents, count = TN._fully_reverse_carbon_free_conjugate(
        TYROSOL_SULFATE,
        "sulfate",
    )
    assert count == 1
    parent = next(value for value in parents if TN._mol_formula(value) == "C8H10O2")
    parent_inchikey = TN._mol_inchikey(parent)
    result = TN.resolve(
        TYROSOL_SULFATE,
        "C8H10O5S",
        {},
        {parent_inchikey.split("-")[0]: {"name": "tyrosol"}},
        name_resolved=False,
    )
    assert result is not None
    assert result.kind == "conjugate_resolved"
    assert result.label == "Tyrosol sulfate"
    assert result.parent_name == "tyrosol"
    assert result.parent_formula == "C8H10O2"
    assert result.conjugate_count == 1
    assert "C8H10O5S" not in result.label


@pytest.mark.skipif(not TN.HAS_RDKIT, reason="RDKit is not installed")
def test_sulfate_uses_parent_formula_when_unresolved() -> None:
    result = TN.resolve(
        TYROSOL_SULFATE,
        "C8H10O5S",
        {},
        {},
        name_resolved=False,
    )
    assert result is not None
    assert result.kind == "conjugate_unresolved"
    assert result.label == "Sulfate conjugate of C8H10O2"
    assert result.parent_formula == "C8H10O2"


@pytest.mark.skipif(not TN.HAS_RDKIT, reason="RDKit is not installed")
def test_sulfate_multiplicity_counts_groups_not_unrelated_sulfur() -> None:
    """X-25954 has one sulfate plus a thioether, not two sulfate conjugates."""

    result = TN.resolve(
        "COCCCC(=NOS(=O)(=O)O)SC1OC(CO)C(O)C(O)C1O",
        "C11H21NO10S2",
        {},
        {},
        name_resolved=False,
    )
    assert result is not None
    assert result.label == "Sulfate conjugate of C11H21NO7S"
    assert result.conjugate_count == 1
    assert "C8H10O5S" not in result.label


@pytest.mark.skipif(not TN.HAS_RDKIT, reason="RDKit is not installed")
@pytest.mark.parametrize(
    ("smiles", "motif", "parent_formula", "parent_name", "expected"),
    [
        (
            O_CRESOL_DISULFATE,
            "sulfate",
            "C7H8O",
            "o-cresol",
            "O-cresol disulfate",
        ),
        (
            MEVALONATE_DIPHOSPHATE,
            "phosphate",
            "C6H12O4",
            "mevalonate",
            "Mevalonate diphosphate",
        ),
    ],
)
def test_multi_conjugates_use_fully_deconjugated_named_parent(
    smiles: str,
    motif: str,
    parent_formula: str,
    parent_name: str,
    expected: str,
) -> None:
    parents, count = TN._fully_reverse_carbon_free_conjugate(smiles, motif)
    parent = next(value for value in parents if TN._mol_formula(value) == parent_formula)
    parent_inchikey = TN._mol_inchikey(parent)
    db_names = {parent_inchikey.split("-")[0]: {"name": parent_name}}
    for candidate in parents:
        if candidate != parent:
            candidate_ik = TN._mol_inchikey(candidate)
            db_names[candidate_ik.split("-")[0]] = {
                "name": f"{parent_name} {motif}"
            }
    result = TN.resolve(smiles, "", {}, db_names, name_resolved=False)
    assert result is not None
    assert result.kind == "conjugate_resolved"
    assert result.label == expected
    assert result.conjugate_count == 2
    assert result.parent_formula == parent_formula


@pytest.mark.skipif(not TN.HAS_RDKIT, reason="RDKit is not installed")
def test_ambiguous_parents_with_one_formula_do_not_claim_an_inchikey() -> None:
    result = TN.resolve(
        AMBIGUOUS_PHOSPHATE,
        "C10H14N5O8P",
        {},
        {},
        name_resolved=False,
    )
    assert result is not None
    assert result.kind == "conjugate_unresolved"
    assert result.label == "Phosphate conjugate of C10H15N5O5"
    assert result.parent_formula == "C10H15N5O5"
    assert result.parent_inchikey is None


@pytest.mark.skipif(not TN.HAS_RDKIT, reason="RDKit is not installed")
def test_phosphate_reversal_prefers_the_larger_parent_scaffold() -> None:
    result = TN.resolve(
        BRIDGED_PHOSPHATE,
        "C31H57O7P",
        {},
        {},
        name_resolved=False,
    )
    assert result is not None
    assert result.kind == "conjugate_unresolved"
    assert result.label == "Phosphate conjugate of C25H44O4"
    assert result.parent_formula == "C25H44O4"


@pytest.mark.skipif(not TN.HAS_RDKIT, reason="RDKit is not installed")
def test_carbon_bearing_conjugate_still_requires_carbon_reduction(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        TN,
        "_reaction_resources",
        lambda: (
            {"glucuronide": ("test",)},
            {"test": "[O:1]>>[O:1]S(=O)(=O)O"},
        ),
    )
    assert TN._reverse_products(TYROSOL_SULFATE, "glucuronide") == []


@pytest.mark.skipif(not TN.HAS_RDKIT, reason="RDKit is not installed")
def test_motif_only_branch(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(TN, "_reverse_products", lambda *_args: [])
    result = TN.resolve("O=S(=O)(O)O", "H2O4S", {}, {}, name_resolved=False)
    assert result is not None
    assert result.kind == "motif_only"
    assert result.label == "Sulfated compound"
    assert result.confidence == "low"


def test_composite_motif_is_readable_and_limited() -> None:
    result = TN._try_motif(
        "",
        ("carnitine", "sulfur_oxygen_rich", "phosphate_like"),
    )
    assert result is not None
    assert result.label == "Carnitine-conjugated and phosphate-like compound"
    assert "+" not in result.label


def test_fallback_and_named_compound_branches() -> None:
    assert TN.resolve("CCO", "C2H6O", {}, {}, name_resolved=False) is None
    assert TN.resolve(
        CAFFEIC_GLUCURONIDE,
        "C15H16O10",
        {},
        {},
        name_resolved=True,
    ) is None


def test_rdkit_absent_rank1_peptide_shim(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(TN, "HAS_RDKIT", False)
    result = TN.resolve(
        "",
        "C8H15N3O4",
        {
            "is_peptide_candidate": True,
            "peptide_composition_label": "tripeptide: Ala, Gly, Ser",
        },
        {},
        name_resolved=False,
    )
    assert result is not None
    assert result.kind == "peptide"
    assert result.label == "Tripeptide (Ala/Gly/Ser composition)"
    assert result.source == "rank1_csv"


def test_rdkit_absent_sulfate_does_not_invent_a_parent(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(TN, "HAS_RDKIT", False)
    assert TN.resolve(
        TYROSOL_SULFATE,
        "C8H10O5S",
        {},
        {},
        name_resolved=False,
    ) is None


def _make_bundle(root: Path, *, peptide: bool = True, name_resolved: bool = False) -> Path:
    report = root / "X-TEST"
    data = report / "data"
    data.mkdir(parents=True)
    formula = "C8H15N3O4" if peptide else "C2H6O"
    smiles = "CC(=O)NCC(=O)NCC(=O)O" if peptide else "CCO"
    fields = {
        "overall": 0.42,
        "identity_confidence": None,
        "identity_status": "unknown",
        "msi_level": 4,
        "tier": "Low",
        "display_name": f"De novo structure ({formula})",
        "name_resolved": name_resolved,
    }
    (data / "x-test_confidence_inputs.json").write_text(json.dumps(fields, indent=2))
    evidence = {
        "target": {
            "control_id": "X-TEST",
            "name": "X-TEST",
            "smiles": smiles,
            "formula": formula,
        }
    }
    (data / "x-test_evidence_summary.json").write_text(json.dumps(evidence, indent=2))
    (report / "manifest.json").write_text(json.dumps({"target": evidence["target"]}, indent=2))
    with (data / "x-test_top_candidates.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "rank",
                "predicted_smiles",
                "is_peptide_candidate",
                "candidate_ordered_sequence",
                "peptide_connectivity_confidence",
            ],
        )
        writer.writeheader()
        writer.writerow(
            {
                "rank": 1,
                "predicted_smiles": smiles,
                "is_peptide_candidate": peptide,
                "candidate_ordered_sequence": "Ala-Gly-Ser" if peptide else "",
                "peptide_connectivity_confidence": 0.9 if peptide else "",
            }
        )
    report_html = f"""<!doctype html><html><head>
<title>X-TEST Structure Resolution (Public)</title></head><body>
<h1>X-TEST</h1><p>Structure-resolution report for <strong>X-TEST</strong>,
formula <code>{formula}</code>.</p>
<span class="badge ok">top-1 candidate</span>
<p><strong>De novo structure ({formula})</strong></p></body></html>"""
    for filename in ("index.html", "public.html", "internal.html", "X_TEST.html"):
        (report / filename).write_text(report_html)
    (report / "guides.html").write_text("<html><h1>X-TEST guides</h1></html>")
    return report


def test_patcher_is_idempotent_and_preserves_confidence(tmp_path: Path) -> None:
    report = _make_bundle(tmp_path)
    first = patcher.patch_report(report, {}, apply=True)
    assert first["changed"] is True
    assert first["summary"]["descriptor_kind"] == "peptide"

    fields = json.loads((report / "data/x-test_confidence_inputs.json").read_text())
    assert fields["display_name"] == "Tripeptide (Ala–Gly–Ser)"
    assert fields["descriptor_parent_formula"] is None
    assert fields["descriptor_conjugate_count"] == 0
    assert fields["overall"] == 0.42
    assert fields["identity_confidence"] is None
    assert fields["msi_level"] == 4
    for filename in ("index.html", "public.html", "internal.html", "X_TEST.html"):
        page = (report / filename).read_text()
        assert "<h1>Tripeptide (Ala–Gly–Ser)</h1>" in page
        assert "<title>Tripeptide (Ala–Gly–Ser) Structure Resolution (Public)</title>" in page
        assert "Tentative descriptor · High" in page
        assert "C8H15N3O4" in page
        assert "— tentative" not in page
    assert (report / "guides.html").read_text() == "<html><h1>X-TEST guides</h1></html>"
    manifest = json.loads((report / "manifest.json").read_text())
    assert manifest["target"]["name"] == "Tripeptide (Ala–Gly–Ser)"

    second = patcher.patch_report(report, {}, apply=True)
    assert second["changed"] is False
    assert second["title_changed"] is False
    assert second["summary"]["new_title"] == "Tripeptide (Ala–Gly–Ser)"


def test_patcher_cleans_fallback_without_badge(tmp_path: Path) -> None:
    report = _make_bundle(tmp_path, peptide=False)
    for path in report.glob("*.html"):
        if path.name != "guides.html":
            path.write_text(
                path.read_text().replace(
                    "De novo structure (C8H15N3O4)", "1-methyladenine"
                )
            )
    outcome = patcher.patch_report(report, {}, apply=True)
    assert outcome["changed"] is True
    fields = json.loads((report / "data/x-test_confidence_inputs.json").read_text())
    assert fields["display_name"] == "Tentative Structure"
    assert fields["descriptor_kind"] is None
    page = (report / "index.html").read_text()
    assert "<title>X-TEST Structure Resolution (Public)</title>" in page
    assert "<h1>X-TEST</h1>" in page
    assert "<p><strong>Tentative Structure</strong></p>" in page
    assert "Tentative descriptor" not in page


def test_patcher_qualifies_named_bundle_as_an_isomer(tmp_path: Path) -> None:
    report = _make_bundle(tmp_path, name_resolved=True)
    data = report / "data"
    confidence_path = data / "x-test_confidence_inputs.json"
    fields = json.loads(confidence_path.read_text())
    fields.update(
        {
            "msi_level": 2,
            "display_name": "ethanol",
            "resolved_name": "ethanol",
        }
    )
    confidence_path.write_text(json.dumps(fields, indent=2))
    for filename in ("x-test_evidence_summary.json",):
        path = data / filename
        payload = json.loads(path.read_text())
        payload["target"]["name"] = "ethanol"
        path.write_text(json.dumps(payload, indent=2))
    manifest_path = report / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["target"]["name"] = "ethanol"
    manifest_path.write_text(json.dumps(manifest, indent=2))
    for path in report.glob("*.html"):
        if path.name != "guides.html":
            path.write_text(
                path.read_text()
                .replace("X-TEST Structure Resolution", "ethanol Structure Resolution")
                .replace("<h1>X-TEST</h1>", "<h1>ethanol</h1>")
                .replace("De novo structure (C8H15N3O4)", "ethanol")
            )

    outcome = patcher.patch_report(report, {}, apply=True)
    assert outcome["status"] == "library_isomer"
    assert outcome["changed"] is True
    assert json.loads(confidence_path.read_text())["display_name"] == "Isomer of ethanol"
    assert json.loads(manifest_path.read_text())["target"]["name"] == "Isomer of ethanol"
    # A named bundle carries no synthesized descriptor, so the page still heads
    # with the compound id and the qualified name appears as the caption.
    page = (report / "public.html").read_text()
    assert "<h1>X-TEST</h1>" in page
    assert "<p><strong>Isomer of ethanol</strong></p>" in page

    second = patcher.patch_report(report, {}, apply=True)
    assert second["changed"] is False


def test_index_humanizes_descriptor_kind() -> None:
    from structural_elucidation.reports.index_render import IndexRow, render_master_index

    page = render_master_index(
        [
            IndexRow(
                control_id="X-TEST",
                report="X-TEST/index.html",
                descriptor_kind="conjugate_unresolved",
                descriptor_confidence="medium",
            )
        ],
        "now",
    )
    assert "Conjugate · Medium" in page


def test_candidate_name_from_rank1() -> None:
    result = TN.resolve(
        "",
        "C19H32O8S2",
        {"candidate_name": "5-androstenetriol 3,17-disulfate"},
        {},
        name_resolved=False,
        control_id="X-21470",
    )
    assert result is not None
    assert result.kind == "candidate_name"
    assert result.label == "Isomer of 5-androstenetriol 3,17-disulfate"
    assert result.source == "rank1_candidate_name"
    fields = TN.descriptor_fields(result)
    assert fields["descriptor_source_row"] == "rank1"


def test_candidate_name_msi1_truth_match_has_no_tentative_suffix() -> None:
    result = TN.resolve(
        "",
        "C19H32O8S2",
        {"candidate_name": "5-androstenetriol 3,17-disulfate"},
        {},
        name_resolved=False,
        control_id="X-21470",
        msi_level="1",
        truth_match=True,
    )
    assert result is not None
    assert result.label == "Isomer of 5-androstenetriol 3,17-disulfate"
    assert "tentative" not in result.label.lower()


def test_candidate_name_equal_to_control_id_falls_through(tmp_path: Path) -> None:
    report = tmp_path / "X-TEST"
    data = report / "data"
    data.mkdir(parents=True)
    (data / "x-test_top_candidates.csv").write_text(
        "rank,candidate_name\n1,X-TEST\n",
        encoding="utf-8",
    )
    result = TN.resolve(
        "CC(=O)NCC(=O)NCC(=O)O",
        "C8H15N3O4",
        {"rank": 1, "candidate_name": "X-TEST", "is_peptide_candidate": True, "candidate_ordered_sequence": "Ala-Gly-Ser"},
        {},
        name_resolved=False,
        bundle_dir=report,
        control_id="X-TEST",
    )
    assert result is not None
    assert result.kind == "peptide"


def test_candidate_name_from_known_x_row(tmp_path: Path) -> None:
    report = tmp_path / "X-TEST"
    data = report / "data"
    data.mkdir(parents=True)
    (data / "x-test_known_x_row.csv").write_text(
        "Level,External Information (Class or modification)\n"
        '1,"KNOWN 5-androstenetriol disulfate (in a Metabolon patent application?)"\n',
        encoding="utf-8",
    )
    result = TN.resolve(
        "",
        "C19H32O8S2",
        {},
        {},
        name_resolved=False,
        bundle_dir=report,
        control_id="X-TEST",
    )
    assert result is not None
    assert result.kind == "candidate_name"
    assert result.label == "5-androstenetriol disulfate"
    assert result.source == "known_x_row"
    assert TN.descriptor_fields(result)["descriptor_source_row"] == "known_x_row"


def test_candidate_name_html_alt_fallback(tmp_path: Path) -> None:
    report = tmp_path / "X-TEST"
    report.mkdir(parents=True)
    (report / "public.html").write_text(
        "<img alt='2&#x27;,3&#x27;-anhydroinosine structure'>",
        encoding="utf-8",
    )
    result = TN.resolve(
        "",
        "C10H10N4O4",
        {},
        {},
        name_resolved=False,
        bundle_dir=report,
        control_id="X-TEST",
        html_alt_name="2',3'-anhydroinosine",
    )
    assert result is not None
    assert result.kind == "candidate_name"
    assert result.source == "html_alt"
    assert "anhydroinosine" in result.label


def test_specificity_guard_preserves_specific_name(tmp_path: Path) -> None:
    report = tmp_path / "X-21470"
    data = report / "data"
    data.mkdir(parents=True)
    formula = "C19H32O8S2"
    fields = {
        "overall": 0.97,
        "identity_confidence": 0.97,
        "identity_status": "known_pending",
        "msi_level": 1,
        "tier": "High",
        "display_name": "5-androstenetriol 3,17-disulfate",
        "name_resolved": False,
        "truth_match": True,
    }
    (data / "x-21470_confidence_inputs.json").write_text(json.dumps(fields, indent=2))
    evidence = {"target": {"control_id": "X-21470", "name": "5-androstenetriol 3,17-disulfate", "formula": formula}}
    (data / "x-21470_evidence_summary.json").write_text(json.dumps(evidence, indent=2))
    (report / "manifest.json").write_text(json.dumps({"target": evidence["target"]}, indent=2))
    (data / "x-21470_top_candidates.csv").write_text(
        "rank,is_lipid_candidate,lipid_class\n1,true,sterol\n",
        encoding="utf-8",
    )
    (report / "public.html").write_text(
        "<html><img alt='X-21470 structure'></html>",
        encoding="utf-8",
    )
    outcome = patcher.patch_report(report, {}, apply=False)
    assert outcome["status"] == "preserved"
    assert outcome["changed"] is False


def test_specificity_guard_allows_patch_from_de_novo(tmp_path: Path) -> None:
    report = _make_bundle(tmp_path)
    data = report / "data"
    fields = json.loads((data / "x-test_confidence_inputs.json").read_text())
    fields["display_name"] = fields["display_name"]
    fields["descriptor_kind"] = None
    (data / "x-test_confidence_inputs.json").write_text(json.dumps(fields, indent=2))
    outcome = patcher.patch_report(report, {}, apply=True)
    assert outcome["status"] == "descriptor"
    assert outcome["changed"] is True
    updated = json.loads((data / "x-test_confidence_inputs.json").read_text())
    assert updated["display_name"] != "De novo structure (C8H15N3O4)"


def test_candidate_name_strips_released_tag() -> None:
    result = TN.resolve(
        "",
        "C16H26N4O6",
        {"candidate_name": "Pyr-leu-gln - RELEASED"},
        {},
        name_resolved=False,
        control_id="X-25744",
        msi_level="1",
        truth_match=True,
    )
    assert result is not None
    assert result.label == "Isomer of Pyr-leu-gln"
    assert result.parent_name == "Pyr-leu-gln"


def test_candidate_name_strips_released_structure_suffix() -> None:
    result = TN.resolve(
        "",
        "C10H18O4",
        {},
        {},
        name_resolved=False,
        control_id="X-22102",
        msi_level="1",
        truth_match=True,
        html_alt_name="lactoyl-isoleucine released structure",
    )
    assert result is not None
    assert result.label == "Isomer of lactoyl-isoleucine"
    assert result.source == "html_alt"


def test_known_x_row_strips_released_tag(tmp_path: Path) -> None:
    report = tmp_path / "X-TEST"
    data = report / "data"
    data.mkdir(parents=True)
    (data / "x-test_known_x_row.csv").write_text(
        "Level,External Information (Class or modification)\n"
        '1,"KNOWN N-lactoylhistidine - RELEASED"\n',
        encoding="utf-8",
    )
    result = TN.resolve(
        "",
        "C9H13N3O4",
        {},
        {},
        name_resolved=False,
        bundle_dir=report,
        control_id="X-TEST",
        msi_level="1",
        truth_match=True,
    )
    assert result is not None
    assert result.label == "N-lactoylhistidine"


def test_index_renders_candidate_name_descriptor() -> None:
    from structural_elucidation.reports.index_render import IndexRow, render_master_index

    page = render_master_index(
        [
            IndexRow(
                control_id="X-21470",
                report="X-21470/index.html",
                descriptor_kind="candidate_name",
                descriptor_confidence="high",
            )
        ],
        "now",
    )
    assert "Candidate name · High" in page


def test_candidate_name_retains_chemical_parenthetical() -> None:
    assert TN._clean_candidate_display_name("acetylcarnitine (C2)") == "acetylcarnitine (C2)"
    assert (
        TN._clean_candidate_display_name("stearoyl sphingomyelin (d18:1/18:0)")
        == "stearoyl sphingomyelin (d18:1/18:0)"
    )
    result = TN.resolve(
        "",
        "C9H17NO4",
        {},
        {},
        name_resolved=False,
        control_id="X-15472",
        html_alt_name="acetylcarnitine (C2)",
    )
    assert result is not None
    assert result.label == "Isomer of acetylcarnitine (C2)"
    assert result.source == "html_alt"


def test_candidate_name_strips_curator_notes_only() -> None:
    assert (
        TN._clean_candidate_display_name("5-androsten-3b,16b,17a-triol sulfate (1) (in patent application?)")
        == "5-androsten-3b,16b,17a-triol sulfate"
    )
    assert TN._clean_candidate_display_name("cortolone glucuronide (2)") == "cortolone glucuronide"


def test_candidate_name_strips_stray_trailing_paren() -> None:
    cleaned = TN._clean_candidate_display_name(
        "N4-acetyl-L-2,4-diaminobutyrate ((2S)-4-(acetylamino)-2-aminobutanoic acid))"
    )
    assert cleaned.endswith(")")
    assert not cleaned.endswith("))")


def test_bare_generic_lipid_class_not_adopted() -> None:
    result = TN.resolve(
        "",
        "C42H82NO8P",
        {},
        {},
        name_resolved=False,
        control_id="X-24893",
        html_alt_name="phosphatidylcholine",
    )
    assert result is None or result.kind != "candidate_name"


def test_qualified_lipid_class_is_adopted() -> None:
    result = TN.resolve(
        "",
        "C42H82NO8P",
        {},
        {},
        name_resolved=False,
        control_id="X-24893",
        html_alt_name="phosphatidylcholine (18:0/20:5, 16:0/22:5n6)",
    )
    assert result is not None
    assert result.kind == "candidate_name"
    assert result.label == "Isomer of phosphatidylcholine (18:0/20:5, 16:0/22:5n6)"


def test_descriptive_known_x_phrase_skipped(tmp_path: Path) -> None:
    report = tmp_path / "X-26062"
    data = report / "data"
    data.mkdir(parents=True)
    (data / "x-26062_known_x_row.csv").write_text(
        "Level,External Information (Class or modification)\n"
        '1,"KNOWN sulfate of small organic alcohol"\n',
        encoding="utf-8",
    )
    result = TN.resolve(
        "",
        "C2H5O4S",
        {},
        {},
        name_resolved=False,
        bundle_dir=report,
        control_id="X-26062",
        html_alt_name="lauryl sulfate",
    )
    assert result is not None
    assert result.source == "html_alt"
    assert result.label == "lauryl sulfate"


def test_specificity_prefers_specific_alt_over_generic_known_x(tmp_path: Path) -> None:
    report = tmp_path / "X-TEST"
    data = report / "data"
    data.mkdir(parents=True)
    (data / "x-test_known_x_row.csv").write_text(
        "Level,External Information (Class or modification)\n"
        '1,"KNOWN phosphatidylcholine"\n',
        encoding="utf-8",
    )
    result = TN.resolve(
        "",
        "C42H82NO8P",
        {},
        {},
        name_resolved=False,
        bundle_dir=report,
        control_id="X-TEST",
        html_alt_name="phosphatidylcholine (18:0/20:5)",
    )
    assert result is not None
    assert result.source == "html_alt"
    assert result.label == "phosphatidylcholine (18:0/20:5)"


def test_conjugate_hint_term_mapping() -> None:
    assert TN._conjugate_hint_term({"transform_name": "glucuronide"}) == "glucuronide"
    assert TN._conjugate_hint_term({"conjugate_name": "hydroxyl_sulfate"}) == "sulfate"
    assert TN._conjugate_hint_term({"transform_name": ""}) == ""


def test_generic_sterol_upgrades_to_glucuronide_hint() -> None:
    result = TN.resolve(
        "",
        "C25H36O9",
        {
            "lipid_class": "sterol",
            "transform_name": "glucuronide",
            "core_formula": "C19H28O3",
        },
        {},
        name_resolved=False,
    )
    assert result is not None
    assert result.kind == "conjugate_unresolved"
    assert result.label == "Glucuronide conjugate of C19H28O3"
    assert result.source == "rank1_transform_hint"


def test_generic_sterol_sulfate_hint_falls_back_to_motif() -> None:
    result = TN.resolve(
        "",
        "C10H12O8S",
        {"lipid_class": "sterol", "conjugate_name": "sulfate"},
        {},
        name_resolved=False,
    )
    assert result is not None
    assert result.kind == "motif_only"
    assert result.label == "Sulfated compound"
    assert result.source == "rank1_transform_hint"


def test_detailed_lipid_shorthand_not_overridden_by_transform_hint() -> None:
    result = TN.resolve(
        "",
        "C42H82NO8P",
        {
            "lipid_class": "PC",
            "lipid_shorthand": "PC (34:1)",
            "transform_name": "glucuronide",
            "core_formula": "C19H28O3",
        },
        {},
        name_resolved=False,
    )
    assert result is not None
    assert result.kind == "lipid"
    assert result.label == "Phosphatidylcholine (34:1)"


def test_bare_phosphatidylcholine_not_relabeled_by_phosphate_hint() -> None:
    result = TN.resolve(
        "",
        "C42H82NO8P",
        {
            "lipid_class": "PC",
            "transform_name": "phosphate",
            "core_formula": "C19H28O3",
        },
        {},
        name_resolved=False,
    )
    assert result is not None
    assert result.kind == "lipid"
    assert result.label == "Phosphatidylcholine"


def _write_known_x_row(
    data_dir: Path,
    *,
    level: str = "3",
    external: str = "",
    tentative: str = "",
    cleaned: str = "",
) -> None:
    data_dir.mkdir(parents=True, exist_ok=True)
    slug = data_dir.parent.name.lower()
    (data_dir / f"{slug}_known_x_row.csv").write_text(
        "Level,External Information (Class or modification),Tentative information,Cleaned External Info\n"
        f'{level},"{external}","{tentative}","{cleaned}"\n',
        encoding="utf-8",
    )


_GLUCURONIDE_SMILES = (
    "C[C@H](CCC(=O)O)[C@H]1CC[C@H]2[C@@H]3[C@@H](O)C[C@@H]4C[C@@H]"
    "(OC5OC(C(=O)O)C(O)C(O)C5O)CC[C@]4(C)[C@@H]3CC[C@]12C"
)
_SULFATE_SMILES = "CC(=O)OS(=O)(=O)O[C@@H]1C=C2CCC3C(CC[C@@]4(C)C3CC[C@@H]4O)[C@@]2(C)CC1"


def test_known_x_conjugate_adopts_specific_glucuronide_name(tmp_path: Path) -> None:
    report = tmp_path / "X-23157"
    _write_known_x_row(
        report / "data",
        tentative="isomer of ursodeoxycholic acid glucuronide",
        cleaned="ursodeoxycholic acid glucuronide",
    )
    result = TN.resolve(
        _GLUCURONIDE_SMILES,
        "C30H48O9",
        {"lipid_class": "sterol"},
        {},
        name_resolved=False,
        bundle_dir=report,
    )
    assert result is not None
    assert result.kind == "known_x_name"
    assert result.label == "ursodeoxycholic acid glucuronide"
    assert result.source == "known_x_row"


def test_known_x_conjugate_class_descriptor_for_steroid_glucuronide(tmp_path: Path) -> None:
    report = tmp_path / "X-11444"
    _write_known_x_row(
        report / "data",
        external="Steroid glucuronide",
        tentative="Possible steroid glucuronide",
        cleaned="Steroid - Glucuronidated",
    )
    result = TN.resolve(
        _GLUCURONIDE_SMILES,
        "C25H36O9",
        {"lipid_class": "sterol"},
        {},
        name_resolved=False,
        bundle_dir=report,
    )
    assert result is not None
    assert result.kind == "chemical_class"
    assert result.label == "Steroid glucuronide"
    assert result.source == "known_x_row"


def test_known_x_conjugate_adopts_specific_sulfate_name(tmp_path: Path) -> None:
    report = tmp_path / "X-21364"
    _write_known_x_row(
        report / "data",
        level="2",
        external="sulfated steroid",
        tentative="sulfate of steroid hormone C21H32O4, 17,21-dihydroxypregnenolone",
        cleaned="Sulfate of 17,21-dihydroxypregnenolone",
    )
    result = TN.resolve(
        _SULFATE_SMILES,
        "C21H32O6S",
        {"lipid_class": "sterol"},
        {},
        name_resolved=False,
        bundle_dir=report,
    )
    assert result is not None
    assert result.kind == "known_x_name"
    assert result.label == "Sulfate of 17,21-dihydroxypregnenolone"
    assert result.confidence == "medium"


def test_known_x_conjugate_disagreement_guard_blocks_upgrade(tmp_path: Path, monkeypatch) -> None:
    report = tmp_path / "X-DISAGREE"
    _write_known_x_row(
        report / "data",
        external="Steroid glucuronide",
        cleaned="Steroid - Glucuronidated",
    )
    monkeypatch.setattr(TN, "_motif_labels", lambda _smiles: ("sulfate",))
    result = TN.resolve(
        _SULFATE_SMILES,
        "C25H36O9",
        {"lipid_class": "sterol"},
        {},
        name_resolved=False,
        bundle_dir=report,
    )
    assert result is not None
    assert result.kind == "lipid"
    assert result.label == "Sterol-like compound"


def test_known_x_conjugate_non_trigger_stays_generic(tmp_path: Path) -> None:
    report = tmp_path / "X-16935"
    _write_known_x_row(
        report / "data",
        external="fatty acid",
        tentative="C26H40O6, likely a fatty acid",
    )
    result = TN.resolve(
        "",
        "C26H40O6",
        {"lipid_class": "sterol"},
        {},
        name_resolved=False,
        bundle_dir=report,
    )
    assert result is not None
    assert result.kind == "lipid"
    assert result.label == "Sterol-like compound"


def test_known_x_conjugate_excludes_phosphate_trigger(tmp_path: Path) -> None:
    report = tmp_path / "X-23126"
    _write_known_x_row(
        report / "data",
        tentative="C24H51O4P, possibly tetracosyl dihydrogen phosphate",
    )
    result = TN.resolve(
        "",
        "C24H51O4P",
        {"lipid_class": "sterol"},
        {},
        name_resolved=False,
        bundle_dir=report,
    )
    assert result is not None
    assert result.kind == "lipid"
    assert result.label == "Sterol-like compound"


def test_known_x_name_specificity_rank_ordering() -> None:
    assert TN.specificity_rank("known_x_name", "ursodeoxycholic acid glucuronide — tentative") == 55
    assert TN.specificity_rank(
        "conjugate_resolved",
        "Testosterone glucuronide",
    ) == 60
    assert TN.specificity_rank(
        "conjugate_unresolved",
        "Glucuronide conjugate of C19H28O3",
    ) == 50
    assert TN.specificity_rank("chemical_class", "Steroid glucuronide — tentative") == 25


def test_fresh_build_known_x_conjugate_on_disk_without_patcher(tmp_path: Path) -> None:
    """Simulates compound.build_report() writing known_x_row.csv before resolve()."""
    report = tmp_path / "X-23157"
    _write_known_x_row(
        report / "data",
        tentative="isomer of ursodeoxycholic acid glucuronide",
        cleaned="ursodeoxycholic acid glucuronide",
    )
    result = TN.resolve(
        _GLUCURONIDE_SMILES,
        "C30H48O9",
        {"lipid_class": "sterol"},
        {},
        name_resolved=False,
        bundle_dir=report,
        control_id="X-23157",
    )
    assert result is not None
    assert result.kind == "known_x_name"
    assert result.label == "ursodeoxycholic acid glucuronide"


def test_known_x_conjugate_resolves_without_lipid_class(tmp_path: Path) -> None:
    report = tmp_path / "X-23157"
    _write_known_x_row(
        report / "data",
        tentative="isomer of ursodeoxycholic acid glucuronide",
        cleaned="ursodeoxycholic acid glucuronide",
    )
    result = TN.resolve(
        _GLUCURONIDE_SMILES,
        "C30H48O9",
        {},
        {},
        name_resolved=False,
        bundle_dir=report,
        control_id="X-23157",
    )
    assert result is not None
    assert result.kind == "known_x_name"
    assert result.label == "ursodeoxycholic acid glucuronide"


def test_patcher_upgrades_legacy_generic_sterol(tmp_path: Path) -> None:
    report = tmp_path / "X-23157"
    data = report / "data"
    data.mkdir(parents=True)
    formula = "C30H48O9"
    fields = {
        "overall": 0.5,
        "identity_confidence": None,
        "identity_status": "unknown",
        "msi_level": 4,
        "tier": "Low",
        "display_name": "Sterol-like compound",
        "descriptor_kind": "lipid",
        "name_resolved": False,
    }
    (data / "x-23157_confidence_inputs.json").write_text(json.dumps(fields, indent=2))
    evidence = {
        "target": {
            "control_id": "X-23157",
            "name": "Sterol-like compound",
            "formula": formula,
            "smiles": _GLUCURONIDE_SMILES,
        }
    }
    (data / "x-23157_evidence_summary.json").write_text(json.dumps(evidence, indent=2))
    (report / "manifest.json").write_text(json.dumps({"target": evidence["target"]}, indent=2))
    (data / "x-23157_top_candidates.csv").write_text(
        f"rank,predicted_smiles,predicted_formula\n1,{_GLUCURONIDE_SMILES},{formula}\n",
        encoding="utf-8",
    )
    _write_known_x_row(
        data,
        tentative="isomer of ursodeoxycholic acid glucuronide",
        cleaned="ursodeoxycholic acid glucuronide",
    )
    (report / "public.html").write_text(
        "<html><h1>Sterol-like compound</h1><p><strong>Sterol-like compound</strong></p></html>",
        encoding="utf-8",
    )
    outcome = patcher.patch_report(report, {}, apply=True)
    assert outcome["status"] == "descriptor"
    assert outcome["changed"] is True
    updated = json.loads((data / "x-23157_confidence_inputs.json").read_text())
    assert updated["display_name"] == "ursodeoxycholic acid glucuronide"
    assert updated["descriptor_kind"] == "known_x_name"
    second = patcher.patch_report(report, {}, apply=True)
    assert second["changed"] is False
    assert second["status"] in {"no_change", "preserved"}


def test_normalize_descriptor_label_strips_legacy_suffix() -> None:
    assert (
        TN.normalize_descriptor_label("Putative analog of salicin (Tanimoto 0.67) — tentative")
        == "Putative analog of salicin (Tanimoto 0.67)"
    )
    assert (
        TN.normalize_descriptor_label("foo — IUPAC (generated, tentative)")
        == "foo — IUPAC (generated)"
    )
    assert TN.normalize_descriptor_label("unchanged name") == "unchanged name"


def test_resolved_labels_never_end_with_tentative_suffix(tmp_path: Path) -> None:
    report = tmp_path / "X-23157"
    _write_known_x_row(
        report / "data",
        tentative="isomer of ursodeoxycholic acid glucuronide",
        cleaned="ursodeoxycholic acid glucuronide",
    )
    cases = [
        TN.resolve(
            "",
            "C19H32O8S2",
            {"candidate_name": "5-androstenetriol 3,17-disulfate"},
            {},
            name_resolved=False,
            control_id="X-21470",
        ),
        TN.resolve(
            _GLUCURONIDE_SMILES,
            "C30H48O9",
            {"lipid_class": "sterol"},
            {},
            name_resolved=False,
            bundle_dir=report,
        ),
    ]
    for result in cases:
        assert result is not None
        assert not result.label.endswith("— tentative")
        assert "(generated, tentative)" not in result.label


def test_patcher_strips_legacy_tentative_suffix_idempotent(tmp_path: Path) -> None:
    report = _make_bundle(tmp_path)
    data = report / "data"
    fields = json.loads((data / "x-test_confidence_inputs.json").read_text())
    suffixed = "Tripeptide (Ala–Gly–Ser) — tentative"
    fields["display_name"] = suffixed
    (data / "x-test_confidence_inputs.json").write_text(json.dumps(fields, indent=2))
    evidence = json.loads((data / "x-test_evidence_summary.json").read_text())
    evidence["target"]["name"] = suffixed
    (data / "x-test_evidence_summary.json").write_text(json.dumps(evidence, indent=2))
    for path in report.glob("*.html"):
        if path.name != "guides.html":
            path.write_text(path.read_text().replace("Tripeptide (Ala–Gly–Ser)", suffixed))

    first = patcher.patch_report(report, {}, apply=True)
    assert first["changed"] is True
    assert first["summary"]["new_title"] == "Tripeptide (Ala–Gly–Ser)"
    page = (report / "public.html").read_text()
    assert "— tentative" not in page
    assert "<h1>Tripeptide (Ala–Gly–Ser)</h1>" in page

    second = patcher.patch_report(report, {}, apply=True)
    assert second["changed"] is False


# --- MSI-consistent naming: what each level lets a report call the compound ---


def test_levels_one_and_two_qualify_a_library_name_as_an_isomer():
    for level in (1, 2):
        assert (
            TN.display_name_for_level(
                "Glycodeoxycholate", msi_level=level, name_is_reference=True
            )
            == "Isomer of Glycodeoxycholate"
        )


def test_level_three_qualifies_a_library_name_without_a_second_hedge():
    hedged = TN.display_name_for_level(
        "Glycodeoxycholate", msi_level=3, name_is_reference=True
    )
    assert hedged == "Isomer of Glycodeoxycholate"


def test_level_three_leaves_a_synthesized_descriptor_alone():
    # "Steroid glucuronide" already states its own uncertainty; hedging it twice
    # would be noise.
    label = TN.display_name_for_level(
        "Steroid glucuronide", msi_level=3, name_is_reference=False
    )
    assert label == "Steroid glucuronide"


def test_level_four_reports_the_formula_instead_of_a_structure():
    label = TN.display_name_for_level(
        "Glycodeoxycholate", msi_level=4, formula="C26H43NO5", name_is_reference=True
    )
    assert "Glycodeoxycholate" not in label
    assert label.startswith("C26H43NO5")
    assert TN.FORMULA_ONLY_NOTE in label
    # With no formula resolved there is still no structure to name.
    assert (
        TN.display_name_for_level("Glycodeoxycholate", msi_level=4, name_is_reference=True)
        == TN.FORMULA_ONLY_LABEL
    )


def test_level_five_reports_only_the_mass():
    label = TN.display_name_for_level(
        "Glycodeoxycholate",
        msi_level=5,
        formula="C26H43NO5",
        observed_mz=464.3007,
        name_is_reference=True,
    )
    assert "Glycodeoxycholate" not in label
    assert "C26H43NO5" not in label
    assert "464.3007" in label
    assert (
        TN.display_name_for_level("Anything", msi_level=5, name_is_reference=True)
        == TN.MASS_ONLY_LABEL
    )


def test_an_unknown_level_leaves_the_name_untouched():
    assert (
        TN.display_name_for_level("Taurine", msi_level="", name_is_reference=True)
        == "Isomer of Taurine"
    )
    assert (
        TN.display_name_for_level("Taurine", msi_level=9, name_is_reference=True)
        == "Isomer of Taurine"
    )


def test_hedging_is_idempotent():
    once = TN.hedged_label("Taurine")
    assert TN.hedged_label(once) == once


def test_sheet_renames_library_matches_as_isomer_of():
    assert (
        TN.sheet_feature_name(
            "valylhistidine",
            name_resolved=True,
            source="exact",
            msi_level=3,
        )
        == "Isomer of valylhistidine"
    )


def test_sheet_leaves_a_synthesized_descriptor_unqualified():
    """"Isomer of Tripeptide (Ala-Gly-Ser)" claims kinship with nothing."""

    for kind, label in (
        ("peptide", "Tripeptide (Ala-Gly-Ser)"),
        ("lipid", "Phosphatidylcholine (34:1)"),
        ("motif_only", "Sulfonated compound"),
    ):
        assert (
            TN.sheet_feature_name(
                label, name_resolved=True, source="exact", descriptor_kind=kind
            )
            == label
        )


def test_sheet_leaves_a_known_x_name_as_the_report_wrote_it():
    """A compound's own known_X identity is not an isomer of itself."""

    assert (
        TN.sheet_feature_name(
            "N-(4-aminobenzoyl)-L-glutamic acid",
            source="exact",
            descriptor_kind="candidate_name",
            msi_level=1,
            from_known_x=True,
        )
        == "N-(4-aminobenzoyl)-L-glutamic acid"
    )


def test_sheet_does_not_qualify_a_withheld_structure_label():
    """Levels 4-5 name no structure, so there is nothing to be an isomer of."""

    for label in (
        f"C4H7N3O; {TN.FORMULA_ONLY_NOTE}",
        TN.FORMULA_ONLY_LABEL,
        TN.MASS_ONLY_LABEL,
        "Exact mass 464.3007; unidentified",
    ):
        assert (
            TN.sheet_feature_name(
                label, name_resolved=True, source="exact", msi_level=4
            )
            == label
        )


def test_sheet_qualifies_confirmed_library_names_too():
    assert (
        TN.sheet_feature_name(
            "valylhistidine",
            name_resolved=True,
            source="exact",
            msi_level=1,
        )
        == "Isomer of valylhistidine"
    )
    assert (
        TN.sheet_feature_name(
            "valylhistidine",
            name_resolved=True,
            source="exact",
            msi_level=3,
            identity_status="known_pending",
        )
        == "Isomer of valylhistidine"
    )


def test_sheet_replaces_de_novo_with_tentative_structure():
    assert TN.sheet_feature_name("De novo structure", name_resolved=False, source="predicted") == (
        "Tentative Structure"
    )
    assert TN.sheet_feature_name("De novo structure (C6H12O6)", source="predicted") == (
        "Tentative Structure"
    )
    assert TN.sheet_feature_name("Tentative Structure (C6H12O6)", source="predicted") == (
        "Tentative Structure"
    )


def test_sheet_leaves_synthesized_descriptors_alone():
    assert (
        TN.sheet_feature_name(
            "Glucuronide conjugate of C17H21NO3",
            name_resolved=False,
            source="predicted",
            msi_level=4,
        )
        == "Glucuronide conjugate of C17H21NO3"
    )


def test_sheet_isomer_prefix_is_idempotent():
    once = TN.sheet_feature_name("isoleucine", name_resolved=True, source="exact", msi_level=3)
    assert once == "Isomer of isoleucine"
    assert TN.sheet_feature_name(once, name_resolved=True, source="exact", msi_level=3) == once


def test_a_borrowed_library_label_is_qualified_as_an_isomer():
    for level, truth_match in (("3", False), ("3", True), ("2", False), ("", False)):
        assert (
            TN._format_candidate_label(
                "Taurine", msi_level=level, truth_match=truth_match
            )
            == "Isomer of Taurine"
        )


def test_a_known_x_label_is_asserted_rather_than_qualified():
    assert (
        TN._format_candidate_label(
            "Taurine", msi_level="3", truth_match=False, from_known_x=True
        )
        == "Taurine"
    )


def test_msi_int_accepts_the_shapes_the_sidecars_use():
    assert TN.msi_int(3) == 3
    assert TN.msi_int("3") == 3
    assert TN.msi_int("3.0") == 3
    assert TN.msi_int("") is None
    assert TN.msi_int(None) is None
    assert TN.msi_int(0) is None
    assert TN.msi_int(6) is None
