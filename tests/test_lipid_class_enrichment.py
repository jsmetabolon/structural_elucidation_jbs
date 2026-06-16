from __future__ import annotations

import pytest

from build_lipid_class_enrichment_v1 import HAS_RDKIT, annotate_from_smiles, parse_metabolon_lipid_name


@pytest.mark.parametrize(
    ("name", "lipid_class", "chains_total", "shorthand"),
    [
        ("LPC(22:3)", "LPC", "22:3", "LPC 22:3"),
        ("DAG(38:7)", "DAG", "38:7", "DAG 38:7"),
        ("PC(16:0/18:1)", "PC", "34:1", "PC 34:1"),
        ("PC(16:1(9Z)/18:0)", "PC", "34:1", "PC 34:1"),
        ("1-palmitoyl-2-oleoyl-GPC", "PC", "34:1", "PC 34:1"),
        ("1-palmitoylglycerophosphocholine", "LPC", "16:0", "LPC 16:0"),
        ("palmitoyl sphingomyelin", "SM", "34:1", "SM d18:1/16:0"),
        ("N-palmitoyl-D-erythrosphingosine", "Cer", "34:1", "Cer d18:1/16:0"),
        ("cholesterol sulfate", "sterol", "", "sterol"),
    ],
)
def test_name_lipid_classification(name: str, lipid_class: str, chains_total: str, shorthand: str) -> None:
    annotation = parse_metabolon_lipid_name(name)
    assert annotation.lipid_class == lipid_class
    assert annotation.lipid_chains_total == chains_total
    assert annotation.lipid_shorthand == shorthand
    assert annotation.lipid_classification_source == "name_parse"


@pytest.mark.skipif(not HAS_RDKIT, reason="RDKit is not installed")
@pytest.mark.parametrize(
    ("smiles", "lipid_class", "chains_total"),
    [
        (
            "CCCCCCCCCCCCCCCC(=O)OCC(COP(=O)(O)OCC[N+](C)(C)C)OC(=O)CCCCCCC/C=C\\CCCCCCCC",
            "PC",
            "34:1",
        ),
        (
            "CCCCCCCCCCCCCCCC(=O)OCC(O)COP(=O)(O)OCC[N+](C)(C)C",
            "LPC",
            "16:0",
        ),
        (
            "CCCCCCCCCCCCCCCC(=O)OCC(COP(=O)(O)OCCN)OC(=O)CCCCCCCCCCCCCCCCC",
            "PE",
            "34:0",
        ),
        (
            "CCCCCCCCCCCCCCCC(=O)NC(CO)C(O)CCCCCCCCCCCCCCOP(=O)(O)OCC[N+](C)(C)C",
            "SM",
            "34:1",
        ),
        (
            "CCCCCCCCCCCCCCCC(=O)NC(CO)C(O)CCCCCCCCCCCCC",
            "Cer",
            "34:1",
        ),
        (
            "CCCCCCCCCCCCCCCC(=O)OCC(O)COC(=O)CCCCCCCCCCCCCCCCC",
            "DAG",
            "34:0",
        ),
        (
            "C[C@H](CCCC(C)C)C1CCC2C3CCC4=CC(O)CCC4(C)C3CCC12C",
            "sterol",
            "",
        ),
    ],
)
def test_smiles_lipid_classification(smiles: str, lipid_class: str, chains_total: str) -> None:
    annotation = annotate_from_smiles(smiles)
    assert annotation.lipid_class == lipid_class
    assert annotation.lipid_chains_total == chains_total
    assert annotation.lipid_classification_source == "smiles_smarts"
