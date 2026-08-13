#!/usr/bin/env python3
"""Export the evidence-score parameters, their bands, and their report language.

Every number here is imported from the code that applies it. That is the whole
point of the file: a parameter table maintained by hand is wrong the first time
somebody edits a weight, and nothing would notice. Importing the constants makes
the table regenerate into a diff instead, and ``tests/test_parameter_table.py``
turns that diff into a failure.

The table covers the five weighted axes in ``confidence.WEIGHTS`` plus the
library-neighbor count, which is a report card rather than an axis and is marked
``scored = false`` so it cannot be read as one.

Usage:
    python3 tools/build_parameter_table.py
    python3 tools/build_parameter_table.py --xlsx evidence_parameters.xlsx
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Any

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from structural_elucidation.inference.mass_windows import (  # noqa: E402
    MASS_SCORE_PLATEAU_DA,
    MASS_SCORE_TOLERANCE_DA,
)
from structural_elucidation.reports.confidence import (  # noqa: E402
    AXIS_LABELS,
    MASS_STRONG,
    MASS_WEAK,
    MS2_PRESENT,
    MS2_PROBABLE_EXACT,
    NEIGHBOR_SIMILARITY_FLOOR,
    RETENTION_IMPLAUSIBLE_Z,
    STRUCT_MODERATE,
    STRUCT_STRONG,
    STRUCT_WEAK,
    WEIGHTS,
)

FIELDNAMES = (
    "parameter",
    "axis_key",
    "scored",
    "weight",
    "parameter_description",
    "value_strength",
    "low_end",
    "high_end",
    "units",
    "report_language",
    "consequence",
)

# The stakeholder sheet's vocabulary, kept so the export drops into the format it
# was asked for. The code says Strong / Partial / Fail for mass and Strong /
# Moderate / Weak on the card; the cutoffs below are the code's regardless of
# which adjective sits on top of them.
STRONG, ACCEPTABLE, BORDERLINE = "Strong", "Acceptable", "Borderline"
UNSCORED = "Not scored"

# Mass thresholds are absolute Da, so their ppm equivalent moves with the mass of
# the compound. Any single ppm figure is therefore only true at one mass, and
# this is the mass the exported ppm equivalents are quoted at.
PPM_REFERENCE_MASS_DA = 300.0


def _ppm(error_da: float) -> float:
    return error_da / PPM_REFERENCE_MASS_DA * 1e6


def _mass_error_for_score(score: float) -> float:
    """Invert the linear arm of ``mass_match_score``: score = 1 - err / tolerance."""

    return (1.0 - score) * MASS_SCORE_TOLERANCE_DA


def _row(
    parameter: str,
    axis_key: str,
    scored: bool,
    weight: float,
    description: str,
    strength: str,
    low: str,
    high: str,
    units: str,
    language: str,
    consequence: str,
) -> dict[str, Any]:
    return {
        "parameter": parameter,
        "axis_key": axis_key,
        "scored": scored,
        "weight": round(weight, 2),
        "parameter_description": description,
        "value_strength": strength,
        "low_end": low,
        "high_end": high,
        "units": units,
        "report_language": language,
        "consequence": consequence,
    }


def structure_rows() -> list[dict[str, Any]]:
    key = "structure"
    description = (
        "Tanimoto similarity between the top candidate and the nearest cataloged reference "
        "metabolite. A candidate taken straight from the catalog pins this axis to 1.00, which "
        "is provenance rather than proof: exact-source candidates are correct less often than "
        "de novo ones on the truth set."
    )
    weight = WEIGHTS[key]
    return [
        _row(
            AXIS_LABELS[key], key, True, weight, description, STRONG,
            f"{STRUCT_STRONG:.2f}", "1.00", "Tanimoto similarity",
            "Top structure is very close to a cataloged metabolite",
            "",
        ),
        _row(
            AXIS_LABELS[key], key, True, weight, description, ACCEPTABLE,
            f"{STRUCT_WEAK:.2f}", f"{STRUCT_STRONG:.2f}", "Tanimoto similarity",
            "Top structure is moderately similar to the nearest cataloged metabolite",
            f"At or above {STRUCT_MODERATE:.2f} the candidate can reach MSI Level 2 given MS2 and "
            f"mass support; at or above {STRUCT_WEAK:.2f} it can reach Level 3.",
        ),
        _row(
            AXIS_LABELS[key], key, True, weight, description, BORDERLINE,
            "0.00", f"{STRUCT_WEAK:.2f}", "Tanimoto similarity",
            "Top structure is distant from any cataloged metabolite (de novo)",
            "Below this the structure cannot support MSI Level 3, so the report withholds a "
            "structure and gives formula or mass only.",
        ),
    ]


def mass_rows() -> list[dict[str, Any]]:
    key = "mass"
    partial_err = _mass_error_for_score(MASS_WEAK)
    description = (
        "Agreement between the candidate's neutral mass and the mass implied by the observed "
        "feature. Called both Mass Match and Mass Formula Match; they are one axis. The score is "
        f"1.00 at or below {MASS_SCORE_PLATEAU_DA:g} Da of error and 0.00 at or above "
        f"{MASS_SCORE_TOLERANCE_DA:g} Da, falling linearly between. Because the windows are "
        f"absolute Da, the ppm equivalents below hold at {PPM_REFERENCE_MASS_DA:.0f} Da and "
        "tighten as the compound gets heavier. "
        "These windows are NOT what scored the delivered reports. That column is filled by "
        "chemistry.conjugates.linear_mass_match_score, a different function running at a 1.0 Da "
        "default, and roughly 80% of delivered scores came from it against 2.6% from the windows "
        "described here. Under a 1.0 Da tolerance the strong-mass cutoff admits an error near "
        "0.1 Da, some hundreds of ppm rather than the few implied above. See "
        "docs/paper/data/thresholds/audit_mass_policy.json."
    )
    weight = WEIGHTS[key]
    return [
        _row(
            AXIS_LABELS[key], key, True, weight, description, STRONG,
            "0.000 Da (0.0 ppm)",
            f"{MASS_SCORE_PLATEAU_DA:.3f} Da ({_ppm(MASS_SCORE_PLATEAU_DA):.1f} ppm)",
            "absolute neutral-mass error",
            "The mass of the proposed structure matches the observed mass, which is a strong "
            "match for this parameter",
            f"Sub-score 1.00. The plateau makes this the only way to clear the {MASS_STRONG:.2f} "
            "strong-mass bar: leaving it drops the score straight to 0.80, so nothing scores "
            "between 0.80 and 1.00.",
        ),
        _row(
            AXIS_LABELS[key], key, True, weight, description, ACCEPTABLE,
            f"{MASS_SCORE_PLATEAU_DA:.3f} Da ({_ppm(MASS_SCORE_PLATEAU_DA):.1f} ppm)",
            f"{partial_err:.4f} Da ({_ppm(partial_err):.1f} ppm)",
            "absolute neutral-mass error",
            "Partial mass agreement",
            f"Sub-score {MASS_WEAK:.2f} to 0.80. Too low for the strong-mass bar, so it cannot "
            "reach MSI Level 2.",
        ),
        _row(
            AXIS_LABELS[key], key, True, weight, description, BORDERLINE,
            f"{partial_err:.4f} Da ({_ppm(partial_err):.1f} ppm)",
            "no limit",
            "absolute neutral-mass error",
            "Candidate mass does NOT match the target neutral mass",
            f"Sub-score below {MASS_WEAK:.2f}, which demotes the report to Tentative whatever the "
            "other axes say, and caps it at MSI Level 4 or below.",
        ),
    ]


def ms2_rows() -> list[dict[str, Any]]:
    key = "ms2"
    description = (
        "How well the candidate explains the observed MS2 spectrum. Scored only where a spectrum "
        "was acquired for the feature; where none was, the axis is reported unavailable rather "
        "than zero, and its weight is redistributed over the axes that are available."
    )
    weight = WEIGHTS[key]
    return [
        _row(
            AXIS_LABELS[key], key, True, weight, description, STRONG,
            f"{MS2_PROBABLE_EXACT:.2f}", "1.00", "spectral fit",
            "Observed MS2 spectrum explained by candidate (fit {fit})",
            "Meets the higher bar an exact catalog hit must clear for MSI Level 2, which exists "
            "because a catalog match supplies no spectral evidence of its own.",
        ),
        _row(
            AXIS_LABELS[key], key, True, weight, description, ACCEPTABLE,
            f"{MS2_PRESENT:.2f}", f"{MS2_PROBABLE_EXACT:.2f}", "spectral fit",
            "Observed MS2 spectrum explained by candidate (fit {fit})",
            "Enough for MSI Level 2 for a de novo candidate, not for one taken from the catalog.",
        ),
        _row(
            AXIS_LABELS[key], key, True, weight, description, BORDERLINE,
            "0.00", f"{MS2_PRESENT:.2f}", "spectral fit",
            "Observed MS2 spectrum poorly explained by candidate (fit {fit})",
            "A spectrum exists and the candidate does not account for it, which is evidence "
            "against the candidate rather than missing evidence.",
        ),
    ]


def pathway_rows() -> list[dict[str, Any]]:
    key = "pathway"
    # The tier fallback is the only place the pathway axis has anything like a
    # band, so it supplies the edges. The live input is a continuous share.
    high, medium, low = 0.90, 0.60, 0.30
    description = (
        "Share of chemically similar library standards that agree on the pathway assignment. The "
        "live input is a continuous vote share; where a report predates the share being recorded, "
        f"the persisted tier stands in as high {high:.2f}, medium {medium:.2f}, low {low:.2f}. The "
        "bands below are those fallback values."
    )
    weight = WEIGHTS[key]
    return [
        _row(
            AXIS_LABELS[key], key, True, weight, description, STRONG,
            f"{high:.2f}", "1.00", "neighbor vote share",
            "Neighbor pathway consensus {share}",
            "",
        ),
        _row(
            AXIS_LABELS[key], key, True, weight, description, ACCEPTABLE,
            f"{medium:.2f}", f"{high:.2f}", "neighbor vote share",
            "Neighbor pathway consensus {share}",
            "",
        ),
        _row(
            AXIS_LABELS[key], key, True, weight, description, BORDERLINE,
            "0.00", f"{medium:.2f}", "neighbor vote share",
            "Neighbor pathway consensus {share}",
            "The pathway is still reported, with the split vote behind it shown, because every "
            "feature receives pathway context whether or not a structure is shown.",
        ),
    ]


def retention_rows() -> list[dict[str, Any]]:
    key = "retention"
    description = (
        "Whether the candidate could plausibly have eluted where the feature did, as standard "
        "deviations between the observed retention index and the one predicted for the candidate. "
        "Available only on libraries whose retention model cleared the feasibility audit, and it "
        "can only lower the evidence score, never raise it."
    )
    weight = WEIGHTS[key]
    return [
        _row(
            AXIS_LABELS[key], key, True, weight, description, STRONG,
            "0.0 sigma", "1.0 sigma", "residual sigma",
            "Retention is consistent with the candidate ({sigma} sigma)",
            "",
        ),
        _row(
            AXIS_LABELS[key], key, True, weight, description, ACCEPTABLE,
            "1.0 sigma", f"{RETENTION_IMPLAUSIBLE_Z:.1f} sigma", "residual sigma",
            "Retention is a loose fit ({sigma} sigma from the predicted index)",
            "",
        ),
        _row(
            AXIS_LABELS[key], key, True, weight, description, BORDERLINE,
            f"{RETENTION_IMPLAUSIBLE_Z:.1f} sigma", "no limit", "residual sigma",
            "Candidate could not have eluted here ({sigma} sigma from its predicted index)",
            "Demotes the report to Tentative on its own and blocks MSI Levels 2 and 3, on the "
            "reasoning that a correct structure does not elute this far from its predicted index.",
        ),
    ]


def library_neighbor_rows() -> list[dict[str, Any]]:
    description = (
        "How many Metabolon library standards are chemically similar to the feature, counted at "
        f"Tanimoto {NEIGHBOR_SIMILARITY_FLOOR:.2f} or better. This is context shown on the report "
        "card, not an evidence axis: it carries no weight and moving it changes no score. A "
        "drafted threshold of five neighbors is not implemented and no count acts as a cutoff."
    )
    return [
        _row(
            "Library neighbors", "", False, 0.0, description, UNSCORED,
            "0", "no limit",
            f"count of library standards at Tanimoto >= {NEIGHBOR_SIMILARITY_FLOOR:.2f}",
            "This metabolite is similar to {count} Metabolon library standards, enabling a higher "
            "confidence structural candidate",
            "Informational. The similarity of the single nearest standard is what the structure "
            "axis scores; the count is not scored anywhere.",
        ),
    ]


def build_rows() -> list[dict[str, Any]]:
    return [
        *structure_rows(),
        *mass_rows(),
        *ms2_rows(),
        *pathway_rows(),
        *retention_rows(),
        *library_neighbor_rows(),
    ]


def build_summary() -> dict[str, Any]:
    rows = build_rows()
    return {
        "method": (
            "Parameters of the identification-confidence evidence score, exported from the "
            "constants that apply them in structural_elucidation/reports/confidence.py and "
            "structural_elucidation/inference/mass_windows.py."
        ),
        "weights": dict(WEIGHTS),
        "weight_total": round(sum(WEIGHTS.values()), 2),
        "ppm_reference_mass_da": PPM_REFERENCE_MASS_DA,
        "parameters": rows,
        "scope_note": (
            "Bands describe how a measurement is read, not how often that reading is right. The "
            "aggregate score is Platt-scaled in-sample on 92 truth-verified controls and is not "
            "externally calibrated, so no band should be read as a probability of correctness."
        ),
    }


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(FIELDNAMES), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_xlsx(path: Path, rows: list[dict[str, Any]]) -> None:
    """Write the stakeholder's six-column layout for hand-off.

    Deliberately not a committed artifact: openpyxl stamps timestamps into the
    zip the same way python-docx does, so the file would never reproduce
    byte-identically and would fail the CI diff that guards everything else here.
    """

    from openpyxl import Workbook
    from openpyxl.styles import Alignment, Font

    workbook = Workbook()
    sheet = workbook.active
    sheet.title = "Evidence Score Ranges"
    headers = (
        "Parameter",
        "Parameter Description",
        "Value Strength",
        "Low End",
        "High End",
        "Report Language",
    )
    sheet.append(list(headers))
    for cell in sheet[1]:
        cell.font = Font(bold=True)
    for row in rows:
        sheet.append(
            [
                row["parameter"],
                row["parameter_description"],
                row["value_strength"],
                row["low_end"],
                row["high_end"],
                row["report_language"],
            ]
        )
    for column, width in zip("ABCDEF", (24, 60, 16, 30, 30, 60)):
        sheet.column_dimensions[column].width = width
    for column in ("B", "F"):
        for cell in sheet[column]:
            cell.alignment = Alignment(wrap_text=True, vertical="top")
    path.parent.mkdir(parents=True, exist_ok=True)
    workbook.save(path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=Path("docs/paper/data/parameters"))
    parser.add_argument(
        "--xlsx",
        type=Path,
        default=None,
        help="also write the stakeholder six-column layout here (not a committed artifact)",
    )
    args = parser.parse_args()

    summary = build_summary()
    rows = summary["parameters"]

    args.out_dir.mkdir(parents=True, exist_ok=True)
    json_path = args.out_dir / "evidence_parameters.json"
    json_path.write_text(json.dumps(summary, indent=2) + "\n")
    csv_path = args.out_dir / "evidence_parameters.csv"
    write_csv(csv_path, rows)

    scored = [row for row in rows if row["scored"]]
    print(
        f"[parameters] {len(rows)} bands across "
        f"{len({row['parameter'] for row in rows})} parameters; "
        f"{len({row['axis_key'] for row in scored})} scored axes totalling "
        f"{summary['weight_total']:.2f} weight"
    )
    print(f"[parameters] wrote {json_path} and {csv_path}")
    if args.xlsx:
        write_xlsx(args.xlsx, rows)
        print(f"[parameters] wrote {args.xlsx}")


if __name__ == "__main__":
    main()
