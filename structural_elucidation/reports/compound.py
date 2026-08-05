#!/usr/bin/env python3
"""Self-contained HTML exploration reports for X-compound structural elucidation.

ARCHITECTURE OVERVIEW
---------------------
The report pipeline has four layers:

1. **Data assembly** — ``build_report`` loads rerank CSVs, ``known_X`` metadata,
   metabolite SQLite, optional peak tables, and WGCNA correlation exports.
   ``predictions_to_candidate_pool`` normalizes heterogeneous prediction columns.

2. **Chemistry & scoring** — RDKit fingerprints (``mol_fp``, ``similarity``),
   Open Babel SVG rendering (``render_smiles_svg``), spectrum alignment
   (``score_smiles_spectrum_alignment``), and MCS consensus heatmaps
   (``consensus_atom_support``).

3. **Network graphs** — Morgan Tanimoto edges (``build_pairwise_similarity_edges``),
   greedy/union clustering, and force-directed layout produce SVG networks for
   library neighbors and extended known/X similarity views.

4. **HTML emission** — ``render_html`` stitches Metro-themed CSS, tabbed candidate
   panels, and public/internal variants; ``render_guides_html`` emits placeholder
   info-card pages. Outputs land under ``output/reports/<control_id>/``.

CLI: ``python -m structural_elucidation.reports.compound --control-id X-25003 ...``
"""


from __future__ import annotations

import argparse
import csv
import html
import json
import math
import re
import shutil
import sqlite3
import subprocess
from pathlib import Path
from typing import Any
from urllib.parse import quote

import pandas as pd

import structural_elucidation.utils.wgcna as wgcna
from structural_elucidation.paths import REPO_ROOT, resolve
from structural_elucidation.reports import confidence as confmod
from structural_elucidation.reports import tentative_naming as tnmod
from structural_elucidation.reports.tentative_naming_ext import default_extended_context
from structural_elucidation.utils.denovo import exact_mass_from_smiles, formula_from_smiles, inchikey_from_smiles

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Draw, rdDepictor, rdFingerprintGenerator, rdFMCS
    from rdkit.DataStructs import TanimotoSimilarity

    HAS_RDKIT = True
    MORGAN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    try:
        # CoordGen gives far better 2D layouts for large/macrocyclic molecules and
        # avoids the atom-position collapse that makes consensus images unreadable.
        rdDepictor.SetPreferCoordGen(True)
    except Exception:  # pragma: no cover - older RDKit without CoordGen toggle
        pass
except ImportError:  # pragma: no cover
    Chem = None
    Draw = None
    MORGAN = None
    HAS_RDKIT = False


DEFAULT_OUT_DIR = REPO_ROOT / "output/reports"
DEFAULT_DB = REPO_ROOT / "vendor/metabolite_synonyms.sqlite"
DEFAULT_KNOWN_X = REPO_ROOT / "data/known_X.csv"
DEFAULT_KNOWN_X_ALL = REPO_ROOT / "data/known_X_all_available_structures_simple.csv"
DEFAULT_CURATION_LEDGER = REPO_ROOT / "data/curation_ledger.csv"
OBABEL_PATH_CANDIDATES = [
    Path.home() / "Downloads" / "ADFRsuite_x86_64Darwin_1.0" / "bin" / "obabel",
]
STRUCTURE_TARGET_BOND_LENGTH = 36.0
STRUCTURE_LINE_WIDTH_RATIO = 0.021 / 0.381
STRUCTURE_BOLD_WIDTH_RATIO = 0.055 / 0.381
STRUCTURE_BOND_SPACING_RATIO = 0.18
STRUCTURE_HASH_SPACING_RATIO = 0.06 / 0.381
STRUCTURE_LABEL_FONT_SIZE = 12.0
STRUCTURE_PADDING = 18.0
STRUCTURE_FONT_FAMILY = "Arial, Helvetica, sans-serif"

# Report default paths resolve via structural_elucidation.paths.REPO_ROOT


# ── IO and SQLite helpers ─────────────────────────────────────────────────────

def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    """Write a list of row dicts to CSV; empty input yields an empty file."""

    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def query_df(db_path: Path, sql: str, params: tuple[Any, ...] = ()) -> pd.DataFrame:
    """Run a parameterized SQL query against the metabolite SQLite DB."""

    with sqlite3.connect(db_path) as conn:
        return pd.read_sql_query(sql, conn, params=params)


_MOL_FP_CACHE: dict[str, tuple[Any, Any]] = {}


def mol_fp(smiles: str) -> tuple[Any, Any] | tuple[None, None]:
    """Parse SMILES and return ``(mol, morgan_fp)`` or ``(None, None)``.

    Results are memoized by SMILES string: the extended similarity network does
    O(N^2) pairwise comparisons and would otherwise re-parse and re-fingerprint
    the same molecules tens of thousands of times per report. The fingerprint is
    deterministic per SMILES, so caching is output-identical.
    """

    if not HAS_RDKIT:
        return None, None
    key = str(smiles or "")
    cached = _MOL_FP_CACHE.get(key)
    if cached is not None:
        return cached
    mol = Chem.MolFromSmiles(key)
    if mol is None:
        result: tuple[Any, Any] = (None, None)
        _MOL_FP_CACHE[key] = result
        return result
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        result = (None, None)
        _MOL_FP_CACHE[key] = result
        return result
    result = (mol, MORGAN.GetFingerprint(mol))
    _MOL_FP_CACHE[key] = result
    return result


def similarity(smiles_a: str, smiles_b: str) -> float:
    """Tanimoto similarity between two SMILES strings."""

    _, fp_a = mol_fp(smiles_a)
    _, fp_b = mol_fp(smiles_b)
    if fp_a is None or fp_b is None:
        return 0.0
    return float(TanimotoSimilarity(fp_a, fp_b))


def compute_2d_coords_robust(mol: Any) -> None:
    """Lay out 2D coordinates, preferring CoordGen but falling back to the default
    depictor when CoordGen collapses atoms onto each other (which makes structure
    and consensus images unreadable). Picks whichever layout best separates atoms.
    """

    if not HAS_RDKIT:
        return

    def min_pair_dist() -> float:
        conf = mol.GetConformer()
        n = mol.GetNumAtoms()
        if n < 2:
            return 9.0
        pts = [(conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y) for i in range(n)]
        return min(
            ((pts[i][0] - pts[j][0]) ** 2 + (pts[i][1] - pts[j][1]) ** 2) ** 0.5
            for i in range(n)
            for j in range(i + 1, n)
        )

    try:
        rdDepictor.SetPreferCoordGen(True)
        rdDepictor.Compute2DCoords(mol)
        coordgen_dist = min_pair_dist()
        if coordgen_dist < 0.30:
            rdDepictor.SetPreferCoordGen(False)
            rdDepictor.Compute2DCoords(mol)
            if min_pair_dist() < coordgen_dist:
                rdDepictor.SetPreferCoordGen(True)
                rdDepictor.Compute2DCoords(mol)
    except Exception:
        try:
            rdDepictor.Compute2DCoords(mol)
        except Exception:
            pass
    finally:
        rdDepictor.SetPreferCoordGen(True)


def structure_image(smiles: str, path: Path, legend: str = "") -> None:
    """Render a 2D structure PNG via RDKit when available."""

    path.parent.mkdir(parents=True, exist_ok=True)
    if not HAS_RDKIT:
        return
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return
    try:
        compute_2d_coords_robust(mol)
        Draw.MolToFile(mol, str(path), size=(520, 360), legend=legend)
    except Exception:
        # Some molecules trip RDKit's drawer (e.g. "zero length vector"); a single
        # un-renderable structure image must not abort the whole report.
        return


# ── Structure rendering (Open Babel / SVG) ────────────────────────────────────

def resolve_obabel_path() -> Path | None:
    """Locate Open Babel binary for SMILES-to-mol-block conversion."""

    which_path = shutil.which("obabel")
    candidates = [Path(which_path).resolve()] if which_path else []
    candidates.extend(candidate.resolve() for candidate in OBABEL_PATH_CANDIDATES)
    seen: set[str] = set()
    for candidate in candidates:
        if not candidate:
            continue
        candidate_str = str(candidate)
        if candidate_str in seen:
            continue
        seen.add(candidate_str)
        if candidate.exists():
            return candidate
    return None


def parse_mol_block(mol_block: str) -> dict[str, Any] | None:
    """Parse an MDL V2000 mol block into atom/bond coordinate dicts."""

    lines = [line.rstrip("\n") for line in mol_block.splitlines()]
    if len(lines) < 3:
        return None
    counts_index = next((index for index, line in enumerate(lines) if "V2000" in line or "V3000" in line), -1)
    if counts_index < 0:
        return None
    counts_line = lines[counts_index]
    if "V3000" in counts_line:
        return None
    try:
        atom_count = int(counts_line[0:3])
        bond_count = int(counts_line[3:6])
    except ValueError:
        return None

    atoms = []
    bonds = []
    atom_lines = lines[counts_index + 1 : counts_index + 1 + atom_count]
    bond_lines = lines[counts_index + 1 + atom_count : counts_index + 1 + atom_count + bond_count]

    for atom_line in atom_lines:
        try:
            atoms.append(
                {
                    "x": float(atom_line[0:10]),
                    "y": float(atom_line[10:20]),
                    "z": float(atom_line[20:30]),
                    "element": atom_line[31:34].strip() or "C",
                }
            )
        except ValueError:
            return None

    for bond_line in bond_lines:
        try:
            bonds.append(
                {
                    "a1": int(bond_line[0:3]) - 1,
                    "a2": int(bond_line[3:6]) - 1,
                    "order": int(bond_line[6:9]),
                    "stereo": int(bond_line[9:12]),
                }
            )
        except ValueError:
            return None

    if not atoms:
        return None
    return {"atoms": atoms, "bonds": bonds}


def generate_smiles_mol_block(smiles: str) -> str | None:
    """Generate a 2D V2000 mol block from SMILES via Open Babel."""

    if not smiles:
        return None
    obabel_path = resolve_obabel_path()
    if obabel_path is None:
        return None
    command = [str(obabel_path), f"-:{smiles}", "-omol", "--gen2d"]
    try:
        result = subprocess.run(command, capture_output=True, text=True, timeout=20)
    except Exception:
        return None
    mol_block = result.stdout.strip()
    if result.returncode != 0 or "V2000" not in mol_block:
        return None
    return mol_block


def _structure_label_for_atom(atom_index: int, structure: dict[str, Any], stereochemical_hydrogens: set[int]) -> str:
    """Choose element label text for SVG atom rendering."""

    atom = structure["atoms"][atom_index]
    element = atom["element"]
    bonds = [bond for bond in structure["bonds"] if bond["a1"] == atom_index or bond["a2"] == atom_index]
    if element == "C":
        return ""
    if element == "H":
        return "H" if atom_index in stereochemical_hydrogens else ""
    if element == "O" and len(bonds) == 1 and bonds[0]["order"] == 1:
        return "OH"
    return element


def _structure_label_radius(label: str) -> float:
    """Estimate label collision radius for bond shortening."""

    if not label:
        return 0.0
    width_estimate = STRUCTURE_LABEL_FONT_SIZE * (0.35 + 0.28 * len(label))
    return max(STRUCTURE_LABEL_FONT_SIZE * 0.62, width_estimate * 0.5)


def _structure_label_color(element: str) -> str:
    """Element-specific fill color for heteroatom labels."""

    palette = {
        "N": "#2E6BE6",
        "O": "#D62728",
        "S": "#C49A00",
        "P": "#E67E22",
        "F": "#2CA02C",
        "CL": "#2CA02C",
        "BR": "#8C564B",
        "I": "#7B52AB",
    }
    return palette.get((element or "").upper(), "#000000")


def _shorten_segment(
    start: tuple[float, float],
    end: tuple[float, float],
    start_shorten: float,
    end_shorten: float,
) -> tuple[tuple[float, float], tuple[float, float]] | None:
    """Trim a bond segment so labels do not overlap endpoints."""

    dx = end[0] - start[0]
    dy = end[1] - start[1]
    length = math.hypot(dx, dy)
    if length <= 1e-6 or start_shorten + end_shorten >= length:
        return None
    ux = dx / length
    uy = dy / length
    return (
        (start[0] + ux * start_shorten, start[1] + uy * start_shorten),
        (end[0] - ux * end_shorten, end[1] - uy * end_shorten),
    )


def _parallel_offset(
    start: tuple[float, float],
    end: tuple[float, float],
    offset: float,
) -> tuple[tuple[float, float], tuple[float, float]]:
    """Offset a bond segment perpendicular to its axis (double bonds)."""

    dx = end[0] - start[0]
    dy = end[1] - start[1]
    length = math.hypot(dx, dy)
    if length <= 1e-6:
        return start, end
    nx = -dy / length
    ny = dx / length
    return ((start[0] + nx * offset, start[1] + ny * offset), (end[0] + nx * offset, end[1] + ny * offset))


def render_mol_block_svg(mol_block: str) -> str | None:
    """Render a V2000 mol block as inline SVG with stereochemistry wedges."""

    structure = parse_mol_block(mol_block)
    if not structure:
        return None
    atoms = structure["atoms"]
    bonds = structure["bonds"]
    if not atoms or not bonds:
        return None

    heavy_bond_lengths = []
    for bond in bonds:
        left = atoms[bond["a1"]]
        right = atoms[bond["a2"]]
        if left["element"] == "H" or right["element"] == "H":
            continue
        heavy_bond_lengths.append(math.hypot(right["x"] - left["x"], right["y"] - left["y"]))
    if not heavy_bond_lengths:
        heavy_bond_lengths = [
            math.hypot(atoms[bond["a2"]]["x"] - atoms[bond["a1"]]["x"], atoms[bond["a2"]]["y"] - atoms[bond["a1"]]["y"])
            for bond in bonds
        ]
    average_bond_length = sum(heavy_bond_lengths) / max(len(heavy_bond_lengths), 1)
    scale = STRUCTURE_TARGET_BOND_LENGTH / max(average_bond_length, 1e-6)

    transformed_points = [(atom["x"] * scale, -atom["y"] * scale) for atom in atoms]
    min_x = min(point[0] for point in transformed_points)
    max_x = max(point[0] for point in transformed_points)
    min_y = min(point[1] for point in transformed_points)
    max_y = max(point[1] for point in transformed_points)
    width = max_x - min_x + STRUCTURE_PADDING * 2
    height = max_y - min_y + STRUCTURE_PADDING * 2
    translated_points = [
        (point[0] - min_x + STRUCTURE_PADDING, point[1] - min_y + STRUCTURE_PADDING) for point in transformed_points
    ]

    stereochemical_hydrogens = {
        bond["a1"] if atoms[bond["a1"]]["element"] == "H" else bond["a2"]
        for bond in bonds
        if bond["stereo"] in {1, 6} and (atoms[bond["a1"]]["element"] == "H" or atoms[bond["a2"]]["element"] == "H")
    }
    labels = {
        index: _structure_label_for_atom(index, structure, stereochemical_hydrogens) for index in range(len(atoms))
    }

    line_width = STRUCTURE_TARGET_BOND_LENGTH * STRUCTURE_LINE_WIDTH_RATIO
    bold_width = STRUCTURE_TARGET_BOND_LENGTH * STRUCTURE_BOLD_WIDTH_RATIO
    double_offset = STRUCTURE_TARGET_BOND_LENGTH * STRUCTURE_BOND_SPACING_RATIO * 0.5
    hash_spacing = STRUCTURE_TARGET_BOND_LENGTH * STRUCTURE_HASH_SPACING_RATIO

    bond_fragments: list[str] = []
    for bond in bonds:
        start = translated_points[bond["a1"]]
        end = translated_points[bond["a2"]]
        shortened = _shorten_segment(
            start,
            end,
            _structure_label_radius(labels[bond["a1"]]),
            _structure_label_radius(labels[bond["a2"]]),
        )
        if shortened is None:
            continue
        start_point, end_point = shortened
        dx = end_point[0] - start_point[0]
        dy = end_point[1] - start_point[1]
        length = math.hypot(dx, dy)
        if length <= 1e-6:
            continue
        nx = -dy / length
        ny = dx / length

        if bond["order"] == 2:
            upper_start, upper_end = _parallel_offset(start_point, end_point, double_offset)
            lower_start, lower_end = _parallel_offset(start_point, end_point, -double_offset)
            bond_fragments.append(
                f'<line x1="{upper_start[0]:.2f}" y1="{upper_start[1]:.2f}" x2="{upper_end[0]:.2f}" y2="{upper_end[1]:.2f}" />'
            )
            bond_fragments.append(
                f'<line x1="{lower_start[0]:.2f}" y1="{lower_start[1]:.2f}" x2="{lower_end[0]:.2f}" y2="{lower_end[1]:.2f}" />'
            )
            continue

        if bond["order"] == 3:
            bond_fragments.append(
                f'<line x1="{start_point[0]:.2f}" y1="{start_point[1]:.2f}" x2="{end_point[0]:.2f}" y2="{end_point[1]:.2f}" />'
            )
            upper_start, upper_end = _parallel_offset(start_point, end_point, double_offset)
            lower_start, lower_end = _parallel_offset(start_point, end_point, -double_offset)
            bond_fragments.append(
                f'<line x1="{upper_start[0]:.2f}" y1="{upper_start[1]:.2f}" x2="{upper_end[0]:.2f}" y2="{upper_end[1]:.2f}" />'
            )
            bond_fragments.append(
                f'<line x1="{lower_start[0]:.2f}" y1="{lower_start[1]:.2f}" x2="{lower_end[0]:.2f}" y2="{lower_end[1]:.2f}" />'
            )
            continue

        if bond["order"] == 1 and bond["stereo"] in {1, 6}:
            half_narrow = line_width * 0.35
            half_wide = bold_width * 0.5
            if bond["stereo"] == 1:
                polygon = [
                    (start_point[0] + nx * half_narrow, start_point[1] + ny * half_narrow),
                    (start_point[0] - nx * half_narrow, start_point[1] - ny * half_narrow),
                    (end_point[0] - nx * half_wide, end_point[1] - ny * half_wide),
                    (end_point[0] + nx * half_wide, end_point[1] + ny * half_wide),
                ]
                points = " ".join(f"{x:.2f},{y:.2f}" for x, y in polygon)
                bond_fragments.append(f'<polygon points="{points}" />')
            else:
                steps = max(6, int(length / max(hash_spacing, 1.0)))
                for step in range(steps):
                    fraction = (step + 1) / (steps + 1)
                    center_x = start_point[0] + dx * fraction
                    center_y = start_point[1] + dy * fraction
                    half_width = half_narrow + (half_wide - half_narrow) * fraction
                    segment = max(hash_spacing * 0.45, 1.2)
                    x1 = center_x - nx * half_width
                    y1 = center_y - ny * half_width
                    x2 = center_x + nx * half_width
                    y2 = center_y + ny * half_width
                    bond_fragments.append(
                        f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" stroke-width="{segment:.2f}" stroke-linecap="round" />'
                    )
            continue

        bond_fragments.append(
            f'<line x1="{start_point[0]:.2f}" y1="{start_point[1]:.2f}" x2="{end_point[0]:.2f}" y2="{end_point[1]:.2f}" />'
        )

    text_fragments: list[str] = []
    for index, label in labels.items():
        if not label:
            continue
        x_pos, y_pos = translated_points[index]
        fill = _structure_label_color(str(atoms[index]["element"]))
        text_fragments.append(f'<text x="{x_pos:.2f}" y="{y_pos:.2f}" fill="{fill}">{html.escape(label)}</text>')

    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width:.0f}px" height="{height:.0f}px" '
        f'viewBox="0 0 {width:.2f} {height:.2f}" role="img" aria-label="2D chemical structure">'
        f"<style>"
        f".bond-layer line{{stroke:#000;stroke-width:{line_width:.2f};fill:none;stroke-linecap:round;}}"
        f".bond-layer polygon{{stroke:#000;stroke-width:{line_width * 0.6:.2f};fill:#000;stroke-linejoin:round;}}"
        f"text{{font-family:{STRUCTURE_FONT_FAMILY};font-size:{STRUCTURE_LABEL_FONT_SIZE:.2f}px;"
        f"font-weight:400;text-anchor:middle;dominant-baseline:middle;"
        f"paint-order:stroke fill;stroke:#fff;stroke-width:3px;stroke-linejoin:round;}}"
        f"</style>"
        f'<g class="bond-layer">{"".join(bond_fragments)}</g>'
        f'<g class="atom-label-layer">{"".join(text_fragments)}</g>'
        f"</svg>"
    )


def render_smiles_svg(smiles: str) -> str | None:
    """Render SMILES as inline SVG via Open Babel mol block fallback chain."""

    if not smiles:
        return None
    mol_block = generate_smiles_mol_block(smiles)
    if mol_block:
        svg_text = render_mol_block_svg(mol_block)
        if svg_text:
            return svg_text
    return None


# ── Scoring, consensus, and candidate cards ───────────────────────────────────

def confidence_bucket(row: pd.Series) -> str:
    """Map neighbor Tanimoto similarity to high/medium/context tier labels."""

    sim = float(row.get("similarity_to_target", 0.0) or 0.0)
    if sim >= 0.75:
        return "high"
    if sim >= 0.50:
        return "medium"
    return "context"


def similarity_color(similarity_value: float, *, target: bool = False) -> str:
    """Map Tanimoto similarity to a scan-friendly warm color ramp."""
    if target:
        return "#001D37"
    sim = max(0.0, min(1.0, float(similarity_value or 0.0)))
    if sim >= 0.60:
        return "#D73A49"
    if sim >= 0.50:
        return "#F28E2B"
    if sim >= 0.40:
        return "#59A14F"
    return "#4E79A7"


def heatmap_color(value: float) -> str:
    """Blue-to-amber-to-red confidence heatmap."""
    v = max(0.0, min(1.0, float(value or 0.0)))
    stops = [
        (0.0, (44, 123, 182)),
        (0.5, (253, 174, 97)),
        (1.0, (215, 25, 28)),
    ]
    for (left_v, left_rgb), (right_v, right_rgb) in zip(stops, stops[1:]):
        if v <= right_v:
            frac = (v - left_v) / max(right_v - left_v, 1e-9)
            rgb = tuple(round(left_rgb[i] + frac * (right_rgb[i] - left_rgb[i])) for i in range(3))
            return f"#{rgb[0]:02X}{rgb[1]:02X}{rgb[2]:02X}"
    return "#D7191C"


def compact_candidate_scores(candidates: pd.DataFrame) -> pd.DataFrame:
    """Select and round report-facing columns from a candidate score frame."""

    candidates = candidates.copy()
    if "candidate_report_score" not in candidates.columns:
        if "lactoyl_rerank_score" in candidates.columns:
            candidates["candidate_report_score"] = pd.to_numeric(candidates["lactoyl_rerank_score"], errors="coerce")
        elif "rerank_score" in candidates.columns:
            candidates["candidate_report_score"] = pd.to_numeric(candidates["rerank_score"], errors="coerce")
        else:
            candidates["candidate_report_score"] = pd.to_numeric(
                candidates.get("mass_match_score", 0.0), errors="coerce"
            )
    columns = [
        "rank",
        "candidate_report_score",
        "predicted_smiles",
        "predicted_formula",
        "abs_mass_error_da",
        "smiles_spectrum_alignment_score",
        "spectrum_fragment_ion_score",
        "spectrum_neutral_loss_score",
        "lactoyl_rerank_score",
        "lactoyl_hypothesis_score",
        "spectrum_matched_peak_count",
        "spectrum_matched_peaks",
        "lactoyl_mass_score",
        "lactoyl_diagnostic_bin_score",
        "lactoyl_aa_specific_score",
        "mass_match_score",
        "evidence_confidence",
        "predicted_pathway",
        "predicted_subpathway",
        "fusion_source_label",
        "rerank_policy",
        "source_comp_id",
    ]
    present = [col for col in columns if col in candidates.columns]
    out = candidates[present].copy()
    for col in out.columns:
        if pd.api.types.is_numeric_dtype(out[col]):
            out[col] = pd.to_numeric(out[col], errors="coerce").round(6)
    return out


_MCS_SMARTS_CACHE: dict[tuple[str, str], str] = {}


def cached_mcs_smarts(smiles_a: str, mol_a: Any, smiles_b: str, mol_b: Any) -> str:
    """Memoized ``rdFMCS`` smartsString for a SMILES pair (symmetric key).

    The consensus heatmap compares every top-10 reference against every top-10
    candidate, i.e. ~100 MCS calls per report over the same ~10 molecules. MCS is
    symmetric, so caching by the unordered SMILES pair roughly halves the calls;
    each avoided call also skips the expensive ``timeout`` budget when a pair
    would otherwise re-hit the MCS time limit.
    """

    key = (smiles_a, smiles_b) if smiles_a <= smiles_b else (smiles_b, smiles_a)
    cached = _MCS_SMARTS_CACHE.get(key)
    if cached is not None:
        return cached
    try:
        mcs = rdFMCS.FindMCS([mol_a, mol_b], timeout=2, matchValences=False, ringMatchesRingOnly=False)
        smarts = mcs.smartsString
    except Exception:
        smarts = ""
    _MCS_SMARTS_CACHE[key] = smarts
    return smarts


def consensus_atom_support(reference_smiles: str, candidates: pd.DataFrame) -> tuple[list[dict[str, Any]], str]:
    """MCS atom-support heatmap SVG and per-atom support CSV rows for one reference."""

    rows: list[dict[str, Any]] = []
    if not HAS_RDKIT or candidates.empty:
        return rows, "<p class='note'>Consensus structure unavailable because RDKit or candidate rows are missing.</p>"
    ref = Chem.MolFromSmiles(reference_smiles)
    if ref is None:
        return (
            rows,
            "<p class='note'>Consensus structure unavailable because the reference SMILES could not be parsed.</p>",
        )
    compute_2d_coords_robust(ref)
    n_atoms = ref.GetNumAtoms()
    support = [0.0] * n_atoms
    valid = 0
    for _, candidate in candidates.iterrows():
        smiles = str(candidate.get("predicted_smiles", "") or "")
        cand = Chem.MolFromSmiles(smiles)
        if cand is None:
            continue
        smarts = cached_mcs_smarts(reference_smiles, ref, smiles, cand)
        if not smarts:
            continue
        patt = Chem.MolFromSmarts(smarts)
        if patt is None:
            continue
        match = ref.GetSubstructMatch(patt)
        if not match:
            continue
        valid += 1
        for atom_idx in match:
            support[int(atom_idx)] += 1.0
    if valid == 0:
        confidence = [0.0] * n_atoms
    else:
        confidence = [value / valid for value in support]

    conf = ref.GetConformer()
    coords = []
    for idx in range(n_atoms):
        pos = conf.GetAtomPosition(idx)
        coords.append((pos.x, pos.y))
    min_x = min(x for x, _ in coords)
    max_x = max(x for x, _ in coords)
    min_y = min(y for _, y in coords)
    max_y = max(y for _, y in coords)

    # Render at a constant pixels-per-bond density so structures of any size read
    # similarly: a small molecule gets a small canvas (scaled up inside its box,
    # large glyphs) and a large molecule gets a proportionally larger canvas that
    # the viewBox scales down, so atom circles never overlap. Because the density
    # is fixed, glyph sizes are constant rather than fixed against a fixed canvas.
    bond_lengths = []
    for bond in ref.GetBonds():
        ax, ay = coords[bond.GetBeginAtomIdx()]
        bx, by = coords[bond.GetEndAtomIdx()]
        dist = ((ax - bx) ** 2 + (ay - by) ** 2) ** 0.5
        if dist > 1e-6:
            bond_lengths.append(dist)
    typical_bond = sorted(bond_lengths)[len(bond_lengths) // 2] if bond_lengths else 1.5
    target_bond_px = 46.0
    scale = target_bond_px / max(typical_bond, 1e-6)
    radius = 16.0
    symbol_font = 13.0
    conf_font = 9.0
    bond_stroke = 6.0
    circle_stroke = 2.5
    margin = radius + 22.0
    canvas_w = max(max_x - min_x, 1e-6) * scale + 2 * margin
    canvas_h = max(max_y - min_y, 1e-6) * scale + 2 * margin

    def xy(atom_idx: int) -> tuple[float, float]:
        """Map RDKit 2D conformer coordinates to SVG canvas space."""
        x, y = coords[atom_idx]
        return margin + (x - min_x) * scale, canvas_h - (margin + (y - min_y) * scale)

    bond_parts = []
    for bond in ref.GetBonds():
        x1, y1 = xy(bond.GetBeginAtomIdx())
        x2, y2 = xy(bond.GetEndAtomIdx())
        bond_parts.append(
            f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" '
            f'stroke="#24445F" stroke-opacity=".42" stroke-width="{bond_stroke:.1f}" stroke-linecap="round" />'
        )

    atom_parts = []
    for atom in ref.GetAtoms():
        idx = atom.GetIdx()
        value = confidence[idx]
        color = heatmap_color(value)
        symbol = html.escape(atom.GetSymbol())
        x, y = xy(idx)
        atom_parts.append(
            f'<g class="consensus-atom" data-confidence="{value:.3f}" transform="translate({x:.1f},{y:.1f})">'
            f'<circle r="{radius:.1f}" fill="{color}" stroke="white" stroke-width="{circle_stroke:.1f}">'
            f'<title>Atom {idx + 1} {symbol}: consensus {value:.2f}</title></circle>'
            f'<text text-anchor="middle" y="{symbol_font * 0.35:.1f}" font-size="{symbol_font:.1f}">{symbol}</text>'
            f'<text class="atom-confidence" text-anchor="middle" y="{radius + conf_font:.1f}" '
            f'font-size="{conf_font:.1f}">{value:.2f}</text></g>'
        )
        rows.append(
            {
                "atom_index": idx + 1,
                "element": atom.GetSymbol(),
                "consensus_support_fraction": round(value, 6),
                "support_count": round(support[idx], 6),
                "valid_candidate_count": valid,
                "heatmap_color": color,
            }
        )

    svg = "\n".join(
        [
            f"<svg class='consensus-svg' width='{canvas_w:.0f}' height='{canvas_h:.0f}' "
            f"viewBox='0 0 {canvas_w:.0f} {canvas_h:.0f}' "
            "role='img' aria-label='Candidate consensus atom heatmap'>",
            *bond_parts,
            *atom_parts,
            "</svg>",
        ]
    )
    return rows, svg


def build_candidate_consensus_maps(
    pool_df: pd.DataFrame,
) -> tuple[dict[int, str], dict[int, dict[str, Any]], int]:
    """Precompute consensus SVGs with each top-10 candidate as the reference structure."""
    top10 = pool_df.head(10).copy()
    maps: dict[int, str] = {}
    summaries: dict[int, dict[str, Any]] = {}
    for _, row in top10.iterrows():
        rank = int(pd.to_numeric(row.get("rank"), errors="coerce") or 0)
        ref_smiles = str(row.get("predicted_smiles", "") or "").strip()
        if not rank or not ref_smiles:
            continue
        _, svg = consensus_atom_support(ref_smiles, top10)
        maps[rank] = svg
        summaries[rank] = {
            "rank": rank,
            "predicted_formula": str(row.get("predicted_formula", "") or ""),
            "predicted_smiles": ref_smiles,
            "candidate_report_score": row.get("candidate_report_score", ""),
            "abs_mass_error_da": row.get("abs_mass_error_da", ""),
            "evidence_confidence": row.get("evidence_confidence", ""),
            "mass_match_score": row.get("mass_match_score", ""),
        }
    default_rank = 1 if 1 in maps else (min(maps) if maps else 0)
    return maps, summaries, default_rank


def candidate_consensus_scale_html() -> str:
    """HTML legend for the consensus atom-support color scale."""

    return """
          <div class="consensus-scale" aria-label="Consensus heatmap scale">
            <div class="consensus-scale-bar"></div>
            <div class="consensus-scale-labels">
              <span>0.0 low support</span>
              <span>0.5 mixed support</span>
              <span>1.0 high support</span>
            </div>
          </div>"""


def _candidate_score_line_html(row: dict[str, Any]) -> str:
    """One-line evidence readout for a candidate.

    Leads with the pipeline's evidence_confidence (library + mass + spectral) and
    mass match - the meaningful identification signals - then the raw rerank score
    explicitly labeled as a relative tie-breaker (it is a small, uncalibrated
    generative-rerank value used to order otherwise-equivalent isomers).
    """

    def fnum(value: Any, nd: int) -> str | None:
        try:
            return f"{float(value):.{nd}f}"
        except (TypeError, ValueError):
            return None

    parts = []
    ev = fnum(row.get("evidence_confidence"), 3)
    mm = fnum(row.get("mass_match_score"), 3)
    rs = fnum(row.get("candidate_report_score"), 4)
    me = fnum(row.get("abs_mass_error_da"), 6)
    if ev is not None:
        parts.append(f"Evidence confidence: <strong>{ev}</strong>")
    if mm is not None:
        parts.append(f"Mass match: <strong>{mm}</strong>")
    if rs is not None:
        parts.append(f"Rank score (relative): <strong>{rs}</strong>")
    if me is not None:
        parts.append(f"Mass error (Da): <strong>{me}</strong>")
    return ("<p>" + " &middot; ".join(parts) + "</p>") if parts else ""


def candidate_consensus_summary_html(summary: dict[str, Any]) -> str:
    """HTML summary block for the selected top candidate."""

    if not summary:
        return "<p class='note'>No candidate summary available.</p>"
    return (
        "<div class='candidate-summary'>"
        f"<p><strong>Rank {html.escape(str(summary.get('rank', '')))}</strong> &middot; "
        f"{html.escape(str(summary.get('predicted_formula', '')))}</p>"
        f"{_candidate_score_line_html(summary)}"
        f"<p><code>{html.escape(str(summary.get('predicted_smiles', '')))}</code></p>"
        "</div>"
    )


def candidate_structure_card_html(
    row: dict[str, Any],
    consensus_svg: str | None,
    *,
    title: str,
) -> str:
    """HTML card combining structure SVG, scores, and optional consensus overlay."""

    smiles = str(row.get("predicted_smiles", "") or "").strip()
    structure_svg = render_smiles_svg(smiles) if smiles else None
    if structure_svg:
        structure_html = f"<div class='structure-svg-frame'>{structure_svg}</div>"
    else:
        structure_html = "<div class='structure-svg-frame'><p class='note'>Structure unavailable.</p></div>"
    consensus_block = ""
    if consensus_svg:
        consensus_block = f"<div class='consensus-box'>{consensus_svg}</div>" f"{candidate_consensus_scale_html()}"
    return (
        "<div class='structure-card'>"
        f"<h3>{html.escape(title)}</h3>"
        f"{structure_html}"
        f"<p><code>{html.escape(smiles)}</code></p>"
        f"<p class='note'>Formula: <strong>{html.escape(str(row.get('predicted_formula', '')))}</strong></p>"
        f"<p class='note'>{_candidate_score_line_html(row)}</p>"
        f"{consensus_block}"
        "</div>"
    )


# ── Correlation rows (legacy X-25003) ───────────────────────────────────────

def load_real_correlation_rows_LEGACY() -> tuple[list[dict[str, Any]], str | None]:
    """Legacy X-25003 WGCNA correlation lookup (superseded by generic loader)."""

    anchor_names = {
        "x-25003",
        "n-lactoylalanine",
        "n-lactoyl alanine",
    }
    searched = []
    for cohort, spec in wgcna.COHORT_INPUTS.items():
        searched.append(str(spec["corr"]))
        searched.append(str(spec["module"]))
        modules = wgcna.load_module_assignments(cohort)
        if modules.empty:
            continue
        matches = modules[
            modules["plot_name"].astype(str).str.strip().str.lower().isin(anchor_names)
            | modules["chem_id"].astype(str).eq("100022476")
            | modules["plot_name"].astype(str).eq("X-25003")
        ].copy()
        if matches.empty:
            continue
        # A direct anchor exists; future callers can implement row-level correlation extraction here.
        # Today the local workspace has no X-25003/N-lactoylalanine anchor rows, so this branch is not reached.
    blocker = (
        "Real sample-abundance correlation rows for X-25003 are not available in the local workspace. "
        "Searched the configured WGCNA cohort exports, but prior notes and current files both indicate "
        "X-25003 did not survive the module-level variance filter and no direct N-lactoylalanine anchor row is present."
    )
    return [], blocker


def build_placeholder_correlation_rows_LEGACY(
    extended_unknown_rows: list[dict[str, Any]],
    *,
    limit: int = 5,
) -> list[dict[str, Any]]:
    """Legacy proxy correlation rows when WGCNA anchors are absent."""

    proxy_rows: list[dict[str, Any]] = []
    for row in extended_unknown_rows:
        x_id = str(row.get("x_compound_id", "") or "")
        if not x_id or x_id == "X-25003":
            continue
        try:
            proxy_score = float(row.get("similarity_to_target", 0.0) or 0.0)
        except (TypeError, ValueError):
            proxy_score = 0.0
        proxy_rows.append(
            {
                "id": x_id,
                "compound": str(row.get("known_name", "") or x_id),
                "correlation_r": f"{proxy_score:.3f}",
                "pathway": "Verified X structure proxy",
            }
        )
        if len(proxy_rows) >= limit:
            break
    return proxy_rows


def candidate_fragment_and_loss_masses(smiles: str) -> tuple[list[float], list[float]]:
    """Enumerate single-bond fragment and neutral-loss masses for spectrum matching."""

    if not HAS_RDKIT:
        return [], []
    mol = Chem.MolFromSmiles(str(smiles or ""))
    if mol is None:
        return [], []
    parent_mass = exact_mass_from_smiles(smiles)
    if not math.isfinite(parent_mass):
        return [], []
    fragments: set[float] = set()
    losses: set[float] = set()
    for bond in mol.GetBonds():
        if bond.IsInRing() or bond.GetBondTypeAsDouble() != 1.0:
            continue
        try:
            fragged = Chem.FragmentOnBonds(mol, [bond.GetIdx()], addDummies=False)
            for frag in Chem.GetMolFrags(fragged, asMols=True, sanitizeFrags=True):
                mass = float(Descriptors.ExactMolWt(frag))
                if 10.0 <= mass <= parent_mass - 1.0:
                    fragments.add(round(mass, 5))
                    losses.add(round(parent_mass - mass, 5))
        except Exception:
            continue
    return sorted(fragments), sorted(losses)


def slugify_label(text: str) -> str:
    """URL-safe slug for guide section anchor IDs."""

    slug = re.sub(r"[^a-z0-9]+", "-", text.lower()).strip("-")
    return slug or "guide"


def render_guides_html_LEGACY(guide_targets: list[tuple[str, str]]) -> str:
    """Legacy X-25003 placeholder guide page HTML."""

    guide_blocks = (
        "".join(
            f"<section class='guide-card' id='guide-{slug}'><h2>{html.escape(label)}</h2><p>{html.escape(label)} guide will go here.</p></section>"
            for slug, label in guide_targets
        )
        or "<section class='guide-card'><h2>No guides yet</h2><p>Guide content will go here.</p></section>"
    )
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>X-25003 Guide Placeholders</title>
  <style>
    :root {{
      --navy:#001D37; --blue:#28496E; --paper:#FFFFFF; --wash:#F1FAFB; --ink:#15304B; --muted:#5C6B7A; --line:rgba(40,73,110,.18);
      --shadow:0 24px 58px rgba(0,29,55,.10);
    }}
    * {{ box-sizing:border-box; }}
    body {{ margin:0; font-family:"Montserrat","Avenir Next","Segoe UI",Arial,sans-serif; color:var(--ink); background:linear-gradient(180deg,#fff 0,#f7fbfc 100%); }}
    main {{ max-width:1080px; margin:0 auto; padding:32px 24px 56px; }}
    h1 {{ margin:0 0 8px; color:var(--navy); }}
    .lede {{ color:var(--muted); margin:0 0 22px; }}
    .guide-grid {{ display:grid; gap:14px; }}
    .guide-card {{ background:var(--paper); border:1px dashed var(--line); border-radius:10px; padding:18px; box-shadow:var(--shadow); scroll-margin-top:18px; }}
    .guide-card h2 {{ margin:0 0 6px; font-size:20px; color:var(--navy); }}
    .guide-card p {{ margin:0; color:var(--muted); }}
    a.back-link {{ color:var(--blue); text-decoration:none; font-weight:700; }}
  </style>
</head>
<body>
  <main>
    <a class="back-link" href="public.html">Back to public report</a>
    <h1>X-25003 Guide Placeholders</h1>
    <p class="lede">Placeholder reference content for mockup info buttons.</p>
    <div class="guide-grid">{guide_blocks}</div>
  </main>
</body>
</html>
"""


# ── Spectrum alignment ─────────────────────────────────────────────────────────

def score_smiles_spectrum_alignment(
    candidates: pd.DataFrame,
    peaks: pd.DataFrame,
    *,
    tolerance_da: float,
) -> pd.DataFrame:
    """Score candidates by matching fragment ions and neutral losses to peak table."""

    out = candidates.copy()
    default_cols = {
        "smiles_spectrum_alignment_score": 0.0,
        "spectrum_fragment_ion_score": 0.0,
        "spectrum_neutral_loss_score": 0.0,
        "spectrum_matched_peak_count": 0,
        "spectrum_matched_peaks": "",
    }
    if out.empty:
        for col, value in default_cols.items():
            out[col] = value
        return out
    if peaks.empty or "predicted_smiles" not in out.columns:
        for col, value in default_cols.items():
            out[col] = value
        return out
    mass_col = "mass" if "mass" in peaks.columns else "mz" if "mz" in peaks.columns else ""
    intensity_col = (
        "intensity" if "intensity" in peaks.columns else "rel_intensity" if "rel_intensity" in peaks.columns else ""
    )
    parent_col = (
        "parent_mz" if "parent_mz" in peaks.columns else "PARENT_MASS" if "PARENT_MASS" in peaks.columns else ""
    )
    if not mass_col:
        for col, value in default_cols.items():
            out[col] = value
        return out
    work = peaks.copy()
    work["_mz"] = pd.to_numeric(work[mass_col], errors="coerce")
    work["_intensity"] = pd.to_numeric(work[intensity_col], errors="coerce") if intensity_col else 1.0
    work = work[work["_mz"].notna() & work["_mz"].gt(0)].copy()
    if work.empty:
        for col, value in default_cols.items():
            out[col] = value
        return out
    max_i = float(pd.to_numeric(work["_intensity"], errors="coerce").fillna(0.0).max() or 1.0)
    work["_weight"] = pd.to_numeric(work["_intensity"], errors="coerce").fillna(1.0).clip(lower=0.0) / max_i
    if parent_col:
        parent_mz = float(pd.to_numeric(work[parent_col], errors="coerce").dropna().iloc[0])
    else:
        parent_mz = float(work["_mz"].max())
    neutral_parent_mass = parent_mz + 1.007276466812
    work["_neutral_loss"] = neutral_parent_mass - (work["_mz"] - 0.000548579909065)
    total_weight = float(work["_weight"].sum() or 1.0)
    rows: list[dict[str, Any]] = []
    for _, candidate in out.iterrows():
        frag_masses, loss_masses = candidate_fragment_and_loss_masses(str(candidate.get("predicted_smiles", "")))
        matched_ion_weight = 0.0
        matched_loss_weight = 0.0
        matched_peaks: list[float] = []
        for _, peak in work.iterrows():
            mz = float(peak["_mz"])
            weight = float(peak["_weight"])
            ion_match = any(abs(mz - pred) <= float(tolerance_da) for pred in frag_masses)
            loss_match = any(abs(float(peak["_neutral_loss"]) - pred) <= float(tolerance_da) for pred in loss_masses)
            if ion_match:
                matched_ion_weight += weight
            if loss_match:
                matched_loss_weight += weight
            if ion_match or loss_match:
                matched_peaks.append(round(mz, 5))
        ion_score = matched_ion_weight / total_weight
        loss_score = matched_loss_weight / total_weight
        alignment_score = 0.55 * ion_score + 0.45 * loss_score
        row = candidate.to_dict()
        row.update(
            {
                "smiles_spectrum_alignment_score": round(float(alignment_score), 6),
                "spectrum_fragment_ion_score": round(float(ion_score), 6),
                "spectrum_neutral_loss_score": round(float(loss_score), 6),
                "spectrum_matched_peak_count": int(len(set(matched_peaks))),
                "spectrum_matched_peaks": "; ".join(f"{value:.4f}" for value in sorted(set(matched_peaks))),
            }
        )
        rows.append(row)
    return pd.DataFrame(rows)


# ── HTML / CSS ───────────────────────────────────────────────────────────────

def css() -> str:
    """Return Metro-themed CSS shared by public and internal report variants."""

    return """
    :root {
      --metro-navy: #001D37;
      --metro-blue: #28496E;
      --metro-coral: #FF4F57;
      --metro-amber: #FCA750;
      --metro-mist: #D3DEE6;
      --metro-text: #15304B;
      --metro-slate: #4D6070;
      --panel: rgba(255,255,255,0.94);
      --line: rgba(40,73,110,0.16);
      --shadow: 0 20px 54px rgba(0,29,55,0.10);
    }
    * { box-sizing: border-box; }
    body {
      margin: 0;
      font-family: "Montserrat", "Avenir Next", "Segoe UI", sans-serif;
      color: var(--metro-text);
      background:
        radial-gradient(circle at top right, rgba(255,79,87,0.16), transparent 25%),
        radial-gradient(circle at 12% 18%, rgba(252,167,80,0.20), transparent 22%),
        linear-gradient(160deg, #f8fbfd 0%, #eef4f8 44%, #dde8ef 100%);
    }
    .brand-bar { height: 12px; background: linear-gradient(90deg, var(--metro-navy), var(--metro-blue), var(--metro-coral), var(--metro-amber)); }
    .page { width: min(1320px, calc(100vw - 32px)); margin: 0 auto 48px; }
    .panel { background: var(--panel); border: 1px solid var(--line); border-radius: 24px; box-shadow: var(--shadow); }
    .hero { margin-top: 18px; padding: 28px; display: grid; grid-template-columns: 1.15fr 0.85fr; gap: 24px; }
    .eyebrow { display: inline-flex; padding: 8px 14px; border-radius: 999px; background: rgba(0,29,55,0.06); color: var(--metro-blue); font-size: 12px; font-weight: 800; letter-spacing: .12em; text-transform: uppercase; }
    h1 { margin: 14px 0 10px; font-size: clamp(36px, 4.6vw, 56px); line-height: 0.98; color: var(--metro-navy); }
    h2 { margin: 0 0 10px; color: var(--metro-navy); }
    h3 { margin: 0 0 8px; color: var(--metro-blue); }
    p { line-height: 1.65; }
    .summary { max-width: 74ch; color: var(--metro-slate); }
    .chips, .filters, .ext-filters { display: flex; flex-wrap: wrap; gap: 10px; align-items: center; margin-top: 14px; }
    .chip, button.filter, button.ext-filter { border: 1px solid var(--line); background: rgba(255,255,255,.78); color: var(--metro-blue); border-radius: 999px; padding: 9px 13px; font-size: 12px; font-weight: 800; }
    button.filter, button.ext-filter { cursor: pointer; }
    button.filter.active, button.ext-filter.active { background: var(--metro-navy); color: white; }
    .grid { display: grid; grid-template-columns: repeat(12, 1fr); gap: 18px; margin-top: 18px; }
    .tabbar { display: flex; flex-wrap: wrap; gap: 10px; margin-top: 18px; }
    .tab-button { border: 1px solid var(--line); background: rgba(255,255,255,.84); color: var(--metro-blue); border-radius: 12px; padding: 11px 16px; font-size: 13px; font-weight: 900; cursor: pointer; }
    .tab-button.active { background: var(--metro-navy); color: white; }
    .tab-panel { display: none; }
    .tab-panel.active { display: grid; }
    .span-4 { grid-column: span 4; }
    .span-5 { grid-column: span 5; }
    .span-7 { grid-column: span 7; }
    .span-8 { grid-column: span 8; }
    .span-12 { grid-column: span 12; }
    .card { padding: 20px; border-radius: 20px; background: rgba(255,255,255,.84); border: 1px solid var(--line); }
    .metric { display: grid; grid-template-columns: 1fr auto; gap: 10px; padding: 10px 0; border-bottom: 1px solid rgba(40,73,110,.10); }
    .metric:last-child { border-bottom: 0; }
    .metric span:first-child { color: var(--metro-slate); }
    .metric strong { color: var(--metro-navy); }
    .structure { width: 100%; max-height: 330px; object-fit: contain; background: white; border-radius: 18px; border: 1px solid var(--line); }
    .structure-svg-frame { background: white; border-radius: 18px; border: 1px solid var(--line); padding: 14px; }
    .structure-svg-frame svg { width: 100%; height: auto; display: block; }
    .structure-note { margin-top: 10px; color: var(--metro-slate); font-size: 12px; line-height: 1.45; }
    table { width: 100%; border-collapse: collapse; font-size: 13px; }
    th, td { text-align: left; border-bottom: 1px solid rgba(40,73,110,.12); padding: 9px 8px; vertical-align: top; }
    th { color: var(--metro-blue); font-size: 12px; text-transform: uppercase; letter-spacing: .05em; }
    .network-wrap { position: relative; overflow: hidden; min-height: 620px; }
    .extended-network-wrap { position: relative; overflow: auto; min-height: 920px; max-height: 1080px; background: linear-gradient(180deg, rgba(255,255,255,.84), rgba(238,244,248,.84)); border: 1px solid var(--line); border-radius: 16px; }
    svg { width: 100%; height: 620px; display: block; }
    .extended-network-svg { width: 100%; min-width: 1480px; height: auto; display: block; }
    .node text { font-size: 11px; fill: var(--metro-navy); font-weight: 700; }
    .edge { stroke-linecap: round; }
    .node circle { stroke: white; stroke-width: 2.5; }
    .node circle, .ext-node circle { cursor: pointer; }
    .node text { pointer-events: none; }
    .node .sim-label { font-size: 10px; fill: var(--metro-slate); font-weight: 800; }
    .node:hover circle { stroke: var(--metro-navy); stroke-width: 3.5; }
    .node.hidden, .edge.hidden, tr.hidden { display: none; }
    .legend { display: flex; flex-wrap: wrap; gap: 12px; align-items: center; margin-top: 12px; color: var(--metro-slate); font-size: 12px; font-weight: 700; }
    .legend-item { display: inline-flex; align-items: center; gap: 6px; }
    .swatch { width: 16px; height: 16px; border-radius: 50%; display: inline-block; border: 2px solid white; box-shadow: 0 0 0 1px rgba(40,73,110,.18); }
    .network-summary { display: grid; grid-template-columns: repeat(4, minmax(0, 1fr)); gap: 10px; margin-top: 12px; }
    .network-stat { border: 1px solid rgba(40,73,110,.12); border-radius: 12px; padding: 10px; background: rgba(255,255,255,.72); }
    .network-stat span { display: block; color: var(--metro-slate); font-size: 11px; font-weight: 800; text-transform: uppercase; letter-spacing: .04em; }
    .network-stat strong { display: block; margin-top: 3px; font-size: 18px; color: var(--metro-navy); }
    .consensus-layout { display: grid; grid-template-columns: minmax(360px, .95fr) 1.05fr; gap: 18px; align-items: start; }
    .consensus-box { background: white; border: 1px solid var(--line); border-radius: 14px; padding: 12px; overflow: hidden; }
    .consensus-svg { max-width: 100%; height: auto; display: block; margin: 0 auto; }
    .consensus-scale { margin-top: 10px; }
    .consensus-scale-bar { height: 16px; border-radius: 999px; background: linear-gradient(90deg, #2C7BB6 0%, #FDAE61 50%, #D7191C 100%); border: 1px solid rgba(40,73,110,.16); }
    .consensus-scale-labels { display: grid; grid-template-columns: repeat(3, 1fr); gap: 8px; margin-top: 6px; color: var(--metro-slate); font-size: 11px; font-weight: 800; }
    .consensus-scale-labels span:nth-child(2) { text-align: center; }
    .consensus-scale-labels span:nth-child(3) { text-align: right; }
    .consensus-atom text { fill: #fff; font-size: 16px; font-weight: 900; pointer-events: none; }
    .consensus-atom .atom-confidence { fill: var(--metro-navy); font-size: 11px; font-weight: 900; }
    .legend-label { fill: var(--metro-slate); font-size: 11px; font-weight: 900; }
    .candidate-table-wrap { overflow-x: auto; }
    .candidate-table td.smiles { font-family: "SFMono-Regular", "Consolas", monospace; font-size: 12px; max-width: 420px; overflow-wrap: anywhere; }
    .ext-node text { font-size: 10px; fill: var(--metro-navy); font-weight: 800; pointer-events: none; }
    .ext-node.hidden, .ext-edge.hidden, tr.ext-hidden { display: none; }
    .hover-tooltip {
      position: fixed;
      z-index: 1000;
      max-width: 320px;
      padding: 10px 12px;
      border-radius: 12px;
      background: rgba(0, 29, 55, 0.94);
      color: white;
      border: 1px solid rgba(255,255,255,.12);
      box-shadow: 0 10px 28px rgba(0, 29, 55, 0.24);
      font-size: 12px;
      line-height: 1.45;
      pointer-events: none;
      opacity: 0;
      transform: translateY(4px);
      transition: opacity 120ms ease, transform 120ms ease;
    }
    .hover-tooltip.visible { opacity: 1; transform: translateY(0); }
    .hover-tooltip strong { display: block; margin-bottom: 2px; font-size: 13px; }
    .note { color: var(--metro-slate); font-size: 13px; }
    .footer { margin-top: 18px; padding: 18px 4px; color: var(--metro-slate); font-size: 13px; }
    @media (max-width: 900px) {
      .hero { grid-template-columns: 1fr; }
      .span-4, .span-5, .span-7, .span-8 { grid-column: span 12; }
      .consensus-layout { grid-template-columns: 1fr; }
    }
    """


# ── Similarity networks ──────────────────────────────────────────────────────

def network_svg(nodes: list[dict[str, Any]], edges: list[dict[str, Any]]) -> str:
    """Render positioned similarity nodes and edges as inline SVG."""

    by_id = {node["id"]: node for node in nodes}
    edge_parts = []
    node_parts = []
    for edge in edges:
        a = by_id.get(edge["source"])
        b = by_id.get(edge["target"])
        if not a or not b:
            continue
        sim = float(edge.get("similarity", 0.0) or 0.0)
        stroke = html.escape(str(edge.get("color", similarity_color(sim))))
        opacity = 0.24 + min(0.70, max(0.0, sim) * 0.66)
        width = 1.0 + max(0.0, sim - 0.25) * 6.0
        pathway = html.escape(str(edge.get("pathway", "")))
        edge_parts.append(
            f'<line class="edge {html.escape(edge["confidence"])}" data-kind="edge" data-confidence="{html.escape(edge["confidence"])}" '
            f'data-similarity="{sim:.3f}" data-pathway="{pathway}" '
            f'data-source-id="{html.escape(str(edge["source"]))}" data-target-id="{html.escape(str(edge["target"]))}" '
            f'stroke="{stroke}" stroke-opacity="{opacity:.3f}" stroke-width="{width:.2f}" '
            f'stroke-linecap="round" x1="{a["x"]:.1f}" y1="{a["y"]:.1f}" x2="{b["x"]:.1f}" y2="{b["y"]:.1f}">'
            f'<title>{html.escape(edge["label"])}</title></line>'
        )
    for node in nodes:
        label = html.escape(str(node["label"]))
        similarity_value = str(node.get("similarity", ""))
        sim_float = float(similarity_value) if similarity_value not in {"", "target"} else 1.0
        fill = html.escape(str(node.get("color", similarity_color(sim_float, target=node.get("kind") == "target"))))
        pathway = html.escape(str(node.get("pathway", "")))
        data_name = html.escape(str(node.get("label", "")))
        raw_label = str(node.get("label", "")).replace("\\n", " ").replace("\n", " ")
        hover_label = html.escape(raw_label)
        tooltip_html = html.escape(f"{raw_label}\nSimilarity: {similarity_value}\n{str(node.get('pathway', ''))}")
        node_parts.append(
            f'<g class="node {html.escape(node["confidence"])} {html.escape(node.get("kind", ""))}" '
            f'data-node-id="{html.escape(str(node.get("id", "")))}" data-kind="{html.escape(str(node.get("kind", "")))}" '
            f'data-confidence="{html.escape(node["confidence"])}" data-similarity="{sim_float:.3f}" '
            f'data-pathway="{pathway}" data-name="{data_name}" data-hover-label="{hover_label}" '
            f'data-tooltip="{tooltip_html}" tabindex="0" transform="translate({node["x"]:.1f},{node["y"]:.1f})">'
            f'<circle r="{node["radius"]:.1f}" fill="{fill}"><title>{label}\\nSimilarity: {similarity_value}\\n{pathway}</title></circle>'
            f'<text x="{node["radius"] + 7:.1f}" y="-2">{label}</text>'
            f'<text class="sim-label" x="{node["radius"] + 7:.1f}" y="13">sim {html.escape(similarity_value)}</text></g>'
        )
    return "\n".join(
        [
            "<svg viewBox='0 0 1000 620' role='img' aria-label='Chemical similarity network'>",
            *edge_parts,
            *node_parts,
            "</svg>",
        ]
    )


def build_pairwise_similarity_edges(
    nodes: list[dict[str, Any]],
    *,
    target_edge_min: float,
    pair_edge_min: float,
    max_edges_per_node: int,
) -> list[dict[str, Any]]:
    """Build capped Morgan Tanimoto edges with target vs pair thresholds."""

    scored: list[tuple[float, str, str, dict[str, Any]]] = []
    for i, left in enumerate(nodes):
        for right in nodes[i + 1 :]:
            sim = similarity(str(left["smiles"]), str(right["smiles"]))
            if "target" in {left["kind"], right["kind"]}:
                if sim < target_edge_min:
                    continue
            elif sim < pair_edge_min:
                continue
            row = {
                "source": str(left["id"]),
                "target": str(right["id"]),
                "similarity": sim,
                "confidence": confidence_bucket(pd.Series({"similarity_to_target": sim})),
                "pathway": f"{left.get('pathway', '')} <-> {right.get('pathway', '')}",
                "color": similarity_color(sim),
                "label": f"{left['label']} <-> {right['label']} · Tanimoto {sim:.3f}",
            }
            scored.append((sim, str(left["id"]), str(right["id"]), row))
    scored.sort(key=lambda item: item[0], reverse=True)

    kept: list[dict[str, Any]] = []
    degree: dict[str, int] = {str(node["id"]): 0 for node in nodes}
    for sim, left_id, right_id, row in scored:
        if "X-25003" in {left_id, right_id}:
            kept.append(row)
            degree[left_id] += 1
            degree[right_id] += 1
            continue
        if degree[left_id] >= max_edges_per_node or degree[right_id] >= max_edges_per_node:
            continue
        kept.append(row)
        degree[left_id] += 1
        degree[right_id] += 1
    return kept


def greedy_similarity_clusters(
    nodes: list[dict[str, Any]],
    *,
    cluster_threshold: float,
) -> list[list[str]]:
    """Group nodes into clades by medoid similarity, avoiding chain-expanded giant components."""
    node_map = {str(node["id"]): node for node in nodes}
    target_id = next((str(node["id"]) for node in nodes if node["kind"] == "target"), None)
    target_cluster = [target_id] if target_id else []
    clusters: list[list[str]] = [target_cluster] if target_cluster else []
    medoids: list[str] = [target_id] if target_id else []

    ordered = sorted(
        [node for node in nodes if str(node["id"]) != target_id],
        key=lambda node: float(node.get("similarity_to_target", 0.0) or 0.0),
        reverse=True,
    )

    for node in ordered:
        node_id = str(node["id"])
        best_idx = None
        best_sim = -1.0
        for idx, medoid_id in enumerate(medoids):
            if medoid_id is None:
                continue
            sim = similarity(str(node["smiles"]), str(node_map[medoid_id]["smiles"]))
            if sim >= cluster_threshold and sim > best_sim:
                best_sim = sim
                best_idx = idx
        if best_idx is None:
            clusters.append([node_id])
            medoids.append(node_id)
        else:
            clusters[best_idx].append(node_id)
    return [cluster for cluster in clusters if cluster]


def pairwise_similarity_clusters(
    nodes: list[dict[str, Any]],
    edges: list[dict[str, Any]],
    *,
    cluster_threshold: float,
) -> list[list[str]]:
    """Cluster visible compounds by non-target pairwise similarity edges."""
    node_ids = [str(node["id"]) for node in nodes]
    target_id = next((str(node["id"]) for node in nodes if node["kind"] == "target"), None)
    parent = {node_id: node_id for node_id in node_ids}

    def find(node_id: str) -> str:
        """Union-find root with path compression."""
        while parent[node_id] != node_id:
            parent[node_id] = parent[parent[node_id]]
            node_id = parent[node_id]
        return node_id

    def union(left: str, right: str) -> None:
        """Merge two node clusters by rank."""
        left_root = find(left)
        right_root = find(right)
        if left_root != right_root:
            parent[right_root] = left_root

    for edge in edges:
        left_id = str(edge.get("source", ""))
        right_id = str(edge.get("target", ""))
        if left_id not in parent or right_id not in parent:
            continue
        if target_id in {left_id, right_id}:
            continue
        if float(edge.get("similarity", 0.0) or 0.0) < cluster_threshold:
            continue
        union(left_id, right_id)

    grouped: dict[str, list[str]] = {}
    for node_id in node_ids:
        if node_id == target_id:
            continue
        grouped.setdefault(find(node_id), []).append(node_id)

    clusters: list[list[str]] = [[target_id]] if target_id else []
    clusters.extend(grouped.values())
    return [cluster for cluster in clusters if cluster]


def cluster_layout_network(
    nodes: list[dict[str, Any]],
    edges: list[dict[str, Any]],
    *,
    width: float = 1680.0,
    height: float = 1180.0,
    cluster_threshold: float = 0.50,
) -> list[dict[str, Any]]:
    """Assign x/y coordinates via cluster-aware force layout for network SVG."""

    node_map = {str(node["id"]): dict(node) for node in nodes}
    target_id = next((str(node["id"]) for node in nodes if node["kind"] == "target"), None)
    clusters = pairwise_similarity_clusters(nodes, edges, cluster_threshold=cluster_threshold)
    cluster_info: list[dict[str, Any]] = []
    for cluster in clusters:
        members = [node_map[node_id] for node_id in cluster]
        max_target_sim = max(float(member.get("similarity_to_target", 0.0) or 0.0) for member in members)
        target_cluster = target_id in cluster if target_id else False
        cluster_info.append(
            {
                "members": cluster,
                "target_cluster": target_cluster,
                "max_target_sim": max_target_sim,
                "size": len(cluster),
            }
        )
    cluster_info.sort(
        key=lambda info: (
            0 if info["target_cluster"] else 1,
            -info["max_target_sim"],
            -info["size"],
        )
    )

    positions: dict[str, tuple[float, float]] = {}
    center_x = width / 2.0
    center_y = height / 2.0

    non_target_clusters = [info for info in cluster_info if not info["target_cluster"]]
    cursor = 0
    ring_idx = 0
    while cursor < len(non_target_clusters):
        capacity = 8 + ring_idx * 8
        ring_clusters = non_target_clusters[cursor : cursor + capacity]
        if not ring_clusters:
            break
        radius = 185.0 + ring_idx * 118.0
        y_scale = 0.76 + min(ring_idx, 2) * 0.04
        for idx, info in enumerate(ring_clusters):
            angle = (2.0 * math.pi * idx) / len(ring_clusters) - math.pi / 2.0
            info["center"] = (
                center_x + math.cos(angle) * radius,
                center_y + math.sin(angle) * radius * y_scale,
            )
        cursor += len(ring_clusters)
        ring_idx += 1

    for info in cluster_info:
        if info["target_cluster"]:
            info["center"] = (center_x, center_y)

    for info in cluster_info:
        members = [node_map[node_id] for node_id in info["members"]]
        cluster_center_x, cluster_center_y = info["center"]
        if any(member["kind"] == "target" for member in members):
            medoid = max(
                members,
                key=lambda member: (
                    1 if member["kind"] == "target" else 0,
                    float(member.get("similarity_to_target", 0.0) or 0.0),
                ),
            )
        else:
            medoid = max(
                members,
                key=lambda member: sum(
                    similarity(str(member["smiles"]), str(other["smiles"]))
                    for other in members
                    if str(other["id"]) != str(member["id"])
                ),
            )
        medoid_id = str(medoid["id"])
        positions[medoid_id] = (cluster_center_x, cluster_center_y)
        others = [member for member in members if str(member["id"]) != medoid_id]
        if not others:
            continue
        ordered = sorted(
            others,
            key=lambda member: similarity(str(member["smiles"]), str(medoid["smiles"])),
            reverse=True,
        )
        for idx, member in enumerate(ordered):
            sim_to_medoid = similarity(str(member["smiles"]), str(medoid["smiles"]))
            angle = (2.0 * math.pi * idx) / max(len(ordered), 1)
            radius = 52.0 + (1.0 - sim_to_medoid) * 72.0 + (idx // 10) * 22.0
            positions[str(member["id"])] = (
                cluster_center_x + math.cos(angle) * radius,
                cluster_center_y + math.sin(angle) * radius * 0.86,
            )

    laid_out: list[dict[str, Any]] = []
    for node in nodes:
        node_id = str(node["id"])
        x, y = positions.get(node_id, (center_x, center_y))
        item = dict(node)
        item["x"] = x
        item["y"] = y
        laid_out.append(item)

    by_id = {str(node["id"]): i for i, node in enumerate(laid_out)}
    edge_pairs: list[tuple[int, int, float, bool]] = []
    for edge in edges:
        left_id = str(edge["source"])
        right_id = str(edge["target"])
        if left_id not in by_id or right_id not in by_id:
            continue
        edge_pairs.append(
            (
                by_id[left_id],
                by_id[right_id],
                float(edge.get("similarity", 0.0) or 0.0),
                target_id in {left_id, right_id},
            )
        )

    margin = 70.0
    damping = 0.84
    repulsion = 2100.0
    for _ in range(180):
        forces = [[0.0, 0.0] for _ in laid_out]
        for i in range(len(laid_out)):
            for j in range(i + 1, len(laid_out)):
                left = laid_out[i]
                right = laid_out[j]
                dx = right["x"] - left["x"]
                dy = right["y"] - left["y"]
                dist_sq = dx * dx + dy * dy + 1.0
                dist = math.sqrt(dist_sq)
                push = repulsion / dist_sq
                fx = push * dx / dist
                fy = push * dy / dist
                forces[i][0] -= fx
                forces[i][1] -= fy
                forces[j][0] += fx
                forces[j][1] += fy

        for left_i, right_i, sim, is_target_edge in edge_pairs:
            left = laid_out[left_i]
            right = laid_out[right_i]
            dx = right["x"] - left["x"]
            dy = right["y"] - left["y"]
            dist = math.sqrt(dx * dx + dy * dy) + 1e-6
            if is_target_edge:
                desired = 150.0 + (1.0 - sim) * 210.0
                strength = 0.0055 * (0.65 + sim)
            else:
                desired = 42.0 + (1.0 - sim) * 175.0
                strength = 0.024 * (0.7 + sim)
            pull = strength * (dist - desired)
            fx = pull * dx / dist
            fy = pull * dy / dist
            forces[left_i][0] += fx
            forces[left_i][1] += fy
            forces[right_i][0] -= fx
            forces[right_i][1] -= fy

        for i, node in enumerate(laid_out):
            if node["kind"] == "target":
                node["x"] = center_x
                node["y"] = center_y
                continue
            cluster_pull_x = (center_x - node["x"]) * 0.0010
            cluster_pull_y = (center_y - node["y"]) * 0.0010
            node["x"] = min(width - margin, max(margin, node["x"] + (forces[i][0] + cluster_pull_x) * damping))
            node["y"] = min(height - margin, max(margin, node["y"] + (forces[i][1] + cluster_pull_y) * damping))

    padding = 12.0
    for _ in range(140):
        moved = False
        for i in range(len(laid_out)):
            for j in range(i + 1, len(laid_out)):
                left = laid_out[i]
                right = laid_out[j]
                dx = right["x"] - left["x"]
                dy = right["y"] - left["y"]
                dist = math.sqrt(dx * dx + dy * dy) + 1e-6
                min_dist = float(left["radius"]) + float(right["radius"]) + padding
                if dist >= min_dist:
                    continue
                overlap = (min_dist - dist) / 2.0
                ux = dx / dist
                uy = dy / dist
                if left["kind"] != "target":
                    left["x"] -= ux * overlap
                    left["y"] -= uy * overlap
                if right["kind"] != "target":
                    right["x"] += ux * overlap
                    right["y"] += uy * overlap
                moved = True
        for node in laid_out:
            if node["kind"] == "target":
                node["x"] = center_x
                node["y"] = center_y
                continue
            node["x"] = min(width - margin, max(margin, node["x"]))
            node["y"] = min(height - margin, max(margin, node["y"]))
        if not moved:
            break

    for node in laid_out:
        node["x"] = round(float(node["x"]), 3)
        node["y"] = round(float(node["y"]), 3)
    return laid_out


def force_layout_network(
    nodes: list[dict[str, Any]],
    edges: list[dict[str, Any]],
    *,
    width: float = 1320.0,
    height: float = 920.0,
    iterations: int = 260,
) -> list[dict[str, Any]]:
    """Simple force-directed layout when cluster layout is not used."""

    center_x = width / 2.0
    center_y = height / 2.0
    ids = [str(node["id"]) for node in nodes]
    idx = {node_id: i for i, node_id in enumerate(ids)}
    positions: list[list[float]] = []
    velocities: list[list[float]] = []
    n = max(len(nodes), 1)
    for i, node in enumerate(nodes):
        if node["kind"] == "target":
            positions.append([center_x, center_y])
            velocities.append([0.0, 0.0])
            continue
        angle = (2.0 * math.pi * i) / n
        if node["kind"] == "unknown":
            radius = 220.0
        else:
            radius = 330.0
        positions.append([center_x + math.cos(angle) * radius, center_y + math.sin(angle) * radius * 0.82])
        velocities.append([0.0, 0.0])

    edge_pairs = []
    for edge in edges:
        if edge["source"] not in idx or edge["target"] not in idx:
            continue
        edge_pairs.append((idx[edge["source"]], idx[edge["target"]], float(edge["similarity"])))

    repulsion = 3200.0
    attraction = 0.020
    damping = 0.88
    margin = 46.0
    for _ in range(iterations):
        forces = [[0.0, 0.0] for _ in nodes]
        for i in range(len(nodes)):
            for j in range(i + 1, len(nodes)):
                dx = positions[j][0] - positions[i][0]
                dy = positions[j][1] - positions[i][1]
                dist_sq = dx * dx + dy * dy + 1.0
                dist = math.sqrt(dist_sq)
                push = repulsion / dist_sq
                fx = push * dx / dist
                fy = push * dy / dist
                forces[i][0] -= fx
                forces[i][1] -= fy
                forces[j][0] += fx
                forces[j][1] += fy

        for left_i, right_i, sim in edge_pairs:
            dx = positions[right_i][0] - positions[left_i][0]
            dy = positions[right_i][1] - positions[left_i][1]
            dist = math.sqrt(dx * dx + dy * dy) + 1e-6
            desired = 120.0 + (1.0 - sim) * 220.0
            pull = attraction * (dist - desired) * (0.7 + sim)
            fx = pull * dx / dist
            fy = pull * dy / dist
            forces[left_i][0] += fx
            forces[left_i][1] += fy
            forces[right_i][0] -= fx
            forces[right_i][1] -= fy

        for i, node in enumerate(nodes):
            if node["kind"] == "target":
                positions[i] = [center_x, center_y]
                velocities[i] = [0.0, 0.0]
                continue
            pull_x = (center_x - positions[i][0]) * 0.0014
            pull_y = (center_y - positions[i][1]) * 0.0014
            velocities[i][0] = (velocities[i][0] + forces[i][0] + pull_x) * damping
            velocities[i][1] = (velocities[i][1] + forces[i][1] + pull_y) * damping
            positions[i][0] = min(width - margin, max(margin, positions[i][0] + velocities[i][0]))
            positions[i][1] = min(height - margin, max(margin, positions[i][1] + velocities[i][1]))

    laid_out = []
    for i, node in enumerate(nodes):
        item = dict(node)
        item["x"] = round(positions[i][0], 3)
        item["y"] = round(positions[i][1], 3)
        laid_out.append(item)
    return laid_out


def extended_network_svg(nodes: list[dict[str, Any]], edges: list[dict[str, Any]]) -> str:
    """Render the expanded known+X similarity network as large-format SVG."""

    by_id = {str(node["id"]): node for node in nodes}
    target_id = next((str(node["id"]) for node in nodes if node.get("kind") == "target"), "target")
    edge_parts = []
    node_parts = []
    for edge in edges:
        left = by_id.get(str(edge["source"]))
        right = by_id.get(str(edge["target"]))
        if not left or not right:
            continue
        sim = float(edge.get("similarity", 0.0) or 0.0)
        edge_parts.append(
            f'<line class="ext-edge" data-kind="edge" data-confidence="{html.escape(str(edge["confidence"]))}" '
            f'data-similarity="{sim:.3f}" data-pathway="{html.escape(str(edge.get("pathway", "")))}" '
            f'data-source-id="{html.escape(str(edge["source"]))}" data-target-id="{html.escape(str(edge["target"]))}" '
            f'data-source-kind="{html.escape(str(left["kind"]))}" data-target-kind="{html.escape(str(right["kind"]))}" '
            f'data-source-name="{html.escape(str(left["label"]))}" data-target-name="{html.escape(str(right["label"]))}" '
            f'x1="{left["x"]:.1f}" y1="{left["y"]:.1f}" x2="{right["x"]:.1f}" y2="{right["y"]:.1f}" '
            f'stroke="#B8C5D0" '
            f'stroke-opacity="{0.08 + min(0.22, sim * 0.18):.3f}" stroke-width="{0.7 + sim * 1.8:.2f}" '
            f'stroke-linecap="round"><title>{html.escape(str(edge["label"]))}</title></line>'
        )
    for node in nodes:
        sim = float(node.get("similarity_to_target", 0.0) or 0.0)
        display_label = str(node.get("display_label", "") or "")
        stroke = "#0B2F4A" if node["kind"] == "unknown" else "white"
        stroke_width = 2.6 if node["kind"] == "unknown" else 1.8
        if node["kind"] == "target":
            stroke = "#001D37"
            stroke_width = 3.4
        hover_label = html.escape(str(node.get("label", "")).replace("\\n", " ").replace("\n", " "))
        label_html = (
            f'<text x="{node["radius"] + 6:.1f}" y="-1">{html.escape(display_label)}</text>' if display_label else ""
        )
        node_parts.append(
            f'<g class="ext-node {html.escape(str(node["kind"]))}" data-node-id="{html.escape(str(node["id"]))}" '
            f'data-kind="{html.escape(str(node["kind"]))}" '
            f'data-confidence="{html.escape(str(node["confidence"]))}" data-similarity="{sim:.3f}" '
            f'data-name="{html.escape(str(node["label"]))}" data-hover-label="{hover_label}" '
            f'data-pathway="{html.escape(str(node.get("pathway", "")))}" tabindex="0" '
            f'transform="translate({node["x"]:.1f},{node["y"]:.1f})">'
            f'<circle r="{node["radius"]:.1f}" fill="{html.escape(str(node["color"]))}" stroke="{stroke}" stroke-width="{stroke_width:.1f}">'
            f'<title>{html.escape(str(node["label"]))}\\nType: {html.escape(str(node["kind"]))}\\nSimilarity to {html.escape(target_id)}: {sim:.3f}\\n{html.escape(str(node.get("pathway", "")))}</title>'
            f"</circle>"
            f"{label_html}</g>"
        )
    return "\n".join(
        [
            "<svg class='extended-network-svg' viewBox='0 0 1680 1180' role='img' aria-label='Expanded known and X similarity network'>",
            *edge_parts,
            *node_parts,
            "</svg>",
        ]
    )


# ── Structural possibility tree (X-25003 template) ─────────────────────────────

def x25003_structural_possibility_summary(
    *,
    candidate_scores: pd.DataFrame,
    neighbors: list[dict[str, Any]],
    target_smiles: str,
    target_formula: str,
) -> list[dict[str, Any]]:
    """Summarize lactoyl-amino-acid structural possibility branches for decision tree."""

    if candidate_scores.empty:
        return []
    smiles = candidate_scores.get("predicted_smiles", pd.Series(dtype=str)).astype(str)
    formulas = candidate_scores.get("predicted_formula", pd.Series(dtype=str)).astype(str)
    exact_mask = smiles.eq(target_smiles)
    formula_mask = formulas.eq(target_formula)
    lactoyl_text = candidate_scores.apply(
        lambda row: " ".join(str(value) for value in row.tolist()), axis=1
    ).str.lower()
    lactoyl_like_mask = lactoyl_text.str.contains("lactoyl|ala|alanine", regex=True, na=False) | exact_mask
    high_known_count = sum(float(row.get("similarity_to_target", 0.0) or 0.0) >= 0.60 for row in neighbors)
    lactoyl_neighbor_count = sum(
        "lactoyl" in f"{row.get('canonical_name', '')} {row.get('metabolon_sub', '')}".lower()
        or "alanine" in f"{row.get('canonical_name', '')} {row.get('metabolon_sub', '')}".lower()
        for row in neighbors
    )

    def top_smiles(mask: pd.Series) -> str:
        """Return rank-1 SMILES among candidates matching ``mask``."""
        rows = candidate_scores.loc[mask] if len(mask) == len(candidate_scores) else pd.DataFrame()
        if rows.empty:
            return ""
        return str(rows.iloc[0].get("predicted_smiles", ""))

    return [
        {
            "structural_possibility_class": "protected_n_lactoylalanine",
            "candidate_count": int(exact_mask.sum()),
            "top_candidate_label": "N-lactoylalanine",
            "top_candidate_smiles": target_smiles,
            "reader_interpretation": "Protected lactoyl-amino-acid rule candidate; current final call.",
            "dominant_gate_status": "rank1 protected candidate; formula and mass supported",
        },
        {
            "structural_possibility_class": "formula_matched_lactoyl_or_alanine_isomers",
            "candidate_count": int((formula_mask & lactoyl_like_mask & ~exact_mask).sum()),
            "top_candidate_label": "Formula-matched lactoyl/alanine-like alternatives",
            "top_candidate_smiles": top_smiles(formula_mask & lactoyl_like_mask & ~exact_mask),
            "reader_interpretation": "Formula-matched candidates retaining lactoyl, alanine, or amino-acid context but not the exact protected call.",
            "dominant_gate_status": "alternative hypotheses below protected rank1",
        },
        {
            "structural_possibility_class": "formula_matched_non_lactoyl_isomers",
            "candidate_count": int((formula_mask & ~lactoyl_like_mask).sum()),
            "top_candidate_label": "Formula-matched non-lactoyl isomers",
            "top_candidate_smiles": top_smiles(formula_mask & ~lactoyl_like_mask),
            "reader_interpretation": "Formula-compatible structural isomers without clear lactoyl/alanine support in local evidence fields.",
            "dominant_gate_status": "lower-ranked structural alternatives",
        },
        {
            "structural_possibility_class": "adduct_or_formula_shift_alternatives",
            "candidate_count": int((~formula_mask).sum()),
            "top_candidate_label": "Lower-formula/adduct alternatives",
            "top_candidate_smiles": top_smiles(~formula_mask),
            "reader_interpretation": "Candidates explained by alternate adduct or formula-shift hypotheses rather than the target neutral formula.",
            "dominant_gate_status": "not primary formula hypothesis",
        },
        {
            "structural_possibility_class": "known_lactoyl_alanine_neighborhood",
            "candidate_count": int(lactoyl_neighbor_count),
            "top_candidate_label": f"{high_known_count} known neighbors at Tanimoto >= 0.60",
            "top_candidate_smiles": "",
            "reader_interpretation": "Known-metabolite neighborhood supports lactoyl/alanine biochemical context around the final call.",
            "dominant_gate_status": "context support, not independent structure proof",
        },
    ]


def x25003_structural_possibility_tree_svg(rows: list[dict[str, Any]]) -> str:
    """Static SVG decision tree annotated with branch candidate counts."""

    by_id = {str(row["structural_possibility_class"]): row for row in rows}

    def count(row_id: str) -> str:
        """Format candidate count for a structural possibility branch leaf."""
        value = by_id.get(row_id, {}).get("candidate_count", 0)
        return str(value)

    return f"""<svg class='decision-tree-svg' viewBox='0 0 1160 640' role='img' aria-label='X-25003 structural possibility decision tree'>
<defs><marker id='x25003-tree-arrow' viewBox='0 0 10 10' refX='9' refY='5' markerWidth='6' markerHeight='6' orient='auto-start-reverse'><path d='M 0 0 L 10 5 L 0 10 z' fill='#6B7C8A' /></marker></defs>
<rect width='1160' height='640' rx='18' fill='#F7F9FB' />
<text class='tree-title' x='60' y='38'>Structural possibility branches from lactoyl-amino-acid evidence</text>
<g class='tree-edges'>
  <path d='M230,320 C250,320 250,170 275,170' fill='none' stroke='#6B7C8A' stroke-width='3' marker-end='url(#x25003-tree-arrow)' /><text class='tree-edge-label' x='250' y='238'>start</text>
  <path d='M455,170 C480,170 480,110 510,110' fill='none' stroke='#4F8F61' stroke-width='3' marker-end='url(#x25003-tree-arrow)' /><text class='tree-edge-label' x='482' y='133'>yes</text>
  <path d='M455,170 C480,170 480,290 510,290' fill='none' stroke='#9AA7B2' stroke-width='3' stroke-dasharray='7 6' marker-end='url(#x25003-tree-arrow)' /><text class='tree-edge-label' x='482' y='228'>no / variant</text>
  <path d='M690,110 C720,110 720,92 755,92' fill='none' stroke='#2E7D32' stroke-width='3' marker-end='url(#x25003-tree-arrow)' /><text class='tree-edge-label' x='722' y='94'>exact</text>
  <path d='M690,110 C720,110 720,202 755,202' fill='none' stroke='#F28E2B' stroke-width='3' marker-end='url(#x25003-tree-arrow)' /><text class='tree-edge-label' x='722' y='153'>isomer</text>
  <path d='M690,290 C720,290 720,312 755,312' fill='none' stroke='#F28E2B' stroke-width='3' marker-end='url(#x25003-tree-arrow)' /><text class='tree-edge-label' x='722' y='294'>formula-only</text>
  <path d='M690,290 C720,290 720,422 755,422' fill='none' stroke='#9AA7B2' stroke-width='3' stroke-dasharray='7 6' marker-end='url(#x25003-tree-arrow)' /><text class='tree-edge-label' x='722' y='352'>formula shift</text>
  <path d='M230,320 C470,320 470,532 755,532' fill='none' stroke='#4F8F61' stroke-width='3' marker-end='url(#x25003-tree-arrow)' /><text class='tree-edge-label' x='486' y='514'>known-neighbor context</text>
</g>
<g class='tree-nodes'>
  <g class='tree-node' data-id='root'><rect x='60' y='288' width='170' height='64' rx='12' fill='#1F4E79' stroke='#24445F' stroke-width='2' /><text class='tree-node-title' fill='#FFFFFF' x='145' y='308'><tspan x='145' dy='0'>X-25003 structural</tspan><tspan x='145' dy='15'>possibilities</tspan></text></g>
  <g class='tree-node' data-id='lactoyl_decision'><rect x='275' y='126' width='180' height='88' rx='8' fill='#FFFFFF' stroke='#C8D3DC' stroke-width='2' /><text class='tree-node-title' fill='#17212B' x='365' y='148'><tspan x='365' dy='0'>Lactoyl amino-acid</tspan><tspan x='365' dy='15'>rule support?</tspan></text><text class='tree-node-metric' x='365' y='194'>rank1 protected call</text></g>
  <g class='tree-node' data-id='formula_decision'><rect x='510' y='246' width='180' height='88' rx='8' fill='#FFFFFF' stroke='#C8D3DC' stroke-width='2' /><text class='tree-node-title' fill='#17212B' x='600' y='268'><tspan x='600' dy='0'>Target formula</tspan><tspan x='600' dy='15'>retained?</tspan></text><text class='tree-node-metric' x='600' y='314'>C6H11NO4 gate</text></g>
  <g class='tree-node' data-id='exact_decision'><rect x='510' y='66' width='180' height='88' rx='8' fill='#FFFFFF' stroke='#C8D3DC' stroke-width='2' /><text class='tree-node-title' fill='#17212B' x='600' y='88'><tspan x='600' dy='0'>Exact final candidate?</tspan></text><text class='tree-node-metric' x='600' y='134'>N-lactoylalanine</text></g>
  <g class='tree-node' data-id='protected_n_lactoylalanine'><rect x='755' y='52' width='300' height='80' rx='12' fill='#2E7D32' stroke='#C8D3DC' stroke-width='2' /><text class='tree-node-title' fill='#FFFFFF' x='905' y='74'><tspan x='905' dy='0'>Protected N-lactoylalanine</tspan></text><text class='tree-leaf-count' x='905' y='104'>n={count("protected_n_lactoylalanine")}</text></g>
  <g class='tree-node' data-id='formula_matched_lactoyl_or_alanine_isomers'><rect x='755' y='162' width='300' height='80' rx='12' fill='#4F8F61' stroke='#C8D3DC' stroke-width='2' /><text class='tree-node-title' fill='#FFFFFF' x='905' y='184'><tspan x='905' dy='0'>Formula-matched lactoyl/alanine</tspan><tspan x='905' dy='15'>isomers</tspan></text><text class='tree-leaf-count' x='905' y='224'>n={count("formula_matched_lactoyl_or_alanine_isomers")}</text></g>
  <g class='tree-node' data-id='formula_matched_non_lactoyl_isomers'><rect x='755' y='272' width='300' height='80' rx='12' fill='#F28E2B' stroke='#C8D3DC' stroke-width='2' /><text class='tree-node-title' fill='#FFFFFF' x='905' y='294'><tspan x='905' dy='0'>Formula-matched non-lactoyl</tspan><tspan x='905' dy='15'>isomers</tspan></text><text class='tree-leaf-count' x='905' y='334'>n={count("formula_matched_non_lactoyl_isomers")}</text></g>
  <g class='tree-node' data-id='adduct_or_formula_shift_alternatives'><rect x='755' y='382' width='300' height='80' rx='12' fill='#7A8691' stroke='#C8D3DC' stroke-width='2' /><text class='tree-node-title' fill='#FFFFFF' x='905' y='404'><tspan x='905' dy='0'>Adduct/formula-shift alternatives</tspan></text><text class='tree-leaf-count' x='905' y='444'>n={count("adduct_or_formula_shift_alternatives")}</text></g>
  <g class='tree-node' data-id='known_lactoyl_alanine_neighborhood'><rect x='755' y='492' width='300' height='80' rx='12' fill='#1F4E79' stroke='#C8D3DC' stroke-width='2' /><text class='tree-node-title' fill='#FFFFFF' x='905' y='514'><tspan x='905' dy='0'>Known lactoyl/alanine</tspan><tspan x='905' dy='15'>neighborhood</tspan></text><text class='tree-leaf-count' x='905' y='554'>n={count("known_lactoyl_alanine_neighborhood")}</text></g>
</g>
</svg>"""


def build_main_similarity_network(
    *,
    target_smiles: str,
    neighbors: list[dict[str, Any]],
    x_neighbors: list[dict[str, Any]],
    target_edge_min: float,
    pair_edge_min: float,
    max_edges_per_node: int,
) -> dict[str, Any]:
    """Build the public bottom network: known + verified X neighbors with pairwise Morgan edges."""
    nodes: list[dict[str, Any]] = [
        {
            "id": "X-25003",
            "label": "X-25003\\nN-lactoylalanine",
            "kind": "target",
            "smiles": target_smiles,
            "similarity_to_target": 1.0,
            "confidence": "target",
            "similarity": "1.000",
            "pathway": "Amino Acid / Lactoyl Amino Acid",
            "color": similarity_color(1.0, target=True),
            "radius": 32.0,
        }
    ]
    for row in neighbors:
        sim = float(row["similarity_to_target"])
        nodes.append(
            {
                "id": str(row["chemical_id"]),
                "label": str(row["canonical_name"]),
                "kind": "known",
                "smiles": str(row["smiles"]),
                "similarity_to_target": sim,
                "confidence": str(row["confidence"]),
                "similarity": f"{sim:.3f}",
                "pathway": f"{row.get('metabolon_super', '')} / {row.get('metabolon_sub', '')}",
                "color": similarity_color(sim),
                "radius": 11.0 + sim * 18.0,
            }
        )
    for row in x_neighbors:
        sim = float(row["similarity_to_target"])
        x_id = str(row.get("x_compound_id", "") or "")
        known_name = str(row.get("known_name", "") or "")
        label = x_id if x_id else known_name
        if known_name:
            label = f"{label}\\n{known_name}"
        nodes.append(
            {
                "id": x_id or str(row.get("known_x_id", "")),
                "label": label,
                "kind": "unknown",
                "smiles": str(row["smiles"]),
                "similarity_to_target": sim,
                "confidence": str(row["confidence"]),
                "similarity": f"{sim:.3f}",
                "pathway": "Verified X structure",
                "color": similarity_color(sim),
                "radius": 10.0 + sim * 14.0,
            }
        )
    edges = build_pairwise_similarity_edges(
        nodes,
        target_edge_min=target_edge_min,
        pair_edge_min=pair_edge_min,
        max_edges_per_node=max_edges_per_node,
    )
    nodes = cluster_layout_network(nodes, edges, width=1000.0, height=620.0)
    return {"render_nodes": nodes, "render_edges": edges}


def build_extended_similarity_network(
    *,
    target_smiles: str,
    all_structured: pd.DataFrame,
    known_x_all: pd.DataFrame,
    max_known_render: int,
    known_min_similarity: float,
    x_min_similarity: float,
    target_edge_min: float,
    pair_edge_min: float,
    max_edges_per_node: int,
    control_id: str = "X-00000",
    target_name: str = "target",
) -> dict[str, Any]:
    """Build expanded known+X similarity network with render and full neighbor tables."""
    known_rows: list[dict[str, Any]] = []
    _, target_fp = mol_fp(target_smiles)
    if target_fp is None:
        return {
            "known_all": pd.DataFrame(),
            "x_all": pd.DataFrame(),
            "render_nodes": [],
            "render_edges": [],
        }

    for _, row in all_structured.iterrows():
        smiles = str(row.get("smiles") or "").split("|")[0].strip()
        mol, fp = mol_fp(smiles)
        if mol is None or fp is None:
            continue
        sim = float(TanimotoSimilarity(target_fp, fp))
        item = row.to_dict()
        item["smiles"] = smiles
        item["similarity_to_target"] = sim
        known_rows.append(item)
    known_all = pd.DataFrame(known_rows)
    if not known_all.empty:
        known_all = (
            known_all.sort_values("similarity_to_target", ascending=False, kind="mergesort")
            .drop_duplicates(subset=["smiles"], keep="first")
            .reset_index(drop=True)
        )
        known_all["node_type"] = "known"
        known_all["confidence"] = known_all.apply(confidence_bucket, axis=1)

    x_rows: list[dict[str, Any]] = []
    for _, row in known_x_all.iterrows():
        smiles = str(row.get("smiles") or "").strip()
        mol, fp = mol_fp(smiles)
        if mol is None or fp is None:
            continue
        sim = float(TanimotoSimilarity(target_fp, fp))
        item = row.to_dict()
        item["similarity_to_target"] = sim
        x_rows.append(item)
    x_all = pd.DataFrame(x_rows)
    if not x_all.empty:
        x_all = x_all.sort_values("similarity_to_target", ascending=False, kind="mergesort").reset_index(drop=True)
        x_all["node_type"] = "unknown"
        x_all["confidence"] = x_all.apply(confidence_bucket, axis=1)

    if known_all.empty:
        known_render = known_all.copy()
    else:
        known_render = (
            known_all[known_all["similarity_to_target"] >= known_min_similarity].head(max_known_render).copy()
        )
    if x_all.empty:
        x_render = x_all.copy()
    else:
        x_render = x_all[x_all["similarity_to_target"] >= x_min_similarity].copy()
        if not x_render.empty and "x_compound_id" in x_render.columns:
            x_render = x_render[x_render["x_compound_id"].astype(str) != control_id].copy()

    nodes: list[dict[str, Any]] = [
        {
            "id": control_id,
            "label": f"{control_id} · {target_name}",
            "display_label": control_id,
            "kind": "target",
            "smiles": target_smiles,
            "similarity_to_target": 1.0,
            "confidence": "target",
            "pathway": "Amino Acid / Lactoyl Amino Acid",
            "radius": 16.0,
            "color": similarity_color(1.0, target=True),
        }
    ]

    for _, row in x_render.iterrows():
        sim = float(row["similarity_to_target"])
        name = str(row.get("known_name", ""))
        label = str(row.get("x_compound_id", ""))
        display_label = label if (sim >= 0.56 or "lactoyl" in name.lower()) else ""
        nodes.append(
            {
                "id": str(row["x_compound_id"]),
                "label": f"{row['x_compound_id']} · {name}",
                "display_label": display_label,
                "kind": "unknown",
                "smiles": str(row["smiles"]),
                "similarity_to_target": sim,
                "confidence": str(row["confidence"]),
                "pathway": "Verified X structure",
                "radius": 6.5 + sim * 6.5,
                "color": similarity_color(sim),
            }
        )

    for _, row in known_render.iterrows():
        sim = float(row["similarity_to_target"])
        name = str(row.get("canonical_name", ""))
        display_label = name if (sim >= 0.60 or "lactoyl" in name.lower()) else ""
        nodes.append(
            {
                "id": f"K-{row['chemical_id']}",
                "label": name,
                "display_label": display_label,
                "kind": "known",
                "smiles": str(row["smiles"]),
                "similarity_to_target": sim,
                "confidence": str(row["confidence"]),
                "pathway": f"{row.get('metabolon_super', '')} / {row.get('metabolon_sub', '')}",
                "radius": 5.5 + sim * 6.0,
                "color": similarity_color(sim),
            }
        )

    edges = build_pairwise_similarity_edges(
        nodes,
        target_edge_min=target_edge_min,
        pair_edge_min=pair_edge_min,
        max_edges_per_node=max_edges_per_node,
    )
    nodes = cluster_layout_network(nodes, edges)
    return {
        "known_all": known_all,
        "x_all": x_all,
        "render_nodes": nodes,
        "render_edges": edges,
        "known_render": known_render,
        "x_render": x_render,
    }


def resolve_pathway_display(db_detail: dict[str, Any], rank1: Any) -> dict[str, Any]:
    """Resolve super/sub pathway for the report.

    Prefers an exact metabolite-DB record; otherwise falls back to the neighbor
    consensus prediction carried on the rank-1 rerank row (``predicted_pathway`` /
    ``predicted_subpathway`` + ``pathway_confidence_tier`` /
    ``pathway_best_reference_similarity``, produced by the pathway enrichment).
    Returns super/sub plus the provenance so the report can label predicted vs
    exact honestly.
    """

    def clean(value: Any) -> str:
        text = str("" if value is None else value).strip()
        return "" if text.lower() in {"", "nan", "none", "null", "<na>"} else text

    exact_super = clean(db_detail.get("super_pathway")) or clean(db_detail.get("metabolon_super"))
    exact_sub = clean(db_detail.get("pathway")) or clean(db_detail.get("metabolon_sub"))
    if exact_super or exact_sub:
        return {
            "super": exact_super,
            "sub": exact_sub,
            "source": "exact",
            "tier": "high",
            "similarity": 1.0,
            "status": "exact_reference_match",
            "note": "From reference DB record",
        }

    pred_super = clean(rank1.get("predicted_pathway"))
    pred_sub = clean(rank1.get("predicted_subpathway"))
    if pred_super or pred_sub:
        tier = clean(rank1.get("pathway_confidence_tier")) or "low"
        status = clean(rank1.get("pathway_assignment_status")) or "neighbor_prediction"
        similarity = 0.0
        sim_raw = rank1.get("pathway_best_reference_similarity")
        if sim_raw is not None and not pd.isna(sim_raw):
            try:
                similarity = float(sim_raw)
            except (TypeError, ValueError):
                similarity = 0.0
        sim_txt = (
            f"nearest reference similarity {similarity:.2f}, {tier} confidence"
            if similarity > 0
            else f"{tier} confidence"
        )
        return {
            "super": pred_super,
            "sub": pred_sub,
            "source": "predicted",
            "tier": tier,
            "similarity": similarity,
            "status": status,
            "note": f"Predicted ({sim_txt})",
        }

    return {
        "super": "",
        "sub": "",
        "source": "none",
        "tier": "",
        "similarity": 0.0,
        "status": clean(rank1.get("pathway_assignment_status")) or "no_reference_match",
        "note": "No reference match above similarity threshold",
    }


def render_html(
    *,
    out_dir: Path,
    target: dict[str, Any],
    evidence: dict[str, Any],
    db_detail: dict[str, Any],
    confidence: Any,
    conf_fields: dict[str, Any],
    pathway_info: dict[str, Any],
    synonyms: list[dict[str, Any]],
    neighbors: list[dict[str, Any]],
    nodes: list[dict[str, Any]],
    edges: list[dict[str, Any]],
    candidate_scores: list[dict[str, Any]],
    fragment_candidate_rows: list[dict[str, Any]] | None = None,
    structure_panel_html: str,
    consensus_svg: str,
    candidate_consensus_maps: dict[int, str],
    candidate_consensus_summaries: dict[int, dict[str, Any]],
    default_consensus_rank: int,
    structural_possibility_tree_svg: str,
    structural_possibility_rows: list[dict[str, Any]],
    feature_total_count: int,
    extended_network_svg_html: str,
    extended_unknown_rows: list[dict[str, Any]],
    extended_known_rows: list[dict[str, Any]],
    extended_stats: dict[str, Any],
    correlation_rows: list[dict[str, Any]],
    correlation_blocker: str | None,
    variant: str,
) -> str:
    """Assemble full HTML report page for public or internal ``variant``."""

    public_mode = variant == "public"
    variant_label = "Public" if variant == "public" else "Internal"
    control_id = str(target.get("control_id", "")).strip() or "compound"
    tentative_descriptor = target.get("tentative_descriptor") or {}
    report_title = str(tentative_descriptor.get("label") or control_id)
    descriptor_confidence = str(tentative_descriptor.get("confidence") or "").strip().title()
    descriptor_badge_html = ""
    if descriptor_confidence:
        descriptor_badge_html = (
            '<!--TENTATIVE_DESCRIPTOR_BADGE-->'
            '<span class="badge" title="Descriptor specificity only; not identification confidence">'
            f"Tentative descriptor \N{MIDDLE DOT} {html.escape(descriptor_confidence)}</span>"
            "<!--/TENTATIVE_DESCRIPTOR_BADGE-->"
        )
    slug = control_slug(control_id)
    pathway_super = str(pathway_info.get("super") or "").strip()
    pathway_badge_html = f"<span class='badge'>{html.escape(pathway_super)}</span>" if pathway_super else ""
    # Top-candidate structural similarity: Tanimoto of the rank-1 predicted structure
    # to its nearest reference metabolite (1.0 when it is an exact reference match).
    # This is the meaningful "how close is our call" number; the old header showed
    # ``score_rank`` (an internal, mostly-empty pool field) mislabeled as Tanimoto.
    try:
        top_similarity = float(pathway_info.get("similarity") or 0.0)
    except (TypeError, ValueError):
        top_similarity = 0.0
    pathway_source = str(pathway_info.get("source") or "").strip()
    if pathway_source == "exact":
        top_similarity_badge = "Exact reference match"
    elif top_similarity > 0:
        top_similarity_badge = f"Nearest-ref Tanimoto {top_similarity:.2f}"
    else:
        top_similarity_badge = "No reference match"
    pubchem_link = f"https://pubchem.ncbi.nlm.nih.gov/#query={quote(str(target['name']))}"
    top_candidate = candidate_scores[0] if candidate_scores else {}
    synonym_lines = [
        html.escape(str(row.get("synonym", ""))) for row in synonyms if str(row.get("synonym", "")).strip()
    ]
    synonym_card_html = (
        "".join(f"<div>{line}</div>" for line in synonym_lines) if synonym_lines else "<div>No synonyms returned.</div>"
    )
    library_neighbors = neighbors[:10]
    parameter_cards_html = confmod.identification_parameter_cards_html(
        confidence,
        candidate_scores=candidate_scores,
        pathway_info=pathway_info,
        top_similarity=top_similarity,
        neighbor_count=len(neighbors),
    )
    top10_candidates = candidate_scores[:10]
    default_rank = (
        default_consensus_rank
        if default_consensus_rank in candidate_consensus_maps
        else (min(candidate_consensus_maps) if candidate_consensus_maps else 0)
    )
    default_consensus_svg = candidate_consensus_maps.get(default_rank, consensus_svg)
    default_summary = candidate_consensus_summaries.get(default_rank, {})
    candidate_select_options = []
    for row in top10_candidates:
        rank = int(pd.to_numeric(row.get("rank"), errors="coerce") or 0)
        if rank not in candidate_consensus_maps:
            continue
        formula = str(row.get("predicted_formula", "") or "n/a")
        ev = row.get("evidence_confidence", "")
        try:
            ev_label = f" &middot; conf {float(ev):.3f}"
        except (TypeError, ValueError):
            ev_label = ""
        selected = " selected" if rank == default_rank else ""
        candidate_select_options.append(
            f'<option value="{rank}"{selected}>Rank {rank} &middot; {html.escape(formula)}{ev_label}</option>'
        )
    candidate_select_html = (
        "\n".join(candidate_select_options)
        if candidate_select_options
        else '<option value="">No candidates available</option>'
    )
    candidate_consensus_data_json = json.dumps(
        {
            "maps": {str(k): v for k, v in candidate_consensus_maps.items()},
            "summaries": {str(k): v for k, v in candidate_consensus_summaries.items()},
        }
    ).replace("</", "<\\/")
    fragment_rows = fragment_candidate_rows if fragment_candidate_rows is not None else candidate_scores
    substructure_rows = confmod.diagnostic_substructure_rows_html(top_candidate, evidence)
    fragment_summary_html = confmod.fragment_summary_section_html(
        conf_fields,
        fragment_rows,
        public=public_mode,
    )
    mass_rows = (
        "<tr>"
        f"<td>{target['exact_mass']:.6f}</td>"
        f"<td>{evidence['parent_mz']:.4f}</td>"
        "<td>Negative (M - H)</td>"
        f"<td>{evidence['mass_error_da']:.6f}</td>"
        "</tr>"
    )
    structural_possibility_table_rows = (
        "\n".join(
            "<tr>"
            f"<td>{html.escape(str(row.get('structural_possibility_class', '')))}</td>"
            f"<td>{html.escape(str(row.get('candidate_count', '')))}</td>"
            f"<td>{html.escape(str(row.get('top_candidate_label', '')))}</td>"
            f"<td>{html.escape(str(row.get('top_candidate_smiles', '')))}</td>"
            f"<td>{html.escape(str(row.get('reader_interpretation', '')))}</td>"
            f"<td>{html.escape(str(row.get('dominant_gate_status', '')))}</td>"
            "</tr>"
            for row in structural_possibility_rows
        )
        or "<tr><td colspan='6'>No structural possibility rows available.</td></tr>"
    )
    ext_unknown_table = (
        "\n".join(
            f"<tr data-kind='unknown' data-similarity='{float(row['similarity_to_target']):.3f}'>"
            f"<td>{html.escape(str(row['x_compound_id']))}</td><td>{html.escape(str(row['known_name']))}</td><td>{float(row['similarity_to_target']):.3f}</td></tr>"
            for row in extended_unknown_rows
        )
        or "<tr><td colspan='3'>No X unknown rows available.</td></tr>"
    )
    correlating_rows = (
        "\n".join(
            f"<tr><td>{html.escape(str(row.get('id', '')))}</td><td>{html.escape(str(row.get('compound', '')))}</td><td>{html.escape(str(row.get('correlation_r', '')))}</td><td>{html.escape(str(row.get('pathway', '')))}</td></tr>"
            for row in correlation_rows
        )
        or "<tr><td colspan='4'>No real correlation rows available.</td></tr>"
    )
    ext_known_table = (
        "\n".join(
            f"<tr data-kind='known' data-similarity='{float(row['similarity_to_target']):.3f}'>"
            f"<td>{html.escape(str(row['chemical_id']))}</td><td>{html.escape(str(row['canonical_name']))}</td><td>{float(row['similarity_to_target']):.3f}</td></tr>"
            for row in extended_known_rows
        )
        or "<tr><td colspan='3'>No known rows available.</td></tr>"
    )
    library_header_html = (
        "<thead><tr><th>ID</th><th>Compound</th><th>Similarity</th><th>Pathway</th></tr></thead>"
        if public_mode
        else "<thead><tr><th>ID</th><th>Compound</th><th>Similarity</th><th>Pathway</th><th>Confidence</th></tr></thead>"
    )
    library_row_html = (
        "\n".join(
            (
                f"<tr data-kind='known' data-confidence='{html.escape(row['confidence'])}' data-similarity='{float(row['similarity_to_target']):.3f}' "
                f"data-name='{html.escape(str(row['canonical_name']))}' data-pathway='{html.escape(str(row.get('metabolon_super', '')))} / {html.escape(str(row.get('metabolon_sub', '')))}' "
                f"style='border-left: 6px solid {similarity_color(float(row['similarity_to_target']))};'>"
                f"<td>{html.escape(row['chemical_id'])}</td><td>{html.escape(row['canonical_name'])}</td><td>{row['similarity_to_target']:.3f}</td>"
                f"<td>{html.escape(row['metabolon_super'])} / {html.escape(row['metabolon_sub'])}</td>"
                + ("" if public_mode else f"<td>{html.escape(row['confidence'])}</td>")
                + "</tr>"
            )
            for row in library_neighbors
        )
        or f"<tr><td colspan='{'4' if public_mode else '5'}'>No library matches available.</td></tr>"
    )
    top10_title = "Top 10 potential structures" if public_mode else "Top 10 Candidate Structures"
    fragment_title = "Fragment Scores" if public_mode else "Substructure Scores"
    ms_fit_title = "MS Fit Scores" if public_mode else "Feature Match Scores"
    params_title = "Parameters Contributing to ID Confidence"
    synonyms_title = "Metabolite Synonyms" if public_mode else "Synonyms Of Top Candidate"
    guide_targets: list[tuple[str, str]] = []

    def info_button(label: str) -> str:
        """Emit public-mode info link and record guide target for guides.html."""
        if not public_mode:
            return ""
        slug = slugify_label(label)
        guide_targets.append((slug, label))
        return (
            f' <a class="info-button" href="guides.html#guide-{slug}" '
            f'aria-label="More information about {html.escape(label)}" '
            f'title="{html.escape(label)} guide will go here">i</a>'
        )

    def heading_with_info(tag: str, label: str, *, include_info: bool = True) -> str:
        """Wrap a section heading with optional public-mode info button."""
        if not include_info:
            return f"<{tag}>{html.escape(label)}</{tag}>"
        return f"<{tag}>{html.escape(label)}{info_button(label)}</{tag}>"

    header_metrics_html = ""
    if not public_mode:
        header_metrics_html = f"""
        <div class="metric-card"><div class="metric-label">Top candidate</div><div class="metric-value">{html.escape(target['name'])}</div><div class="metric-note">Rank-1 displayed candidate</div></div>
        <div class="metric-card"><div class="metric-label">Candidate rows</div><div class="metric-value">{len(candidate_scores)}</div><div class="metric-note">Displayed SMILES candidates</div></div>
        <div class="metric-card"><div class="metric-label">Known neighbors</div><div class="metric-value">{len(neighbors)}</div><div class="metric-note">Structured DB neighborhood</div></div>
        <div class="metric-card"><div class="metric-label">Mass error</div><div class="metric-value">{evidence['mass_error_da']:.6f}</div><div class="metric-note">Neutral mass Da</div></div>"""
    top10_section = f"""
    <section class="section">
      {heading_with_info("h2", top10_title)}
      <p class="note">Select a candidate to view its consensus atom-position heatmap. The selected structure is the reference; atom circles show maximum-common-substructure retention across the top-10 candidate pool.</p>
      <p class="note"><strong>Reading the scores:</strong> "Evidence confidence" (library + mass + spectral corroboration) and "Mass match" are the identification signals. The "Rank score (relative)" is a small, uncalibrated value the reranker uses only to order otherwise-equivalent candidates (e.g. mass-consistent isomers), so a low absolute value does not indicate a poor match.</p>
      <div class="candidate-picker">
        <label for="candidate-consensus-select">Candidate</label>
        <select id="candidate-consensus-select" aria-label="Select candidate structure for consensus heatmap">
          {candidate_select_html}
        </select>
      </div>
      <div id="candidate-consensus-summary">{candidate_consensus_summary_html(default_summary)}</div>
      <div class="consensus-box" id="candidate-consensus-panel">{default_consensus_svg}</div>
      {candidate_consensus_scale_html()}
      <script type="application/json" id="candidate-consensus-data">{candidate_consensus_data_json}</script>
    </section>"""
    correlation_note = (
        html.escape(correlation_blocker)
        if correlation_blocker
        else "In addition to spectral information, we also evaluate the correlation of metabolite abundance values across samples."
    )
    similarity_section = f"""
    <section class="section">
      {heading_with_info("h2", "Similarity Network")}
      <p class="note">Map of similar molecules with distances proportionate to the degree of similarity.</p>
      {heading_with_info("h3", "Predicted Candidate Similarity Network", include_info=False)}
      <p class="note">Nodes include both structured known metabolites and verified X-compound structures from the same expanded local similarity graph used in the internal view.</p>
      <div class="ext-filters">
        <button class="ext-filter active" data-filter="all">All</button>
        <button class="ext-filter" data-filter="known">Known only</button>
        <button class="ext-filter" data-filter="unknown">X only</button>
        <button class="ext-filter" data-filter="ge060">Similarity &gt;= 0.60</button>
        <button class="ext-filter" data-filter="lactoyl">Lactoyl / amino acid context</button>
      </div>
      <div class="network-summary">
        <div class="network-stat"><span>Known in full graph</span><strong>{int(extended_stats.get('n_known_all', 0))}</strong></div>
        <div class="network-stat"><span>X compounds in full graph</span><strong>{int(extended_stats.get('n_x_all', 0))}</strong></div>
        <div class="network-stat"><span>Known rendered</span><strong>{int(extended_stats.get('n_known_render', 0))}</strong></div>
        <div class="network-stat"><span>X rendered</span><strong>{int(extended_stats.get('n_x_render', 0))}</strong></div>
      </div>
      <div class="extended-network-wrap">{extended_network_svg_html}</div>
      <p class="note">Sidecars: <a href="data/extended_similarity_network_nodes.csv">nodes CSV</a>, <a href="data/extended_similarity_network_edges.csv">edges CSV</a>, <a href="data/all_known_similarity_to_{slug}.csv">all known similarities</a>, <a href="data/all_unknown_similarity_to_{slug}.csv">all X similarities</a>.</p>
    </section>"""
    internal_sections = ""
    if variant == "internal":
        internal_sections = f"""

    <section class="section">
      <h2>Structure Consensus Heat Map</h2>
      <p class="note">
        The final consensus atom-position confidence map is shown beside the most likely structure at the top of the report.
        This compound currently has one candidate-consensus support map, so that same map is used as the final consensus-guided confidence visual.
      </p>
      <p class="note">Sidecars: <a href="data/{slug}_consensus_atom_support.csv">{slug}_consensus_atom_support.csv</a> and <a href="data/{slug}_evidence_summary.csv">{slug}_evidence_summary.csv</a>.</p>
    </section>

    <section class="section">
      <h2>Structural Possibility Decision Tree</h2>
      <p class="note">Expanded internal view showing the main structure classes retained after mass, fragment, and neighborhood gating.</p>
      <div class="decision-tree-box">{structural_possibility_tree_svg}</div>
      <p class="note">Sidecar: <a href="data/{slug}_structural_possibility_summary.csv">{slug}_structural_possibility_summary.csv</a>.</p>
      <div class="table-wrap">
        <table>
          <thead><tr><th>Class</th><th>Candidate Count</th><th>Top Candidate</th><th>Top Candidate SMILES</th><th>Interpretation</th><th>Dominant Gate</th></tr></thead>
          <tbody>{structural_possibility_table_rows}</tbody>
        </table>
      </div>
    </section>

    <section class="section">
      <h2>Expanded Known + X Similarity Network</h2>
      <p class="note">Internal context view including both structured known metabolites and neighboring X-compounds from the broader local similarity graph.</p>
      <div class="ext-filters">
        <button class="ext-filter active" data-filter="all">All</button>
        <button class="ext-filter" data-filter="known">Known only</button>
        <button class="ext-filter" data-filter="unknown">X only</button>
        <button class="ext-filter" data-filter="ge060">Similarity &gt;= 0.60</button>
        <button class="ext-filter" data-filter="lactoyl">Lactoyl / amino acid context</button>
      </div>
      <div class="network-summary">
        <div class="network-stat"><span>Known in full graph</span><strong>{int(extended_stats.get('n_known_all', 0))}</strong></div>
        <div class="network-stat"><span>X compounds in full graph</span><strong>{int(extended_stats.get('n_x_all', 0))}</strong></div>
        <div class="network-stat"><span>Known rendered</span><strong>{int(extended_stats.get('n_known_render', 0))}</strong></div>
        <div class="network-stat"><span>X rendered</span><strong>{int(extended_stats.get('n_x_render', 0))}</strong></div>
      </div>
      <div class="extended-network-wrap">{extended_network_svg_html}</div>
      <p class="note">Sidecars: <a href="data/extended_similarity_network_nodes.csv">nodes CSV</a>, <a href="data/extended_similarity_network_edges.csv">edges CSV</a>, <a href="data/all_known_similarity_to_{slug}.csv">all known similarities</a>, <a href="data/all_unknown_similarity_to_{slug}.csv">all X similarities</a>.</p>
    </section>

    <section class="section">
      <h2>Closest X And Known Tables</h2>
      <div class="call-grid">
        <div>
          <h3>Closest X Compounds</h3>
          <div class="table-wrap">
            <table>
              <thead><tr><th>X ID</th><th>Closest Known</th><th>Similarity</th></tr></thead>
              <tbody>{ext_unknown_table}</tbody>
            </table>
          </div>
        </div>
        <div>
          <h3>Closest Known Metabolites</h3>
          <div class="table-wrap">
            <table>
              <thead><tr><th>ID</th><th>Compound</th><th>Similarity</th></tr></thead>
              <tbody>{ext_known_table}</tbody>
            </table>
          </div>
        </div>
      </div>
    </section>"""
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>{html.escape(report_title)} Structure Resolution ({variant_label})</title>
  <style>
    :root {{
      color-scheme: light;
      --navy:#001D37; --blue:#28496E; --coral:#FF4F57; --purple:#855DA3; --amber:#FCA750; --mist:#D3DEE6;
      --paper:#FFFFFF; --wash:#F1FAFB; --ink:#15304B; --muted:#5C6B7A; --line:rgba(40,73,110,.18);
      --green:#0E8A68; --shadow:0 24px 58px rgba(0,29,55,.14);
    }}
    * {{ box-sizing:border-box; }}
    body {{ margin:0; font-family:"Montserrat","Avenir Next","Segoe UI",Arial,sans-serif; color:var(--ink); background:linear-gradient(180deg,#fff 0,#f7fbfc 22%,#eef6f9 100%); line-height:1.5; }}
    header {{ position:relative; padding:54px 54px 34px; background:linear-gradient(135deg,rgba(0,29,55,.98),rgba(40,73,110,.92) 56%,rgba(255,79,87,.88)); color:#fff; border-bottom:0; overflow:hidden; }}
    header::after {{ content:""; position:absolute; inset:auto 0 0 0; height:10px; background:linear-gradient(90deg,var(--amber),var(--coral) 36%,var(--purple) 68%,var(--mist)); }}
    main {{ max-width:1380px; margin:0 auto; padding:28px 28px 64px; }}
    h1 {{ margin:0; font-size:46px; line-height:1; letter-spacing:0; }}
    h2 {{ margin:0 0 16px; font-size:26px; line-height:1.08; letter-spacing:0; color:var(--navy); }}
    h3 {{ margin:20px 0 10px; font-size:18px; letter-spacing:0; color:var(--navy); }}
    p {{ margin:8px 0; }}
    .eyebrow {{ margin:0 0 10px; color:rgba(255,255,255,.74); font-size:12px; font-weight:800; letter-spacing:.12em; text-transform:uppercase; }}
    .subtitle {{ margin-top:10px; max-width:1040px; color:rgba(255,255,255,.88); font-size:17px; }}
    .hero-grid {{ display:grid; grid-template-columns:minmax(0,1.2fr) minmax(280px,.8fr); gap:24px; align-items:start; }}
    .section {{ position:relative; background:var(--paper); border:1px solid var(--line); border-radius:8px; padding:24px; margin:20px 0; box-shadow:var(--shadow); overflow:hidden; }}
    .section::before {{ content:""; display:block; width:88px; height:5px; border-radius:999px; background:linear-gradient(90deg,var(--navy),var(--coral) 58%,var(--amber)); margin-bottom:16px; }}
    .metrics {{ display:grid; grid-template-columns:repeat(4,minmax(0,1fr)); gap:12px; margin-top:18px; }}
    .metric-card {{ background:linear-gradient(180deg,rgba(241,250,251,.98),#fff); border:1px solid var(--line); border-radius:8px; padding:13px 14px; box-shadow:0 10px 24px rgba(0,29,55,.06); }}
    .copyable-card {{ position:relative; cursor:pointer; transition:border-color 120ms ease, box-shadow 120ms ease, transform 120ms ease; }}
    .copyable-card:hover,.copyable-card:focus {{ border-color:rgba(40,73,110,.45); box-shadow:0 14px 28px rgba(0,29,55,.1); transform:translateY(-1px); outline:none; }}
    .copy-hint {{ position:absolute; top:10px; right:10px; opacity:0; pointer-events:none; font-size:11px; font-weight:800; color:var(--blue); background:#fff; border:1px solid var(--line); border-radius:999px; padding:3px 8px; transition:opacity 120ms ease; }}
    .copyable-card:hover .copy-hint,.copyable-card:focus .copy-hint {{ opacity:1; }}
    .sub-block {{ margin-top:18px; padding-top:14px; border-top:1px solid var(--line); }}
    .sub-block h4 {{ margin:0 0 10px; font-size:15px; color:var(--navy); }}
    .info-button {{ display:inline-flex; align-items:center; justify-content:center; width:20px; height:20px; margin-left:8px; border-radius:50%; border:1px solid var(--line); background:#f7fbfc; color:var(--blue); font-size:12px; font-weight:800; text-decoration:none; vertical-align:middle; }}
    .info-button:hover,.info-button:focus {{ background:var(--blue); color:#fff; border-color:var(--blue); outline:none; }}
    .metric-list {{ margin-top:6px; color:var(--ink); font-size:14px; line-height:1.45; }}
    .metric-list div + div {{ margin-top:4px; }}
    .metric-label {{ color:var(--blue); font-size:12px; text-transform:uppercase; letter-spacing:.08em; font-weight:800; }}
    .metric-value {{ margin-top:4px; font-size:23px; font-weight:700; overflow-wrap:anywhere; color:var(--navy); }}
    .metric-note {{ color:var(--muted); font-size:12px; margin-top:2px; }}
    .call-grid {{ display:grid; grid-template-columns:repeat(2,minmax(0,1fr)); gap:18px; }}
    .structure-card {{ border:1px solid var(--line); border-radius:8px; padding:16px; background:linear-gradient(180deg,rgba(241,250,251,.9),#fff 38%); }}
    .structure-card img {{ width:100%; max-height:260px; object-fit:contain; background:#fff; border:1px solid #e8edf1; border-radius:6px; }}
    .structure-svg-frame {{ background:#fff; border:1px solid #e8edf1; border-radius:6px; padding:12px; }}
    .structure-svg-frame svg {{ width:100%; height:auto; max-height:300px; display:block; }}
    .structure-note {{ margin-top:10px; color:var(--muted); font-size:12px; line-height:1.45; }}
    .badge {{ display:inline-block; border-radius:999px; padding:5px 10px; font-size:12px; font-weight:800; border:1px solid rgba(255,255,255,.16); background:rgba(255,255,255,.09); color:#fff; margin-right:6px; }}
    .badge.ok {{ color:#fff; background:rgba(14,138,104,.24); border-color:rgba(14,138,104,.42); }}
    .badge.warn {{ color:#fff; background:rgba(252,167,80,.2); border-color:rgba(252,167,80,.4); }}
    .badge.fail {{ color:#fff; background:rgba(255,79,87,.2); border-color:rgba(255,79,87,.4); }}
    .note {{ color:var(--muted); }}
    .table-wrap,.candidate-table-wrap {{ overflow:auto; border:1px solid var(--line); border-radius:8px; }}
    table {{ border-collapse:collapse; width:100%; min-width:860px; background:#fff; }}
    th,td {{ padding:8px 10px; border-bottom:1px solid #e8edf1; text-align:left; vertical-align:top; font-size:13px; }}
    th {{ position:sticky; top:0; background:linear-gradient(180deg,#f7fbfc,#eef6f9); z-index:1; font-size:12px; color:var(--blue); }}
    code {{ background:#eef6f9; border-radius:4px; padding:1px 4px; overflow-wrap:anywhere; }}
    .network-controls,.filters,.ext-filters,.legend {{ display:flex; gap:10px; flex-wrap:wrap; align-items:center; margin:12px 0; }}
    .filters button,.ext-filters button {{ border:1px solid var(--line); background:#fff; border-radius:6px; padding:7px 10px; cursor:pointer; }}
    .filters button.active,.ext-filters button.active {{ background:#17324d; color:#fff; border-color:#17324d; }}
    .legend {{ color:var(--muted); font-size:12px; font-weight:700; }}
    .legend-item {{ display:inline-flex; align-items:center; gap:6px; }}
    .swatch {{ width:16px; height:16px; border-radius:50%; display:inline-block; border:2px solid white; box-shadow:0 0 0 1px rgba(40,73,110,.18); }}
    .network-summary {{ display:grid; grid-template-columns:repeat(4,minmax(0,1fr)); gap:10px; margin:12px 0; }}
    .network-stat {{ border:1px solid var(--line); border-radius:8px; padding:10px; background:#fbfcfd; }}
    .network-stat span {{ display:block; color:var(--muted); font-size:11px; font-weight:800; text-transform:uppercase; letter-spacing:.04em; }}
    .network-stat strong {{ display:block; margin-top:3px; font-size:18px; color:var(--ink); overflow-wrap:anywhere; }}
    .network-wrap {{ position:relative; overflow:hidden; min-height:620px; }}
    .network-wrap.labels-hidden .node text,
    .network-wrap.labels-hidden .node .sim-label {{ display:none; }}
    .network-wrap svg {{ width:100%; height:620px; display:block; border:1px solid var(--line); border-radius:8px; background:#f7f9fb; }}
    .candidate-picker {{ display:flex; flex-wrap:wrap; gap:12px; align-items:center; margin:12px 0; }}
    .candidate-picker label {{ font-size:13px; font-weight:700; color:var(--blue); }}
    .candidate-picker select {{ min-width:min(520px,100%); padding:8px 10px; border:1px solid var(--line); border-radius:6px; font:inherit; }}
    .candidate-summary {{ margin:12px 0; padding:12px 14px; border:1px solid var(--line); border-radius:8px; background:#fbfcfd; }}
    .candidate-summary p {{ margin:4px 0; }}
    .extended-network-wrap {{ position:relative; overflow:auto; min-height:760px; max-height:1040px; background:#f7f9fb; border:1px solid var(--line); border-radius:8px; }}
    .extended-network-svg {{ width:100%; min-width:1480px; height:auto; display:block; }}
    .node text,.ext-node text {{ font-size:11px; fill:var(--ink); font-weight:700; pointer-events:none; paint-order:stroke; stroke:#fff; stroke-width:4px; stroke-linejoin:round; }}
    .node .sim-label {{ font-size:10px; fill:var(--muted); font-weight:800; }}
    .node circle,.ext-node circle {{ cursor:pointer; }}
    .node.hidden,.edge.hidden,tr.hidden,.ext-node.hidden,.ext-edge.hidden,tr.ext-hidden {{ display:none; }}
    .consensus-layout {{ display:grid; grid-template-columns:minmax(360px,.9fr) minmax(0,1.1fr); gap:18px; align-items:start; }}
    .consensus-box {{ background:#fff; border:1px solid var(--line); border-radius:8px; padding:12px; overflow:hidden; }}
    .consensus-svg {{ max-width:100%; height:auto; display:block; margin:0 auto; }}
    .consensus-scale {{ margin-top:10px; }}
    .consensus-scale-bar {{ height:16px; border-radius:999px; background:linear-gradient(90deg,#2C7BB6 0%,#FDAE61 50%,#D7191C 100%); border:1px solid var(--line); }}
    .consensus-scale-labels {{ display:grid; grid-template-columns:repeat(3,1fr); gap:8px; margin-top:5px; color:var(--muted); font-size:11px; font-weight:700; }}
    .consensus-scale-labels span:nth-child(2) {{ text-align:center; }}
    .consensus-scale-labels span:nth-child(3) {{ text-align:right; }}
    .consensus-atom text {{ fill:#fff; font-size:16px; font-weight:800; pointer-events:none; }}
    .consensus-atom .atom-confidence {{ fill:var(--ink); font-size:11px; font-weight:800; }}
    .decision-tree-box {{ background:#fff; border:1px solid var(--line); border-radius:8px; padding:12px; overflow:auto; }}
    .decision-tree-svg {{ width:100%; min-width:1040px; height:auto; display:block; }}
    .tree-title {{ font-size:20px; font-weight:800; fill:var(--ink); }}
    .tree-node-title {{ font-size:13px; font-weight:800; text-anchor:middle; }}
    .tree-node-metric {{ font-size:12px; font-weight:800; text-anchor:middle; fill:#394956; }}
    .tree-node-note,.tree-leaf-count,.tree-leaf-top {{ font-size:11px; text-anchor:middle; fill:#394956; }}
    .tree-leaf-count,.tree-leaf-top {{ fill:#fff; font-weight:800; }}
    .tree-edge-label {{ font-size:11px; font-weight:800; fill:#394956; paint-order:stroke; stroke:#f7f9fb; stroke-width:4px; }}
    .candidate-table td.smiles {{ font-family:"SFMono-Regular","Consolas",monospace; font-size:12px; max-width:420px; overflow-wrap:anywhere; }}
    .hover-tooltip {{ position:fixed; z-index:1000; max-width:320px; padding:10px 12px; border-radius:8px; background:rgba(23,33,43,.94); color:#fff; border:1px solid rgba(255,255,255,.12); box-shadow:0 10px 28px rgba(0,29,55,.24); font-size:12px; line-height:1.45; pointer-events:none; opacity:0; transform:translateY(4px); transition:opacity 120ms ease, transform 120ms ease; }}
    .hover-tooltip.visible {{ opacity:1; transform:translateY(0); }}
    .hover-tooltip strong {{ display:block; margin-bottom:2px; font-size:13px; }}
    footer {{ color:var(--muted); font-size:12px; padding:0 28px 34px; max-width:1320px; margin:0 auto; }}
    @media (max-width:980px) {{ .hero-grid,.call-grid,.metrics,.consensus-layout,.network-summary {{ grid-template-columns:1fr; }} header {{ padding:34px 24px 26px; }} h1 {{ font-size:38px; }} }}
  </style>
</head>
<body>
  <header>
    <div class="hero-grid">
      <div>
        <p class="eyebrow">X-Compound Elucidation Dashboard • {variant_label}</p>
        <h1>{html.escape(report_title)}</h1>
        <p class="subtitle">
          Structure-resolution report for <strong>{html.escape(control_id)}</strong>, formula <code>{html.escape(target['formula'])}</code>.
          The most likely candidate structure is shown below.
        </p>
        <p>
          <span class="badge ok">top-1 candidate</span>
          {descriptor_badge_html}
          {confmod.confidence_badge_html(confidence)}
          {confmod.msi_badge_html(confidence)}
          {confmod.identity_badge_html(confidence)}
          <span class="badge ok">{html.escape(top_similarity_badge)}</span>
          {pathway_badge_html}
        </p>
      </div>
      <div><div class="metrics">{header_metrics_html}</div></div>
    </div>
  </header>
  <main>
    <section class="section">
      {heading_with_info("h2", "Top Candidate Compound")}
      <p class="note">
        On average, 80% of compounds identified in Metabolon data sets meet Level 1 criteria through fingerprint matching to an authentic reference standard analyzed with the same experimental methods. The remaining compounds are identified through Metabolon’s structural elucidation workflow; the evidence stack below follows that dashboard structure.
      </p>
      <div class="call-grid">
        <div class="structure-card">
          {heading_with_info("h3", "Most Likely Structure", include_info=False)}
          {structure_panel_html}
          <p><strong>{html.escape(target['name'])}</strong></p>
          <p class="note">This is the top ranked candidate structure for this feature. Note that exact stereochemistry is not known from available evidence.</p>
        </div>
        <div class="structure-card">
          {heading_with_info("h3", "Atom-position confidence mapped to top candidate", include_info=False)}
          <div class="consensus-box">{consensus_svg}</div>
          <div class="consensus-scale" aria-label="Consensus heatmap scale">
            <div class="consensus-scale-bar"></div>
            <div class="consensus-scale-labels">
              <span>0.0 low support</span>
              <span>0.5 mixed support</span>
              <span>1.0 high support</span>
            </div>
          </div>
          <p class="note">Sidecar: <a href="data/{slug}_consensus_atom_support.csv">{slug}_consensus_atom_support.csv</a>.</p>
        </div>
      </div>
      <div class="metrics">
        {confmod.confidence_metric_card_html(confidence)}
        {confmod.msi_metric_card_html(confidence)}
        {confmod.identity_metric_card_html(confidence)}
      </div>
      <div class="structure-card" style="margin-top:18px;">
        {heading_with_info("h3", "Identification Evidence", include_info=False)}
        <div class="metrics">
          <div class="metric-card copyable-card" tabindex="0" data-copy-value="{evidence['parent_mz']:.4f}"><span class="copy-hint">Click to copy</span><div class="metric-label">Observed parent m/z</div><div class="metric-value">{evidence['parent_mz']:.4f}</div><div class="metric-note">{html.escape(control_id)} feature</div></div>
          <div class="metric-card copyable-card" tabindex="0" data-copy-value="{target['exact_mass']:.6f}"><span class="copy-hint">Click to copy</span><div class="metric-label">Exact mass</div><div class="metric-value">{target['exact_mass']:.6f}</div><div class="metric-note">Candidate neutral mass</div></div>
          <div class="metric-card copyable-card" tabindex="0" data-copy-value="{html.escape(target['smiles'], quote=True)}"><span class="copy-hint">Click to copy</span><div class="metric-label">Top Candidate SMILES</div><div class="metric-value"><code>{html.escape(target['smiles'])}</code></div><div class="metric-note"><a href="{pubchem_link}">Other DB Links (PubChem of Top Candidate)</a></div></div>
          <div class="metric-card copyable-card" tabindex="0" data-copy-value="{html.escape(target['inchikey'], quote=True)}"><span class="copy-hint">Click to copy</span><div class="metric-label">InChIKey</div><div class="metric-value">{html.escape(target['inchikey'])}</div><div class="metric-note">From candidate SMILES</div></div>
          <div class="metric-card"><div class="metric-label">{synonyms_title}</div><div class="metric-list">{synonym_card_html}</div></div>
        </div>
      </div>
    </section>

    <section class="section">
      {heading_with_info("h2", params_title)}
      {parameter_cards_html}
      {confmod.scored_parameters_section_html(confidence, public=public_mode)}
    </section>

    {top10_section}

    <section class="section">
      {heading_with_info("h2", "Top 10 Metabolon Library Matches")}
      <p class="note">Structural elucidation is most effective when the unknown compound has features that resemble library metabolites. Scores above 0.5 typically support higher-confidence candidate structures.</p>
      <div class="table-wrap">
        <table>
          {library_header_html}
          <tbody>{library_row_html}</tbody>
        </table>
      </div>
    </section>

    <section class="section">
      {heading_with_info("h2", fragment_title)}
      <p class="note">Metabolon’s elucidation pipeline uses diagnostic spectral information including adducts and biotransformations to help with elucidation and confidence in the identification.</p>
      <div class="table-wrap">
        <table>
          <thead><tr><th>Fragment</th><th>Match Score</th></tr></thead>
          <tbody>{substructure_rows}</tbody>
        </table>
      </div>
    </section>

    <section class="section">
      {heading_with_info("h2", "Mass Match Scores")}
      <p class="note">Mass matching is a key parameter for structural elucidation. With high-resolution accurate-mass instruments, the expected mass should closely match the experimental value.</p>
      <div class="table-wrap">
        <table>
          <thead><tr><th>Top Candidate Predicted Mass</th><th>Observed m/z</th><th>Ionization Mode / ion</th><th>Mass Error (Predicted - Observed)</th></tr></thead>
          <tbody>{mass_rows}</tbody>
        </table>
      </div>
    </section>

    <section class="section">
      {heading_with_info("h2", ms_fit_title)}
      <p class="note">Metabolon’s elucidation pipeline uses mass fragment (MS2) matching as a key parameter for structural elucidation. The tables below summarize observed fragment counts and how well each candidate explains the spectrum.</p>
      {fragment_summary_html}
    </section>

    <section class="section">
      {heading_with_info("h2", "Top Correlating Metabolites")}
      <p class="note">{correlation_note}</p>
      <div class="table-wrap">
        <table>
          <thead><tr><th>ID</th><th>Compound</th><th>Correlation (R)</th><th>Pathway</th></tr></thead>
          <tbody>{correlating_rows}</tbody>
        </table>
      </div>
    </section>

    {similarity_section}
{internal_sections}
  </main>
  <footer>
    Generated by the structural elucidation report builder from Spec2Mol outputs plus the local metabolite synonym DB.
  </footer>
  <div id="hover-tooltip" class="hover-tooltip" aria-hidden="true"></div>
  <script>
    const buttons = document.querySelectorAll('button.filter');
    function isLactoylAlanineContext(el) {{
      const text = ((el.dataset.name || '') + ' ' + (el.dataset.pathway || '')).toLowerCase();
      return /lactoyl|alanine|alanyl|acetylalanine|propionylalanine|carbamoylalanine/.test(text);
    }}
    function showMainNode(mode, el) {{
      const confidence = el.dataset.confidence || '';
      const similarity = Number(el.dataset.similarity || '0');
      const kind = el.dataset.kind || '';
      if (confidence === 'target') return true;
      if (mode === 'all') return true;
      if (mode === 'known') return kind === 'known';
      if (mode === 'unknown') return kind === 'unknown';
      if (mode === 'ge060') return similarity >= 0.60;
      if (mode === 'ge050') return similarity >= 0.50;
      return true;
    }}
    function applyMainFilter(mode) {{
      const visibleNodeIds = new Set();
      document.querySelectorAll('.node').forEach(el => {{
        const visible = showMainNode(mode, el);
        el.classList.toggle('hidden', !visible);
        if (visible) visibleNodeIds.add(el.dataset.nodeId || '');
      }});
      document.querySelectorAll('.edge').forEach(el => {{
        const visible = visibleNodeIds.has(el.dataset.sourceId || '') && visibleNodeIds.has(el.dataset.targetId || '');
        el.classList.toggle('hidden', !visible);
      }});
      document.querySelectorAll('tbody tr[data-confidence]').forEach(el => {{
        el.classList.toggle('hidden', !showMainNode(mode, el));
      }});
    }}
    buttons.forEach(btn => btn.addEventListener('click', () => {{
      if (btn.id === 'network-label-toggle') return;
      buttons.forEach(b => b.classList.remove('active'));
      btn.classList.add('active');
      applyMainFilter(btn.dataset.filter);
    }}));

    const networkWrap = document.querySelector('.network-wrap');
    const networkLabelToggle = document.getElementById('network-label-toggle');
    if (networkWrap && networkLabelToggle) {{
      networkLabelToggle.addEventListener('click', () => {{
        const hidden = networkWrap.classList.toggle('labels-hidden');
        networkLabelToggle.textContent = hidden ? 'Show labels' : 'Hide labels';
        networkLabelToggle.setAttribute('aria-pressed', hidden ? 'true' : 'false');
      }});
    }}

    const candidateConsensusDataEl = document.getElementById('candidate-consensus-data');
    const candidateConsensusSelect = document.getElementById('candidate-consensus-select');
    const candidateConsensusPanel = document.getElementById('candidate-consensus-panel');
    const candidateConsensusSummary = document.getElementById('candidate-consensus-summary');
    if (candidateConsensusDataEl && candidateConsensusSelect && candidateConsensusPanel && candidateConsensusSummary) {{
      const candidateConsensusData = JSON.parse(candidateConsensusDataEl.textContent || '{{}}');
      const consensusMaps = candidateConsensusData.maps || {{}};
      const consensusSummaries = candidateConsensusData.summaries || {{}};
      function renderCandidateSummary(rank) {{
        const summary = consensusSummaries[String(rank)] || {{}};
        const bits = [];
        if (summary.rank) bits.push('<p><strong>Rank ' + summary.rank + '</strong> · ' + (summary.predicted_formula || '') + '</p>');
        const fmt = (v, nd) => (v !== undefined && v !== '' && Number.isFinite(Number(v))) ? Number(v).toFixed(nd) : null;
        const seg = [];
        const ev = fmt(summary.evidence_confidence, 3);
        const mm = fmt(summary.mass_match_score, 3);
        const rs = fmt(summary.candidate_report_score, 4);
        const me = fmt(summary.abs_mass_error_da, 6);
        if (ev !== null) seg.push('Evidence confidence: <strong>' + ev + '</strong>');
        if (mm !== null) seg.push('Mass match: <strong>' + mm + '</strong>');
        if (rs !== null) seg.push('Rank score (relative): <strong>' + rs + '</strong>');
        if (me !== null) seg.push('Mass error (Da): <strong>' + me + '</strong>');
        if (seg.length) bits.push('<p>' + seg.join(' · ') + '</p>');
        if (summary.predicted_smiles) bits.push('<p><code>' + summary.predicted_smiles + '</code></p>');
        candidateConsensusSummary.innerHTML = bits.join('') || '<p class="note">No candidate summary available.</p>';
      }}
      candidateConsensusSelect.addEventListener('change', () => {{
        const rank = candidateConsensusSelect.value;
        candidateConsensusPanel.innerHTML = consensusMaps[String(rank)] || '<p class="note">Consensus map unavailable for this candidate.</p>';
        renderCandidateSummary(rank);
      }});
    }}

    const extButtons = document.querySelectorAll('button.ext-filter');
    function showExtendedNode(mode, el) {{
      const kind = el.dataset.kind || '';
      const similarity = Number(el.dataset.similarity || '0');
      const text = ((el.dataset.name || '') + ' ' + (el.dataset.pathway || '')).toLowerCase();
      if (kind === 'target') return true;
      if (mode === 'all') return true;
      if (mode === 'unknown') return kind === 'unknown';
      if (mode === 'known') return kind === 'known';
      if (mode === 'ge060') return similarity >= 0.60;
      if (mode === 'lactoyl') return /lactoyl|alanine|glycine|amino acid/.test(text);
      return true;
    }}
    function applyExtendedFilter(mode) {{
      const visibleNodeIds = new Set();
      document.querySelectorAll('.ext-node').forEach(el => {{
        const visible = showExtendedNode(mode, el);
        el.classList.toggle('hidden', !visible);
        if (visible) visibleNodeIds.add(el.dataset.nodeId || '');
      }});
      document.querySelectorAll('.ext-edge').forEach(el => {{
        const visible = visibleNodeIds.has(el.dataset.sourceId || '') && visibleNodeIds.has(el.dataset.targetId || '');
        el.classList.toggle('hidden', !visible);
      }});
    }}
    extButtons.forEach(btn => btn.addEventListener('click', () => {{
      extButtons.forEach(b => b.classList.remove('active'));
      btn.classList.add('active');
      applyExtendedFilter(btn.dataset.filter);
    }}));

    const tooltip = document.getElementById('hover-tooltip');
    function tooltipHtml(el) {{
      const label = (el.dataset.hoverLabel || el.dataset.name || '').replace(/\\n/g, ' ').trim();
      const similarity = el.dataset.similarity || '';
      const pathway = (el.dataset.pathway || '').trim();
      const bits = [];
      if (similarity) bits.push('Similarity ' + similarity);
      if (pathway) bits.push(pathway);
      return '<strong>' + label + '</strong>' + bits.map(x => '<div>' + x + '</div>').join('');
    }}
    function showTooltip(el, evt) {{
      tooltip.innerHTML = tooltipHtml(el);
      tooltip.classList.add('visible');
      tooltip.setAttribute('aria-hidden', 'false');
      if (evt && Number.isFinite(evt.clientX) && Number.isFinite(evt.clientY)) {{
        moveTooltip(evt);
        return;
      }}
      const rect = el.getBoundingClientRect();
      moveTooltip({{ clientX: rect.left + rect.width / 2, clientY: rect.top + rect.height / 2 }});
    }}
    function moveTooltip(evt) {{
      if (!tooltip.classList.contains('visible')) return;
      const offset = 16;
      const x = Math.min(window.innerWidth - tooltip.offsetWidth - 12, evt.clientX + offset);
      const y = Math.min(window.innerHeight - tooltip.offsetHeight - 12, evt.clientY + offset);
      tooltip.style.left = x + 'px';
      tooltip.style.top = y + 'px';
    }}
    function hideTooltip() {{
      tooltip.classList.remove('visible');
      tooltip.setAttribute('aria-hidden', 'true');
    }}
    async function copyValue(text, el) {{
      if (!text) return;
      try {{
        await navigator.clipboard.writeText(text);
        const hint = el.querySelector('.copy-hint');
        if (hint) {{
          const original = hint.textContent;
          hint.textContent = 'Copied';
          setTimeout(() => {{
            hint.textContent = original;
          }}, 1200);
        }}
      }} catch (err) {{
        console.warn('Clipboard write failed', err);
      }}
    }}
    document.querySelectorAll('.node, .ext-node').forEach(el => {{
      el.addEventListener('mouseenter', evt => showTooltip(el, evt));
      el.addEventListener('mousemove', moveTooltip);
      el.addEventListener('mouseleave', hideTooltip);
      el.addEventListener('focus', evt => showTooltip(el, evt));
      el.addEventListener('blur', hideTooltip);
    }});
    document.querySelectorAll('.copyable-card').forEach(el => {{
      el.addEventListener('click', () => copyValue(el.dataset.copyValue || '', el));
      el.addEventListener('keydown', evt => {{
        if (evt.key === 'Enter' || evt.key === ' ') {{
          evt.preventDefault();
          copyValue(el.dataset.copyValue || '', el);
        }}
      }});
    }});
  </script>
</body>
</html>
"""


# ── Report assembly ────────────────────────────────────────────────────────────

def control_slug(control_id: str) -> str:
    """Filesystem-safe slug from a control id (e.g. ``X-25003`` -> ``x-25003``)."""

    return re.sub(r"[^a-z0-9]+", "-", control_id.lower()).strip("-") or "compound"


def predictions_to_candidate_pool(predictions: pd.DataFrame, control_id: str) -> pd.DataFrame:
    """Filter and normalize rerank/prediction rows for one ``control_id``."""

    work = predictions.copy()
    id_col = "control_id" if "control_id" in work.columns else "x_compound_id"
    work = work[work[id_col].astype(str).eq(control_id)].copy()
    if work.empty:
        return work
    smiles_col = next(
        (c for c in ("predicted_smiles", "smiles", "candidate_smiles", "canonical_2d_smiles") if c in work.columns),
        None,
    )
    if smiles_col is None:
        return pd.DataFrame()
    rank_col = "rank" if "rank" in work.columns else "_rank" if "_rank" in work.columns else None
    if rank_col:
        work = work.sort_values(rank_col, kind="mergesort")
        work["rank"] = pd.to_numeric(work[rank_col], errors="coerce")
    else:
        work["rank"] = range(1, len(work) + 1)
    work["control_id"] = control_id
    work["predicted_smiles"] = work[smiles_col].astype(str)
    if "predicted_formula" not in work.columns and "candidate_formula" in work.columns:
        work["predicted_formula"] = work["candidate_formula"]
    if "candidate_name" not in work.columns and "base_name" in work.columns:
        work["candidate_name"] = work["base_name"]
    score_col = next(
        (
            c
            for c in ("unified_final_score", "rerank_score_calibrated", "rerank_score", "score_rank")
            if c in work.columns
        ),
        None,
    )
    if score_col:
        work["candidate_report_score"] = pd.to_numeric(work[score_col], errors="coerce")
    return work


def lookup_known_x_row(known_x: pd.DataFrame, control_id: str) -> dict[str, Any]:
    """Return the first ``known_X`` row matching ``control_id`` in Unknown column."""

    unk = known_x.get("Unknown", pd.Series("", index=known_x.index)).astype(str)
    rows = known_x[unk.eq(control_id)]
    return rows.iloc[0].to_dict() if not rows.empty else {}


def resolve_chemical_id(db_path: Path, known_x_row: dict[str, Any], target_name: str) -> str:
    """Resolve metabolite DB chemical_id from known_X or canonical name lookup."""

    if known_x_row.get("KNOWN"):
        return str(known_x_row["KNOWN"]).strip()
    if not db_path.exists():
        return ""
    df = query_df(
        db_path,
        "SELECT chemical_id FROM metabolites WHERE lower(canonical_name) = lower(?) LIMIT 1",
        (target_name,),
    )
    if not df.empty:
        return str(df.iloc[0]["chemical_id"])
    return ""


def is_real_chemical_id(db_path: Path, chemical_id: str) -> bool:
    """True when ``chemical_id`` is an actual ``metabolites.chemical_id`` row."""

    if not db_path.exists() or not chemical_id:
        return False
    try:
        df = query_df(db_path, "SELECT 1 FROM metabolites WHERE chemical_id = ? LIMIT 1", (chemical_id,))
    except Exception:
        return False
    return not df.empty


def resolve_db_identity_by_inchikey(db_path: Path, inchikey: str) -> tuple[str, str] | None:
    """Match a target InChIKey skeleton to a DB metabolite.

    Returns ``(chemical_id, canonical_name)`` when the 2D structure (first
    InChIKey block) matches a metabolite row, else ``None``. This lets reports
    populate the library-match/synonyms sections and show a real display name
    even when ``known_X`` only carries a free-text ``KNOWN`` label.
    """

    if not db_path.exists() or not inchikey:
        return None
    skeleton = inchikey.split("-")[0]
    if len(skeleton) < 14:
        return None
    try:
        df = query_df(
            db_path,
            "SELECT chemical_id, canonical_name FROM metabolites "
            "WHERE inchikeys LIKE ? AND canonical_name IS NOT NULL AND trim(canonical_name) <> '' LIMIT 1",
            (f"%{skeleton}%",),
        )
    except Exception:
        return None
    if df.empty:
        return None
    row = df.iloc[0]
    return str(row["chemical_id"]), str(row["canonical_name"])


def load_real_correlation_rows(
    control_id: str,
    chemical_id: str,
    target_name: str,
) -> tuple[list[dict[str, Any]], str | None]:
    """Attempt WGCNA module-anchor lookup for real abundance correlation rows."""

    anchor_names = {control_id.lower(), target_name.lower()}
    if target_name:
        anchor_names.add(target_name.strip().lower())
    for cohort, spec in wgcna.COHORT_INPUTS.items():
        modules = wgcna.load_module_assignments(cohort)
        if modules.empty:
            continue
        matches = modules[
            modules["plot_name"].astype(str).str.strip().str.lower().isin(anchor_names)
            | modules["plot_name"].astype(str).eq(control_id)
        ]
        if chemical_id:
            matches = matches | modules["chem_id"].astype(str).eq(str(chemical_id))
        if not matches.empty:
            break
    blocker = (
        f"Real sample-abundance correlation rows for {control_id} are not available in the configured WGCNA exports."
    )
    return [], blocker


def build_placeholder_correlation_rows(
    extended_unknown_rows: list[dict[str, Any]],
    *,
    exclude_control_id: str,
    limit: int = 5,
) -> list[dict[str, Any]]:
    """Build proxy correlation rows from nearest verified X-structure neighbors."""

    proxy_rows: list[dict[str, Any]] = []
    for row in extended_unknown_rows:
        x_id = str(row.get("x_compound_id", "") or "")
        if not x_id or x_id == exclude_control_id:
            continue
        try:
            proxy_score = float(row.get("similarity_to_target", 0.0) or 0.0)
        except (TypeError, ValueError):
            proxy_score = 0.0
        proxy_rows.append(
            {
                "id": x_id,
                "compound": str(row.get("known_name", "") or x_id),
                "correlation_r": f"{proxy_score:.3f}",
                "pathway": "Verified X structure proxy",
            }
        )
        if len(proxy_rows) >= limit:
            break
    return proxy_rows


def render_guides_html(control_id: str, guide_targets: list[tuple[str, str]]) -> str:
    """Render placeholder guide cards linked from public report info buttons."""

    guide_blocks = (
        "".join(
            f"<section class='guide-card' id='guide-{slug}'><h2>{html.escape(label)}</h2><p>{html.escape(label)} guide will go here.</p></section>"
            for slug, label in guide_targets
        )
        or "<section class='guide-card'><h2>No guides yet</h2><p>Guide content will go here.</p></section>"
    )
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>{html.escape(control_id)} Guide Placeholders</title>
  <style>
    :root {{
      --navy:#001D37; --blue:#28496E; --paper:#FFFFFF; --wash:#F1FAFB; --ink:#15304B; --muted:#5C6B7A; --line:rgba(40,73,110,.18);
      --shadow:0 24px 58px rgba(0,29,55,.10);
    }}
    * {{ box-sizing:border-box; }}
    body {{ margin:0; font-family:"Montserrat","Avenir Next","Segoe UI",Arial,sans-serif; color:var(--ink); background:linear-gradient(180deg,#fff 0,#f7fbfc 100%); }}
    main {{ max-width:1080px; margin:0 auto; padding:32px 24px 56px; }}
    h1 {{ margin:0 0 8px; color:var(--navy); }}
    .lede {{ color:var(--muted); margin:0 0 22px; }}
    .guide-grid {{ display:grid; gap:14px; }}
    .guide-card {{ background:var(--paper); border:1px dashed var(--line); border-radius:10px; padding:18px; box-shadow:var(--shadow); scroll-margin-top:18px; }}
    .guide-card h2 {{ margin:0 0 6px; font-size:20px; color:var(--navy); }}
    .guide-card p {{ margin:0; color:var(--muted); }}
  </style>
</head>
<body>
  <main>
    <h1>{html.escape(control_id)} guides</h1>
    <p class="lede">Placeholder guide cards for dashboard info buttons.</p>
    <div class="guide-grid">{guide_blocks}</div>
  </main>
</body>
</html>"""


def build_report(args: argparse.Namespace) -> dict[str, Any]:
    """Orchestrate data load, scoring, network build, and HTML/CSV artifact emission."""

    repo_root = REPO_ROOT
    control_id = str(args.control_id).strip()
    slug = control_slug(control_id)
    out_dir = resolve(args.out_dir, repo_root)
    assets = out_dir / "assets"
    data_dir = out_dir / "data"
    out_dir.mkdir(parents=True, exist_ok=True)
    assets.mkdir(parents=True, exist_ok=True)
    data_dir.mkdir(parents=True, exist_ok=True)

    db_path = resolve(args.db, repo_root)
    predictions_path = resolve(args.predictions, repo_root)
    predictions = pd.read_csv(predictions_path, low_memory=False)
    if args.candidate_pool and Path(resolve(args.candidate_pool, repo_root)).exists():
        candidate_pool = pd.read_csv(resolve(args.candidate_pool, repo_root), low_memory=False)
        top = candidate_pool[
            candidate_pool.get("control_id", pd.Series("", index=candidate_pool.index)).astype(str).eq(control_id)
        ].copy()
    else:
        top = predictions_to_candidate_pool(predictions, control_id)
    if top.empty:
        raise SystemExit(f"No prediction rows found for control_id={control_id} in {predictions_path}")

    known_x_path = resolve(args.known_x, repo_root)
    known_x = pd.read_csv(known_x_path, low_memory=False)
    known_x_all_path = resolve(args.known_x_all, repo_root)
    known_x_all = pd.read_csv(known_x_all_path, low_memory=False) if known_x_all_path.exists() else pd.DataFrame()

    peak_csv = resolve(args.peak_csv, repo_root) if args.peak_csv else None
    peaks = pd.read_csv(peak_csv, low_memory=False) if peak_csv and peak_csv.exists() else pd.DataFrame()

    rank1 = top.sort_values("rank", kind="mergesort").iloc[0]
    smiles_col = "predicted_smiles" if "predicted_smiles" in rank1.index else "smiles"
    target_smiles = str(args.target_smiles or rank1.get(smiles_col) or "").strip()
    if not target_smiles:
        raise SystemExit(f"Could not resolve target SMILES for {control_id}")

    target_name = str(args.target_name or rank1.get("candidate_name") or rank1.get("base_name") or control_id).strip()
    if target_name.lower() in {"", "nan", "none"}:
        target_name = control_id
    kx_row = lookup_known_x_row(known_x, control_id)
    chemical_id = str(args.chemical_id or resolve_chemical_id(db_path, kx_row, target_name) or "").strip()

    target_inchikey = inchikey_from_smiles(target_smiles)
    db_identity = resolve_db_identity_by_inchikey(db_path, target_inchikey)
    name_resolved = bool(db_identity)
    if not args.chemical_id and db_identity:
        resolved_chem_id, resolved_name = db_identity
        if not is_real_chemical_id(db_path, chemical_id):
            chemical_id = resolved_chem_id
        if target_name == control_id and resolved_name:
            target_name = resolved_name

    target_formula = formula_from_smiles(target_smiles)
    pd.DataFrame([kx_row]).to_csv(data_dir / f"{slug}_known_x_row.csv", index=False)
    tentative = tnmod.resolve(
        target_smiles,
        target_formula,
        rank1,
        tnmod.SqliteNameLookup(db_path),
        name_resolved=name_resolved,
        bundle_dir=data_dir.parent,
        extended_ctx=default_extended_context(reports_root=data_dir.parent.parent),
        inchikey=target_inchikey,
        control_id=control_id,
    )
    if tentative:
        target_name = tentative.label
    elif not name_resolved:
        target_name = "De novo structure"
    target = {
        "control_id": control_id,
        "name": target_name,
        "smiles": target_smiles,
        "formula": target_formula,
        "exact_mass": exact_mass_from_smiles(target_smiles),
        "inchikey": target_inchikey,
    }
    if tentative:
        target["tentative_descriptor"] = tentative.as_dict()
    target_structure_svg = render_smiles_svg(target_smiles)
    asset_base = f"{slug}_structure"
    if target_structure_svg:
        (assets / f"{asset_base}.svg").write_text(target_structure_svg, encoding="utf-8")
        structure_panel_html = (
            "<div class='structure-svg-frame'>"
            f"{target_structure_svg}"
            "<p class='structure-note'>2D structure rendered from SMILES.</p>"
            "</div>"
        )
    else:
        structure_image(target_smiles, assets / f"{asset_base}.png", f"{control_id}: {target_name}")
        structure_panel_html = (
            f"<img class='structure' src='assets/{asset_base}.png' alt='{html.escape(target_name)} structure'>"
        )

    top = top.sort_values("rank", kind="mergesort")
    if not peaks.empty:
        feature_col = next((c for c in ("feature_key", "FEATURE_KEY") if c in peaks.columns), "")
        peak_mask = pd.Series(False, index=peaks.index)
        if feature_col:
            peak_mask = peak_mask | peaks[feature_col].astype(str).str.contains(
                re.escape(control_id), regex=True, na=False
            )
        if peak_mask.any():
            peaks = peaks[peak_mask].copy()

    top.head(int(args.candidate_tab_n)).copy().to_csv(data_dir / f"{slug}_top_candidates.csv", index=False)
    candidate_score_input = score_smiles_spectrum_alignment(
        top.head(int(args.candidate_tab_n)).copy(),
        peaks,
        tolerance_da=float(args.spectrum_tolerance_da),
    )
    candidate_score_df = compact_candidate_scores(candidate_score_input)
    candidate_score_df.to_csv(data_dir / f"{slug}_candidate_score_summary.csv", index=False)
    top10_candidates = candidate_score_df.head(10).copy()
    candidate_consensus_maps, candidate_consensus_summaries, default_consensus_rank = build_candidate_consensus_maps(
        top10_candidates
    )
    (data_dir / f"{slug}_candidate_consensus_maps.json").write_text(
        json.dumps(
            {
                "maps": {str(k): v for k, v in candidate_consensus_maps.items()},
                "summaries": {str(k): v for k, v in candidate_consensus_summaries.items()},
                "default_rank": default_consensus_rank,
            },
            indent=2,
        )
    )
    consensus_rows, consensus_svg = consensus_atom_support(target_smiles, top.head(int(args.consensus_top_n)).copy())
    write_csv(data_dir / f"{slug}_consensus_atom_support.csv", consensus_rows)

    db_detail = {}
    syn_df = pd.DataFrame()
    if db_path.exists() and chemical_id:
        db_detail_df = query_df(
            db_path,
            """
            SELECT chemical_id, canonical_name, smiles, inchikeys, super_pathway, pathway,
                   metabolon_super, metabolon_sub, resolution_bucket, resolution_status,
                   resolution_method, resolution_sources, resolution_query,
                   resolution_match_name, resolution_score, resolution_candidate_count
            FROM metabolites
            WHERE chemical_id = ?
            """,
            (chemical_id,),
        )
        db_detail = db_detail_df.iloc[0].to_dict() if not db_detail_df.empty else {}
        db_detail_df.to_csv(data_dir / "metabolite_db_detail.csv", index=False)
        syn_df = query_df(
            db_path,
            """
            SELECT entry_value AS synonym, normalized_value AS normalized, sources
            FROM synonym_index
            WHERE chemical_id = ? AND entry_type = 'synonym'
            ORDER BY lower(entry_value), entry_value
            """,
            (chemical_id,),
        )
        syn_df.to_csv(data_dir / "metabolite_db_synonyms.csv", index=False)

    pathway_info = resolve_pathway_display(db_detail, rank1)
    pd.DataFrame(
        [
            {
                "control_id": control_id,
                "predicted_pathway": pathway_info["super"],
                "predicted_subpathway": pathway_info["sub"],
                "pathway_source": pathway_info["source"],
                "pathway_confidence_tier": pathway_info["tier"],
                "pathway_best_reference_similarity": pathway_info["similarity"],
                "pathway_assignment_status": pathway_info["status"],
            }
        ]
    ).to_csv(data_dir / f"{slug}_pathway_prediction.csv", index=False)

    correlation_rows, correlation_blocker = load_real_correlation_rows(control_id, chemical_id, target_name)

    all_structured = pd.DataFrame()
    if db_path.exists():
        all_structured = query_df(
            db_path,
            """
            SELECT chemical_id, canonical_name, smiles, inchikeys, super_pathway, pathway,
                   metabolon_super, metabolon_sub, resolution_status
            FROM metabolites
            WHERE smiles IS NOT NULL AND trim(smiles) <> ''
            """,
        )

    extended = build_extended_similarity_network(
        target_smiles=target_smiles,
        all_structured=all_structured,
        known_x_all=known_x_all,
        max_known_render=int(args.extended_known_render_limit),
        known_min_similarity=float(args.extended_known_min_similarity),
        x_min_similarity=float(args.extended_x_min_similarity),
        target_edge_min=float(args.extended_target_edge_min),
        pair_edge_min=float(args.extended_pair_edge_min),
        max_edges_per_node=int(args.extended_max_edges_per_node),
        control_id=control_id,
        target_name=target_name,
    )
    if not extended["known_all"].empty:
        extended["known_all"].to_csv(data_dir / f"all_known_similarity_to_{slug}.csv", index=False)
    if not extended["x_all"].empty:
        extended["x_all"].to_csv(data_dir / f"all_unknown_similarity_to_{slug}.csv", index=False)
    write_csv(data_dir / "extended_similarity_network_nodes.csv", extended["render_nodes"])
    write_csv(data_dir / "extended_similarity_network_edges.csv", extended["render_edges"])
    (data_dir / "extended_similarity_network.json").write_text(
        json.dumps({"nodes": extended["render_nodes"], "edges": extended["render_edges"]}, indent=2)
    )

    rows: list[dict[str, Any]] = []
    _, target_fp = mol_fp(target_smiles)
    for _, row in all_structured.iterrows():
        smiles = str(row.get("smiles") or "").split("|")[0].strip()
        mol, fp = mol_fp(smiles)
        if mol is None or fp is None or target_fp is None:
            continue
        sim = float(TanimotoSimilarity(target_fp, fp))
        if sim < float(args.min_similarity):
            continue
        item = row.to_dict()
        item["smiles"] = smiles
        item["similarity_to_target"] = sim
        rows.append(item)
    sim_df = pd.DataFrame(rows)
    if not sim_df.empty:
        sim_df = sim_df.sort_values("similarity_to_target", ascending=False, kind="mergesort").head(
            int(args.max_neighbors)
        )
        sim_df["confidence"] = sim_df.apply(confidence_bucket, axis=1)
    neighbors = sim_df.to_dict(orient="records") if not sim_df.empty else []
    write_csv(data_dir / "similarity_neighbors.csv", neighbors)

    structural_possibility_rows: list[dict[str, Any]] = []
    structural_possibility_tree_svg = ""
    if args.hypotheses and Path(resolve(args.hypotheses, repo_root)).exists():
        structural_possibility_rows = x25003_structural_possibility_summary(
            candidate_scores=candidate_score_df,
            neighbors=neighbors,
            target_smiles=target_smiles,
            target_formula=target["formula"],
        )
        structural_possibility_tree_svg = x25003_structural_possibility_tree_svg(structural_possibility_rows)
        write_csv(data_dir / f"{slug}_structural_possibility_summary.csv", structural_possibility_rows)

    for row in neighbors[:12]:
        safe = "".join(ch if ch.isalnum() else "-" for ch in str(row.get("canonical_name", "")).lower()).strip("-")[:64]
        structure_image(
            str(row["smiles"]), assets / f"neighbor-{row['chemical_id']}-{safe}.png", str(row.get("canonical_name", ""))
        )

    nodes = extended["render_nodes"]
    edges = extended["render_edges"]
    write_csv(data_dir / "similarity_network_nodes.csv", nodes)
    write_csv(data_dir / "similarity_network_edges.csv", edges)
    (data_dir / "similarity_network.json").write_text(json.dumps({"nodes": nodes, "edges": edges}, indent=2))

    parent_mz = float(pd.to_numeric(rank1.get("parent_mz"), errors="coerce") or 0.0)
    evidence = {
        "nearest_reference_tanimoto": float(
            pd.to_numeric(pathway_info.get("similarity"), errors="coerce") or 0.0
        ),
        "nearest_reference_source": str(pathway_info.get("source") or ""),
        "rerank_pool_score": float(
            pd.to_numeric(rank1.get("unified_final_score"), errors="coerce") or 0.0
        ),
        "parent_mz": parent_mz,
        "mass_error_da": float(target["exact_mass"] - parent_mz) if parent_mz else 0.0,
        "lactoyl_hypothesis_score": "",
        "diagnostic_bin_score": "",
        "matched_diagnostic_bins_json": "",
    }
    write_csv(data_dir / f"{slug}_evidence_summary.csv", [{"field": k, "value": v} for k, v in evidence.items()])
    (data_dir / f"{slug}_evidence_summary.json").write_text(
        json.dumps({"target": target, "evidence": evidence, "db_detail": db_detail}, indent=2)
    )

    extended_x_all = extended["x_all"]
    if not extended_x_all.empty and "x_compound_id" in extended_x_all.columns:
        extended_unknown_rows = (
            extended_x_all[extended_x_all["x_compound_id"].astype(str) != control_id]
            .head(20)
            .fillna("")
            .to_dict(orient="records")
        )
    else:
        extended_unknown_rows = []
    if not correlation_rows:
        correlation_rows = build_placeholder_correlation_rows(
            extended_unknown_rows, exclude_control_id=control_id, limit=5
        )
        correlation_blocker = (
            "Proxy correlation rows use nearest verified X-structure neighbors from the expanded similarity network."
        )
    write_csv(data_dir / f"{slug}_top_correlating_metabolites.csv", correlation_rows)
    (data_dir / f"{slug}_top_correlating_metabolites_status.json").write_text(
        json.dumps({"n_rows": len(correlation_rows), "blocker": correlation_blocker}, indent=2)
    )
    extended_known_rows = (
        extended["known_all"].head(20).fillna("").to_dict(orient="records") if not extended["known_all"].empty else []
    )
    known_render = extended.get("known_render", pd.DataFrame())
    x_render = extended.get("x_render", pd.DataFrame())

    # Hybrid identification confidence (shared model; see reports/confidence.py).
    _css_records = candidate_score_df.fillna("").to_dict(orient="records")

    def _n(value: Any) -> float | None:
        v = pd.to_numeric(value, errors="coerce")
        return None if pd.isna(v) else float(v)

    _ms2_candidates = [
        _n(rank1.get("predicted_spectrum_score_v2")),
        _n(rank1.get("forward_spectrum_score")),
        _n(rank1.get("ms2_match_score")),
    ]
    _ms2_candidates = [c for c in _ms2_candidates if c is not None]
    _ms2_score = max(_ms2_candidates) if _ms2_candidates else None
    _observed_mz = _n(rank1.get("parent_mz"))
    _conf_fields = confmod.fields_from_parts(
        _css_records[0] if _css_records else {},
        _css_records[1] if len(_css_records) > 1 else None,
        {
            "pathway_source": pathway_info.get("source", ""),
            "pathway_best_reference_similarity": pathway_info.get("similarity"),
            "pathway_vote_share": pathway_info.get("vote_share"),
            "pathway_confidence_tier": pathway_info.get("tier", ""),
        },
        observed_mz_available=_observed_mz is not None and _observed_mz > 0,
        ms2_spectrum_present=_ms2_score is not None and _ms2_score > 0,
        ms2_score=_ms2_score,
        library_neighbors_ge_050=sum(
            1 for nb in neighbors if float(nb.get("similarity_to_target", 0.0) or 0.0) >= 0.50
        ),
        correlations_real=not correlation_blocker,
    )
    _conf_fields["observed_mz"] = _observed_mz
    _conf_fields["observed_neutral_mass"] = (
        _n(rank1.get("target_neutral_mass")) or _n(rank1.get("neutral_mass"))
    )
    _conf_fields["candidate_exact_mass"] = target["exact_mass"]
    # Parity with the enrich -> patch flow: Metabolon Level, mass-consistent
    # alternative, and observed-MS2 detail so the from-source render matches.
    _level = str(kx_row.get("Level", "") or "").strip()
    _conf_fields["metabolon_level"] = "" if _level.lower() in {"nan", "artifact", "exclude", "4?"} else _level
    _conf_fields["mass_alt"] = confmod.mass_consistent_alt(_css_records)
    _n_peaks = int(len(peaks)) if hasattr(peaks, "empty") and not peaks.empty else 0
    _conf_fields["ms2_detail"] = {
        # The observed-MS2 panel must reflect an *observed* spectrum (from --peak-csv /
        # the observed lookup), not the pipeline's predicted spectral score.
        "present": _n_peaks > 0,
        "n_pairs": _n_peaks,
        "n_ion_bins": 0,
        "n_loss_bins": 0,
        "forward": _n(rank1.get("forward_spectrum_score")),
        "predicted_v2": _n(rank1.get("predicted_spectrum_score_v2")),
        "match": _n(rank1.get("ms2_match_score")),
    }
    _conf_fields["ms2_n_pairs"] = _n_peaks
    # Reviewed verdict, matching the enrich -> patch flow so both delivery paths
    # honor the same curation ledger.
    confmod.apply_curation_entry(
        _conf_fields,
        confmod.load_curation_ledger(
            resolve(getattr(args, "curation_ledger", DEFAULT_CURATION_LEDGER), repo_root)
        ).get(control_id),
    )
    confidence_obj = confmod.compute_confidence(_conf_fields)
    # The MSI level decides what the compound may be called. Levels 4-5 withhold
    # the structure, so they cannot carry a structural name however well it
    # resolved; Level 3 must read as a hypothesis rather than an assertion.
    unqualified_name = target_name
    target_name = tnmod.display_name_for_level(
        unqualified_name,
        msi_level=confidence_obj.msi_level,
        formula=target_formula,
        observed_mz=_observed_mz,
        name_is_reference=name_resolved,
    )
    target["name"] = target_name
    _conf_fields["smiles"] = target_smiles
    _conf_fields["resolved_name"] = unqualified_name
    _conf_fields["display_name"] = target_name
    _conf_fields["synonyms"] = (
        []
        if confidence_obj.msi_level >= 4
        else [
            str(value).strip()
            for value in syn_df.get("synonym", pd.Series(dtype=str)).tolist()
            if str(value).strip()
        ][:8]
    )
    _conf_fields["name_resolved"] = name_resolved
    _conf_fields.update(tnmod.descriptor_fields(tentative))
    _conf_fields["msi_level"] = confidence_obj.msi_level
    _conf_fields["msi_label"] = confidence_obj.msi_label
    _conf_fields["msi_rationale"] = confidence_obj.msi_rationale
    _conf_fields["identity_status"] = confidence_obj.identity_status
    _conf_fields["identity_confidence"] = confidence_obj.identity_confidence
    _conf_fields["identity_basis"] = confidence_obj.identity_basis
    (data_dir / f"{slug}_confidence_inputs.json").write_text(
        json.dumps(_conf_fields, indent=2)
    )

    fragment_candidate_rows = top.head(int(args.candidate_tab_n)).fillna("").to_dict(orient="records")
    render_kwargs = dict(
        out_dir=out_dir,
        target=target,
        evidence=evidence,
        db_detail=db_detail,
        confidence=confidence_obj,
        conf_fields=_conf_fields,
        pathway_info=pathway_info,
        # Synonyms name a structure, so they go with it: withheld at Levels 4-5.
        synonyms=(
            [] if confidence_obj.msi_level >= 4 else syn_df.fillna("").to_dict(orient="records")
        ),
        neighbors=neighbors,
        nodes=nodes,
        edges=edges,
        candidate_scores=candidate_score_df.fillna("").to_dict(orient="records"),
        fragment_candidate_rows=fragment_candidate_rows,
        structure_panel_html=structure_panel_html,
        consensus_svg=consensus_svg,
        candidate_consensus_maps=candidate_consensus_maps,
        candidate_consensus_summaries=candidate_consensus_summaries,
        default_consensus_rank=default_consensus_rank,
        structural_possibility_tree_svg=structural_possibility_tree_svg,
        structural_possibility_rows=structural_possibility_rows,
        feature_total_count=int(len(peaks)),
        extended_network_svg_html=extended_network_svg(extended["render_nodes"], extended["render_edges"]),
        extended_unknown_rows=extended_unknown_rows,
        extended_known_rows=extended_known_rows,
        extended_stats={
            "n_known_all": int(len(extended["known_all"])),
            "n_x_all": int(len(extended["x_all"])),
            "n_known_render": int(len(known_render)),
            "n_x_render": int(len(x_render)),
        },
        correlation_rows=correlation_rows,
        correlation_blocker=correlation_blocker,
    )
    public_html = render_html(**render_kwargs, variant="public")
    internal_html = render_html(**render_kwargs, variant="internal")
    public_guide_targets = [
        (slugify_label(label), label)
        for label in [
            "Top Candidate Compound",
            "Parameters Contributing to ID Confidence",
            "Top 10 potential structures",
            "Top 10 Metabolon Library Matches",
            "Fragment Scores",
            "Mass Match Scores",
            "MS Fit Scores",
            "Top Correlating Metabolites",
            "Similarity Network",
        ]
    ]
    guides_html = render_guides_html(control_id, public_guide_targets)
    (out_dir / "public.html").write_text(public_html)
    (out_dir / "guides.html").write_text(guides_html)
    (out_dir / "internal.html").write_text(internal_html)
    (out_dir / "index.html").write_text(public_html)
    compound_html = out_dir / f"{control_id.replace('-', '_')}.html"
    compound_html.write_text(public_html)
    payload = {
        "out_dir": str(out_dir),
        "html": str(out_dir / "index.html"),
        "compound_html": str(compound_html),
        "public_html": str(out_dir / "public.html"),
        "guides_html": str(out_dir / "guides.html"),
        "internal_html": str(out_dir / "internal.html"),
        "n_neighbors": int(len(neighbors)),
        "target": target,
        "db_chemical_id": chemical_id,
        "control_id": control_id,
    }
    (out_dir / "manifest.json").write_text(json.dumps(payload, indent=2))
    print(json.dumps(payload, indent=2, sort_keys=True))
    return payload


# ── CLI entrypoint ────────────────────────────────────────────────────────────

def main() -> None:
    """CLI entrypoint: parse args and invoke ``build_report``."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--control-id", required=True, help="X compound control id, e.g. X-25003")
    parser.add_argument("--predictions", type=Path, required=True, help="Rerank or candidate pool CSV")
    parser.add_argument("--out-dir", type=Path, default=None)
    parser.add_argument("--db", type=Path, default=DEFAULT_DB)
    parser.add_argument("--known-x", type=Path, default=DEFAULT_KNOWN_X)
    parser.add_argument("--known-x-all", type=Path, default=DEFAULT_KNOWN_X_ALL)
    parser.add_argument(
        "--curation-ledger",
        type=Path,
        default=DEFAULT_CURATION_LEDGER,
        help="Reviewed human/laboratory verdicts; see data/curation_ledger.csv.",
    )
    parser.add_argument("--candidate-pool", type=Path, default=None)
    parser.add_argument("--hypotheses", type=Path, default=None)
    parser.add_argument("--peak-csv", type=Path, default=None)
    parser.add_argument("--target-smiles", type=str, default="")
    parser.add_argument("--target-name", type=str, default="")
    parser.add_argument("--chemical-id", type=str, default="")
    parser.add_argument("--spectrum-tolerance-da", type=float, default=1.0)
    parser.add_argument("--min-similarity", type=float, default=0.28)
    parser.add_argument("--max-neighbors", type=int, default=10)
    parser.add_argument("--candidate-tab-n", type=int, default=20)
    parser.add_argument("--consensus-top-n", type=int, default=20)
    parser.add_argument("--extended-known-render-limit", type=int, default=120)
    parser.add_argument("--extended-known-min-similarity", type=float, default=0.30)
    parser.add_argument("--extended-x-min-similarity", type=float, default=0.0)
    parser.add_argument("--extended-target-edge-min", type=float, default=0.32)
    parser.add_argument("--extended-pair-edge-min", type=float, default=0.42)
    parser.add_argument("--extended-max-edges-per-node", type=int, default=6)
    args = parser.parse_args()
    if args.out_dir is None:
        args.out_dir = Path("output/reports") / args.control_id
    build_report(args)


if __name__ == "__main__":
    main()
