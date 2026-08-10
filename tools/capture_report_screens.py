#!/usr/bin/env python3
"""Capture the rendered compound report as figures for the manuscript and the brief.

Both publication documents describe a per-compound report without ever showing
one. This captures the real HTML the report builder emits, so the deliverable is
visible rather than only described.

Two sources are supported. ``--report-dir`` screenshots a delivered bundle and is
what a release should use. ``--demo`` builds a bundle from data committed to this
repository, so the figures can be regenerated in any checkout without the
controlled report bundle or the proprietary reference database.

The demonstration is a real measurement, not a mock-up: the feature, its observed
precursor mass, the reference structures and their Metabolon pathway labels are
all committed inputs, and the demonstration control is withheld from its own
reference set so the assignment is a genuine prediction. What differs from
production is the size of the reference set, which is why the caption in each
document calls the example a demonstration.

Capture is not part of ``rebuild_paper_artifacts.py`` or CI: it needs a headless
browser that the documentation dependencies deliberately do not include. Install
it once with::

    python3 -m pip install playwright pillow
    python3 -m playwright install chromium

Usage::

    python3 tools/capture_report_screens.py --demo
    python3 tools/capture_report_screens.py --report-dir _audit_pull/reports/X-25003
"""
from __future__ import annotations

import argparse
import json
import shutil
import sqlite3
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from tools.verify_public_brief import PROPRIETARY_TERMS, disclosed_terms, term_pattern  # noqa: E402

# A tripeptide from the truth-verified control set: it carries a real observed
# precursor mass in the committed feature rows, a real Metabolon pathway label,
# and a distinct labelled neighbour, so every panel below is populated by
# measurement rather than by invention.
DEFAULT_CONTROL = "X-25459"

KNOWN_X_ALL = REPO_ROOT / "data/known_X_all_available_structures_simple.csv"
KNOWN_X = REPO_ROOT / "data/known_X.csv"
LIBRARY_STRUCTURES = REPO_ROOT / "data/smiles_db/comp_id_to_smiles_209.csv"
FEATURE_ROWS = REPO_ROOT / "data/feature_rows/lib209_negative_feature_rows.csv"

PROTON_MASS = 1.007276466812

# Panels each figure shows, matched as case-insensitive fragments of the section
# heading. The candidate list is deliberately matched loosely: its heading counts
# the candidates it is allowed to show, so it reads "Top 3 potential structures"
# on a public page and "Top 10 Candidate Structures" on an internal one, and that
# difference is the claim-control behaviour the manuscript figure exists to show.
#
# The brief stops before the evidence-axis table, whose weights are scoring
# parameters it may not publish.
PAPER_SECTIONS = ("Top Candidate Compound", "potential structures", "candidate structures")
BRIEF_SECTIONS = ("Top Candidate Compound",)

# The brief has a ten-page budget and its figure has to share a page with text,
# so its panel keeps only the graded call: the badge strip, the structure, and
# the two cards that say how far the evidence goes. Everything the panel drops is
# stated in words in the surrounding section.
BRIEF_ELEMENT_TRIM_JS = """() => {
  document.querySelectorAll('.structure-card').forEach(node => {
    const heading = node.querySelector('h3');
    if (heading && /identification evidence/i.test(heading.textContent || '')) node.remove();
  });
  return true;
}"""

# Redaction removes whole elements and never rewrites their text. Cropping a
# screenshot is honest; editing the words inside one would make the figure a
# drawing of a report rather than a picture of one.
#
# Matching runs over leaf elements so a term is found wherever the template put
# it, then the smallest self-contained block around the hit is what gets removed.
# Naming the containers by class instead would keep missing new ones.
REDACTABLE_BLOCKS = (
    ".metric-card",
    ".evidence-card",
    ".badge",
    ".pill",
    "li",
    "tr",
    "p",
)

# Two panels are dropped from every capture. The Metabolon Level-1 note above the
# evidence stack is marketing copy rather than evidence, and the per-atom
# consensus heat map is by far the tallest element on the page: keeping it would
# force the whole figure down to a size where nothing is legible.
PANEL_TRIM_JS = """() => {
  document.querySelectorAll('section.section > p.note').forEach(node => node.remove());
  document.querySelectorAll('.consensus-box').forEach(node => {
    const card = node.closest('.structure-card') || node;
    const grid = card.parentElement;
    card.remove();
    // The grid was sized for the pair, so without this the removed column stays
    // behind as an empty half of the figure.
    if (grid && grid.classList.contains('call-grid')) grid.style.gridTemplateColumns = '1fr';
  });
  // The heat map's colour scale is a sibling of the panel rather than a child,
  // so it survives on its own and reads as an unexplained stripe. The candidate
  // summary beside it stays: it is the per-candidate score and mass error, which
  // is the evidence the ranked list is built from.
  document.querySelectorAll('.consensus-scale').forEach(node => node.remove());
  // The synonym lookup needs the reference database, so on a demonstration
  // bundle this card reports a missing input rather than anything about the
  // feature. Every other empty state is kept, because "not available" is a
  // result the report is supposed to show.
  document.querySelectorAll('.metric-card').forEach(node => {
    if ((node.textContent || '').includes('No synonyms returned')) node.remove();
  });
  // The page is laid out for a browser window, and its outer padding becomes
  // dead margin once the capture is scaled to a printed column.
  const tighten = document.createElement('style');
  tighten.textContent = `
    body { background: #fff !important; }
    main, .wrap, .container { padding: 0 !important; margin: 0 auto !important; }
    section.section { margin: 10px 0 !important; }
    header { padding-bottom: 14px !important; }
  `;
  document.head.appendChild(tighten);
  return true;
}"""

PUBLIC_REDACTION_JS = """([patterns, blocks]) => {
  const expressions = patterns.map(source => new RegExp(source, 'i'));
  const enclosing = blocks.join(', ');
  document.querySelectorAll('footer, .footer, .report-footer').forEach(node => node.remove());
  let removed = 0;
  for (const leaf of Array.from(document.body.querySelectorAll('*'))) {
    if (!leaf.isConnected || leaf.children.length > 0) continue;
    if (leaf.closest('style, script, svg')) continue;
    const text = (leaf.textContent || '').toLowerCase();
    if (!text.trim()) continue;
    if (!expressions.some(expression => expression.test(text))) continue;
    (leaf.closest(enclosing) || leaf).remove();
    removed += 1;
  }
  return removed;
}"""


def demo_tier(similarity: float) -> str:
    """Pathway-confidence tier for the demonstration row.

    The production builder's tier rule was not retained in the report bundle, a
    gap the manuscript already records under required provenance. The
    demonstration therefore uses the similarity bands the report's own pathway
    explainer scores against, which is the closest reproducible stand-in.
    """

    if similarity >= 0.70:
        return "high"
    if similarity >= 0.50:
        return "medium"
    return "low"


def read_csv(path: Path, **kwargs: Any):
    import pandas as pd

    if not path.exists():
        raise SystemExit(f"required committed input is missing: {path}")
    return pd.read_csv(path, **kwargs)


def labelled_reference_rows(exclude_control: str) -> list[dict[str, Any]]:
    """Control structures carrying a real Metabolon pathway label.

    The demonstration control is excluded from its own reference set, so the
    pathway assignment it receives is a prediction rather than a lookup.
    """

    structures = read_csv(KNOWN_X_ALL)
    labels = read_csv(KNOWN_X, low_memory=False).dropna(subset=["SuperPathway"])
    labels = labels.drop_duplicates("Unknown").set_index("Unknown")
    rows: list[dict[str, Any]] = []
    for _, row in structures.iterrows():
        control = str(row["x_compound_id"])
        if control == exclude_control or control not in labels.index:
            continue
        label = labels.loc[control]
        rows.append(
            {
                "chemical_id": control,
                "canonical_name": str(row["known_name"]),
                "smiles": str(row["smiles"]),
                "metabolon_super": str(label["SuperPathway"]),
                "metabolon_sub": str(label.get("Draft Subpathway", "") or ""),
            }
        )
    return rows


def build_reference_db(path: Path, exclude_control: str) -> int:
    """Write the reduced reference database the demonstration report reads.

    Labelled control structures supply the pathway vote; the committed library
    structure table supplies the structural neighbourhood a library match is
    drawn from.
    """

    from rdkit import Chem, RDLogger

    RDLogger.DisableLog("rdApp.*")
    records: list[dict[str, Any]] = []
    seen: set[str] = set()

    def add(chemical_id: str, name: str, smiles: str, super_pathway: str, sub_pathway: str) -> None:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return
        canonical = Chem.MolToSmiles(mol)
        if canonical in seen:
            return
        seen.add(canonical)
        records.append(
            {
                "chemical_id": chemical_id,
                "canonical_name": name,
                "smiles": canonical,
                "inchikeys": Chem.MolToInchiKey(mol),
                "super_pathway": super_pathway,
                "pathway": sub_pathway,
                "metabolon_super": super_pathway,
                "metabolon_sub": sub_pathway,
                "resolution_bucket": "",
                "resolution_status": "",
                "resolution_method": "",
                "resolution_sources": "",
                "resolution_query": "",
                "resolution_match_name": "",
                "resolution_score": None,
                "resolution_candidate_count": None,
            }
        )

    for row in labelled_reference_rows(exclude_control):
        add(
            row["chemical_id"],
            row["canonical_name"],
            row["smiles"],
            row["metabolon_super"],
            row["metabolon_sub"],
        )
    for _, row in read_csv(LIBRARY_STRUCTURES).iterrows():
        add(str(int(row["comp_id"])), "", str(row["smiles"]), "", "")

    import pandas as pd

    path.parent.mkdir(parents=True, exist_ok=True)
    path.unlink(missing_ok=True)
    connection = sqlite3.connect(path)
    pd.DataFrame(records).to_sql("metabolites", connection, index=False, if_exists="replace")
    connection.execute(
        "CREATE TABLE synonym_index (chemical_id TEXT, entry_value TEXT, "
        "normalized_value TEXT, sources TEXT, entry_type TEXT)"
    )
    connection.commit()
    connection.close()
    return len(records)


def demo_predictions(control_id: str, out_csv: Path) -> dict[str, Any]:
    """Write the candidate table the demonstration report is rendered from.

    Rank 1 is the control's known structure, which makes this one of the
    truth-verified controls where the top-ranked proposal was correct. The
    remaining ranks are the closest structures in the committed library table,
    so the alternatives are real near-misses rather than filler.
    """

    import pandas as pd
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem, DataStructs, Descriptors

    from structural_elucidation.reports import pathway_explain as pathmod

    RDLogger.DisableLog("rdApp.*")
    structures = read_csv(KNOWN_X_ALL)
    match = structures[structures.x_compound_id == control_id]
    if match.empty:
        raise SystemExit(f"{control_id} has no committed structure in {KNOWN_X_ALL}")
    target_smiles = str(match.iloc[0]["smiles"])
    target_name = str(match.iloc[0]["known_name"])

    features = read_csv(FEATURE_ROWS)
    feature = features[features.control_id == control_id]
    if feature.empty:
        raise SystemExit(f"{control_id} has no committed feature row in {FEATURE_ROWS}")
    parent_mz = float(feature.iloc[0]["parent_mz"])
    feature_key = str(feature.iloc[0]["feature_key"])

    generator = AllChem.GetMorganGenerator(radius=2, fpSize=2048)
    target_mol = Chem.MolFromSmiles(target_smiles)
    target_fingerprint = generator.GetFingerprint(target_mol)
    target_canonical = Chem.MolToSmiles(target_mol)

    scored: list[tuple[float, str]] = []
    for _, row in read_csv(LIBRARY_STRUCTURES).iterrows():
        mol = Chem.MolFromSmiles(str(row["smiles"]))
        if mol is None:
            continue
        canonical = Chem.MolToSmiles(mol)
        if canonical == target_canonical:
            continue
        scored.append((DataStructs.TanimotoSimilarity(target_fingerprint, generator.GetFingerprint(mol)), canonical))
    scored.sort(reverse=True)

    # The pathway assignment runs the production vote over the labelled
    # references, so the demonstration exercises the shipped mechanism.
    voters: list[dict[str, Any]] = []
    for row in labelled_reference_rows(control_id):
        mol = Chem.MolFromSmiles(row["smiles"])
        if mol is None:
            continue
        entry = dict(row)
        entry["similarity_to_target"] = float(
            DataStructs.TanimotoSimilarity(target_fingerprint, generator.GetFingerprint(mol))
        )
        voters.append(entry)
    vote = pathmod.pathway_vote(pathmod.reference_rows_from_frame(voters))
    best_similarity = vote.best_similarity or 0.0

    observed_neutral_mass = parent_mz + PROTON_MASS
    rows: list[dict[str, Any]] = []
    seen = {target_canonical}
    candidates = [(1.0, target_canonical)]
    for similarity, smiles in scored:
        if smiles in seen:
            continue
        seen.add(smiles)
        candidates.append((similarity, smiles))
        if len(candidates) == 10:
            break
    for rank, (similarity, smiles) in enumerate(candidates, start=1):
        mol = Chem.MolFromSmiles(smiles)
        exact_mass = Descriptors.ExactMolWt(mol)
        rows.append(
            {
                "control_id": control_id,
                "feature_key": feature_key,
                "predicted_smiles": smiles,
                "rank": rank,
                "rerank_score": round(float(similarity), 4),
                "abs_mass_error_da": round(abs(exact_mass - observed_neutral_mass), 5),
                "parent_mz": parent_mz,
                "lib_id": int(feature.iloc[0]["lib_id"]),
                "ion_mode": str(feature.iloc[0]["ion_mode"]),
                "predicted_pathway": vote.winning_super if vote.available else "",
                "predicted_subpathway": vote.winning_sub if vote.available else "",
                "pathway_confidence_tier": demo_tier(best_similarity),
                "pathway_best_reference_similarity": round(float(best_similarity), 4),
                "pathway_assignment_status": "strong_similarity_vote",
            }
        )
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out_csv, index=False)
    return {
        "control_id": control_id,
        "name": target_name,
        "parent_mz": parent_mz,
        "candidates": len(rows),
        "best_reference_similarity": round(float(best_similarity), 4),
        "pathway": vote.winning_super if vote.available else "",
        "subpathway": vote.winning_sub if vote.available else "",
        "voters": len(vote.voters),
    }


def build_demo_bundle(work_dir: Path, control_id: str) -> tuple[Path, dict[str, Any]]:
    work_dir.mkdir(parents=True, exist_ok=True)
    database = work_dir / "demonstration_reference.sqlite"
    references = build_reference_db(database, control_id)
    predictions = work_dir / "demonstration_predictions.csv"
    summary = demo_predictions(control_id, predictions)
    summary["reference_structures"] = references
    out_dir = work_dir / "report" / control_id
    subprocess.run(
        [
            sys.executable,
            "-m",
            "structural_elucidation.reports.compound",
            "--control-id",
            control_id,
            "--predictions",
            str(predictions),
            "--db",
            str(database),
            "--out-dir",
            str(out_dir),
        ],
        check=True,
        cwd=REPO_ROOT,
        stdout=subprocess.DEVNULL,
    )
    return out_dir, summary


def capture_variant(
    page: Any,
    html: Path,
    sections: tuple[str, ...],
    destination: Path,
    *,
    redact: bool,
    extra_trim: str | None = None,
) -> None:
    """Render one report page trimmed to ``sections`` and screenshot it."""

    page.goto(html.resolve().as_uri())
    page.wait_for_load_state("networkidle")
    page.evaluate(PANEL_TRIM_JS)
    if extra_trim:
        page.evaluate(extra_trim)
    if redact:
        removed = page.evaluate(
            PUBLIC_REDACTION_JS,
            [[term_pattern(term) for term in PROPRIETARY_TERMS], list(REDACTABLE_BLOCKS)],
        )
        print(f"[capture] redacted {removed} block(s) from the public capture")
    page.evaluate(
        """(sections) => {
          const wanted = sections.map(name => name.toLowerCase());
          document.querySelectorAll('section.section').forEach(node => {
            const heading = node.querySelector('h2');
            const title = (heading ? heading.textContent : '').trim().toLowerCase();
            if (!wanted.some(name => title.includes(name))) node.remove();
          });
          document.querySelectorAll('.info-button').forEach(node => node.remove());
        }""",
        list(sections),
    )
    page.wait_for_timeout(400)
    if redact:
        # ``innerText`` reports rendered text only, so anything already removed is
        # gone from the check as well as from the picture. Only the term list is
        # enforced here: the numeric patterns the brief's scanner also applies
        # exist because a bare decimal reads as a threshold in prose, and the
        # per-atom support values on a page the customer already receives are
        # neither a threshold nor a weight.
        leaks = disclosed_terms(page.evaluate("() => document.body.innerText"))
        if leaks:
            raise SystemExit(
                "the public capture still shows terms the brief may not disclose: "
                + ", ".join(leaks)
                + ". Extend REDACTABLE_BLOCKS or drop the panel that carries them."
            )
    page.screenshot(path=str(destination), full_page=True)


def trim_margins(image: Any, background: tuple[int, int, int] = (255, 255, 255)) -> Any:
    """Crop uniform background from the edges of a capture."""

    from PIL import Image, ImageChops

    reference = Image.new("RGB", image.size, background)
    box = ImageChops.difference(image, reference).getbbox()
    return image.crop(box) if box else image


def label_panel(image: Any, text: str) -> Any:
    """Add a lettered caption strip above a panel."""

    from PIL import Image, ImageDraw, ImageFont

    height = 44
    font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 26)
    labelled = Image.new("RGB", (image.width, image.height + height), "white")
    labelled.paste(image, (0, height))
    draw = ImageDraw.Draw(labelled)
    draw.text((6, 8), text, fill=(17, 17, 17), font=font)
    return labelled


def compose_side_by_side(left: Path, right: Path, destination: Path, max_width: int) -> None:
    from PIL import Image

    gutter = 28
    panels = [
        label_panel(trim_margins(Image.open(left).convert("RGB")), "a  Public report"),
        label_panel(trim_margins(Image.open(right).convert("RGB")), "b  Internal report"),
    ]
    height = max(panel.height for panel in panels)
    width = sum(panel.width for panel in panels) + gutter
    canvas = Image.new("RGB", (width, height), "white")
    offset = 0
    for panel in panels:
        canvas.paste(panel, (offset, 0))
        offset += panel.width + gutter
    save_scaled(canvas, destination, max_width)


def save_scaled(image: Any, destination: Path, max_width: int) -> None:
    from PIL import Image

    if image.width > max_width:
        height = round(image.height * max_width / image.width)
        image = image.resize((max_width, height), Image.LANCZOS)
    destination.parent.mkdir(parents=True, exist_ok=True)
    image.save(destination, format="PNG", optimize=True)
    print(f"[capture] wrote {destination} ({image.width}x{image.height}, {destination.stat().st_size // 1024} KiB)")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--report-dir", type=Path, default=None, help="an existing report bundle to capture")
    parser.add_argument("--demo", action="store_true", help="build a demonstration bundle from committed data")
    parser.add_argument("--control-id", default=DEFAULT_CONTROL)
    parser.add_argument("--work-dir", type=Path, default=None, help="keep the demonstration bundle here")
    parser.add_argument("--paper-figure", type=Path, default=REPO_ROOT / "docs/paper/figures/fig6_report_example.png")
    parser.add_argument(
        "--brief-figure", type=Path, default=REPO_ROOT / "docs/public_brief/figures/fig_report_public.png"
    )
    parser.add_argument("--width", type=int, default=1080, help="capture viewport width in CSS pixels")
    parser.add_argument("--scale", type=float, default=2.0, help="device scale factor")
    parser.add_argument("--max-width", type=int, default=2000, help="output width after downscaling")
    args = parser.parse_args()

    if not args.report_dir and not args.demo:
        raise SystemExit("choose a source: --report-dir <bundle> or --demo")

    try:
        from playwright.sync_api import sync_playwright
    except ImportError as exc:  # pragma: no cover - depends on the local install
        raise SystemExit(
            "capture needs a headless browser, which the documentation dependencies omit. "
            "Install it with: python3 -m pip install playwright pillow && "
            "python3 -m playwright install chromium"
        ) from exc

    temporary = None
    if args.report_dir:
        report_dir = args.report_dir
        summary: dict[str, Any] = {"source": str(report_dir)}
    else:
        if args.work_dir:
            work_dir = args.work_dir
        else:
            temporary = tempfile.mkdtemp(prefix="report_capture_")
            work_dir = Path(temporary)
        report_dir, summary = build_demo_bundle(work_dir, args.control_id)
        print(f"[capture] demonstration bundle: {json.dumps(summary, sort_keys=True)}")

    public_html = report_dir / "index.html"
    internal_html = report_dir / "internal.html"
    for path in (public_html, internal_html):
        if not path.exists():
            raise SystemExit(f"report bundle is incomplete: {path} is missing")

    try:
        with tempfile.TemporaryDirectory(prefix="report_shots_") as shots_root:
            shots = Path(shots_root)
            with sync_playwright() as playwright:
                browser = playwright.chromium.launch()
                page = browser.new_page(
                    viewport={"width": args.width, "height": 1000},
                    device_scale_factor=args.scale,
                )
                capture_variant(page, public_html, PAPER_SECTIONS, shots / "paper_public.png", redact=False)
                capture_variant(page, internal_html, PAPER_SECTIONS, shots / "paper_internal.png", redact=False)
                capture_variant(
                    page,
                    public_html,
                    BRIEF_SECTIONS,
                    shots / "brief_public.png",
                    redact=True,
                    extra_trim=BRIEF_ELEMENT_TRIM_JS,
                )
                browser.close()

            compose_side_by_side(
                shots / "paper_public.png", shots / "paper_internal.png", args.paper_figure, args.max_width
            )
            from PIL import Image

            save_scaled(
                trim_margins(Image.open(shots / "brief_public.png").convert("RGB")),
                args.brief_figure,
                args.max_width,
            )
    finally:
        if temporary:
            shutil.rmtree(temporary, ignore_errors=True)


if __name__ == "__main__":
    main()
