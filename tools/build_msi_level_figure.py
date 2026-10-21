#!/usr/bin/env python3
"""Draw the L1-L5 identification-level distribution, one column per level.

Same panel treatment as the left half of fig4: one bar per level in the brand
series order, count and share above the bar, square-root height scale so L5
stays visible next to L3. Two cohorts are supported, because the two questions
people ask of this chart have different denominators:

``reports``
    The 3,897 delivered X-compound reports, from the exec-summary export. This
    is the cohort fig4 and the presentation audit use.

``inventory``
    The 8,992 levelled metabolites of the all-metabolites workbook, known
    library and X compounds together.

Two readings of the L2/L3 boundary are supported:

``delivered``
    The level as shipped.

``generous``
    Anything at L3 is promoted to L2 when it carries a specific compound name.
    Only rows whose "name" is a placeholder -- a class label, a formula-based
    conjugate descriptor, a peptide composition, or a machine-generated IUPAC
    string -- stay at L3. A "putative analog of X (Tanimoto 0.65)" note names
    its compound by the neighbour it resembles, and counts as named. L1, L4 and
    L5 are untouched: they turn on the reference standard, the mass mismatch,
    and the missing observed mass, none of which a name can change.

Usage:
    python3 tools/build_msi_level_figure.py --png
    python3 tools/build_msi_level_figure.py --titled   # add title and caption
"""
from __future__ import annotations

import argparse
import csv
import math
import re
import subprocess
import sys
from collections import Counter
from pathlib import Path
from typing import Any, Iterable

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from openpyxl import load_workbook  # noqa: E402

from tools.brand import DATA_LABEL, MUTED, SERIES  # noqa: E402
from tools.paper_figures import axis_text, canvas, rule, text, value_text  # noqa: E402

LEVELS = ("1", "2", "3", "4", "5")
KNOWN = "Known"
UNKNOWN = "X"

# Names that describe a family, a composition, or a derivation rather than
# naming a compound. A row holding one of these has not been named, so the
# generous rule leaves it at L3.
PLACEHOLDER_NAMES = frozenset(
    {
        "acylcarnitine",
        "alkaloid",
        "amino acid",
        "amino acid derivative",
        "aromatic molecule",
        "artifact",
        "ask dts",
        "bile acid",
        "can add something here?",
        "can add this",
        "carbohydrate",
        "carnitine conjugate",
        "collection tube artifact",
        "de novo structure",
        "dipeptide",
        "do not include",
        "exclude",
        "fatty acid",
        "fatty acid derivative",
        "flavonoid",
        "glucoside",
        "hydroxy acid derivative",
        "lipid",
        "no further information",
        "nucleoside",
        "nucleotide",
        "organic acid",
        "partially characterized molecule",
        "pentapeptide",
        "peptide",
        "phenol",
        "polyphenol",
        "small peptide",
        "steroid",
        "sterol",
        "sugar",
        "terpenoid",
        "tetrapeptide",
        "tripeptide",
        "unknown",
    }
)

PLACEHOLDER_PATTERNS = (
    # A generated IUPAC string is a reading of the structure, not an identification.
    re.compile(r"\biupac \(generated\)", re.I),
    # "Sulfate conjugate of C9H10O3" / "glucuronide of C12H20O3 (5)" carry a formula, not an identity.
    re.compile(r"\bconjugate of c\d", re.I),
    re.compile(r"^[a-z]+(?: [a-z]+)? of c\d", re.I),
    re.compile(r"\b(derivative|isomer) c\d", re.I),
    # Peptide composition or ambiguous-residue sequence placeholders.
    re.compile(r"^(?:di|tri|tetra|penta)?peptide\b.*(?:composition|with |[\u2013\u2014/-])", re.I),
    re.compile(r"^small peptide", re.I),
    # Bare class labels, with or without a conjugation or derivative suffix.
    re.compile(
        r"^(?:amino acid|fatty acid|hydroxy acid|dicarboxylic acid|organic acid|bile acid|steroid|"
        r"sterol|lipid|peptide|sugar|carbohydrate|nucleotide|nucleoside|acylcarnitine|flavonoid|"
        r"phenol|alkaloid|terpenoid|glucoside|aromatic molecule|mono|di|oligo)saccharide?"
        r"[ \-]*(?:derivative|conjugate|glucuronidated|sulfated)?$",
        re.I,
    ),
    re.compile(r"-like compound$", re.I),
    # Bookkeeping strings that leaked into the name column.
    re.compile(r"^confidence \d", re.I),
    re.compile(r"^x-\d+$", re.I),
)


def is_named(value: Any) -> bool:
    """True when the name column identifies a compound rather than a family."""

    name = str(value or "").strip()
    if not name:
        return False
    if name.lower() in PLACEHOLDER_NAMES:
        return False
    return not any(pattern.search(name) for pattern in PLACEHOLDER_PATTERNS)


def read_inventory(workbook: Path, sheet: str) -> list[dict[str, str]]:
    """Level, cohort, and name per metabolite from the all-metabolites workbook."""

    book = load_workbook(workbook, read_only=True, data_only=True)
    if sheet not in book.sheetnames:
        raise SystemExit(f"{workbook} has no sheet {sheet!r}; found {book.sheetnames}")
    stream = book[sheet].iter_rows(values_only=True)
    header = [str(cell) for cell in next(stream)]
    rows = []
    for values in stream:
        row = dict(zip(header, values))
        rows.append(
            {
                "level": str(row.get("MSI Level") or "").strip(),
                "cohort": KNOWN if str(row.get("Entity Type") or "").strip() == KNOWN else UNKNOWN,
                "name": str(row.get("Name") or "").strip(),
            }
        )
    return rows


def read_reports(export: Path) -> list[dict[str, str]]:
    """Level and display name per delivered X-compound report."""

    with export.open(newline="", encoding="utf-8") as handle:
        return [
            {
                "level": (row.get("msi_level") or "").strip(),
                "cohort": UNKNOWN,
                "name": (row.get("feature_name") or "").strip(),
            }
            for row in csv.DictReader(handle)
        ]


def tally(rows: Iterable[dict[str, str]], *, generous: bool) -> tuple[dict[str, Counter[str]], Counter[str]]:
    """Counts per level and cohort, plus the rows that carry no L1-L5 level."""

    counts: dict[str, Counter[str]] = {level: Counter() for level in LEVELS}
    unlevelled: Counter[str] = Counter()
    for row in rows:
        level = row["level"]
        if level not in LEVELS:
            unlevelled[level or "(blank)"] += 1
            continue
        if generous and level == "3" and is_named(row["name"]):
            level = "2"
        counts[level][row["cohort"]] += 1
    return counts, unlevelled


def wrap(value: str, width: int) -> list[str]:
    """Break a caption on spaces at roughly `width` characters per line."""

    lines: list[str] = []
    line = ""
    for word in value.split():
        if line and len(line) + 1 + len(word) > width:
            lines.append(line)
            line = word
        else:
            line = f"{line} {word}".strip()
    if line:
        lines.append(line)
    return lines


def level_chart(
    path: Path,
    counts: dict[str, Counter[str]],
    *,
    title: str = "",
    subtitle: str = "",
) -> None:
    """One column per level: brand series colour, count and share above the bar.

    Heights are on a square-root scale, as fig4 draws them, so a level holding a
    handful of compounds still reads as a bar rather than a line. That makes the
    columns comparable in rank but not in area, which is why every bar carries
    its own count.
    """

    width, height = 1024, 736
    left, right = 60, 40
    caption = wrap(subtitle, 128) if subtitle else []
    header = (54 if title else 0) + 15 * len(caption)
    top, bottom = header + 46, 100
    plot_w, plot_h = width - left - right, height - top - bottom
    totals = [sum(counts[level].values()) for level in LEVELS]
    grand_total = sum(totals)
    max_root = math.sqrt(max(totals) or 1)
    slot = plot_w / len(LEVELS)
    bar_w = slot * 0.56

    svg = canvas(width, height)
    if title:
        svg.append(text(left, 34, title, size=20, anchor="start", bold=True))
    for index, line in enumerate(caption):
        svg.append(text(left, 58 + 15 * index, line, size=12, anchor="start", color=MUTED))

    baseline = top + plot_h
    for index, level in enumerate(LEVELS):
        total = totals[index]
        bar_h = math.sqrt(total) / max_root * plot_h if total else 0.0
        center = left + slot * (index + 0.5)
        x = center - bar_w / 2
        y = baseline - bar_h
        if bar_h > 0:
            svg.append(
                f"<rect x='{x:.1f}' y='{y:.1f}' width='{bar_w:.1f}' "
                f"height='{bar_h:.1f}' fill='{SERIES[index]}'/>"
            )
        share = 100 * total / grand_total if grand_total else 0.0
        svg.append(value_text(center, y - 12, f"{total} ({share:.1f}%)", size=17, color=DATA_LABEL))
        svg.append(axis_text(center, baseline + 34, f"L{level}", size=15))
    svg.append(rule(left, baseline, left + plot_w, baseline))
    svg.append("</svg>")
    path.write_text("\n".join(svg) + "\n")


def write_counts(path: Path, tallies: dict[tuple[str, str], dict[str, Counter[str]]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["cohort", "rule", "msi_level", "known", "x", "total"])
        for (source, name), counts in tallies.items():
            for level in LEVELS:
                known = counts[level][KNOWN]
                unknown = counts[level][UNKNOWN]
                writer.writerow([source, name, f"L{level}", known, unknown, known + unknown])


def render_png(svg_path: Path) -> Path | None:
    """Rasterize for slides, in the vendored typeface, if a renderer is present."""

    from tools import brand

    brand.install_fonts()
    png = svg_path.with_suffix(".png")
    try:
        import cairosvg
    except ImportError:
        pass
    else:
        cairosvg.svg2png(url=str(svg_path), write_to=str(png), scale=2.0)
        return png
    try:
        subprocess.run(
            ["rsvg-convert", "-z", "2", str(svg_path), "-o", str(png)],
            check=True,
            capture_output=True,
        )
    except (FileNotFoundError, subprocess.CalledProcessError):
        return None
    return png


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--workbook",
        type=Path,
        default=Path("docs/paper/data/all_metabolites_matrix_detection.xlsx"),
    )
    parser.add_argument("--sheet", default="Latest_Unique_X")
    parser.add_argument(
        "--reports",
        type=Path,
        default=Path("docs/paper/data/exec_summary/exec_summary_by_compound.csv"),
    )
    parser.add_argument("--out-dir", type=Path, default=Path("docs/paper/figures"))
    parser.add_argument("--counts-csv", type=Path, default=Path("docs/paper/data/msi_level_counts_by_rule.csv"))
    parser.add_argument("--titled", action="store_true", help="draw the title and provenance caption")
    parser.add_argument("--png", action="store_true", help="also render PNG copies for slides")
    args = parser.parse_args()

    if not args.reports.exists():
        raise SystemExit(f"missing {args.reports}")
    if not args.workbook.exists():
        raise SystemExit(f"missing {args.workbook}")
    sources = {
        "reports": (read_reports(args.reports), f"{args.reports.name}, one row per delivered X report"),
        "inventory": (
            read_inventory(args.workbook, args.sheet),
            f"{args.sheet} of {args.workbook.name}, one row per metabolite",
        ),
    }
    rules = {
        "delivered": "levels exactly as shipped",
        "generous": (
            "every named L3 promoted to L2; only class labels, formula-based conjugate "
            "descriptors, peptide compositions and generated IUPAC strings stay at L3"
        ),
    }

    args.out_dir.mkdir(parents=True, exist_ok=True)
    tallies: dict[tuple[str, str], dict[str, Counter[str]]] = {}
    for source, (rows, provenance) in sources.items():
        for name, rule_text in rules.items():
            counts, unlevelled = tally(rows, generous=name == "generous")
            tallies[(source, name)] = counts
            path = args.out_dir / f"msi_level_{source}_{name}.svg"
            dropped = sum(unlevelled.values())
            level_chart(
                path,
                counts,
                title=f"Identification level, {name} L2/L3 boundary" if args.titled else "",
                subtitle=(
                    f"{provenance}: {rule_text}. Square-root height scale; "
                    f"{dropped:,} rows without an L1-L5 level excluded."
                )
                if args.titled
                else "",
            )
            line = "  ".join(f"L{level}={sum(counts[level].values()):,}" for level in LEVELS)
            print(f"[msi-levels] {path}  {line}  (excluded {dropped:,})")
            if args.png:
                png = render_png(path)
                print(f"[msi-levels] {png}" if png else "[msi-levels] no SVG renderer, skipped PNG")

    write_counts(args.counts_csv, tallies)
    print(f"[msi-levels] wrote {args.counts_csv}")


if __name__ == "__main__":
    main()
