#!/usr/bin/env python3
"""The Metabolon brand palette and typography, as one source of truth.

Values are taken from the brand-guide pages and the chart formatting carried in
``figure_examples.xlsx``: the guide fixes the colours and their names, and the
workbook's own charts fix how they are applied. Both documents' figures are
generated from here so the manuscript and the public brief cannot drift apart
by one of them being edited and the other not.

The ordered ``SERIES`` sequence is not invented. It is the order of the LEGEND
swatch block on the workbook's reference sheet, and the same order the workbook
walks when a chart has more categories than roles.
"""

from __future__ import annotations

import atexit
import os
import shutil
import tempfile
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
FONT_DIR = PROJECT_ROOT / "assets" / "fonts" / "montserrat"

# Primary.
DARK_MODE_BLUE = "#001D38"
ACTIVATED_RED = "#FF4F58"
PURE_WHITE = "#FFFFFF"

# Accent.
MOLECULAR_BLUE = "#284A6F"
LIVING_PURPLE = "#855DA3"
INSIGHTFUL_ORANGE = "#FDA851"
METHOD_GREY = "#B2BBC3"
PRECISION_GREY = "#D3DFE7"
METABOLITE_WHITE = "#F2FAFC"

# Application colours, used here only as additional categorical steps.
AZURE_BLUE = "#357DED"
MINT_GREEN = "#52B788"
COTTON_CANDY_PINK = "#F4B6D1"

SERIES = (
    ACTIVATED_RED,
    DARK_MODE_BLUE,
    INSIGHTFUL_ORANGE,
    METHOD_GREY,
    LIVING_PURPLE,
    MINT_GREEN,
    COTTON_CANDY_PINK,
    AZURE_BLUE,
)

# Roles rather than colours, so a chart asks for what a mark means. The workbook
# is consistent about these: one measured series is orange, a second is the dark
# accent blue, and a trend drawn over bars is red so it reads as a different
# kind of thing rather than a third category.
BAR_PRIMARY = INSIGHTFUL_ORANGE
BAR_SECONDARY = MOLECULAR_BLUE
ACCENT_LINE = ACTIVATED_RED

# Chart furniture. The workbook draws gridlines and the value-axis rule in the
# same pale blue, frames the canvas in a near-white grey, and keeps every axis
# label well below the weight of the data.
GRID = "#CADFF4"
FRAME = "#D9D9D9"
AXIS_LABEL = "#7F7F7F"
DATA_LABEL = "#404040"
ERROR_BAR = "#808080"

INK = DARK_MODE_BLUE
MUTED = "#5C6B7A"
WHITE = PURE_WHITE

HAIRLINE = 0.75

# Axis labels are small, bold, and set in caps; values are the same size but
# regular weight, so the number reads before its label does.
AXIS_LABEL_SIZE = 9
DATA_LABEL_SIZE = 9

FAMILY = "Montserrat"

# Deliberately unquoted. These strings are interpolated into SVG attributes that
# are sometimes single- and sometimes double-quoted, and CSS accepts multi-word
# family names without quotes, so leaving them bare keeps one stack usable in
# both places. Mirrors the stack the delivered HTML reports already use.
FONT_STACK = "Montserrat, Avenir Next, Segoe UI, Arial, sans-serif"

_FONTCONFIG_DIR: Path | None = None


def on_fill(fill: str) -> str:
    """The text colour that stays legible on `fill`.

    Most of the palette is mid-toned. Insightful Orange and Method Grey in
    particular are light enough that white type on them all but disappears,
    while Dark Mode Blue needs white, so a label inside a bar cannot pick one
    colour and hold it across figures.
    """

    red, green, blue = (int(fill.lstrip("#")[index : index + 2], 16) / 255 for index in (0, 2, 4))
    luminance = 0.2126 * red + 0.7152 * green + 0.0722 * blue
    return WHITE if luminance < 0.45 else INK


def _fontconfig_document(font_dir: Path, cache_dir: Path) -> str:
    return (
        '<?xml version="1.0"?>\n'
        '<!DOCTYPE fontconfig SYSTEM "urn:fontconfig:fonts.dtd">\n'
        "<fontconfig>\n"
        f"  <dir>{font_dir}</dir>\n"
        f"  <cachedir>{cache_dir}</cachedir>\n"
        '  <include ignore_missing="yes">/etc/fonts/fonts.conf</include>\n'
        "</fontconfig>\n"
    )


def install_fonts() -> bool:
    """Make the vendored fonts visible to anything that renders through cairo.

    Call this before importing ``cairosvg``. Fontconfig reads ``FONTCONFIG_FILE``
    once, when it first initialises, and cairo initialises it on the first piece
    of text it draws, so setting the variable later is a race that usually loses
    silently: the render succeeds and the type is wrong.

    Returns False when the vendored directory is missing, leaving the ambient
    font configuration alone. A caller that needs the font to be certain should
    use ``font_is_resolvable``; this only reports that the configuration was
    written, not that the family resolved.
    """

    global _FONTCONFIG_DIR
    if not FONT_DIR.is_dir():
        return False
    if _FONTCONFIG_DIR is not None:
        return True

    directory = Path(tempfile.mkdtemp(prefix="metabolon_fontconfig_"))
    cache = directory / "cache"
    cache.mkdir()
    config = directory / "fonts.conf"
    config.write_text(_fontconfig_document(FONT_DIR, cache))
    os.environ["FONTCONFIG_FILE"] = str(config)
    _FONTCONFIG_DIR = directory
    atexit.register(shutil.rmtree, directory, True)
    return True


def font_is_resolvable(family: str = FAMILY) -> bool:
    """Report whether `family` actually renders, rather than falling back.

    An unresolved family is the one failure in this pipeline that nothing else
    catches: every gate still passes and the figures ship in the wrong typeface.
    So draw the same glyph twice, once in `family` and once in a family that
    cannot exist, and compare the pixels. Identical output means both requests
    landed on the same fallback and `family` was never found.
    """

    install_fonts()
    import cairosvg

    def render(requested: str) -> bytes:
        svg = (
            "<svg xmlns='http://www.w3.org/2000/svg' width='240' height='60'>"
            f"<text x='4' y='44' font-size='40' font-family='{requested}'>Rgy</text>"
            "</svg>"
        )
        return cairosvg.svg2png(bytestring=svg.encode())

    return render(family) != render("Metabolon No Such Family")


if __name__ == "__main__":
    installed = install_fonts()
    print(f"[brand] vendored fonts at {FONT_DIR}: {'configured' if installed else 'MISSING'}")
    print(f"[brand] {FAMILY} resolves: {font_is_resolvable()}")
