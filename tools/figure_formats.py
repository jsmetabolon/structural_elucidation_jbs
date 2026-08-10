#!/usr/bin/env python3
"""Validate the figure formats the publication documents are allowed to link.

Both verifiers used to parse every linked figure as XML, which was correct while
every figure was a generated SVG. Report screenshots are raster by nature, so the
format check now branches: vector figures still have to parse, and raster figures
are validated from their own header rather than being waved through on the file
extension alone.
"""
from __future__ import annotations

import struct
import xml.etree.ElementTree as ET
from pathlib import Path

PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"
RASTER_SUFFIXES = frozenset({".png"})

# A screenshot narrower than this cannot carry legible interface text once it is
# scaled to the printed column width, so it is a defect rather than a style choice.
MIN_RASTER_WIDTH = 900


def png_dimensions(path: Path) -> tuple[int, int]:
    """Pixel width and height read from the PNG IHDR chunk."""

    header = path.read_bytes()[:24]
    if len(header) < 24 or not header.startswith(PNG_SIGNATURE):
        raise ValueError(f"not a PNG file: {path}")
    if header[12:16] != b"IHDR":
        raise ValueError(f"PNG is missing its IHDR chunk: {path}")
    width, height = struct.unpack(">II", header[16:24])
    return width, height


def figure_error(path: Path) -> str | None:
    """Return why a linked figure is unusable, or ``None`` when it is fine."""

    suffix = path.suffix.lower()
    if suffix == ".svg":
        try:
            ET.parse(path)
        except ET.ParseError as exc:
            return f"invalid SVG: {exc}"
        return None
    if suffix in RASTER_SUFFIXES:
        try:
            width, height = png_dimensions(path)
        except (OSError, ValueError) as exc:
            return f"invalid PNG: {exc}"
        if width < MIN_RASTER_WIDTH:
            return f"raster figure is only {width}px wide; {MIN_RASTER_WIDTH}px is the minimum for legible text"
        if height <= 0:
            return f"raster figure has a non-positive height: {height}"
        return None
    return f"unsupported figure format '{suffix}'; use .svg for charts or .png for screenshots"
