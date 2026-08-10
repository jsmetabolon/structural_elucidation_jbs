"""Tests for the figure-format gate shared by both publication verifiers.

The gate was relaxed so report screenshots can be linked as PNG. That relaxation
is only safe if the raster branch actually validates something, so these tests
assert that a mislabelled or undersized file still fails.
"""

import struct
import zlib
from pathlib import Path

import pytest

from tools.figure_formats import MIN_RASTER_WIDTH, figure_error, png_dimensions

SVG = '<svg xmlns="http://www.w3.org/2000/svg" width="10" height="10"></svg>'


def write_png(path: Path, width: int, height: int) -> Path:
    """Write the smallest PNG that carries a valid signature and IHDR chunk."""

    ihdr = struct.pack(">IIBBBBB", width, height, 8, 6, 0, 0, 0)
    chunk = struct.pack(">I", len(ihdr)) + b"IHDR" + ihdr
    chunk += struct.pack(">I", zlib.crc32(b"IHDR" + ihdr))
    path.write_bytes(b"\x89PNG\r\n\x1a\n" + chunk)
    return path


def test_valid_svg_passes(tmp_path: Path) -> None:
    path = tmp_path / "chart.svg"
    path.write_text(SVG)
    assert figure_error(path) is None


def test_malformed_svg_fails(tmp_path: Path) -> None:
    path = tmp_path / "chart.svg"
    path.write_text("<svg><unclosed>")
    assert "invalid SVG" in figure_error(path)


def test_valid_png_passes(tmp_path: Path) -> None:
    path = write_png(tmp_path / "screen.png", MIN_RASTER_WIDTH, 1200)
    assert figure_error(path) is None
    assert png_dimensions(path) == (MIN_RASTER_WIDTH, 1200)


def test_png_extension_on_a_non_png_fails(tmp_path: Path) -> None:
    path = tmp_path / "screen.png"
    path.write_text(SVG)
    assert "invalid PNG" in figure_error(path)


def test_undersized_png_fails(tmp_path: Path) -> None:
    path = write_png(tmp_path / "screen.png", MIN_RASTER_WIDTH - 1, 400)
    assert "minimum for legible text" in figure_error(path)


def test_unsupported_format_fails(tmp_path: Path) -> None:
    path = tmp_path / "screen.jpg"
    path.write_bytes(b"\xff\xd8\xff")
    assert "unsupported figure format" in figure_error(path)


def test_png_dimensions_rejects_a_non_png(tmp_path: Path) -> None:
    path = tmp_path / "chart.svg"
    path.write_text(SVG)
    with pytest.raises(ValueError):
        png_dimensions(path)
