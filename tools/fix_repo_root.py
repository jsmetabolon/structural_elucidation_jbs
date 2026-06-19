#!/usr/bin/env python3
"""Fix REPO_ROOT references after package migration.

One-shot codemod: replaces per-module ``Path(__file__).parents[N]`` repo-root
guesses with ``from structural_elucidation.paths import REPO_ROOT`` so all
modules resolve paths consistently after the flat ``src/`` refactor.
"""

from __future__ import annotations

import re
from pathlib import Path

# ── Paths ─────────────────────────────────────────────────────────────────────

PKG = Path(__file__).resolve().parents[1] / "structural_elucidation"

# ── Regex patterns ────────────────────────────────────────────────────────────

REPO_ROOT_LINE = re.compile(r"^REPO_ROOT = Path\(__file__\)\.resolve\(\)\.parents\[1\]\s*\n", re.MULTILINE)
REPO_ROOT_MAIN = re.compile(r"repo_root = Path\(__file__\)\.resolve\(\)\.parent")
SRC_ROOT_BLOCK = re.compile(
    r"SRC_ROOT = Path\(__file__\)\.resolve\(\)\.parent\n"
    r"if str\(SRC_ROOT\) not in sys\.path:\n"
    r"    sys\.path\.insert\(0, str\(SRC_ROOT\)\)\n\n",
    re.MULTILINE,
)


def ensure_repo_import(text: str) -> str:
    """Insert ``REPO_ROOT`` import when the file references it but lacks the import."""
    if "from structural_elucidation.paths import REPO_ROOT" in text:
        return text
    if "from structural_elucidation.paths import REPO_ROOT, resolve" in text:
        return text
    if RESOLVE_IMPORT := "from structural_elucidation.paths import resolve" in text:
        return text.replace(
            "from structural_elucidation.paths import resolve",
            "from structural_elucidation.paths import REPO_ROOT, resolve",
        )
    anchor = "from __future__ import annotations\n\n"
    if anchor in text:
        return text.replace(anchor, anchor + "from structural_elucidation.paths import REPO_ROOT\n\n", 1)
    return "from structural_elucidation.paths import REPO_ROOT\n\n" + text


def fix_file(path: Path) -> bool:
    """Rewrite a single Python file; return True when content changed."""
    text = path.read_text(encoding="utf-8")
    original = text
    if path.name == "paths.py":
        return False
    text = REPO_ROOT_LINE.sub("", text)
    text = SRC_ROOT_BLOCK.sub("", text)
    text = REPO_ROOT_MAIN.sub("repo_root = REPO_ROOT", text)
    text = text.replace(
        'vendor_root = Path(__file__).resolve().parents[1] / "vendor"',
        'vendor_root = REPO_ROOT / "vendor"',
    )
    text = text.replace(
        "repo_root=Path(__file__).resolve().parent",
        "repo_root=REPO_ROOT",
    )
    text = text.replace(
        'repo_python = Path(__file__).resolve().parent / ".venv" / "bin" / "python"',
        'repo_python = REPO_ROOT / ".venv" / "bin" / "python"',
    )
    if "REPO_ROOT" in text and path.name != "paths.py":
        text = ensure_repo_import(text)
    if text != original:
        path.write_text(text, encoding="utf-8")
        return True
    return False


def main() -> None:
    """Walk ``structural_elucidation/`` and fix every ``*.py`` file."""
    changed = 0
    for path in PKG.rglob("*.py"):
        if fix_file(path):
            changed += 1
            print(f"fixed {path.relative_to(PKG.parent)}")
    print(f"Updated {changed} files")


if __name__ == "__main__":
    main()
