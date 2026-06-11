"""Repository path resolution helpers."""

from __future__ import annotations

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]


def resolve(path: Path, repo_root: Path | None = None) -> Path:
    """Resolve a path relative to the repository root."""
    root = REPO_ROOT if repo_root is None else repo_root
    path = path.expanduser()
    return path if path.is_absolute() else root / path
