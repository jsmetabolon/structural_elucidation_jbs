"""Repository path resolution helpers.

Anchors all relative CLI paths to the ``spec2mol_v21_inference`` package root
so scripts run correctly regardless of the caller's working directory.

Pipeline context:
    Imported by every pipeline CLI (candidates, verify, pool, rerank,
    enrichment) to resolve ``--output-root`` and input CSV paths.

Key I/O:
    ``REPO_ROOT``: absolute path to the inference package root.
    ``resolve()``: joins relative paths against ``REPO_ROOT``.
"""

from __future__ import annotations

from pathlib import Path

# ── Repository anchor ────────────────────────────────────────────────
REPO_ROOT: Path = Path(__file__).resolve().parents[1]  # spec2mol_v21_inference/


def resolve(path: Path, repo_root: Path | None = None) -> Path:
    """Resolve a path relative to the repository root.

    Args:
        path: User-supplied path (may be relative, absolute, or ``~``).
        repo_root: Override anchor; defaults to ``REPO_ROOT``.

    Returns:
        Absolute path. Relative inputs are joined under ``repo_root``.
    """
    root = REPO_ROOT if repo_root is None else repo_root
    path = path.expanduser()
    return path if path.is_absolute() else root / path
