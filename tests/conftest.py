"""Pytest configuration for structural_elucidation package tests.

Ensures the repository root is on ``sys.path`` so tests can import the package
without an editable install — mirrors how smoke scripts run from the repo root.
"""

from __future__ import annotations

import sys
from pathlib import Path

# ── Path bootstrap ────────────────────────────────────────────────────────────

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
