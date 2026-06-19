#!/usr/bin/env python3
"""Feature-table provenance audit helpers for RI / source-html lineage.

Reads companion ``*_feature_source_audit.json`` artifacts emitted during
feature-table construction so downstream prediction runs can attach
source-audit metadata without re-parsing HTML inputs.
"""

from __future__ import annotations

import json
from pathlib import Path


def derive_feature_artifact_prefix(feature_table_path: Path) -> str:
    """Strip encoder/meta/peaks suffixes from a feature-table stem.

    Args:
        feature_table_path: Path to ``*_encoder.csv``, ``*_meta.pkl``, etc.

    Returns:
        Shared artifact prefix used for sibling audit JSON files.
    """
    stem = feature_table_path.stem
    for suffix in ("_encoder", "_meta", "_peaks"):
        if stem.endswith(suffix):
            return stem[: -len(suffix)]
    return stem


def load_feature_source_audit_summary(
    feature_table_path: Path,
) -> tuple[Path | None, dict[str, object] | None]:
    """Load optional feature-source audit JSON beside a feature table.

    Args:
        feature_table_path: Encoder-ready feature table path.

    Returns:
        Tuple of ``(audit_json_path, payload)``. Both are ``None`` when no
        audit file exists. ``payload`` is ``None`` when the file is unreadable.
    """
    prefix = derive_feature_artifact_prefix(feature_table_path)
    candidate = feature_table_path.parent / f"{prefix}_feature_source_audit.json"
    if not candidate.exists():
        return None, None

    try:
        payload = json.loads(candidate.read_text())
    except Exception:
        return candidate, None
    if not isinstance(payload, dict):
        return candidate, None
    return candidate, payload
