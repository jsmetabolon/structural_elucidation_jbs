from __future__ import annotations

import json
from pathlib import Path


def derive_feature_artifact_prefix(feature_table_path: Path) -> str:
    stem = feature_table_path.stem
    for suffix in ("_encoder", "_meta", "_peaks"):
        if stem.endswith(suffix):
            return stem[: -len(suffix)]
    return stem


def load_feature_source_audit_summary(feature_table_path: Path) -> tuple[Path | None, dict[str, object] | None]:
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
