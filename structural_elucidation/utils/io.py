"""Shared table and file I/O helpers.

Pipeline context:
    Used by evaluation, domain holdout, enrichment, and report builders to
    read heterogeneous CSV exports and write audit artifacts (CSV, JSON,
    markdown tables) with consistent column normalization.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


def read_table(path: Path) -> pd.DataFrame:
    """Read a delimiter-sniffed CSV with normalized column headers.

    Args:
        path: CSV path (comma or other separator).

    Returns:
        DataFrame with stripped, BOM-free column names.
    """
    try:
        df = pd.read_csv(path, sep=None, engine="python")
    except Exception:
        df = pd.read_csv(path)
    df.columns = [str(col).strip().strip('"').lstrip("\ufeff") for col in df.columns]
    return df


def pick_col(df: pd.DataFrame, options: list[str]) -> str | None:
    """Return the first matching column name from ``options`` (case-insensitive).

    Args:
        df: Input table.
        options: Candidate column names in priority order.

    Returns:
        Actual column name from ``df``, or ``None`` when no match exists.
    """
    lookup = {str(col).strip().lower(): str(col) for col in df.columns}
    for option in options:
        found = lookup.get(str(option).strip().lower())
        if found is not None:
            return found
    return None


def dataframe_to_markdown(df: pd.DataFrame) -> str:
    """Render a small DataFrame as a GitHub-flavored markdown table.

    Args:
        df: Table to render.

    Returns:
        Markdown string, or ``_No rows._`` when empty.
    """
    if df.empty:
        return "_No rows._"
    cols = [str(col) for col in df.columns]
    lines = ["| " + " | ".join(cols) + " |", "| " + " | ".join(["---"] * len(cols)) + " |"]
    for _, row in df.iterrows():
        values = []
        for col in df.columns:
            value = row.get(col, "")
            if isinstance(value, float):
                value = round(value, 4) if np.isfinite(value) else ""
            values.append(str(value).replace("|", "\\|"))
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)


def write_csv(path: Path, df: pd.DataFrame) -> None:
    """Write a CSV, creating parent directories as needed.

    Args:
        path: Output CSV path.
        df: Table to persist.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)


def write_json(path: Path, payload: dict[str, Any]) -> None:
    """Write sorted, indented JSON, creating parent directories as needed.

    Args:
        path: Output JSON path.
        payload: Serializable mapping.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
