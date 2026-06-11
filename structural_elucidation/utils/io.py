"""Shared table and file IO helpers."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


def read_table(path: Path) -> pd.DataFrame:
    try:
        df = pd.read_csv(path, sep=None, engine="python")
    except Exception:
        df = pd.read_csv(path)
    df.columns = [str(col).strip().strip('"').lstrip("\ufeff") for col in df.columns]
    return df


def pick_col(df: pd.DataFrame, options: list[str]) -> str | None:
    lookup = {str(col).strip().lower(): str(col) for col in df.columns}
    for option in options:
        found = lookup.get(str(option).strip().lower())
        if found is not None:
            return found
    return None


def dataframe_to_markdown(df: pd.DataFrame) -> str:
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
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
