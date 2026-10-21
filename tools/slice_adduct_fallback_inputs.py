#!/usr/bin/env python3
"""Slice feature rows and hypotheses for the adduct-derived fallback IDs."""

from __future__ import annotations

import json
import sys
from pathlib import Path

from prepare_ppm10_recovery_inputs import assign_inputs, load_input_ids, write_inputs


def main() -> int:
    ids_path, feature_dir, recovered_dir, out_dir = map(Path, sys.argv[1:5])
    ids = [line.strip() for line in ids_path.read_text(encoding="utf-8").split() if line.strip()]
    assignments = assign_inputs(ids, *load_input_ids(feature_dir, recovered_dir))
    out_dir.mkdir(parents=True, exist_ok=True)
    used = write_inputs(assignments, feature_dir, recovered_dir, out_dir)
    payload = {
        "run_libs": ",".join(used),
        "selected": len(assignments),
        "no_inputs": len(ids) - len(assignments),
        "assignments": assignments,
    }
    (out_dir / "selection.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print("[slice]", json.dumps({k: payload[k] for k in ("run_libs", "selected", "no_inputs")}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
