#!/usr/bin/env python3
"""Write the disjoint reasons eligible compounds are missing from the summary."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from prepare_ppm10_recovery_inputs import catalog_targets

GROUP_SOURCES = (
    ("mass_only", "mass_only"),
    ("no_inputs", "no_inputs"),
    ("nohit_decoded", "nohit10"),
    ("searched10_over", "searched10_over"),
    ("library10_still", "library10"),
)


def load_ids(path: Path) -> set[str]:
    if not path.exists():
        return set()
    return {line.strip() for line in path.read_text(encoding="utf-8").split() if line.strip()}


def disjoint_missing(missing: set[str], lists: dict[str, set[str]]) -> dict[str, list[str]]:
    """Split ``missing`` by the historical ID lists, earliest reason first.

    A compound that appears in more than one source list is kept in the first
    group only, so the five files do not overlap.
    """

    used: set[str] = set()
    groups: dict[str, list[str]] = {}
    for name, source in GROUP_SOURCES:
        group = (missing & lists.get(source, set())) - used
        used |= group
        groups[name] = sorted(group)
    groups["unlisted"] = sorted(missing - used)
    return groups


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--mass-only", type=Path, required=True)
    parser.add_argument("--library10", type=Path, required=True)
    parser.add_argument("--searched10-over", type=Path, required=True)
    parser.add_argument("--nohit10", type=Path, required=True)
    parser.add_argument("--no-inputs", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args(argv)

    kept, _excluded = catalog_targets(pd.read_csv(args.catalog, low_memory=False))
    published = set(
        pd.read_csv(args.summary, usecols=["x_compound_id"], low_memory=False).x_compound_id.astype(str)
    )
    eligible = set(kept)
    missing = eligible - published
    groups = disjoint_missing(
        missing,
        {
            "mass_only": load_ids(args.mass_only),
            "library10": load_ids(args.library10),
            "searched10_over": load_ids(args.searched10_over),
            "nohit10": load_ids(args.nohit10),
            "no_inputs": load_ids(args.no_inputs),
        },
    )
    args.out_dir.mkdir(parents=True, exist_ok=True)
    for name, ids in groups.items():
        (args.out_dir / f"{name}_ids.txt").write_text("\n".join(ids) + ("\n" if ids else ""), encoding="utf-8")
    decode = groups["library10_still"] + groups["searched10_over"]
    (args.out_dir / "decode_ids.txt").write_text("\n".join(decode) + ("\n" if decode else ""), encoding="utf-8")
    payload = {
        "eligible": len(eligible),
        "published_eligible": len(eligible & published),
        "missing": len(missing),
        "counts": {name: len(ids) for name, ids in groups.items()},
        "decode": len(decode),
    }
    (args.out_dir / "gap.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print("[gap]", json.dumps(payload))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
