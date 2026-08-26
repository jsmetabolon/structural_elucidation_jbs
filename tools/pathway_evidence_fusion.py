#!/usr/bin/env python3
"""Evidence fusion for pathway assignment (T4).

Fuses cheap, independent evidence channels with the structure-based pathway
vote, targeted at its weak segment (low-confidence tier: 36.9% super accuracy
over 40% of the evaluable cohort) and at the dominant Lipid <-> Peptide <->
Amino Acid confusions:

- **Peptide detector**: ``peptide_composition_confidence`` from the peptide
  enrichment stage; a confident residue-composition fit routes the compound to
  the Peptide super-pathway regardless of the fingerprint vote.
- **Lipid detector**: a non-empty ``lipid_class`` from the lipid enrichment
  stage routes to Lipid.
- **Neutral-mass prior**: reference-DB metabolites within a ppm window of the
  feature's neutral mass vote a super-pathway; used only to backfill/replace
  weak structure votes (never confident ones).

Precedence: peptide detector > lipid detector > mass prior (weak baselines
only) > structure vote. Every output row records which channel decided.

Inputs and gating
-----------------
- Baseline assignments: the committed audit cohort
  ``docs/paper/data/pathway/audit_pathway_accuracy.csv`` (or any variant CSV
  with the same columns).
- Mass prior (``mass-prior`` subcommand): committed feature m/z
  (``data/feature_rows/``) + the metabolite synonym SQLite (fetched asset).
- Detector columns (optional): fetched
  ``models/pipeline_equalization/rerank_lib*/
  unified_rerank_all_top10_with_pathways_lipid_class_peptide_composition.csv``.

Usage::

    python tools/pathway_evidence_fusion.py mass-prior \
        --out docs/paper/data/pathway/mass_prior_predictions.csv
    python tools/pathway_evidence_fusion.py fuse \
        --mass-prior docs/paper/data/pathway/mass_prior_predictions.csv \
        --out docs/paper/data/pathway/fusion_predictions.csv
    python tools/pathway_variant_benchmark.py score \
        --pred docs/paper/data/pathway/fusion_predictions.csv --label fusion
"""

from __future__ import annotations

import argparse
import bisect
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
BASELINE_AUDIT = REPO_ROOT / "docs" / "paper" / "data" / "pathway" / "audit_pathway_accuracy.csv"
FEATURE_DIR = REPO_ROOT / "data" / "feature_rows"
REFERENCE_DB_CANDIDATES = (
    REPO_ROOT / "Metabolite Synonym Database" / "metabolite_synonyms_merged_with_mtx_giga.sqlite",
    REPO_ROOT / "vendor" / "metabolite_synonyms.sqlite",
)
DEFAULT_LIBS = ("209", "305", "400", "402")
DETECTOR_INPUT_NAME = "unified_rerank_all_top10_with_pathways_lipid_class_peptide_composition.csv"

PPM_DEFAULT = 10.0
MASS_SHARE_THRESHOLD = 0.6
MASS_MIN_CANDIDATES = 2
PEPTIDE_CONFIDENCE_THRESHOLD = 0.6
WEAK_SIMILARITY_THRESHOLD = 0.5


def clean(value: Any) -> str:
    text = str("" if value is None else value).strip()
    return "" if text.lower() in {"nan", "none", "null", "<na>"} else text


def to_float(value: Any) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return float("nan")
    return number


def find_reference_db() -> Path | None:
    for candidate in REFERENCE_DB_CANDIDATES:
        if candidate.exists():
            return candidate
    return None


# ---------------------------------------------------------------------------
# Mass prior channel
# ---------------------------------------------------------------------------

def load_feature_masses(feature_dir: Path = FEATURE_DIR, libs: tuple[str, ...] = DEFAULT_LIBS) -> pd.DataFrame:
    """control_id -> neutral mass, first-lib-wins (audit convention)."""

    frames = []
    for lib in libs:
        for path in sorted(feature_dir.glob(f"lib{lib}_*_feature_rows.csv")):
            frames.append(pd.read_csv(path, low_memory=False)[["control_id", "neutral_mass_fallback"]])
    if not frames:
        raise FileNotFoundError(f"no feature_rows CSVs under {feature_dir}")
    feat = pd.concat(frames, ignore_index=True)
    return feat.drop_duplicates(subset=["control_id"], keep="first").reset_index(drop=True)


def build_mass_table(db_path: Path) -> tuple[list[float], list[tuple[str, str]]]:
    """Sorted (monoisotopic mass, (super, sub)) table from reference SMILES."""

    import sqlite3

    from rdkit import Chem
    from rdkit.Chem.Descriptors import ExactMolWt

    con = sqlite3.connect(str(db_path))
    rows = con.execute(
        "SELECT smiles, metabolon_super, metabolon_sub FROM metabolites "
        "WHERE TRIM(COALESCE(smiles,'')) <> '' AND TRIM(COALESCE(metabolon_super,'')) <> ''"
    ).fetchall()
    con.close()
    entries: list[tuple[float, tuple[str, str]]] = []
    seen: set[str] = set()
    for smiles, super_, sub in rows:
        smiles = str(smiles).split("|")[0].strip()
        if not smiles or smiles in seen:
            continue
        seen.add(smiles)
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue
        entries.append((ExactMolWt(mol), (clean(super_), clean(sub))))
    entries.sort(key=lambda e: e[0])
    return [e[0] for e in entries], [e[1] for e in entries]


def mass_prior_vote(
    neutral_mass: float,
    masses: list[float],
    labels: list[tuple[str, str]],
    ppm: float = PPM_DEFAULT,
) -> dict[str, Any]:
    """Uniform vote of reference metabolites within the ppm window."""

    tol = 1e-6 * ppm * neutral_mass
    lo = bisect.bisect_left(masses, neutral_mass - tol)
    hi = bisect.bisect_right(masses, neutral_mass + tol)
    hits = labels[lo:hi]
    if not hits:
        return {"mass_super": "", "mass_sub": "", "mass_share": 0.0, "mass_n_candidates": 0}
    tally = Counter(super_ for super_, _ in hits)
    winner = sorted(tally.items(), key=lambda kv: (-kv[1], kv[0]))[0][0]
    sub_tally = Counter(sub for super_, sub in hits if super_ == winner and sub)
    return {
        "mass_super": winner,
        "mass_sub": sub_tally.most_common(1)[0][0] if sub_tally else "",
        "mass_share": round(tally[winner] / len(hits), 4),
        "mass_n_candidates": len(hits),
    }


def run_mass_prior(args: argparse.Namespace) -> int:
    db_path = args.reference_db or find_reference_db()
    if db_path is None or not Path(db_path).exists():
        print(
            "[fusion] SKIPPED: metabolite synonym SQLite not present "
            f"(looked for {[str(p) for p in REFERENCE_DB_CANDIDATES]}). "
            "Fetch it with ./fetch_assets.sh and re-run."
        )
        return 0
    features = load_feature_masses()
    masses, labels = build_mass_table(Path(db_path))
    rows = []
    for _, feat in features.iterrows():
        neutral = to_float(feat["neutral_mass_fallback"])
        vote = (
            mass_prior_vote(neutral, masses, labels, args.ppm)
            if neutral == neutral
            else {"mass_super": "", "mass_sub": "", "mass_share": 0.0, "mass_n_candidates": 0}
        )
        rows.append({"control_id": clean(feat["control_id"]), "neutral_mass": neutral, **vote})
    out = pd.DataFrame(rows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False)
    called = int(out["mass_super"].map(clean).ne("").sum())
    print(f"[fusion] mass prior: {called}/{len(out)} controls with a vote -> {args.out}")
    return 0


# ---------------------------------------------------------------------------
# Detector channel (rank-1 rows of the fetched enrichment chain)
# ---------------------------------------------------------------------------

def load_detector_signals(root: Path, libs: tuple[str, ...], input_name: str) -> pd.DataFrame:
    """Rank-1 peptide/lipid detector columns per control, first-lib-wins."""

    frames = []
    seen: set[str] = set()
    for lib in libs:
        path = root / f"rerank_lib{lib}" / input_name
        if not path.exists():
            print(f"[fusion] WARN missing {path}")
            continue
        df = pd.read_csv(path, low_memory=False)
        id_col = "control_id" if "control_id" in df.columns else "x_compound_id"
        df = df.rename(columns={id_col: "control_id"})
        df = df.sort_values("rank", kind="mergesort").drop_duplicates(subset=["control_id"], keep="first")
        keep = [
            c
            for c in (
                "control_id",
                "peptide_composition_label",
                "peptide_composition_confidence",
                "lipid_class",
            )
            if c in df.columns
        ]
        df = df.loc[~df["control_id"].map(clean).isin(seen), keep]
        seen |= set(df["control_id"].map(clean))
        frames.append(df)
    if not frames:
        return pd.DataFrame(columns=["control_id"])
    return pd.concat(frames, ignore_index=True)


# ---------------------------------------------------------------------------
# Fusion policy
# ---------------------------------------------------------------------------

def fuse_one(
    baseline_super: str,
    baseline_sub: str,
    baseline_tier: str,
    baseline_similarity: float,
    *,
    peptide_label: str = "",
    peptide_confidence: float = float("nan"),
    lipid_class: str = "",
    mass_super: str = "",
    mass_sub: str = "",
    mass_share: float = 0.0,
    mass_n_candidates: int = 0,
    peptide_threshold: float = PEPTIDE_CONFIDENCE_THRESHOLD,
    mass_share_threshold: float = MASS_SHARE_THRESHOLD,
    mass_min_candidates: int = MASS_MIN_CANDIDATES,
    weak_similarity: float = WEAK_SIMILARITY_THRESHOLD,
) -> dict[str, str]:
    """One control's fused super/sub and the channel that decided it.

    Precedence: peptide detector > lipid detector > mass prior (weak baseline
    only) > structure vote.
    """

    if clean(peptide_label) and peptide_confidence == peptide_confidence and peptide_confidence >= peptide_threshold:
        sub = baseline_sub if clean(baseline_super) == "Peptide" else ""
        return {"predicted_super": "Peptide", "predicted_sub": sub, "fusion_source": "peptide_detector"}
    if clean(lipid_class):
        sub = baseline_sub if clean(baseline_super) == "Lipid" else ""
        return {"predicted_super": "Lipid", "predicted_sub": sub, "fusion_source": "lipid_detector"}
    baseline_is_weak = clean(baseline_tier).lower() == "low" or (
        baseline_similarity == baseline_similarity and baseline_similarity < weak_similarity
    )
    if (
        baseline_is_weak
        and clean(mass_super)
        and mass_share >= mass_share_threshold
        and mass_n_candidates >= mass_min_candidates
    ):
        return {"predicted_super": clean(mass_super), "predicted_sub": clean(mass_sub), "fusion_source": "mass_prior"}
    return {
        "predicted_super": clean(baseline_super),
        "predicted_sub": clean(baseline_sub),
        "fusion_source": "structure_vote",
    }


def run_fuse(args: argparse.Namespace) -> int:
    baseline = pd.read_csv(args.baseline)
    mass = pd.read_csv(args.mass_prior) if args.mass_prior and args.mass_prior.exists() else None
    if args.mass_prior and mass is None:
        print(f"[fusion] WARN no mass prior at {args.mass_prior}; channel disabled")
    detector_root = (
        args.equalization_root if args.equalization_root.is_absolute() else REPO_ROOT / args.equalization_root
    )
    detectors = (
        load_detector_signals(detector_root, tuple(str(lib) for lib in args.libs), args.detector_input_name)
        if detector_root.exists()
        else None
    )
    if detectors is None:
        print(f"[fusion] detector outputs not present under {detector_root}; detector channels disabled")

    mass_by_control = (
        {clean(r["control_id"]): r for _, r in mass.iterrows()} if mass is not None else {}
    )
    det_by_control = (
        {clean(r["control_id"]): r for _, r in detectors.iterrows()}
        if detectors is not None and not detectors.empty
        else {}
    )

    rows = []
    source_counts: dict[str, int] = defaultdict(int)
    for _, r in baseline.iterrows():
        cid = clean(r["control_id"])
        mass_row = mass_by_control.get(cid, {})
        det_row = det_by_control.get(cid, {})
        fused = fuse_one(
            clean(r.get("predicted_pathway")),
            clean(r.get("predicted_subpathway")),
            clean(r.get("tier")),
            to_float(r.get("best_similarity")),
            peptide_label=clean(det_row.get("peptide_composition_label")),
            peptide_confidence=to_float(det_row.get("peptide_composition_confidence")),
            lipid_class=clean(det_row.get("lipid_class")),
            mass_super=clean(mass_row.get("mass_super")),
            mass_sub=clean(mass_row.get("mass_sub")),
            mass_share=to_float(mass_row.get("mass_share")) if mass_row else 0.0,
            mass_n_candidates=int(to_float(mass_row.get("mass_n_candidates")) or 0) if mass_row else 0,
            peptide_threshold=args.peptide_threshold,
            mass_share_threshold=args.mass_share_threshold,
            mass_min_candidates=args.mass_min_candidates,
            weak_similarity=args.weak_similarity,
        )
        source_counts[fused["fusion_source"]] += 1
        rows.append({"control_id": cid, **fused})
    out = pd.DataFrame(rows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False)
    print(f"[fusion] wrote {len(out)} fused predictions -> {args.out}")
    print(f"[fusion] decision sources: {dict(source_counts)}")
    print(
        "[fusion] score with: python tools/pathway_variant_benchmark.py score "
        f"--pred {args.out} --pred-sub-col predicted_sub --label fusion"
    )
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)

    pm = sub.add_parser("mass-prior", help="neutral-mass reference-DB vote per control (needs DB)")
    pm.add_argument("--reference-db", type=Path, default=None)
    pm.add_argument("--ppm", type=float, default=PPM_DEFAULT)
    pm.add_argument("--out", type=Path, default=REPO_ROOT / "docs/paper/data/pathway/mass_prior_predictions.csv")

    pf = sub.add_parser("fuse", help="fuse detector + mass channels with the structure vote")
    pf.add_argument("--baseline", type=Path, default=BASELINE_AUDIT)
    pf.add_argument("--mass-prior", type=Path, default=REPO_ROOT / "docs/paper/data/pathway/mass_prior_predictions.csv")
    pf.add_argument("--equalization-root", type=Path, default=Path("models/pipeline_equalization"))
    pf.add_argument("--libs", nargs="+", type=int, default=[209, 305, 400, 402])
    pf.add_argument("--detector-input-name", default=DETECTOR_INPUT_NAME)
    pf.add_argument("--peptide-threshold", type=float, default=PEPTIDE_CONFIDENCE_THRESHOLD)
    pf.add_argument("--mass-share-threshold", type=float, default=MASS_SHARE_THRESHOLD)
    pf.add_argument("--mass-min-candidates", type=int, default=MASS_MIN_CANDIDATES)
    pf.add_argument("--weak-similarity", type=float, default=WEAK_SIMILARITY_THRESHOLD)
    pf.add_argument("--out", type=Path, default=REPO_ROOT / "docs/paper/data/pathway/fusion_predictions.csv")

    args = parser.parse_args()
    if args.command == "mass-prior":
        return run_mass_prior(args)
    if args.command == "fuse":
        return run_fuse(args)
    raise AssertionError(args.command)


if __name__ == "__main__":
    raise SystemExit(main())
