"""Run the full elucidation stack on brand-new, un-curated untargeted features.

This driver is the single entry point invoked (in the elucidation venv) by the
``untargeted_workflow`` structural-elucidation step.  Given a neutral handoff
file it:

1. Synthesizes per-library input tables via
   :mod:`structural_elucidation.adapters.untargeted_inputs`.
2. Optionally runs encoder de-novo/retrieval inference
   (:mod:`structural_elucidation.inference.merged`).
3. Runs library-anchored candidates -> verifier -> unified pool -> rerank
   for every ``(lib, ion)`` scope present in the handoff.
4. Optionally builds per-compound HTML reports
   (:mod:`structural_elucidation.reports.batch`).
5. Writes ``novel_elucidation_summary.csv`` mapping each rerank result back to
   the originating untargeted ``feature_id``.

Unlike ``run_inference.sh`` (which reuses pre-built ``data/base_pools``), this
driver is designed for features that have never been through Metabolon curation:
it runs encoder inference itself, tolerates the verifier's non-zero acceptance
exit code, and seeds a synthetic ``all_x_catalog`` so the pool's non-S0 split
surfaces every novel control in ``unified_rerank_all_top10.csv``.

Example::

    python -m structural_elucidation.pipeline.run_novel \\
        --handoff /path/handoff.json \\
        --out-root output/novel_run \\
        --jobs 8
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pandas as pd

from structural_elucidation.adapters.untargeted_inputs import (
    AdapterManifest,
    ScopeInputs,
    build_inputs_from_handoff,
    feature_id_from_control_id,
)
from structural_elucidation.paths import REPO_ROOT

LOGGER = logging.getLogger("structural_elucidation.run_novel")

ALL_STAGES = ("candidates", "verify", "pool", "rerank", "reports")
# The verifier exits 2 when its acceptance diagnostic finds no known winners,
# which is expected for novel-only runs; its scored output is still written.
VERIFY_OK_RETURNCODES = frozenset({0, 2})


@dataclass
class ScopeResult:
    """Per-``(lib, ion)`` execution status."""

    lib_id: int
    ion_mode: str
    n_features: int
    stages_ok: dict[str, bool] = field(default_factory=dict)
    rerank_top10: Path | None = None
    error: str = ""


def _candidate_metrics_path(la_root: Path, lib_id: int, ion_mode: str) -> Path:
    """Canonical ``candidate_metrics_all.csv`` location produced by candidates."""
    return (
        la_root
        / "runs"
        / "library_anchored_v1"
        / "track_control_blind"
        / f"lib_{lib_id}_{ion_mode}"
        / "selfies_profiles"
        / "strict_1da"
        / "candidate_metrics_all.csv"
    )


class StageRunner:
    """Runs elucidation CLI modules as subprocesses in the elucidation venv."""

    def __init__(self, repo_root: Path, dry_run: bool = False) -> None:
        self.repo_root = Path(repo_root)
        self.dry_run = dry_run
        self.env = dict(os.environ)
        existing = self.env.get("PYTHONPATH", "")
        self.env["PYTHONPATH"] = (
            f"{self.repo_root}{os.pathsep}{existing}" if existing else str(self.repo_root)
        )

    def run(self, module: str, args: list[str], *, ok_returncodes=frozenset({0})) -> int:
        """Run ``python -m module args``; return the return code.

        Raises:
            RuntimeError: When the return code is not in ``ok_returncodes``.
        """
        cmd = [sys.executable, "-m", module, *[str(a) for a in args]]
        LOGGER.info("+ %s", " ".join(cmd))
        if self.dry_run:
            return 0
        proc = subprocess.run(
            cmd,
            cwd=str(self.repo_root),
            env=self.env,
            capture_output=True,
            text=True,
        )
        if proc.stdout:
            LOGGER.debug("[%s stdout]\n%s", module, proc.stdout[-4000:])
        if proc.returncode not in ok_returncodes:
            tail = (proc.stderr or proc.stdout or "").strip().splitlines()
            raise RuntimeError(
                f"{module} failed (rc={proc.returncode}): " + "\n".join(tail[-8:])
            )
        return proc.returncode


def _asset(assets_root: Path, rel: str) -> Path | None:
    """Return ``assets_root/rel`` if it exists, else ``None``."""
    path = Path(assets_root) / rel
    return path if path.exists() else None


def run_denovo(
    runner: StageRunner,
    scope: ScopeInputs,
    out_root: Path,
    assets_root: Path,
    extra: list[str],
) -> Path | None:
    """Optionally run encoder de-novo/retrieval inference for one scope.

    Returns the ``*_merged_predictions.csv`` path when the run succeeds and the
    required encoder checkpoint is available, else ``None``.
    """
    encoder = _asset(assets_root, "models/encoder/best_encoder.pt")
    if encoder is None:
        LOGGER.warning(
            "[lib%s] --enable-denovo set but %s is missing; skipping de-novo inference.",
            scope.lib_id,
            Path(assets_root) / "models/encoder/best_encoder.pt",
        )
        return None
    output_dir = out_root / f"denovo_lib{scope.lib_id}"
    prefix = f"novel_lib{scope.lib_id}"
    args: list[str] = [
        "--feature_table", scope.encoder_table,
        "--lib_id", scope.lib_id,
        "--ion_mode", scope.ion_mode,
        "--encoder_ckpt", encoder,
        "--output_dir", output_dir,
        "--artifact_prefix", prefix,
    ]
    for rel, flag in (
        ("models/decoder_weights", "--decoder_dir"),
        ("SMILES_DB/comp_id_to_smiles_expanded_950.csv", "--retrieval_smiles_csv"),
        ("SMILES_DB/name_dict_runtime_combined.csv", "--retrieval_name_csv"),
    ):
        found = _asset(assets_root, rel)
        if found is not None:
            args.extend([flag, str(found)])
    args.extend(extra)
    try:
        runner.run("structural_elucidation.inference.merged", args)
    except RuntimeError as exc:
        LOGGER.warning("[lib%s] de-novo inference failed: %s", scope.lib_id, exc)
        return None
    merged = output_dir / f"{prefix}_merged_predictions.csv"
    return merged if merged.exists() else None


def run_scope(
    runner: StageRunner,
    scope: ScopeInputs,
    manifest: AdapterManifest,
    out_root: Path,
    assets_root: Path,
    *,
    stages: tuple[str, ...],
    enable_denovo: bool,
    denovo_extra: list[str],
    top_k: int,
) -> ScopeResult:
    """Run candidates -> verify -> pool -> rerank for a single ``(lib, ion)``."""
    lib_id, ion = scope.lib_id, scope.ion_mode
    result = ScopeResult(lib_id=lib_id, ion_mode=ion, n_features=scope.n_features)

    la_root = out_root / f"library_anchored_lib{lib_id}"
    verifier_root = out_root / f"verifier_lib{lib_id}"
    pool_root = out_root / f"unified_pool_lib{lib_id}"
    rerank_root = out_root / f"rerank_lib{lib_id}"

    candidate_inputs: list[str] = []

    if enable_denovo:
        merged = run_denovo(runner, scope, out_root, assets_root, denovo_extra)
        if merged is not None:
            candidate_inputs.append(str(merged))

    try:
        if "candidates" in stages:
            runner.run(
                "structural_elucidation.pipeline.candidates",
                [
                    "--output-root", la_root,
                    "--diagnosis-csv", scope.gap_manifest,
                    "--feature-rows", scope.feature_rows,
                    "--neutral-mass-hypotheses", scope.neutral_mass_hypotheses,
                    "--lib-id", lib_id,
                    "--ion-mode", ion,
                ],
            )
            result.stages_ok["candidates"] = True
        la_metrics = _candidate_metrics_path(la_root, lib_id, ion)
        if la_metrics.exists():
            candidate_inputs.append(str(la_metrics))

        if not candidate_inputs:
            result.error = "no candidate inputs produced (candidates + de-novo empty)"
            LOGGER.warning("[lib%s] %s; skipping verify/pool/rerank", lib_id, result.error)
            return result

        if "verify" in stages:
            rc = runner.run(
                "structural_elucidation.pipeline.verify",
                [
                    "--output-root", verifier_root,
                    "--observed-lookup", scope.observed_ms2_lookup,
                    "--all-x-observed-lookup", scope.observed_ms2_lookup,
                    "--current-integrated-best", manifest.integrated_best,
                    "--input-candidate-metrics", ",".join(candidate_inputs),
                    "--all-x-catalog", manifest.all_x_catalog,
                    "--lib-id", lib_id,
                    "--ion-mode", ion,
                    "--include-all-x",
                ],
                ok_returncodes=VERIFY_OK_RETURNCODES,
            )
            scored = verifier_root / "scored_candidates.csv"
            if not scored.exists():
                result.error = f"verify produced no scored_candidates.csv (rc={rc})"
                LOGGER.warning("[lib%s] %s", lib_id, result.error)
                return result
            result.stages_ok["verify"] = True

        if "pool" in stages:
            runner.run(
                "structural_elucidation.pipeline.pool",
                [
                    "--output-root", pool_root,
                    "--scored-candidates", verifier_root / "scored_candidates.csv",
                    "--extra-inputs", "",
                    "--current-integrated-best", manifest.integrated_best,
                    "--all-x-catalog", manifest.all_x_catalog,
                    "--lib-id", lib_id,
                    "--ion-mode", ion,
                    "--include-all-x",
                ],
            )
            result.stages_ok["pool"] = True

        if "rerank" in stages:
            runner.run(
                "structural_elucidation.pipeline.rerank",
                [
                    "--output-root", rerank_root,
                    "--unified-pool", pool_root / "unified_candidates_per_control.csv",
                    "--non-s0-pool", pool_root / "unified_candidates_per_control_non_s0.csv",
                    "--lib-id", lib_id,
                    "--ion-mode", ion,
                    "--top-k", top_k,
                    "--model-selection-mode", "auto",
                ],
            )
            top10 = rerank_root / "unified_rerank_all_top10.csv"
            if top10.exists():
                result.rerank_top10 = top10
                result.stages_ok["rerank"] = True
    except RuntimeError as exc:
        result.error = str(exc)
        LOGGER.error("[lib%s] stage failure: %s", lib_id, exc)
    return result


def build_summary(
    manifest: AdapterManifest,
    results: list[ScopeResult],
    reports_dir: Path,
    out_path: Path,
    *,
    top_k: int,
) -> pd.DataFrame:
    """Join rerank top-N predictions back to originating ``feature_id``s."""
    feature_map = pd.read_csv(manifest.novel_feature_map) if manifest.novel_feature_map else pd.DataFrame()
    map_by_control: dict[str, dict[str, Any]] = {}
    if not feature_map.empty:
        for record in feature_map.to_dict(orient="records"):
            map_by_control[str(record["control_id"])] = record

    rows: list[dict[str, Any]] = []
    for result in results:
        if result.rerank_top10 is None or not Path(result.rerank_top10).exists():
            continue
        try:
            top10 = pd.read_csv(result.rerank_top10, low_memory=False)
        except (pd.errors.EmptyDataError, FileNotFoundError):
            continue
        if top10.empty or "control_id" not in top10.columns:
            continue
        smiles_col = next(
            (c for c in ("predicted_smiles", "canonical_2d_smiles", "smiles") if c in top10.columns),
            None,
        )
        score_col = next(
            (c for c in ("rerank_score", "unified_final_score", "unified_oof_score") if c in top10.columns),
            None,
        )
        if "rank" in top10.columns:
            top10 = top10.sort_values(["control_id", "rank"])
        for control_id, group in top10.groupby("control_id"):
            control_id = str(control_id)
            feature_id = feature_id_from_control_id(control_id)
            if feature_id is None:
                continue
            meta = map_by_control.get(control_id, {})
            head = group.head(top_k)
            top_smiles = [str(s) for s in head[smiles_col].tolist()] if smiles_col else []
            report_index = reports_dir / control_id / "index.html"
            rows.append(
                {
                    "feature_id": feature_id,
                    "control_id": control_id,
                    "lib_id": result.lib_id,
                    "ion_mode": result.ion_mode,
                    "se_n_candidates": int(len(group)),
                    "se_rank1_smiles": top_smiles[0] if top_smiles else "",
                    "se_rank1_score": (
                        float(head.iloc[0][score_col]) if score_col and len(head) else float("nan")
                    ),
                    "se_topk_smiles_json": json.dumps(top_smiles),
                    "se_report_path": str(report_index) if report_index.exists() else "",
                    "prior_match_level": meta.get("prior_match_level", ""),
                    "prior_matched_name": meta.get("prior_matched_name", ""),
                }
            )
    summary = pd.DataFrame(
        rows,
        columns=[
            "feature_id",
            "control_id",
            "lib_id",
            "ion_mode",
            "se_n_candidates",
            "se_rank1_smiles",
            "se_rank1_score",
            "se_topk_smiles_json",
            "se_report_path",
            "prior_match_level",
            "prior_matched_name",
        ],
    )
    summary.to_csv(out_path, index=False)
    return summary


def run(args: argparse.Namespace) -> dict[str, Any]:
    """Execute the novel-feature elucidation pipeline end to end."""
    out_root = Path(args.out_root).resolve()
    out_root.mkdir(parents=True, exist_ok=True)
    inputs_dir = Path(args.inputs_dir).resolve() if args.inputs_dir else (out_root / "inputs")
    assets_root = Path(args.assets_root).resolve() if args.assets_root else REPO_ROOT
    reports_dir = out_root / "reports"
    stages = tuple(s.strip() for s in str(args.stages).split(",") if s.strip())
    denovo_extra = [tok for tok in str(args.denovo_extra or "").split() if tok]

    LOGGER.info("Building elucidation inputs from handoff %s -> %s", args.handoff, inputs_dir)
    manifest = build_inputs_from_handoff(
        args.handoff, inputs_dir, write_encoder_table=bool(args.enable_denovo)
    )
    if not manifest.scopes:
        LOGGER.warning("No supported (lib, ion) scopes found in handoff; nothing to do.")
        summary_path = out_root / "novel_elucidation_summary.csv"
        build_summary(manifest, [], reports_dir, summary_path, top_k=int(args.top_k))
        return {"scopes": [], "summary": str(summary_path)}

    runner = StageRunner(REPO_ROOT, dry_run=bool(args.dry_run))
    wanted_libs = (
        {int(x) for x in str(args.libs).split(",") if x.strip()} if args.libs else None
    )

    results: list[ScopeResult] = []
    for scope in manifest.scopes:
        if wanted_libs is not None and scope.lib_id not in wanted_libs:
            continue
        LOGGER.info(
            "=== lib %s (%s): %s novel feature(s) ===", scope.lib_id, scope.ion_mode, scope.n_features
        )
        results.append(
            run_scope(
                runner,
                scope,
                manifest,
                out_root,
                assets_root,
                stages=stages,
                enable_denovo=bool(args.enable_denovo),
                denovo_extra=denovo_extra,
                top_k=int(args.top_k),
            )
        )

    if "reports" in stages and not args.dry_run and any(r.rerank_top10 for r in results):
        try:
            runner.run(
                "structural_elucidation.reports.batch",
                [
                    "--equalization-root", out_root,
                    "--out-dir", reports_dir,
                    "--observed-lookup-dir", inputs_dir,
                    "--jobs", args.jobs,
                ],
            )
        except RuntimeError as exc:
            LOGGER.warning("report generation failed (non-fatal): %s", exc)

    summary_path = out_root / "novel_elucidation_summary.csv"
    summary = build_summary(manifest, results, reports_dir, summary_path, top_k=int(args.top_k))
    LOGGER.info(
        "Wrote %s (%s feature summaries across %s scope(s))",
        summary_path,
        len(summary),
        len(results),
    )

    payload = {
        "out_root": str(out_root),
        "inputs_dir": str(inputs_dir),
        "summary": str(summary_path),
        "n_feature_summaries": int(len(summary)),
        "scopes": [
            {
                "lib_id": r.lib_id,
                "ion_mode": r.ion_mode,
                "n_features": r.n_features,
                "stages_ok": r.stages_ok,
                "rerank_top10": str(r.rerank_top10) if r.rerank_top10 else None,
                "error": r.error,
            }
            for r in results
        ],
    }
    (out_root / "run_novel_summary.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return payload


def main() -> None:
    """CLI entry point for the novel-feature elucidation driver."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--handoff", type=Path, required=True, help="Handoff JSON from untargeted_workflow.")
    parser.add_argument("--out-root", type=Path, required=True, help="Root directory for all outputs.")
    parser.add_argument(
        "--inputs-dir", type=Path, default=None, help="Where to write synthesized inputs (default: <out-root>/inputs)."
    )
    parser.add_argument(
        "--assets-root",
        type=Path,
        default=None,
        help="Root holding fetched model/reference assets (default: repo root).",
    )
    parser.add_argument(
        "--stages",
        default=",".join(ALL_STAGES),
        help=f"Comma-separated subset of stages to run (default all: {','.join(ALL_STAGES)}).",
    )
    parser.add_argument("--libs", default="", help="Optional comma-separated lib id filter (e.g. 400,402).")
    parser.add_argument("--enable-denovo", action="store_true", help="Run encoder de-novo/retrieval inference.")
    parser.add_argument("--denovo-extra", default="", help="Extra whitespace-separated flags for inference.merged.")
    parser.add_argument("--top-k", type=int, default=10, help="Top-N predictions to keep per feature.")
    parser.add_argument("--jobs", type=int, default=4, help="Parallel report jobs.")
    parser.add_argument("--dry-run", action="store_true", help="Build inputs and print planned commands only.")
    parser.add_argument("--verbose", action="store_true", help="Enable debug logging.")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )
    payload = run(args)
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
