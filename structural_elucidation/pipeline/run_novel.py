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
import hashlib
import json
import logging
import os
import shlex
import subprocess
import sys
from contextlib import contextmanager
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

import pandas as pd

from structural_elucidation.adapters.untargeted_inputs import (
    AdapterManifest,
    ScopeInputs,
    build_inputs_from_handoff,
    feature_id_from_control_id,
    load_handoff,
)
from structural_elucidation.paths import REPO_ROOT

LOGGER = logging.getLogger("structural_elucidation.run_novel")

ALL_STAGES = ("candidates", "verify", "pool", "rerank", "reports")
# The verifier exits 2 when its acceptance diagnostic finds no known winners,
# which is expected for novel-only runs; its scored output is still written.
VERIFY_OK_RETURNCODES = frozenset({0, 2})
SCORE_COLUMNS = (
    "se_msi_level", "se_msi_label", "se_msi_rationale",
    "se_confidence", "se_confidence_tier", "se_evidence_score",
    "se_reported_level", "se_reported_identification",
    "se_candidate_identity_confidence", "se_confidence_flags_json",
)


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

    def __init__(self, repo_root: Path, dry_run: bool = False, assets_root: Path | None = None) -> None:
        self.repo_root = Path(repo_root)
        self.dry_run = dry_run
        self.env = dict(os.environ)
        existing = self.env.get("PYTHONPATH", "")
        self.env["PYTHONPATH"] = (
            f"{self.repo_root}{os.pathsep}{existing}" if existing else str(self.repo_root)
        )
        if assets_root is not None:
            self.env["NCE_REFERENCE_DIR"] = str(assets_root / "vendor" / "modelfiles")

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


def _clear_outputs(runner: StageRunner, *paths: Path) -> None:
    """A successful stage must produce its own outputs, not reuse stale files."""
    if not runner.dry_run:
        for path in paths:
            path.unlink(missing_ok=True)


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
    merged = output_dir / f"{prefix}_merged_predictions.csv"
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
        _clear_outputs(runner, merged)
        runner.run("structural_elucidation.inference.merged", args)
    except (RuntimeError, OSError) as exc:
        LOGGER.warning("[lib%s] de-novo inference failed: %s", scope.lib_id, exc)
        return None
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
    la_metrics = _candidate_metrics_path(la_root, lib_id, ion)
    scored = verifier_root / "scored_candidates.csv"
    pool_path = pool_root / "unified_candidates_per_control.csv"
    non_s0_path = pool_root / "unified_candidates_per_control_non_s0.csv"
    top10 = rerank_root / "unified_rerank_all_top10.csv"

    if enable_denovo:
        merged = run_denovo(runner, scope, out_root, assets_root, denovo_extra)
        if merged is not None:
            candidate_inputs.append(str(merged))

    try:
        if "candidates" in stages:
            _clear_outputs(runner, la_metrics)
            runner.run(
                "structural_elucidation.pipeline.candidates",
                [
                    "--output-root", la_root,
                    "--diagnosis-csv", scope.gap_manifest,
                    "--feature-rows", scope.feature_rows,
                    "--neutral-mass-hypotheses", scope.neutral_mass_hypotheses,
                    "--hmdb-csv", assets_root / "vendor/nce/elucidation/table_builders/hmdb_minimal.csv",
                    "--hmdb-extended-csv", assets_root / "vendor/nce/elucidation/table_builders/HMDB_36_classyfire_21_annotations.csv",
                    "--nce-full-smiles", assets_root / "data/nce_full_smiles_candidates.csv",
                    "--lib-id", lib_id,
                    "--ion-mode", ion,
                ],
            )
            result.stages_ok["candidates"] = True
        if la_metrics.exists() or runner.dry_run:
            candidate_inputs.append(str(la_metrics))

        if "verify" in stages and not candidate_inputs:
            result.error = "no candidate inputs produced (candidates + de-novo empty)"
            LOGGER.warning("[lib%s] %s; skipping verify/pool/rerank", lib_id, result.error)
            return result

        if "verify" in stages:
            _clear_outputs(runner, scored)
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
            if not scored.exists() and not runner.dry_run:
                result.error = f"verify produced no scored_candidates.csv (rc={rc})"
                LOGGER.warning("[lib%s] %s", lib_id, result.error)
                return result
            result.stages_ok["verify"] = True

        if "pool" in stages:
            _clear_outputs(runner, pool_path, non_s0_path)
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
            _clear_outputs(runner, top10)
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
            if top10.exists() or runner.dry_run:
                result.stages_ok["rerank"] = True
            else:
                raise RuntimeError("rerank produced no unified_rerank_all_top10.csv")
        if top10.exists() and not runner.dry_run:
            result.rerank_top10 = top10
    except (RuntimeError, OSError) as exc:
        result.error = str(exc)
        LOGGER.error("[lib%s] stage failure: %s", lib_id, exc)
    return result


def _prediction_groups(manifest: AdapterManifest, results: list[ScopeResult], top_k: int):
    """Yield only current handoff features with real ranked structures."""
    if not results:
        return
    feature_map = pd.read_csv(manifest.novel_feature_map, keep_default_na=False)
    by_control = {str(row["control_id"]): row for row in feature_map.to_dict(orient="records")}
    for result in results:
        if result.rerank_top10 is None:
            continue
        try:
            frame = pd.read_csv(result.rerank_top10, keep_default_na=False, low_memory=False)
        except (pd.errors.EmptyDataError, FileNotFoundError):
            continue
        if frame.empty or "control_id" not in frame:
            continue
        smiles_col = next((col for col in ("predicted_smiles", "canonical_2d_smiles", "smiles") if col in frame), None)
        if smiles_col is None:
            continue
        present = ~frame[smiles_col].astype(str).str.strip().str.lower().isin({"", "nan", "none", "null", "<na>"})
        frame = frame[present].copy()
        if "rank" in frame:
            frame["rank"] = pd.to_numeric(frame["rank"], errors="coerce")
            frame = frame.sort_values(["control_id", "rank"], kind="stable")
        for control_id, group in frame.groupby("control_id", sort=True):
            control_id = str(control_id)
            meta = by_control.get(control_id)
            if meta is None or int(meta["lib_id"]) != result.lib_id:
                continue
            yield result, control_id, meta, group.head(top_k), smiles_col


def version_metadata() -> dict[str, Any]:
    """Record the source checkout and the exact shipped calibration."""
    calibration = REPO_ROOT / "tools" / "identity_confidence_calibration.json"
    metadata: dict[str, Any] = {
        "summary_schema_version": 2,
        "companion_commit": None,
        "companion_dirty": None,
        "calibration_path": str(calibration),
        "calibration_checksum": None,
    }
    if calibration.is_file():
        try:
            metadata["calibration_checksum"] = hashlib.sha256(calibration.read_bytes()).hexdigest()
        except OSError as exc:
            metadata["calibration_error"] = f"Could not read calibration checksum: {exc}"
    try:
        metadata["companion_commit"] = subprocess.check_output(
            ["git", "-C", str(REPO_ROOT), "rev-parse", "HEAD"], text=True, stderr=subprocess.DEVNULL
        ).strip()
        metadata["companion_dirty"] = bool(subprocess.check_output(
            ["git", "-C", str(REPO_ROOT), "status", "--porcelain"], text=True, stderr=subprocess.DEVNULL
        ).strip())
    except (OSError, subprocess.CalledProcessError):
        pass
    return metadata


@contextmanager
def _nce_assets(reference_dir: str):
    """Use the same asset root for in-process scores and stage subprocesses."""
    previous = os.environ.get("NCE_REFERENCE_DIR")
    os.environ["NCE_REFERENCE_DIR"] = reference_dir
    matcher_available = None
    try:
        from tools.rescore_pure_nce import matcher_available

        matcher_available.cache_clear()
        yield
    finally:
        if matcher_available is not None:
            matcher_available.cache_clear()
        if previous is None:
            os.environ.pop("NCE_REFERENCE_DIR", None)
        else:
            os.environ["NCE_REFERENCE_DIR"] = previous


def score_predictions(
    manifest: AdapterManifest,
    results: list[ScopeResult],
    handoff_path: Path,
    scores_dir: Path,
    assets_root: Path,
    *,
    top_k: int,
    persisted_ids: set[str] | None = None,
) -> tuple[dict[str, dict[str, Any]], list[dict[str, Any]]]:
    """Persist scores before optional reports; isolate failures to each feature."""
    from structural_elucidation.reports.novel_scoring import score_novel_feature

    metadata = {}
    for feature in load_handoff(handoff_path):
        metadata[feature.control_id] = {
            **asdict(feature), "control_id": feature.control_id,
            "neutral_mass": feature.effective_neutral_mass, "feature_key": feature.feature_key,
        }
    scores_dir.mkdir(parents=True, exist_ok=True)
    scores: dict[str, dict[str, Any]] = {}
    errors: list[dict[str, Any]] = []
    db = _asset(assets_root, "vendor/metabolite_synonyms.sqlite")
    ledger = _asset(assets_root, "data/curation_ledger.csv")
    for result, control_id, _, group, _ in _prediction_groups(manifest, results, top_k):
        try:
            payload = score_novel_feature(
                group.to_dict(orient="records"), metadata[control_id], db_path=db, curation_ledger=ledger,
            )
            serialized = json.dumps(payload, indent=2, allow_nan=False)
        except Exception as exc:
            message = f"{type(exc).__name__}: {exc}"
            LOGGER.warning("Scoring %s failed: %s", control_id, message)
            payload = {
                "schema_version": 1, "control_id": control_id, "status": "error",
                "summary": {col: None for col in SCORE_COLUMNS}, "errors": [message],
            }
            payload["summary"]["se_confidence_flags_json"] = json.dumps([f"Scoring unavailable: {message}"])
            serialized = json.dumps(payload, indent=2)
        scores[control_id] = payload
        try:
            (scores_dir / f"{control_id}.json").write_text(serialized, encoding="utf-8")
            if persisted_ids is not None and payload.get("status") == "scored":
                persisted_ids.add(control_id)
        except OSError as exc:
            message = f"Could not persist scoring JSON: {exc}"
            payload.setdefault("errors", []).append(message)
            flags = json.loads(payload["summary"].get("se_confidence_flags_json") or "[]")
            payload["summary"]["se_confidence_flags_json"] = json.dumps([*flags, message])
        for error in payload.get("errors", []):
            errors.append({"control_id": control_id, "lib_id": result.lib_id, "error": error})
    return scores, errors


def build_summary(
    manifest: AdapterManifest,
    results: list[ScopeResult],
    reports_dir: Path,
    out_path: Path,
    *,
    top_k: int,
    scores: dict[str, dict[str, Any]] | None = None,
    report_paths: dict[str, str] | None = None,
) -> pd.DataFrame:
    """Join current predictions and scores back to originating feature ids."""
    scores = scores or {}
    report_paths = report_paths or {}
    rows = []
    for result, control_id, meta, group, smiles_col in _prediction_groups(manifest, results, top_k):
        feature_id = feature_id_from_control_id(control_id)
        if feature_id is None:
            continue
        score_col = next((col for col in ("rerank_score", "unified_final_score", "unified_oof_score") if col in group), None)
        score = pd.to_numeric(group.iloc[0][score_col], errors="coerce") if score_col else None
        score_summary = scores.get(control_id, {}).get("summary", {})
        rows.append({
            "feature_id": feature_id, "control_id": control_id, "lib_id": result.lib_id,
            "ion_mode": result.ion_mode, "se_n_candidates": len(group),
            "se_rank1_smiles": str(group.iloc[0][smiles_col]), "se_rank1_score": score,
            "se_topk_smiles_json": json.dumps(group[smiles_col].astype(str).tolist()),
            "se_report_path": report_paths.get(control_id, ""),
            "prior_match_level": meta.get("prior_match_level", ""),
            "prior_matched_name": meta.get("prior_matched_name", ""),
            **{col: score_summary.get(col) for col in SCORE_COLUMNS},
        })
    summary = pd.DataFrame(rows, columns=[
        "feature_id", "control_id", "lib_id", "ion_mode", "se_n_candidates",
        "se_rank1_smiles", "se_rank1_score", "se_topk_smiles_json", "se_report_path",
        "prior_match_level", "prior_matched_name", *SCORE_COLUMNS,
    ])
    summary["se_msi_level"] = pd.to_numeric(summary["se_msi_level"], errors="coerce").astype("Int64")
    summary["se_reported_level"] = summary["se_reported_level"].astype("string")
    summary.to_csv(out_path, index=False)
    return summary


def run(args: argparse.Namespace) -> dict[str, Any]:
    """Execute novel-feature elucidation and publish scores before rendering."""
    stages = tuple(s.strip() for s in str(args.stages).split(",") if s.strip())
    if set(stages) - set(ALL_STAGES):
        raise ValueError(f"Unknown stages: {sorted(set(stages) - set(ALL_STAGES))}")
    if int(args.top_k) < 1 or int(args.jobs) < 1:
        raise ValueError("--top-k and --jobs must be positive")
    out_root = Path(args.out_root).expanduser().resolve()
    out_root.mkdir(parents=True, exist_ok=True)
    inputs_dir = Path(args.inputs_dir).expanduser().resolve() if args.inputs_dir else out_root / "inputs"
    assets_root = Path(args.assets_root).expanduser().resolve() if args.assets_root else REPO_ROOT
    reports_dir, scores_dir = out_root / "reports", out_root / "scores"
    manifest = build_inputs_from_handoff(args.handoff, inputs_dir, write_encoder_table=bool(args.enable_denovo))
    runner = StageRunner(REPO_ROOT, dry_run=bool(args.dry_run), assets_root=assets_root)
    wanted_libs = {int(x) for x in str(args.libs).split(",") if x.strip()} if args.libs else None
    results = []
    for scope in manifest.scopes:
        if wanted_libs is not None and scope.lib_id not in wanted_libs:
            continue
        results.append(run_scope(
            runner, scope, manifest, out_root, assets_root, stages=stages,
            enable_denovo=bool(args.enable_denovo), denovo_extra=shlex.split(args.denovo_extra or ""),
            top_k=int(args.top_k),
        ))

    if args.dry_run and "reports" in stages and results:
        runner.run("structural_elucidation.reports.batch", [
            "--equalization-root", out_root, "--out-dir", reports_dir,
            "--observed-lookup-dir", inputs_dir, "--scoring-dir", scores_dir, "--jobs", args.jobs,
        ])

    scores: dict[str, dict[str, Any]] = {}
    persisted_ids: set[str] = set()
    scoring_errors: list[dict[str, Any]] = []
    report_errors: list[dict[str, Any]] = []
    report_paths: dict[str, str] = {}
    if not args.dry_run and any(r.rerank_top10 for r in results):
        try:
            with _nce_assets(runner.env["NCE_REFERENCE_DIR"]):
                scores, scoring_errors = score_predictions(
                    manifest, results, Path(args.handoff), scores_dir, assets_root,
                    top_k=int(args.top_k), persisted_ids=persisted_ids,
                )
        except Exception as exc:
            scoring_errors.append({"error": f"{type(exc).__name__}: {exc}"})
            LOGGER.warning("Scoring unavailable: %s", exc)

    if "reports" in stages and not args.dry_run and persisted_ids:
        batch_summary = reports_dir / "batch_summary.json"
        _clear_outputs(runner, batch_summary)
        try:
            report_args = [
                "--predictions", *[str(r.rerank_top10) for r in results if r.rerank_top10],
                "--out-dir", reports_dir, "--observed-lookup-dir", inputs_dir,
                "--scoring-dir", scores_dir, "--control-ids", ",".join(sorted(persisted_ids)),
                "--overwrite", "--jobs", args.jobs,
            ]
            for relative, flag in (("vendor/metabolite_synonyms.sqlite", "--db"), ("data/known_X.csv", "--known-x"),
                                   ("data/known_X_all_available_structures_simple.csv", "--known-x-all")):
                found = _asset(assets_root, relative)
                if found:
                    report_args.extend([flag, found])
            runner.run("structural_elucidation.reports.batch", report_args)
            if not batch_summary.exists():
                raise RuntimeError("reports produced no batch_summary.json")
            report_errors.extend(json.loads(batch_summary.read_text()).get("failed", []))
            failed = {str(error.get("control_id")) for error in report_errors}
            for control_id in persisted_ids:
                index = reports_dir / control_id / "index.html"
                if index.exists() and control_id not in failed:
                    report_paths[control_id] = str(index)
        except (RuntimeError, OSError, ValueError) as exc:
            report_errors.append({"error": str(exc)})
            LOGGER.warning("Report generation failed (scores retained): %s", exc)

    summary_path = out_root / "novel_elucidation_summary.csv"
    summary = build_summary(manifest, results, reports_dir, summary_path, top_k=int(args.top_k),
                            scores=scores, report_paths=report_paths)
    provenance = version_metadata()
    if provenance.get("calibration_error"):
        scoring_errors.append({"stage": "calibration", "error": provenance["calibration_error"]})
    payload = {
        **provenance, "out_root": str(out_root), "inputs_dir": str(inputs_dir),
        "assets_root": str(assets_root), "summary": str(summary_path), "scores_dir": str(scores_dir),
        "n_feature_summaries": len(summary), "scoring_errors": scoring_errors, "report_errors": report_errors,
        "scopes": [{**asdict(result), "rerank_top10": str(result.rerank_top10) if result.rerank_top10 else None}
                   for result in results],
    }
    (out_root / "run_novel_summary.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    LOGGER.info("Wrote %s (%s features)", summary_path, len(summary))
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
