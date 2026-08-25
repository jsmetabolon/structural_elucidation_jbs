# Retention screen for reported top-10 candidates

Status: implemented demote-only. A preliminary measurement on committed data
reads GO (below); the deciding top-10 audit is pending pipeline outputs (see
"Running the audit").

## What exists now

The calibrated retention axis (fitted by `tools/audit_rt_feasibility.py`,
constants in `docs/paper/data/retention/ri_calibration.json`) has until now
been applied to the rank-1 candidate only, inside the confidence model. This
change extends the same standardized residual to every reported top-10 row:

- `structural_elucidation/reports/compound.py::candidate_retention_flags`
  attaches `retention_z` and `retention_implausible` (|z| >= 3.0, the
  `RETENTION_IMPLAUSIBLE_Z` bar) to each candidate row. The columns ride into
  the `*_top_candidates.csv` sidecar and the candidate score summary.
- The internal report lists flagged candidates last (candidate picker and
  per-candidate spectral fit table), marks them "retention-implausible
  (z=...)", and says why in a note. Ranks, scores, list membership, the MSI
  level and the public variant are untouched.
- `tools/audit_topk_retention.py` measures, on real reranked outputs, how often
  the flag fires, whom it hits, and whether it ever hits a `label_ge_0p7`
  candidate (one the truth set marks essentially correct).

## Why demote-only, and why these constants

The feasibility audit (`docs/paper/data/retention/audit_rt_feasibility.json`)
bounds what retention can do here: recalibrated decoy AUC 0.68 (lib 209) and
0.62 (lib 400) with residual RMSE ~520-600 index units, no usable axis on lib
305 (the fitted slope is negative) or lib 402 (two labeled controls). That is
enough to notice a candidate that could not have eluted where the feature did,
and not enough to argue that one did - so the screen withdraws display
priority and never supplies support, deletes nothing, and stays out of the
reranker (`plans/isomer_discrimination.plan.md` already evaluated and rejected
retention as a ranking signal).

The 3.0-sigma bar is not adjustable downward: the worst *true* structure in
the feasibility cohort sits at z ~ 2.9 (X-25273, lib 400), so any lower bar
starts rejecting correct answers.

## The go / no-go gate

`tools/audit_topk_retention.py` writes
`docs/paper/data/retention/audit_topk_retention.json` with a verdict:

- **GO** - keep the internal-report flag - only if the screen fires on >= 5% of
  scoreable controls on some library **and** flags zero labeled-correct rows.
- **NO_GO** otherwise: a screen that fires almost never is not worth a
  rendering surface, and one that fires on known-correct answers is worse than
  none. Record the negative result here and remove the render wiring.

`tests/test_topk_retention.py` guards the committed audit's internal
consistency once it exists (it skips until then).

## Preliminary measurement on committed data

`tools/audit_integrated_best_retention.py` runs the same screen on what a
plain checkout does carry: the delivered rank-1 and best-in-top-10 candidates
per truth-verified control (`data/integrated_best/`), joined to observed
retention indices from `data/feature_rows/`. Committed result
(`docs/paper/data/retention/audit_integrated_best_retention.json`):

- **Library 209** (51 scoreable controls): 4 rank-1 candidates flag (7.8%),
  and all four are wrong against the known structure (Tanimoto 0.0-0.14:
  X-15136 z=3.7, X-26094 z=3.5, X-26218 z=8.0, X-26517 z=7.2). No correct
  candidate flags anywhere. One full improvement case: X-26094's wrong rank-1
  is demoted while the correct candidate sitting at rank 7 (Tanimoto 0.71,
  z=1.9) passes the screen, so the internal list surfaces it.
- **Library 400** (30 scoreable controls): the screen never fires
  (max |z| = 2.5), consistent with that library's tighter calibration
  headroom.

Reading: preliminary **GO** - already past the 5% bar on library 209 with zero
false kills while scoring only two candidates per control, so the full-list
audit can only raise the fire rate. It remains a lower bound on the
truth-verified cohort, not the deciding measurement.

## Running the audit

The audit needs the reranked outputs, which this checkout does not carry.
On a machine with the fetched assets (see `fetch_assets.sh`):

```bash
# Regenerate outputs for the calibrated libraries (long-running):
START_AT=verifier RUN_LIBS=209,400 ./run_inference.sh

# Score the reported top-10s and write the audit + verdict:
python3 tools/audit_topk_retention.py \
  --rerank-root output/pipeline_equalization --libs 209,400
```

Commit the resulting `audit_topk_retention.json` / `.md` beside the other
retention artifacts. If the verdict is NO_GO, revert the render wiring and
keep the audit as the written negative result.

## Interpretation notes

- The observed retention index belongs to the feature, not the candidate: it is
  identical across a control's ten rows (parsed from the `feature_key`), so all
  within-list discrimination comes from the per-candidate proxy
  (`ri_proxy_from_smiles`: exact mass, MolLogP, ring count, TPSA).
- Rows without evidence - no index in the feature key (`_ri0_` or a
  `unified_rerank_*` fallback key), no computable proxy, or an unqualified
  library - are reported as unavailable, never as a zero score.
- `label_ge_0p7` exists only for truth-verified controls, so the false-kill
  count is bounded by that cohort, not the full X set.
