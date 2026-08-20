---
name: Isomer discrimination
overview: 'Correct the failure where the right structure sits below rank 1 because every candidate fits the same mass equally well. Separates the reporting fixes that are safe now from the ranking fix that needs a new evidence axis and a validation set larger than 92 controls.'
todos:
  - id: render-unavailable-blank
    content: 'Carry evidence_confidence_basis into the report sidecars and render an unavailable score as blank with a footnote, never as 0.0 (the open display half of pipeline_refinement_handoff item 7)'
    status: pending
  - id: recovery-diagnosis-column
    content: 'Add a per-compound recovery_diagnosis to build_exec_summary.py so matched_truth TRUE beside a blank evidence cell reads as a stated finding rather than a contradiction'
    status: pending
  - id: surface-ties
    content: 'Present evidentially tied candidates at equal rank instead of an arbitrary order, and stop rank_frame breaking a true tie on predicted_smiles'
    status: pending
  - id: propagate-feature-mass
    content: 'Propagate the feature-level neutral mass onto every candidate row of a control before the evidence backfill, so retrieval rows stop carrying 0.0 on mass_match_score and evidence_confidence'
    status: pending
  - id: formula-for-every-smiles
    content: 'Assign an RDKit formula to every predicted SMILES so mass_match_score stops defaulting to 0 for rows with a blank formula'
    status: pending
  - id: fix-spectral-axis
    content: 'Pool the predicted-fragment and predicted-loss axes against both observed lists, fix the hydrogen bookkeeping, and remove the 0.5 empty-set default in overlap_score'
    status: pending
  - id: occurrence-prior
    content: 'Prototype a structure-occurrence prior as an isomer-specific reranker feature, measured against the truth controls before it is allowed to move any published number'
    status: pending
  - id: real-fragmenter
    content: 'Scope replacing candidate_fragment_bins with a combinatorial fragmenter that breaks rings and multiple bonds, models hydrogen rearrangement, and weights by intensity'
    status: pending
  - id: holdout-validation
    content: 'Establish the library hold-out set before fitting any new axis; 61 controls with a usable spectrum is not enough to tell a real gain from a fitted one'
    status: pending
isProject: false
---
# Correcting the isomer failure

Status: proposal. The audit that motivates it is implemented and committed
(`tools/audit_isomer_discrimination.py`); none of the corrections below are.

## The question this answers

`docs/paper/data/exec_summary/exec_summary_top10.csv` shows X-23080 with
`matched_truth` **TRUE** at rank 2 and an **empty** `evidence_confidence`. That
looks like a contradiction and is not one. The two columns come from different
places and neither depends on the other:

- `matched_truth` is an audit column, computed against the blinded reference
  structure. X-23080 is a truth control (Metabolon Level 1, reference answer
  4-methylbenzenesulfonate) and every candidate carries a precomputed
  `truth_tanimoto` against that hidden answer. The rank-2 candidate
  `Cc1ccc(cc1)S(=O)(=O)O` *is* 4-methylbenzenesulfonate, Tanimoto 1.0, well over
  the 0.70 threshold, so it is TRUE and the compound rolls up as recovered in
  the top 3.
- `evidence_confidence` is what the pipeline knew at prediction time. That row
  arrived from `ec2_bundle_retrieval_selfies_v1` with no mass fields, so no
  score could be computed and the cell is honestly blank.

So the TRUE is the finding, not the bug: the correct structure was in the list
at rank 2, and a wrong isomer took rank 1. The report is MSI 4 and withholds a
structure, which is the right call on the evidence it had.

This plan is about correcting the underlying failure. The first thing it does is
establish that the obvious correction is the wrong one.

## The finding that reframes it

Run `python3 tools/audit_isomer_discrimination.py`. Measured on the delivered
3,897-report bundle:

- **24 of 92 truth controls have the truth in the top 10 but below rank 1.** (A
  further 24 never produce the truth at all, and 44 already have it at rank 1.
  The two 24s are a coincidence; only the first is this plan's subject.)
- **20 of those 24 (83%) share the rank-1 candidate's molecular formula.** Every
  candidate in the group fits the same measured neutral mass equally well. The
  mass axis, the formula prior and `evidence_confidence` are all degenerate
  across them — frequently identical to six decimal places. No mass gate,
  tolerance change or formula prior can separate them. This is a different
  population from the rank-1 mass mismatches in
  `docs/pipeline_refinement_handoff.md` item 2, and none of that item's fixes
  touch it.
- **X-23080 is the only control whose truth row is unscored while rank 1 is
  scored.** The mechanism the blank cell exposes is real, but it is the operative
  cause for one compound out of 92.
- X-23662 is the mirror image and is worth keeping in mind before attributing
  much power to these features: its rank-1 candidate carries `mass_match_score`
  0.0 and no evidence score at all, and still beat a rank-4 candidate carrying
  0.9997 and 0.9981. The generative score dominates the ordering.

For X-23080 concretely, ranks 1 through 4 are all `C7H8O3S` against an observed
neutral mass of 172.019226, and all fit it to 1.9e-4 Da:

| Rank | SMILES | What it is | Tanimoto |
|---|---|---|---|
| 1 | `COCc1ccc(C(=O)O)s1` | 5-(methoxymethyl)thiophene-2-carboxylic acid | 0.158 |
| 2 | `Cc1ccc(cc1)S(=O)(=O)O` | 4-methylbenzenesulfonate — **the truth** | 1.000 |
| 3 | `Cc1ccccc1S(=O)(=O)O` | 2-methylbenzenesulfonate | 0.407 |
| 4 | `O=S(=O)(O)C1=CCC=CC=C1` | cycloheptatriene sulfonate | 0.242 |

Rank 1 leads rank 2 by 0.0106 of `candidate_confidence`.

## The axis that should have decided it is inert

MS2 is the only shipped axis that can distinguish structures of identical
formula. The spectra exist: all 20 isomer-degenerate controls have an observed
MS2 spectrum, and so do 85 of the 85 truth controls that have a local feature
row. The scorer is not using them.

`nce/evidence.py::forward_spectrum_score` compares predicted fragments against
observed **ions**, and predicted losses against observed **losses**. On **56 of
the 61** truth controls with a local spectrum, that same-axis overlap is zero
for every candidate in the top 10. The feature contributes nothing to those
rankings. Where it is non-zero it is 1 or 2 bins out of 7 to 105 observed peaks.

Three compounding defects, in the order they matter:

1. **The two axes are never crossed.** A fragment of mass `m` and a loss of
   `parent - m` are the same event. For X-23080 the truth's predicted **loss** of
   SO3 (79.9573) lands in exactly the observed **ion** bin 47374 (79.95729, the
   textbook aryl-sulfonate diagnostic ion), and its predicted **fragment**
   toluene (92.0630) lands in exactly the observed **loss** bin 61472
   (92.06302). Two exact agreements, both invisible to the scorer. The rank-1
   thiophene explains neither.
2. **Hydrogen bookkeeping.** `candidate_fragment_bins` calls
   `FragmentOnBonds(addDummies=False)` then sanitizes, so both pieces of a split
   gain an implicit hydrogen. Predicted fragment masses sit about 1 Da off the
   deprotonated species actually observed in negative mode, and the two pieces do
   not sum to the parent.
3. **`overlap_score` returns 0.5 when either set is empty**, so a candidate that
   predicts no fragments outscores one that predicts fragments and misses.

## What was tried and did not work

Do not repeat this. Pooling the axes and allowing a hydrogen shift revives all
56 dead controls — at least one observed peak becomes explainable. It is still
not a ranking fix. Re-ordering the delivered top 10 by five normalizations of
the pooled score, against a delivered baseline of 37 at rank 1 and 52 in the top
3 over the 61 assessable controls:

| Variant | Rank 1 | Top 3 | Δ rank 1 |
|---|---|---|---|
| `raw_count` | 18 | 37 | -19 |
| `recall` | 18 | 37 | -19 |
| `f1_recall_precision` | 18 | 37 | -19 |
| `geometric_mean` | 15 | 36 | -22 |
| `recall_over_log_npredicted` | 15 | 35 | -22 |

Every one loses, and confining the score to a tie-break within a 0.10
`candidate_confidence` gap gains at most one control while moving five the wrong
way. The fragmenter breaks only acyclic single bonds and is blind to
regiochemistry, so it cannot rank isomers even when read on the right axis. For
X-23080 it produces byte-identical predicted bins for the para and the ortho
toluenesulfonate: pooling moves the truth from rank 2 into a three-way tie for
rank 1, which is a class assignment, not an identification.

This is also a warning about method. Five plausible hand-built scorers, tuned
against the only 61 labelled controls available, all looked reasonable and all
lost 19 to 22 rank-1 hits. Anything fitted on this set needs held-out
validation before it is believed.

## What is ruled out

- **A different mass window.** Every candidate in a degenerate group fits the
  same mass. See `plans/mass_policy_impact_measurement.plan.md`, which reaches
  the same conclusion from the mass-mismatch side.
- **Backfilling `evidence_confidence` alone.** Worth doing for correctness, but
  it hands the unscored row the *same* score as its competitors, not a
  distinguishing one. Expect X-23080's rank-2 row to go from blank to ~0.9987,
  matching rank 1 exactly. It will not promote it.
- **Retention index.** `docs/paper/data/retention/audit_rt_feasibility.json`
  measures decoy AUC 0.651 overall and 0.590 for lib 305. For X-23080's own
  library the proxy correlation is *negative* (r = -0.765, n = 13). Not a route.

## The work

### Tier 1 — make the report honest. No retrain, no pipeline run.

This is what actually answers the question the blank cell raises, and none of it
is blocked on anything.

1. **Render an uncomputed score as blank.** Carry `evidence_confidence_basis`
   into the report sidecars and have `tools/patch_candidate_scores.py`,
   `reports/compound.py` and `tools/build_exec_summary.py` show `unavailable`
   rows as blank with a footnote. This is the open display half of
   `docs/pipeline_refinement_handoff.md` item 7. Note `build_exec_summary.py`
   already does this; the report bundle does not.
   - Verify: no delivered sidecar shows `0.0` where the basis is `unavailable`.
2. **State the diagnosis per compound.** Add a `recovery_diagnosis` column to
   `exec_summary_by_compound.csv` with the categories the audit already computes:
   `truth_at_rank_1`, `evidence_degenerate_same_formula`, `truth_row_unscored`,
   `different_formula`, `truth_absent_from_top10`. Reuse
   `tools/audit_isomer_discrimination.py::recovery_breakdown` rather than
   reimplementing it. A reader then sees "the correct structure was present at
   rank 2 and the available evidence could not distinguish it from rank 1"
   instead of a TRUE beside a blank.
   - Verify: the column's counts match `audit_isomer_discrimination.json`.
3. **Stop inventing an order among tied candidates.** On 152 of 3,897 compounds
   (3.9%) rank 1 is tied with at least one lower-ranked candidate on every
   displayed evidence column; X-23583 is a three-way tie in which the truth is
   second. `rerank.py::rank_frame` currently breaks a full tie on
   `predicted_smiles` ascending, which is arbitrary with respect to evidence.
   Emit equal ranks, or a `tied_with_rank` marker, and have the report present
   them as a set. This costs nothing in claims — those compounds are already
   withheld at MSI 3-4 — and removes a false implication.
   - Verify: `n_rank1_tied` in the audit is matched by the number of reports
     showing a tie marker; no tie is broken on a SMILES string.

### Tier 2 — close the plumbing gaps. Small code changes, needs a pool re-run.

Correctness work. It fixes the artefact that prompted the question and will not
move the headline.

4. **Propagate the feature-level neutral mass to every candidate row.** The
   observed neutral mass is a property of the feature, not of the candidate. For
   X-23080 it is 172.019226, sitting in
   `data/feature_rows/lib305_negative_feature_rows.csv::neutral_mass_fallback`
   and again in `lib305_negative_neutral_mass_hypotheses.csv::neutral_mass`.
   `pipeline/candidates.py::mass_hypotheses` already unions these into a
   per-control target table. Join it onto the pool by `control_id` before
   `chemistry/evidence.py::backfill_evidence_confidence` runs, so retrieval rows
   stop reaching the reranker with `mass_match_score` 0.0 and no evidence.
   - **Do not remove the guard in `evidence.py`.** It exists because
     `formula_prior_score` falls back to element-ratio plausibility and returns
     0.65 for any well-formed formula when the mass is missing, which would
     manufacture evidence. Using the control's *own measured* mass is not
     manufacturing; inventing one for a control that has none still is. Keep
     `unavailable` for controls with no mass anywhere.
   - Watch for controls with several distinct hypothesis masses (X-23080 has
     three adduct rows agreeing to 3e-4 Da, but that will not hold generally).
     Prefer the observed consensus, and record which was used.
   - Verify: no pool row has `evidence_confidence` 0.0 alongside a non-zero
     `formula_prior_score`; every remaining 0.0 carries basis `unavailable`;
     `n_evidence_backfilled` rises in `unified_candidate_pool_summary.json`.
5. **Give every predicted SMILES a formula.** RDKit, at pool time. This is fix 2
   of `pipeline_refinement_handoff.md` item 2 and closes the 248 blank-formula
   rank-1 rows that default to `mass_match = 0`.

Because `evidence_confidence` and `mass_match_score` are trained reranker
features, both changes require the model to be refit, not merely re-applied.
That means the EC2 path: LOOCV, the G3 gate and the per-source isotonic
calibrators all regenerate. Steps 4 to 7 of `pipeline_refinement_handoff.md`
item 7 already describe that sequence; follow it rather than duplicating it.

### Tier 3 — the actual ranking fix. Needs a new evidence axis.

The mass and formula axes are degenerate for 83% of the failures and the
spectral axis is inert. Nothing already in the reranker can fix this. Two
candidate axes, both isomer-specific, in cost order.

6. **Repair the spectral axis first, because it is cheap and it is wrong
   regardless.** Pool the predicted fragment and loss sets and match them against
   both observed lists; allow an explicit hydrogen shift instead of relying on
   the sanitizer's accidental one; drop the 0.5 empty-set default. This does not
   improve ranking on its own — that is measured, see above — but it is a
   precondition for any better fragmenter, and it makes the axis honest.
   - Verify: `n_same_axis_dead` in the audit falls from 56 toward 0; the sweep is
     re-run and reported, not assumed.
7. **A structure-occurrence prior.** The cheapest genuinely isomer-specific
   signal, and standard practice: SIRIUS/CSI:FingerID's database bonus and
   MetFrag's reference-count score both do this. p-Toluenesulfonate is a common
   xenobiotic; 2-methylbenzenesulfonate is rare; 5-(methoxymethyl)thiophene-2-
   carboxylic acid is not a plausible human metabolite. Unlike the formula prior,
   this term differs *between* isomers of one formula, which is exactly the
   discrimination that is missing. Starting points in-tree are
   `data/smiles_db/comp_id_to_smiles_209.csv` and
   `vendor/metabolite_synonyms.sqlite`; a real implementation wants HMDB, ChEBI
   or PubChem occurrence counts.
   - **The risk is structural, not technical.** This pipeline exists to propose
     structures for compounds that are *not* in the databases. An occurrence
     prior biases against exactly that. It must be a modest additive reranker
     feature whose weight is learned, never a gate, and its effect on novel
     compounds has to be reported separately from its effect on the catalogued
     truth controls.
8. **A real fragmenter.** Replace `candidate_fragment_bins` rather than patching
   it: multi-bond and ring cleavage, bond-energy weighting, explicit hydrogen
   rearrangement, intensity weighting. This is the MetFrag / CFM-ID / SIRIUS
   class of tool. Note the observed side currently stores 10 ppm log-grid bins
   with no intensity, so this also needs the raw peaks propagated — item 3 of
   `pipeline_refinement_handoff.md`, which is where the "MS2 present for 44%"
   report-layer figure comes from. Feature-level availability is much higher
   (85 of 85 truth controls here), so the gap is propagation, not acquisition.

### Tier 4 — validate before shipping any of Tier 3.

9. **Establish the library hold-out set.** 61 truth controls with a usable
   spectrum cannot distinguish a real gain from a fitted one; the failed sweep
   above is the demonstration. `docs/pipeline_refinement_handoff.md` item 4 and
   `docs/library_holdout_validation_scope.md` describe running the reference
   library against itself with each compound suppressed from retrieval, giving
   thousands of paired comparisons. The one number to establish before scoping
   anything further is how many library standards have a usable acquired MS2
   spectrum. Sequence this **before** Tier 3 is fitted, not after.
   - Report any resulting rate as measured on catalogued chemistry, the way
     `tools/audit_pathway_predictions.py --loo` already reports the pathway LOO.
     This work measures the gap at 38.1% exact-source rank-1 recovery against
     56.0% predicted-source, so the two populations are not interchangeable.

## Expected outputs

- Tier 1: a delivered bundle in which no uncomputed score renders as `0.0`, each
  compound states why the truth was or was not recovered, and tied candidates
  are shown as tied. No published rate changes.
- Tier 2: a pool in which no row carries a zero on a feature its identical-formula
  neighbours score 0.99 on. Published recovery figures must be re-measured on
  EC2 and may move in either direction; the change is correctness-motivated.
- Tier 3: the first reranker feature that can separate two structures of the
  same formula. This is where the 20 degenerate controls can actually move.
- Tier 4: a validation set big enough for the Tier 3 question to be answerable.

## Verification

Re-run `python3 tools/audit_isomer_discrimination.py` after each tier and diff
`docs/paper/data/isomer/audit_isomer_discrimination.json`. The numbers that
should move, and the tier that should move them:

| Metric | Now | Moved by |
|---|---|---|
| `recovery.counts.truth_row_unscored` | 1 | Tier 2 |
| `evidence_basis.all_compounds.all_rows.unavailable` | 2,236 | Tier 2 |
| `spectral_axis.n_same_axis_dead` | 56 of 61 | Tier 3 step 6 |
| `recovery.n_same_formula_as_rank1` | 20 of 24 | Tier 3 steps 7-8 only |
| `ties_all_compounds.n_rank1_tied` | 152 of 3,897 | Tier 1 step 3 (presentation) |

Nothing in Tier 1 or Tier 2 should move `n_same_formula_as_rank1`. If a change
claims to, check it has not simply reordered a degenerate group by something
arbitrary.

## Open questions

- Does the occurrence prior's bias against novel structures disqualify it for
  this pipeline's stated purpose, or can a learned low weight contain it? This
  needs a decision before step 7 is built, not after.
- `evidence_confidence` still has three definitions in-tree (`pipeline/candidates.py`
  max-based, `chemistry/conjugates.py` `0.35 + 0.35 * prior`, and a third meaning
  in `chemistry/diagnostics.py`). Unifying them changes existing non-zero values
  and needs the same retrain as Tier 2, so it should ride along with it. Carried
  over from `pipeline_refinement_handoff.md` item 7 open question (a).
- The 24 controls whose truth never reaches the top 10 are a separate failure —
  candidate generation, not ranking — and are out of scope here.
- Should a compound whose top candidates are evidentially tied be capped at MSI 4
  by rule? X-23080 already lands there, but by a different route.
