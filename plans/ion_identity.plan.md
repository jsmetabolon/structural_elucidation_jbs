---
name: Ion identity before structure (cluster co-eluting ion forms, then score spectra under each adduct hypothesis)
overview: 'Treat the ion form of a precursor as an inferred quantity rather than an assumption. Cluster the co-eluting ion forms already recorded per feature, score each adduct hypothesis against the observed MS2 (diagnostic neutral loss, then fragments of the implied neutral), carry the surviving hypotheses into candidate generation as a set, and report the winner as "compound, observed as [adduct]". Fixes the class of failure where an ammonium adduct is read as [M+H]+ and the true structure is withheld as a 17.0265 Da mass mismatch.'
todos:
  - id: audit-tool
    content: 'Write tools/audit_ion_identity.py (stdlib + numpy) to count the population at risk, parse the per-feature ion-form clusters out of consensus_observed_adducts, and report how often the diagnostic-loss test is decidable; emit to docs/paper/data/ion_identity/ and docs/ion_identity_findings.md'
    status: pending
  - id: ion-form-module
    content: 'Add ion_form_candidates() and the charge-explicit forward/inverse arithmetic to structural_elucidation/chemistry/adducts.py, reusing the grammar promoted out of tools/audit_identification_rules.py; no behaviour change'
    status: pending
  - id: fragment-convention
    content: 'Fix the one-proton offset between predicted fragment bins (neutral masses) and observed ion bins (ion m/z), reconcile the two disagreeing observed-bin producers, and make the charge carrier an explicit argument so a hypothesis can enter the fragment axis'
    status: pending
  - id: ion-form-clusters
    content: 'Build the per-feature co-eluting ion-form cluster from the committed consensus_observed_adducts list plus RI, and emit a cluster table keyed by implied neutral mass'
    status: pending
  - id: hypothesis-scoring
    content: 'Score each ion-form hypothesis per feature from cluster support, diagnostic neutral loss, and fragment agreement with the implied neutral; emit a posterior and an explicit undecidable class'
    status: pending
  - id: candidates-consume
    content: 'Have pipeline/candidates.py::mass_hypotheses emit one target per surviving ion form carrying ion_form and ion_form_posterior, behind a flag, with the affected-feature count reported before it is enabled'
    status: pending
  - id: report-surface
    content: 'Name the ion form in the report identity line, and stop scoring an adduct-explained mass as a mass mismatch in reports/confidence.py; cap any ion-form claim that has neither cluster nor loss evidence'
    status: pending
  - id: upstream-export
    content: 'In untargeted_workflow: carry adduct_type, neutral_mass and group_id into consensus_ms2.mgf instead of PEPMASS alone, and implement the correlation gating at the three existing TODOs so the cluster is trustworthy (separate repo, separate PR)'
    status: pending
  - id: ci
    content: 'Add tools/audit_ion_identity.py to the regenerate-and-diff step and the ruff list of .github/workflows/publication-checks.yml so a regression in any of the above shows up as a diff'
    status: pending
isProject: false
---
# Ion identity before structure

Status: proposal, not started. Nothing below has been implemented. Every count
quoted here was measured from the committed feature rows in this checkout and is
reproducible by the stage 0 tool; until that tool exists the numbers live only in
this document.

## Goal

The pipeline reads every precursor as `[M+H]+` in positive mode and `[M-H]-` in
negative mode (CHG-2 in [docs/identification_rules.md](../docs/identification_rules.md#chg-2--the-neutral-mass-assumes-a-singly-charged-monomer)).
When the ion is really `[M+NH4]+`, the neutral mass it computes is **M + NH₃**,
17.0265 Da too high. That is not a random error in a harmless direction. It is a
shift by an intact nitrogen-containing group, so the target it hands to candidate
generation is the formula of the truth *plus one nitrogen and three hydrogens*.

Worked through on glucose:

| Quantity | Value |
|---|---:|
| glucose, neutral monoisotopic (C₆H₁₂O₆) | 180.063388 |
| `[glucose+NH4]+` precursor | 198.097211 |
| neutral mass the pipeline derives (`m/z − proton`) | 197.089935 |
| offset from the truth | **+17.026547** (NH₃) |

At the 5 ppm window `pipeline/candidates.py` uses for library retrieval, the only
formula within reach of 197.0899 across C≤20 H≤40 N≤4 O≤12 is **C₆H₁₅NO₆** —
glucose plus ammonia. No nitrogen-free formula fits. Widen to 0.005 Da and four of
the five formulas that fit carry nitrogen. So the search does not merely lose
glucose, it is steered onto nitrogenous analogues, and glucose itself is
unreachable at every mass window in force (5 ppm, 0.005 Da, and 1.0 Da all need
17.03 Da to reach it). A reader then sees an amino sugar, or sees the structure
withheld as a mass mismatch, and neither is a statement about the evidence.

**This is an ion-identity problem, not a candidate-structure problem.** No amount
of better candidate generation, reranking, or isomer discrimination fixes it,
because the correct answer is never in the candidate set.

### The population at risk, measured

Ammonium is a positive-mode phenomenon and the delivered ladder confirms it:
`data/feature_rows/*_neutral_mass_hypotheses.csv` carries **41 NH-bearing
hypothesis rows across 36 features in lib400, and 25 rows across 18 features in
lib402**, with **zero** in lib209 and lib305 — 54 controls in all.
Every one of those 66 rows stores the `[M+H]`-derived neutral mass; the median
offset from the correct `[M+NH4]+` inversion is **17.0265 Da**, i.e. exactly the
error above, on every row. The `2m+H+NH3` rows are worse still, because the dimer
error and the ammonium error compound: `X-25217` at m/z 146.06342 stores 145.056
where the label implies 64.015.

The three rows labelled `m+NH4` outright are all in lib402: `X-26332` (m/z
467.373), `X-26324` (539.430), `X-26300` (595.493). The rest use the CHRO_LIB
spelling `m+H+NH3`, which is the same ion.

### The evidence the fix needs already exists, unparsed

The `CHRO_LIB:observed_consensus` row of each feature carries, in its free-text
`note` column, a `consensus_observed_adducts=` list — the ion forms observed for
that compound. This is the co-eluting ion-form cluster, already computed upstream
and already committed:

| | lib400 | lib402 | lib209 |
|---|---:|---:|---:|
| features carrying a consensus ion-form list | 547 / 572 | 344 / 353 | 1689 / 3381 |
| features listing 2 or more ion forms | 456 | 256 | 1418 |
| features listing both a plain (de)protonation **and** an NH-bearing ion | **54** | **24** | 0 |

`X-23416` is the shape of the thing: its list contains `m+H`, `m+H+NH3`,
`m+H+H2O`, `m+H+CH3OH`, `m+H[13C1]` and `2m+H` — a protonation, its ammonium and
solvent clusters, its carbon isotopologue, and its dimer. Read as a cluster, that
is one compound seen six ways, and five of the six ions confirm the sixth. Read one
row at a time, as the pipeline does today, it is six unrelated mass targets, five
of which are wrong.

Nothing in the package parses that column. `pipeline/candidates.py::mass_hypotheses`
(line 280) takes `target_neutral_mass` from the `neutral_mass` column and the
`note` never leaves the CSV.

## What this does not plan

**The adduct arithmetic itself.** Deriving a neutral mass from a label that names
an ion is stage 0 and stage 1 of
[plans/identification_rules_hardening.plan.md](identification_rules_hardening.plan.md),
which also owns `chemistry/adducts.py`, the charge/multiplier fallback, and the
integer-multiple invariant. **This plan depends on that module and must not write
a second copy of the grammar.** The division of labour: that plan makes the
arithmetic follow the label; this plan decides *which label is true*.

That distinction is load-bearing, not cosmetic — see the first unresolved
question below. If the ladder's NH labels name co-observed buddy peaks rather than
the feature's own ion, then deriving a neutral mass from the label makes these 66
rows worse, and only the cluster-plus-spectrum route in stages 3 and 4 gets them
right.

**Mass tolerance.** Four windows disagree by three orders of magnitude. Measured
and owned by [plans/mass_policy_impact_measurement.plan.md](mass_policy_impact_measurement.plan.md)
and [docs/mass_axis_decisions_needed.md](../docs/mass_axis_decisions_needed.md).
This plan changes *which mass a candidate is compared against*, never *how close
it has to be*.

**Isotope and halogen evidence.** ISO-4 through ISO-8 and the `iso_count` key
mismatch are stages 3 and 4 of the hardening plan. Isotopologues appear here only
as members of an ion-form cluster.

**Isomer ranking.** Once the neutral mass is right, choosing among structures at
that mass is [plans/isomer_discrimination.plan.md](isomer_discrimination.plan.md).

## Key decisions already taken

Recorded so the next agent does not relitigate them. Reverse any deliberately,
not by accident.

- **Ion identity is inferred, with a posterior, not read off a label.** The ladder
  records ion forms that were observed near a feature; which one *is* the
  fragmented precursor is a separate question that the label does not answer.
  Every stage below treats the label as evidence, weighted, never as ground truth.
- **`[M±H]` stays the default and the fallback.** It is the right answer for the
  large majority of features (1,839 of 1,839 `protonated` rows already follow it
  exactly). The change is that it becomes a hypothesis that can lose, not the only
  hypothesis there is. A feature with no usable evidence keeps today's answer.
- **Downstream receives a set, not a scalar.** A feature yields one target per
  surviving ion form, each carrying `ion_form`, `neutral_mass`, `ion_form_posterior`
  and the evidence that produced it. Collapsing to a single winner before candidate
  generation throws away the thing that makes the claim auditable.
- **Undecidable is a verdict, and it is not the same as `[M+H]+`.** Where the MS2
  cannot distinguish hypotheses, say so and carry both. Silently defaulting is what
  produces today's confident wrong answer.
- **Fragment conventions get fixed before hypotheses are scored.** Stage 4 asks
  "do the fragments match protonated glucose"; that question cannot be answered on
  a fragment axis that is off by a proton (stage 2). Scoring hypotheses on the
  current axis would rank noise.
- **Do not reimplement clustering upstream.** `untargeted_workflow` already runs
  OpenMS `MetaboliteFeatureDeconvolution` with `[M+NH4]+` in its adduct table and
  already computes `group_id`, `adduct_type` and `neutral_mass`. The upstream
  defect is that it drops them at the export boundary, not that it lacks them.
- **Nothing widens ranking before it is measured.** Every stage that moves a
  neutral mass ships behind a flag and reports its affected-feature count first.

## Constraints

- **The raw MS2 is not in this checkout.** `run_inference.sh:71` reads
  `data/observed_ms2_lookups/`, which is a fetched asset (`fetch_assets.sh:179,205`)
  and absent here. That is where `ms2_pairs_json`, the only source of fragment m/z
  at full precision, lives. Stages 2 through 5 can be built and unit-tested here;
  their effect on delivered results can only be measured where the bundle exists
  (`docs/ec2_export_handoff.md`).
- **The committed MS2 is 10 ppm bins, and that is enough for the feasibility
  counts but not for the scoring.** `observed_ion_bins_json` and
  `observed_loss_bins_json` in the feature rows invert to about 0.002 Da at 200 Da,
  which separates an NH₃ loss (17.0265) from an OH loss (17.0027) by 13 bins. Good
  enough to count how often the test is decidable; not the input the scorer should
  take when the raw peaks are available.
- **No rdkit, pandas, scipy or sklearn in this checkout** (numpy only). The stage 0
  tool must be stdlib plus numpy, following `tools/audit_identification_rules.py`,
  which is in the CI documentation job precisely because it does not need the
  inference stack.
- **Positive-mode MS2 here is depleted near the precursor**, which caps how much
  the diagnostic-loss test can carry:

  | | lib400 | lib402 |
  |---|---:|---:|
  | median gap from precursor down to the nearest fragment below it | 54.5 Da | 159.1 Da |
  | features with any fragment within 20 Da below the precursor | 171 (29.9%) | 59 (16.7%) |
  | features with any fragment within 5 Da below the precursor | 34 (5.9%) | 12 (3.4%) |

  An NH₃ loss lands 17.03 Da below the precursor, so the test is **decidable for at
  most 30% of lib400 and 17% of lib402 features** on this data. Of the 25
  NH-labelled lib402 rows, **none** shows a fragment at precursor − 17.0265. The
  cluster (stage 3) has to carry the rest, which is why stage 3 is not optional.
- **Two of the axes the request names are not available in this checkout.**
  `prevalence` is 1.000 on every row of all four libraries (min = median = max), so
  reproducibility across samples and injections cannot be sourced from the
  delivered ladder. And there is one bin set per feature, with no per-collision-energy
  breakdown, so consistency across collision energies cannot be measured either.
  Both axes exist upstream in `untargeted_workflow` (per-sample intensities in
  `feature_df`, per-scan MS2 in `peak_ms2_df`) and reaching them is stage 7.
- **Moving the mass axis invalidates the Platt calibration** (`PLATT_A`, `PLATT_B`
  in `reports/confidence.py`, refit by `tools/calibrate_confidence.py`) and every
  audited number under `docs/paper/`. Stage 5 moves it. Budget the refit and the
  re-audit into stage 5, once, not per stage.
- **CI runs a regenerate-and-diff.** Anything written under `docs/paper/data/` must
  reproduce byte-identically from committed inputs
  (`.github/workflows/publication-checks.yml`).
- **Compliance.** Observed MS2 ion detail is internal-only. Stage 6 adds an
  identity line naming an ion form; keep fragment m/z out of the public bundle,
  following `tests/test_report_disclosure.py`.

## Stage 0 — make the failure countable

Measurement only, no behaviour change, and the gate for everything after it.
Write `tools/audit_ion_identity.py`, stdlib plus numpy, read-only, reproducible
from committed inputs. It emits, per library:

1. **Population at risk.** Hypothesis rows by ion family (ammonium, sodium,
   potassium, dimer, multiply charged, solvent cluster, isotopologue, numeric
   offset), with the neutral-mass shift each implies against the stored value.
2. **Ion-form clusters.** `consensus_observed_adducts` parsed out of the `note`
   column into a per-feature ion-form set, with the multiplicity distribution and
   the count of features carrying both a plain (de)protonation and a non-trivial
   adduct.
3. **Decidability of the diagnostic-loss test.** Per feature, the gap from the
   precursor to the nearest fragment below it, and whether the window a given
   diagnostic loss would occupy is populated at all. This is the number that says
   how much of the problem MS2 can solve and how much needs the cluster.
4. **The fragment-convention check of stage 2**, as a standing regression: the
   position of common neutral losses in `observed_loss_bins_json`.

Write to `docs/paper/data/ion_identity/` and summarise in
`docs/ion_identity_findings.md`, following the shape of
`docs/paper/data/identification/`. Every number in this plan must come back out of
this tool; where it does not, the tool is right and this document is stale.

## Stage 1 — one ion-form arithmetic, charge carrier explicit

Depends on `structural_elucidation/chemistry/adducts.py` from stage 0 of the
hardening plan. If that has not landed, land it first — do not fork the grammar.

Add to that module:

- `ion_form_candidates(parent_mz, ion_mode, *, allow_multimer=True) -> list[IonForm]`,
  enumerating the hypotheses worth testing for an observed precursor, each with its
  implied neutral mass, charge, multiplier, and the neutral gain or loss that
  distinguishes it. Seed it from the phase-0 ladder
  (`utils/phase0._adduct_neutral_mass_hypotheses`, line 163), whose arithmetic is
  already correct for `[M+NH4]+`, `[M+Na]+`, `[M+K]+`, `[M+2H]2+` and `[2M+H]+`,
  and which is inert today only because nothing on the delivered path calls it.
- `charge_carrier_mass(ion_form)` and a fragment-m/z predictor that takes the ion
  form as an argument, for stage 2.
- `diagnostic_losses(ion_form)`: the neutral loss that would corroborate the
  hypothesis. `[M+NH4]+` → NH₃ at 17.026549 yielding `[M+H]+`; `[M+Na]+` → no
  cheap diagnostic in MS2 (say so, do not invent one); `[2M+H]+` → the monomer
  `[M+H]+`; solvent clusters → the solvent.

Then retire the duplication that currently makes this ambiguous. There are three
proton constants and two neutralization conventions in the package
(`nce/evidence.py:45`, `utils/phase0.py:39`, `vendor/nce/elucidation/chemFacts.py:25`),
and `phase0.closest_neutral_mass` (line 150) takes a `predicted_mass` argument,
discards it, and returns the `[M±H]` value regardless. Leave the constants alone
for now — they agree to 1.4e-8 Da — but mark `closest_neutral_mass` as the
misnamed function it is, because a later stage will need it to mean what it says.

**Expected output:** no delivered number changes. Unit tests on the enumeration
per ion family, including that `[M+NH4]+` inverts 198.0972 to 180.0634 and that
`[2M+H+NH3]+` inverts 146.06342 to 64.0148.

## Stage 2 — the fragment axis is off by a proton

Blocking for stage 4, and worth doing on its own because it silently weakens every
MS2 score in the pipeline today.

**What is wrong.** `nce/evidence.candidate_fragment_bins` (line 387) predicts bins
of *neutral* fragment masses (`Descriptors.ExactMolWt` of each single-bond
fragment). The observed bins it is compared against are bins of fragment *ion m/z*:
`evidence.neutral_ion_mass` (line 136) adjusts by the electron mass, which at 10 ppm
bins is indistinguishable from the raw m/z. Confirmed empirically on the features
where the precursor survives into its own spectrum: 144 of 572 lib400 features have
an observed bin at `mass_bin(m/z)`, and only 3 have one at `mass_bin(m/z − proton)`.
So predicted and observed ion bins sit a proton apart — about 5,600 ppm at 180 Da,
which is hundreds of 10 ppm bins. They cannot match except by accident.

The loss axis carries the same offset, and negative mode proves it decisively.
`ion_and_loss_bins` (line 177) computes the loss as
`parent_neutral(proton) − ion_neutral(electron)`, so a true loss `L` is recorded at
`L + 1.00783` in negative mode. Counting features whose loss bins are populated at
each position:

| Neutral loss | at its true mass | at true mass + 1.0078 |
|---|---:|---:|
| H₂O (18.0106) | 0 | **71** |
| CO₂ (43.9898) | 2 | **173** |
| HCOOH (46.0055) | 3 | **103** |
| SO₃ (79.9568) | 0 | **97** |
| hexose (162.0528) | 2 | **82** |
| glucuronide (176.0321) | 4 | **76** |

Every common loss is in the wrong place, by exactly a proton plus an electron.
Note that `chemistry/diagnostics.py` specifies its diagnostic losses at their
true masses (sulfate 79.9568, glucuronide 176.0321, glycine 75.0320) with a
0.010 Da window, so those rules are looking 1.0078 Da away from where the data
sits.

Worse, the two producers of observed bins disagree with each other.
`verifier_v2._positive_adduct_bins_from_pairs` (line 51) deprotonates each
fragment before binning (`fragment_neutral = mz - PROTON_MASS`), which *is*
compatible with `candidate_fragment_bins`, while the committed bins are not — and
`annotate_predicted_spectrum_scores_v2` (line 200) takes `max()` over scores from
both, i.e. over two different conventions.

**What to do.**

1. Pick one convention and state it in `docs/identification_rules.md`. The
   recommendation: **compare ion m/z to predicted ion m/z**, with the charge
   carrier of the hypothesised ion form applied explicitly. That is the convention
   that lets an ion-form hypothesis change the fragment axis, which is what stage 4
   needs. A neutral-to-neutral convention hides the charge carrier inside a
   constant and cannot express "these fragments are protonated glucose".
2. Define the neutral loss as `precursor m/z − fragment m/z`, ion minus ion, with
   no neutralization on either side. This quantity is invariant to the adduct, so
   it is the one axis that is usable *before* ion identity is known — which is
   exactly why the NH₃-loss test works.
3. Fix the producers to agree, and add a contract test that the two paths give the
   same bins for the same peaks. The class of bug here is a convention
   disagreement, so test the agreement, not just the behaviour.
4. Re-check `reports/compound.score_smiles_spectrum_alignment` (line 1111), which
   hardcodes `neutral_parent_mass = parent_mz + 1.007276466812` — the negative-mode
   form, applied regardless of polarity, and so wrong by 2.01 Da on every positive
   library.

**Verification:** the stage 0 loss-position table flips, with common losses landing
at their true masses; a fabricated spectrum of protonated glucose fragments scores
near 1.0 against glucose and near 0 against the C₆H₁₅NO₆ analogue; the two
observed-bin producers agree bin for bin.

**Expect MS2 scores to move for every feature**, in both directions. This is the
one stage here whose blast radius is the whole matrix rather than a labelled
subset, so land it alone, and re-run `tools/audit_presentation.py` after it.

## Stage 3 — cluster the co-eluting ion forms into one identity

Buildable and testable entirely from committed data, and the load-bearing stage
given how rarely the MS2 test is decidable.

Parse `consensus_observed_adducts` from the `note` column into an ion-form set per
feature, then group features into compound-level clusters keyed by
`(lib_id, retention window, implied neutral mass)` under the stage 1 arithmetic, so
that `[M+H]+`, `[M+NH4]+`, `[M+Na]+`, the dimer and the isotopologues of one
compound land together. Emit a cluster table: cluster id, member features, the ion
form each member is hypothesised to be, the consensus neutral mass, and the
evidence that joined them.

Two details that will bite:

- **Retention is partly missing.** The RI is embedded in `feature_key` (for example
  `feature_X_400_mz105.0914_ri0_X-22996`) and a substantial share are `ri0`, which
  means absent rather than zero. Co-elution cannot be asserted for those; join them
  on mass relationship alone and record the weaker evidence rather than pretending
  to a retention match. `tools/audit_rt_feasibility.py` already deals with this
  column and is the precedent.
- **The cluster is evidence about a compound, not proof about a feature.** A
  cluster containing `m+H` and `m+H+NH3` says the compound forms both ions. It does
  not by itself say which one *this* feature's precursor was. That is stage 4.

**Verification:** the 54 lib400 and 24 lib402 features that list both a plain
protonation and an NH-bearing ion each land in a cluster with a consistent neutral
mass; clusters do not merge two compounds whose neutral masses differ by more than
the mass window; a feature whose list has one entry stays a singleton.

## Stage 4 — score each ion-form hypothesis

For each feature, for each hypothesis from `ion_form_candidates`, combine:

1. **Cluster support.** Does the feature's cluster contain a plain protonation at
   the m/z this hypothesis implies? For `[glucose+NH4]+` at 198.0972 the supporting
   observation is an `m+H` ion at 181.0707 in the same cluster. Co-occurring
   `[M+Na]+`, `[M+K]+` and dimers of the same neutral all count.
2. **Diagnostic neutral loss.** From `diagnostic_losses(ion_form)`: for an ammonium
   adduct, a fragment at `precursor − 17.026549`, which *is* `[M+H]+`. Record three
   outcomes, not two — present, absent, or **out of range**, the last when no
   fragment is observed anywhere in that window (70% of lib400 and 83% of lib402
   features). Absent is evidence against; out of range is no evidence at all, and
   conflating them is how a correct hypothesis gets rejected for want of data.
3. **Fragment agreement with the implied neutral.** Under the stage 2 convention,
   score the fragments below the diagnostic loss against the neutral this
   hypothesis implies. This is the "subsequent fragments match protonated glucose"
   test, and it is the one that distinguishes a genuine `[M+NH4]+` from a real
   nitrogenous compound that happens to lose NH₃.
4. **Plausibility of the implied neutral.** Formula prior at the implied mass
   (`chemistry/formula.formula_prior_score`), so a hypothesis implying an
   implausible neutral is penalised.

Emit per feature a ranked hypothesis table with a posterior and an explicit
`ion_form_status` of `supported`, `undecidable`, or `refuted`. Do not collapse to a
winner here.

**The asymmetry to preserve.** The failure this plan exists to fix is a false
nitrogen. The mirror failure — calling a genuine aminated metabolite an ammonium
adduct of its deaminated parent — is equally wrong and equally invisible. The
discriminator is real: a true `[M+NH4]+` loses NH₃ readily and then fragments as
the protonated parent, whereas an aminated compound retains the nitrogen in its
fragments. That discriminator only works on the 30% of lib400 and 17% of lib402
features with fragments in the window, so on the rest the honest output is
`undecidable` plus the cluster prior, not a coin flip.

**Verification:** on the truth controls in
`data/known_X_all_available_structures_simple.csv` (92 with structures), no control
whose truth structure matches the `[M+H]+` reading may be moved to an ammonium
reading; a synthetic `[M+NH4]+` spectrum built from glucose fragments must select
the ammonium hypothesis; a synthetic spectrum of a genuine C₆H₁₅NO₆ compound,
whose fragments keep the nitrogen, must not.

## Stage 5 — feed the hypotheses into candidate generation

`pipeline/candidates.py::mass_hypotheses` (line 280) currently emits one target per
row from the `neutral_mass` column plus one from `neutral_mass_fallback`. Change it
to emit one target per surviving ion form, carrying `ion_form`,
`ion_form_posterior`, `neutral_mass_source` and `neutral_mass_shift_da`.

Ship behind `--derive-neutral-mass-from-ion-form`, default **off** for the first
run, so both policies can be compared on one input.

**The gate, before it is turned on.** Report how many features gain a target that
moves by more than the mass tolerance, how many controls that touches, and how many
candidate rows change mass-consistency class. From the ladder the expected scope is
the 66 NH-bearing rows plus the sodium, dimer and multiply-charged families the
hardening plan's audit counted at roughly 4,100 hypotheses — large enough that the
effect must be measured, not assumed benign.

**The risk to measure rather than argue.** A correct target set finds structures a
`[M±H]`-only set cannot, and also admits candidates at masses that were previously
unreachable, on a mass axis worth 0.30 of the confidence score. Measure paired, on
the truth controls, the way `tools/audit_mass_policy.py --rerank` already does it:
same candidates, both policies, exact McNemar on rank-1 and rank-3.

Downstream gates that will need to see the new column rather than recompute the old
mass: `utils/merge.rank_selfies_pool` (line 83, drops on `abs_mass_error_da`),
`phase0.passes_mass_filter` (line 250) and `soft_mass_score` (line 262), and
`inference/predict.compute_mass_gated_retrieval`, whose `MASS_GATES_DA` of
`[0.05, 0.5, 1.0]` compares an observed `parent_mz` against a bank of theoretical
`[M+H]+` parents and so is 17 Da outside every gate for an ammonium precursor.

## Stage 6 — say which ion it was

Two changes, and the second matters as much as the first.

**Name the ion form.** The report identity line should read *"glucose, observed as
an ammonium adduct `[M+NH4]+`"*, with the supporting evidence stated: the
co-eluting protonated form, the ammonia loss, the fragments of the protonated
parent. `reports/confidence.py` already has somewhere to put this — the evidence
table at line 1558 surfaces `adduct_label` and `adduct_score` when present — but an
observed ion form is a claim about identity, not a score in a list, and belongs on
the identity line.

**Stop scoring an adduct-explained mass as a mismatch.** `assign_msi_level`
(line 647) demotes to Level 4 when the mass axis falls below `MASS_WEAK = 0.50`
(line 744). Under the `[M+H]+` reading an ammonium adduct carries a 17 Da error, so
any correct structure that reaches the report by another route — de novo
generation, a curated overlay — scores zero on the mass axis and is demoted for the
one reason that is not its fault. (Structures that come through library retrieval
never get that far; they are filtered out upstream, which is what stage 5 fixes.)
Once the ion form is inferred, the mass axis must be scored against the neutral
mass that ion form implies. Withholding a structure for a mass mismatch that the
ion form explains is the specific behaviour this plan removes.

**And a cap, so the fix cannot be abused.** An ion-form claim resting on neither
cluster support nor loss evidence — `ion_form_status = undecidable` — must not
raise the MSI level above what the `[M±H]` reading would have earned. Put it where
the refuting-curation cap sits. This is the guard that stops "observed as an
ammonium adduct" from becoming a way to explain away any inconvenient 17 Da.

## Stage 7 — stop discarding ion identity upstream

Separate repository (`untargeted_workflow`), separate PR, and independent of
stages 1 through 6, so it can run in parallel.

The upstream pipeline already does most of this work. `OpenMSGroupRelatedPeaksStep`
runs `MetaboliteFeatureDeconvolution` with `[M+NH4]+` in its adduct table
(`src/goo/pipeline/openms_peak_grouper.py:46-51`) and writes `group_id`,
`adduct_type`, `is_base_peak` and `neutral_mass` per feature (lines 285-295);
`IsotopeLinkingStep` then links C13 isotopologues. MS1 library matching already
honours it, computing `target_mz` from `neutral_mass` rather than from the raw m/z
(`src/goo/pipeline/spectral_matcher.py:217-224`).

Then it throws it away. The MGF export writes `TITLE`, `FEATURE_ID`, `SCANS`,
`N_SPECTRA` and `PEPMASS` and nothing else (`spectral_matcher.py:601-616`), so
`consensus_ms2.mgf` — the natural hand-off to elucidation — carries a bare
precursor m/z and no ion identity at all. Downstream then has no choice but to
assume `[M+H]+`, which is the bug this whole plan is working around.

1. **Carry the annotation into the export.** Add `ADDUCT`, `NEUTRAL_MASS`,
   `GROUP_ID` and `CHARGE` to each MGF block from the grouper's metadata. Small,
   self-contained, and it is the single highest-value change in this plan per line
   of code.
2. **Implement the correlation gating.** Three TODOs already mark the right places:
   `openms_peak_grouper.py:358` (split group members whose across-sample intensity
   vectors do not correlate with the base peak) and `isotope_linker.py:85,168`
   (require correlation before accepting an isotope link). Until these exist the
   cluster stage 3 depends on can over-group, and an over-grouped cluster is
   affirmative evidence for a wrong ion form.
3. **Supply the two axes the elucidation repo cannot see.** Reproducibility across
   samples and injections is computable from `feature_df` rows; consistency across
   collision energies is computable from `peak_ms2_df` before consensus. Both are
   in the request and neither is available downstream, where `prevalence` is
   constant at 1.000. Export them per feature.

## What must not change while measuring

- The mass tolerances. Stage 5 changes the target, not the window; changing both at
  once makes neither attributable.
- The `[M±H]` default for features with no usable evidence.
- `linear_mass_match_score` and `predict.mass_match_score`, owned by the
  mass-policy plan.
- The Platt constants, until stage 5 is decided. Refit once, at the end.
- The isotope and halogen work of the hardening plan. Isotopologues are cluster
  members here and nothing more.

## Verification, end to end

```bash
python3 tools/audit_ion_identity.py                      # must reproduce byte-identically
python3 tools/audit_identification_rules.py              # must be unaffected by stages 0-4
python3 -m pytest tests/test_ion_identity.py tests/test_identification_rules.py -q
python3 -m ruff check structural_elucidation/chemistry/adducts.py tools/ tests/
git diff --exit-code docs/paper/data/
```

Then, where the fetched bundle exists: `tools/audit_mass_policy.py --rerank` for the
paired policy comparison, `tools/prototype_mass_gate_rerank.py` re-run because
stage 5 changes its inputs, and `tools/audit_presentation.py` for the MSI and tier
distributions, which stage 2 moves for every feature and stage 6 moves by design.

## Unresolved questions

- **Does an NH label name the feature's own ion, or a buddy peak observed beside
  it?** The highest-value single measurement in this plan, because the two readings
  demand opposite fixes. Evidence for the buddy reading: none of the 25 NH-labelled
  lib402 rows (18 distinct features) shows a fragment at `precursor − 17.0265`, and
  each such feature's consensus list also contains `m+H`, which cannot be the same
  m/z as `m+H+NH3`.
  Evidence against concluding it: only about 17% of lib402 features have any
  fragment in that window, so the absence is weakly informative, and the stored
  neutral mass is `[M+H]`-derived either way. Resolve it on the raw MS2 where the
  bundle exists, on the 78 features that carry both forms, before stage 5 changes a
  single delivered mass. This is the same question the hardening plan raises for
  `m-H+62.9959`, and the two should be answered together.
- **What decides ion identity for the 70-83% of positive features with no fragment
  near the precursor?** Stage 3's cluster is the only in-repo answer, and it is
  weaker than a spectrum. Upstream co-elution and correlation (stage 7) is the real
  answer, which is an argument for doing stage 7 early rather than last.
- **Should an undecidable ion form widen the candidate set or hold the structure
  back?** Widening admits nitrogenous candidates at the `[M+H]`-derived mass, which
  is precisely the failure being fixed, in the other direction. Holding back costs
  coverage. Decide on the measured rate of each error on the 92 truth controls, and
  note that 92 controls resolve very little — the same argument for the
  library-holdout validation in `docs/mass_axis_decisions_needed.md`.
- **Is `prevalence` genuinely 1.000, or unpopulated?** It is constant across all
  10,568 rows of all four libraries, which reads like a default rather than a
  measurement. Answer before any stage relies on reproducibility.
- **Is the CHRO_LIB exporter fixable at source?** If it is, parsing ion forms out of
  a free-text `note` column is a stopgap and should be marked as one. If it is not,
  the parser is permanent and belongs behind a stable interface rather than inline.
- **Which fragment convention?** Stage 2 recommends ion-to-ion with an explicit
  charge carrier. The alternative, neutral-to-neutral, is a smaller change to the
  existing code and cannot express an adduct hypothesis on the fragment axis.
  Getting this wrong makes stage 4 unbuildable.

## How to dispatch this

**Stage 0 first and alone.** It is measurement only, it is the gate for every
later stage, and it converts this document's numbers into a regenerable artifact.

Then two tracks that do not touch the same files and answer different questions:

- **Track A, the fragment axis:** stage 2, then stage 4. Stage 2 is independently
  valuable — it fixes a defect that weakens every MS2 score today, labelled subset
  or not — and stage 4 is unbuildable without it.
- **Track B, the cluster:** stage 3, which needs only stage 1 and committed data.

Stage 1 is a prerequisite for both and should follow the hardening plan's stage 0
rather than duplicating it. Stages 5 and 6 need A and B to have landed and
measured. Do not merge either track into the other before both are measured: stage 2
moves every MS2 score and stage 3 changes which features share an identity, and a
combined result is attributable to neither.

**Stage 7 is a different repository and can start immediately.** Its first item —
four extra fields in the MGF writer — is the smallest change in this plan and
unblocks every downstream consumer of ion identity, including consumers outside
this repository.

If only one thing gets done, do stage 0 plus the MGF fields in stage 7: the first
makes the failure countable, the second stops the information being destroyed while
the rest is being built.

## Task list

1. `tools/audit_ion_identity.py`: population at risk by ion family, per-feature
   ion-form clusters from `consensus_observed_adducts`, decidability of the
   diagnostic-loss test, and the loss-position convention check. Emit to
   `docs/paper/data/ion_identity/` and `docs/ion_identity_findings.md`.
2. `chemistry/adducts.py`: `ion_form_candidates`, `charge_carrier_mass`,
   `diagnostic_losses`, seeded from the phase-0 ladder; unit tests per ion family.
3. Fix the one-proton offset between predicted and observed fragment bins, make the
   charge carrier an explicit argument, reconcile the two observed-bin producers,
   and fix the hardcoded negative-mode parent in
   `reports/compound.score_smiles_spectrum_alignment`.
4. Build the per-feature ion-form cluster table, keyed by implied neutral mass and
   retention window, handling absent RI honestly.
5. Score ion-form hypotheses per feature from cluster support, diagnostic loss,
   fragment agreement and formula plausibility; emit a posterior and an explicit
   `undecidable` status.
6. `pipeline/candidates.py::mass_hypotheses`: emit one target per surviving ion
   form behind `--derive-neutral-mass-from-ion-form`, with the affected-feature
   count and the paired McNemar comparison reported before it is enabled.
7. Report the ion form on the identity line; score the mass axis against the
   implied neutral mass; cap undecidable ion-form claims.
8. `untargeted_workflow`: add `ADDUCT`, `NEUTRAL_MASS`, `GROUP_ID` and `CHARGE` to
   the MGF writer; implement correlation gating at the three existing TODOs; export
   the per-sample and per-collision-energy reproducibility axes.
9. Add `python3 tools/audit_ion_identity.py` to the "Confirm derived artifacts are
   current" step of `.github/workflows/publication-checks.yml` (line 44, beside the
   existing `audit_identification_rules.py`), and to the ruff list in the "Lint the
   publication tooling" step (line 109).
