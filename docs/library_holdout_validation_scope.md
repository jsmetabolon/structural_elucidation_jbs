# Scope: validating against the library, held out against itself

A costing, not a design. It exists so the second route under handoff item 4 can be
funded or declined on facts rather than on the idea sounding reasonable.

## The proposition

Every question of the form "does this reranking change improve results" currently
dies on the same number: **9 of 598** rank-1 mass-mismatch reports have a known
correct answer, and rank-1 recovery on 92 controls carries a twenty-point
confidence interval. Nothing smaller than a ten-point shift is measurable.

The reference library is already thousands of compounds holding both a known
structure and an acquired spectrum. Running those spectra through the pipeline with
each compound suppressed from retrieval makes every one of them a held-out control,
at the cost of compute rather than acquisition.

`tools/audit_pathway_predictions.py --loo` establishes the pattern in this
repository, over 5,095 reference metabolites, and reports the result as mechanism
validation rather than as accuracy.

## What has to be established first

One unknown decides whether this is worth costing further, and it is not answerable
from this checkout.

**How many library standards have an acquired MS2 spectrum in a form the pipeline
can consume?** The pathway leave-one-out needed only structures and labels, both of
which sit in `vendor/metabolite_synonyms.sqlite`. A structural version needs
spectra as well. The relevant store is `data/observed_ms2_lookups/lib{lib}_{ion}_observed_ms2_lookup.csv`,
described in [docs/ec2_export_handoff.md](ec2_export_handoff.md).

If the answer is in the thousands, this is the cheapest validation available. If it
is in the low hundreds, it is still several times the current truth set but the
argument for it weakens. **Answer this before scoping anything else.**

## What it needs to run

Not runnable here. `models/` in a clean checkout is 32 KB of decoder configuration;
the weights and banks are fetched.

- Assets from `fetch_assets.sh`: `models/encoder/best_encoder.pt`,
  `models/decoder_weights/*.pt`, `models/embedding_banks/embeddings_clean/`,
  `vendor/modelfiles/`, `vendor/metabolite_synonyms.sqlite`.
- The library standard spectra, plus their structures for scoring.
- GPU time on the EC2 host, since candidate generation is an encoder-decoder over
  spectra. Unlike the pathway version this cannot be replayed offline from the
  database.

## The suppression mechanism, and where it bites

Suppression has to happen everywhere a compound could leak into its own answer, and
there are three places rather than one.

1. **Library-anchored generation**, [structural_elucidation/pipeline/candidates.py](../structural_elucidation/pipeline/candidates.py).
   The largest candidate source in the delivered set at 13,127 candidates, so this
   is the one that matters. The held-out compound must not be available as an
   anchor for its own query.
2. **Retrieval against the embedding bank**, which would otherwise return the
   compound's own embedding as its nearest neighbour.
3. **Pathway assignment**, which already knows how to do this and can be reused
   directly.

The failure mode to guard against is a partial hold-out that leaks the answer
through one of the three and reports an accuracy that is really a lookup. The
pathway audit's `j != i` exclusion is the shape of the check; the structural
version needs the same guarantee at three sites, and a test that a held-out
compound never appears among its own candidates.

## What comes out, and what it is worth

Each held-out compound yields a candidate list with a known correct answer, which
is the same shape the existing harness already consumes. `tools/audit_mass_policy.py --rerank`
would take it unchanged, and the paired McNemar comparison it already runs would
then have thousands of pairs rather than 92.

Concretely, that is the difference between resolving a ten-point shift and
resolving a one-point shift, which is the range every proposed rerank change so far
has fallen into.

## What it is not

**Not a substitute for truth on unknowns.** Catalogued compounds are a different
population, and this work already measures the gap: exact-source rank-1 recovery is
38.1% against 56.0% for predicted-source. Any absolute rate from this must be
labelled as measured on catalogued chemistry, exactly as the pathway leave-one-out
is labelled.

**Not a replacement for handoff item 4.** Authentic-standard expansion is still the
only thing that grows truth on real unknowns. This changes the sequencing rather
than the destination: it is cheaper, it arrives sooner, and it tells you whether a
proposed change is worth the standards budget before that budget is spent.

## Suggested sequencing

1. Count the library standards with usable spectra. Stop here if the number is
   small.
2. Build the hold-out with the three suppression sites and a leakage test.
3. Run it on a few hundred compounds first and confirm recovery is plausible
   against the exact-source figure of 38.1%, which is the closest existing
   comparator.
4. Scale up, then re-run the mass-policy comparison and handoff item 2's gate
   against it, which is the thing that could not be validated before.
