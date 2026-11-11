---
name: x-compound-report
description: "Build a client X-compound deliverable from a list of X compound ids: by default a client folder with a summary CSV and one PowerPoint (one slide per compound, the UPMC slide process), generated from the latest naming-v2 10 ppm results. Also builds the older per-compound public HTML reports, master index and executive_summary_top10.xlsx. Use when someone supplies X-compound ids (X-12345) and asks for a client deliverable, summary CSV, slides, X compound reports, or a refresh of past client reports."
---

# X-compound client deliverables

Two deliveries live in this skill. **Default to the slide deliverable (A)** unless the user asks for HTML reports.

**A. Summary CSV + PowerPoint (UPMC slide process, current).** Given X ids, produce a client folder:

```
<Client>_X-compound_summary.csv          one row per compound
<Client>_X-compound_slides_<yyyymmdd>.pptx   one deck, one slide per compound
README.md
_build/                                  provenance, slide-input xlsx, screen; internal, do not send
```

**B. Public HTML reports (earlier clients).** Per-compound `public.html`, master `index.html` and
`reports_index.csv`, and `executive_summary_top10.xlsx`. See "Delivery B" below.

Everything runs in the elucidation virtual environment. The helper scripts are in
`.cursor/skills/x-compound-report/scripts/`.

## Delivery A: summary CSV and slides

The process is the newest UPMC slide run on the machine where it was recorded:
`x_compound_info_summary_20261014/build_x_compound_slides.py` (2026-10-14 12:04), which wrote
`UPMC_X-compound_slides_20261014.pptx`. It opens the template deck `20260919 external slide deck.pptx`
(Montserrat, Metabolon colours, footer), drops its slides, and adds one slide per compound on the
template's **"2 Columns"** layout: a grey title band (`X-#####: Preliminary Structural Information`),
the candidate super/sub-pathway line, three cards (rank-1 predicted structure + SMILES, feature summary,
most similar named metabolite + Tanimoto) and a takeaway band. `scripts/build_slides.py` is that script,
vendored. Its only changes: paths are arguments, and blank fields print "Not assigned".

The UPMC summary sheet was typed by hand. Here `scripts/client_deck.py` generates it from the results, with
the UPMC columns (Compound ID, Observed m/z, Putative formula, Candidate superpathway/subpathway,
Rank-1 predicted SMILES, Rank-1 additional Information, Most similar named metabolite, Tanimoto similarity)
plus MSI level, MSI label, pathway basis and Metabolon level.

### Which results to use

Use only the naming-v2 10 ppm delivery:

- `output/reports_inventory_naming_v2/` (the `X-#####` folders link to `output/reports_ppm10_full_naming_v2/`)
  and its `exec_summary/exec_summary_by_compound.csv`. Each compound's pathway, neighbours, m/z and formula
  come from `X-#####/data/`.
- There is no `scores/<id>.json` or `run_novel` output on the reference machine, so the deck does not use them.
  If a later machine has them, stop and ask which is authoritative; do not mix sources silently.
- **Never use the June 2026 bundle (`handoff/se_report_inputs.tgz`) or `reports 4/`.** `client_deck.py`
  refuses those paths.
- If a requested compound is not in the naming-v2 inventory, it is listed as `missing` and left out. Say so in
  the report-back. Do not back-fill it from older data. A compound with no structure, neighbour or pathway
  in naming-v2 stays on the slide with "Not assigned".
- `output/` is not in git. Get it from the EC2 workspace (`fetch_assets.sh`, `docs/ec2_export_handoff.md`)
  or point `--reports-root` at the copy on the machine you are on.
- A newer tree (for example `output/reports_v12_ibp`) changes calls for many compounds. Do not switch to it
  unless the user says it is approved for delivery.

### Run it

```bash
.venv/bin/python .cursor/skills/x-compound-report/scripts/client_deck.py \
  --compounds ids.txt --client "NIH" --out NIH_deliverable \
  [--reports-root output/reports_inventory_naming_v2] \
  [--template "/path/to/20260919 external slide deck.pptx"] [--dry-run]
```

- `--compounds` takes literal ids and `.xlsx`/`.csv`/`.tsv`/`.txt` files in any mix; ids are normalized
  (`x 1911` becomes `X-01911`), de-duplicated and screened.
- **Screen:** compounds that `data/known_X.csv` marks `Exclude`, `Artifact` or `DO NOT INCLUDE` are other
  clients' proprietary or spiked compounds and are dropped (`--allow-excluded` is for internal QA only).
- The template is not in git (22 MB). It defaults to `../x_compound_info_summary_20261014/` next to the repo;
  otherwise pass `--template`.
- Python needs python-pptx, openpyxl, lxml and Pillow. RDKit is optional and gives crisp drawings; without it
  the slides fall back to the dashboard PNGs. Run `--dry-run` first to see built, excluded and missing ids.
- It takes seconds, so run it directly.

### Verify, then report back

1. `_build/deck_manifest.json`: `built`, `excluded`, `missing`, and each source's path and mtime.
2. Slide count equals the number of built compounds; the CSV has the same rows; only the requested ids appear.
3. Open the deck once (Quick Look or PowerPoint) and check one slide with a structure and one "Not assigned".
4. Report the folder path, file list, sources with mtimes and why they are the latest, missing and excluded ids,
   and that slides are one deck.

## Delivery B: public HTML reports

Given a list of X compounds, produce the same delivery that earlier clients received:

- one self-contained report per compound (`public.html`, served as `<id>/index.html`, with the
  similarity CSVs it links),
- a sortable master `index.html` and `reports_index.csv`,
- `executive_summary_top10.xlsx` (top-10 candidates, by-compound rollup, naming conventions, summary).

These helpers only orchestrate code that already ships in the repo (`structural_elucidation.reports.batch`,
`tools/build_exec_summary.py`, `tools/build_exec_summary_workbook.py`, `tools/audit_reports.py`).

### Rules that are not negotiable

Rules 1 and 3 apply to both deliveries; 2, 4 and 5 are for Delivery B.

1. **Never deliver a compound that `data/known_X.csv` marks `Exclude`, `Artifact` or
   `DO NOT INCLUDE`.** These are other clients' proprietary compounds or spiked standards (Merck
   spiked compounds, a Novo Nordisk drug, a MARS compound, "client drug B-199"). The screen
   drops them; `--allow-excluded` exists for internal QA only.
2. **Clients get the public package only.** `internal.html`, scoring JSONs, build logs and
   observed MS2 ion detail stay in `work/`. The package step copies `public.html` and the files
   it links, then fails the run if anything internal slipped in.
3. **Use the latest results, and say which ones you used.** Run the resolver first and read its
   warnings. If it reports the June 2026 bundle or no `unified_rerank_all_top10_ppm10.csv`, ask the
   user for the current delivery root before sending anything to a client.
4. Do not rerun inference (`run_inference.sh`, `run_full_pipeline_ec2.sh`) for a report request.
   Reports are rendered from existing rerank CSVs, scoring JSONs or rendered bundles.
5. Regeneration takes about 2 seconds per compound at `--jobs 4` with the synonym DB. Past
   `--max-compounds` (default 60, roughly 2 minutes) the script refuses and writes the exact command into
   `delivery_manifest.json` (`regenerate_command`). Hand that command to the user instead of running it.

### Steps

#### 0. Environment (once per machine)

```bash
python3.10 -m venv .venv            # or: uv venv -p 3.10 .venv
.venv/bin/pip install -r requirements.txt   # reports need only pandas<2, numpy<2, rdkit, openpyxl, scipy
```

Report inputs are not in git (`output/`, `vendor/metabolite_synonyms.sqlite` are ignored). Get them from
the EC2 workspace with `./fetch_assets.sh`, `fetch_get_assets.sh`, or the handoff in
`docs/ec2_export_handoff.md`. The June 26 2026 bundle is on branch `handoff/se-report-inputs`
(`handoff/se_report_inputs.tgz`, unpack at the repo root). It works for testing, but it is stale for
deliveries (see step 1).

#### 1. Find the latest results

```bash
.venv/bin/python .cursor/skills/x-compound-report/scripts/resolve_latest.py \
  --search output models /path/to/other/deliveries --compounds ids.txt --json-out latest.json
```

The resolver scans for three kinds of source and picks one of each:

| Kind | What it looks for | Preference |
|---|---|---|
| predictions | `<root>/rerank_lib{209,305,400,402}/unified_rerank_all_top10*.csv` | `_ppm10` refill first, then not the June bundle, then newest |
| scores | directories of `X-*.json` with `schema_version: 1` (from `run_novel` / `reports.novel_scoring`) | calibration checksum equal to `tools/identity_confidence_calibration.json`, real ids, newest |
| reports | `<root>/X-*/manifest.json` | naming policy v2, newest |

It reports coverage of the requested ids and flags these problems: stale calibration, synthetic-only ids
(`X-9xxxxxxxx`, untargeted features), producing commits that are not in this checkout's history, dirty
checkouts, and requested compounds whose rank-1 structure differs between prediction roots.

Known lineage, newest first, recorded in git history:

1. `codex/structural-elucidation-latest` (merged to `main` as PR #46, Oct 21 2026). Canonical
   `scores/<control_id>.json` are written before rendering, and `reports.batch --scoring-dir` consumes them.
2. October 15 2026 delivery: `output/reports_ppm10_full_naming_v2` from `scripts/generate_ppm10_naming_v2.sh`
   on the 10 ppm-refilled CSVs (`unified_rerank_all_top10_ppm10.csv`), merged with
   `tools/merge_inventory_delivery.py` (stored mass, then adduct fallback, then earlier adduct-derive). The
   reviewed workbook is `docs/paper/data/exec_summary/executive_summary_top10.xlsx` (3,197 compounds).
3. June 26 2026 bundle: 3,897 compounds, no ppm gate, naming v1 era. Superseded.

#### 2. Build the delivery

```bash
.venv/bin/python .cursor/skills/x-compound-report/scripts/client_report.py \
  --compounds client_list.xlsx \
  --client "UPMC" --out deliveries/upmc_2026-11 --jobs 4 --zip \
  [--reports-root output/reports_ppm10_full_naming_v2] \
  [--predictions-root <delivery_root> --predictions-name unified_rerank_all_top10_ppm10.csv] \
  [--scores-dir output/novel_run/scores] [--observed-lookup-dir data/observed_ms2_lookups]
```

- `--compounds` accepts `.xlsx`, `.csv`, `.tsv` or `.txt` files and literal ids, in any mix. Ids are
  normalized (`x 1911` and `X1911` both become `X-01911`), de-duplicated, and screened.
- `--mode auto` (the default) reuses a rendered v2 bundle when one covers every requested id; this is an
  exact copy of what was delivered. Otherwise it regenerates from predictions. Force a choice with
  `--mode reuse` or `--mode regenerate`.
- Run `--dry-run` first on large lists. It screens the ids, resolves sources, estimates runtime and
  writes the manifest without building anything.

Output under `--out`:

```
compounds_normalized.txt      ids that will be built
compound_screen.csv           every input id: status, known_X level/sublevel/libs, exclusion reason, report found
latest_results.json           resolver output (sources, coverage, warnings)
delivery_manifest.json        provenance: repo HEAD, calibration sha256, sources, counts, audit, warnings
work/reports/                 full bundles (public + internal), master index    <- internal only
work/exec_summary/            exec_summary_*.csv, exec_summary.json, executive_summary_top10.xlsx
work/audit/                   audit_findings.csv, audit_summary.json, audit_gallery.html
<Client>_x_compound_reports/  the client package (and .zip with --zip)
```

#### 3. Verify before handing over

1. `delivery_manifest.json`: `status` is `ok`, `package.problems` is empty, the `warnings` are understood,
   `excluded` and `no_report` are explained to the user, and `audit` shows zero template leaks, wrong
   titles and broken links.
2. Compare with the most recent past delivery:

```bash
.venv/bin/python .cursor/skills/x-compound-report/scripts/compare_reference.py \
  --new deliveries/upmc_2026-11/work/exec_summary/exec_summary_by_compound.csv \
  --reference docs/paper/data/exec_summary/executive_summary_top10.xlsx \
  --out deliveries/upmc_2026-11/compare
```

   Differences in `n_candidates`, `identity_candidate_rank` or rank-1 structure usually mean the inputs
   are not the October 10 ppm delivery. A changed `msi_level` or `feature_name` with the same inputs means
   the report code changed, and needs a reason.
3. Open two or three package pages (`<id>/index.html`) in a browser: a Level 1, a Level 3 and a Level 4.

#### 4. Report back

Give the package path, the counts by MSI level, which sources were used (from `delivery_manifest.json`),
the excluded and missing ids, and every resolver warning.

## What the report may claim (enforced by the renderer, not by this skill)

| MSI level | Structure | Name shown |
|---|---|---|
| 1 | shown | bare name (authentic standard) |
| 2 | shown | `Probable <name>` |
| 3 | shown as hypothesis, up to 3 public candidates | `Possible <name>` |
| 4 | withheld | `<formula>; structure unresolved` (the workbook may say `Isomer of <library metabolite>`) |
| 5 | withheld | `Exact mass <m>; unidentified` |

Human and laboratory verdicts enter only through `data/curation_ledger.csv`.

## Caveats

- Scoring JSONs from `run_novel` score untargeted features and deliberately ignore `known_X` levels. If
  they are applied to registered X compounds with a curated Level 1/2, those compounds lose the curated
  level. The script warns when this happens. Confirm the JSONs were produced for registered X compounds
  before using them.
- Co-abundance (WGCNA) sections use reference cohorts (HELIOS, Twins UK, PRIMM feces/serum) or a
  structural-neighbor proxy. Real client-specific co-abundance would need a per-client
  `parent_neighbors` table (`docs/pipeline_refinement_handoff.md` section 1).
- The client "Curated" and "Match conf" index columns (commit `9e6e359`) are on the unmerged branch
  `cursor/candidate-confidence-rework-f8a2` and do not appear in the master index built from `main`.

## Smoke test (Delivery B)

`examples/sample_client_list.csv` holds 12 compounds from the October delivery (truth controls, Level 2C,
other), one duplicate, two badly formatted ids, one excluded proprietary compound and one unregistered id.
On the June bundle it builds 12 reports in about 25 seconds, excludes `X-26724`, and lists `X-99999` under
`no_report`:

```bash
.venv/bin/python .cursor/skills/x-compound-report/scripts/client_report.py \
  --compounds .cursor/skills/x-compound-report/examples/sample_client_list.csv \
  --client "Sample Client" --out /tmp/x_report_smoke --jobs 4
.venv/bin/python -m pytest -q tests/test_x_compound_report_skill.py
```
