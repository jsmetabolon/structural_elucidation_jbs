# X-compound report audit - findings

Audit of all **3897** delivered report bundles (`tools/audit_reports.py` over
`output/reports`; raw output in `output/_audit/audit_findings.csv` +
`audit_summary.json` + `audit_gallery.html`). All issues trace to the report
renderer [structural_elucidation/reports/compound.py](../structural_elucidation/reports/compound.py),
which was authored for the single compound X-25003 and never fully parameterized.

## Findings (by prevalence)

| Issue | Reports affected | Cause (compound.py) |
|-------|------------------|---------------------|
| Stale `X-25003` text in page | 3897 / 3897 | hardcoded literals throughout `render_html` |
| Broken intra-report links | 3897 / 3897 | `x25003_*` and `_to_target` link literals vs slug-named files |
| Fabricated confidence ("99.9% / Level 2C", "10/10") | 3897 / 3897 | hardcoded metric cards (lines ~2542-2545) |
| Hardcoded identity "N-lactoylalanine" | 3897 / 3897 | subtitle/badges/narrative literals |
| Wrong `<h1>` / `<title>` (shows `X-25003`) | 3896 / 3897 | lines 2391, 2500 |
| Display name is just the control id | 3892 / 3897 | `target_name` falls back to `control_id`; no DB name resolution |
| Empty synonyms section | 3897 / 3897 | `resolve_chemical_id` returns `KNOWN` name, not a DB `chemical_id` |
| Consensus image overlapping (overlap_ratio > 1) | 931 / 3897 | fixed atom radius `r=24` in `consensus_atom_support` (line 687); max 85 atoms |

## Specific code sites
- Identity literals: `<title>` 2391, `<h1>` 2500, subtitle/identity 2502-2503, pathway badge 2508, narrative 2518, metric-note 2550, internal narrative 2325, extended-network tooltip 1703.
- Fabricated metric cards: 2542 (99.9% / Level 2C), 2543 (Lactoyl Amino Acid pathway), 2545 (10/10 high), 2544 subpathway fallback.
- Broken links: `data/x25003_consensus_atom_support.csv` (2538, 2327), `data/x25003_evidence_summary.csv` (2327), `data/all_known_similarity_to_target.csv` + `_unknown_` (2315, 2360). Actual files are slug-named (`{slug}_...`, `all_known_similarity_to_{slug}.csv`) at 3090/3151/3153/3219.
- Consensus crowding: fixed `r=24`, `stroke-width=4`, fixed label offsets and always-on numeric labels in `consensus_atom_support` (660-709); `scale` shrinks for large molecules so circles overlap.

## Planned fixes
1. Parameterize `render_html` from `target['control_id']` + `control_slug(...)`: per-compound title/h1/subtitle; remove fabricated confidence/level cards; drive pathway/subpathway from `db_detail` (honest fallback); parameterize narrative + metric-note + extended-network tooltip.
2. Fix all sidecar links to slug-named files.
3. Adaptive `consensus_atom_support` sizing: atom radius, font, bond width and label offsets scale with the per-molecule scale; drop per-atom numeric labels when too small; keep color + hover tooltip.
4. (P2) resolve a real display name / DB `chemical_id` to populate name + synonyms where possible.

## Verification
Re-run `tools/audit_reports.py` after fixes: target zero stale `X-25003`, zero
broken links, zero fabricated-confidence cards, `h1`/`title` == control id, and
consensus overlap_ratio under threshold on the worst-case sample.
