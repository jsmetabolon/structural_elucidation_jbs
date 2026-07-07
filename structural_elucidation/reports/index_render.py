"""Master-index row model + HTML renderer (pure standard library).

Extracted from ``batch.py`` so the renderer can be reused by stdlib-only tooling
(e.g. ``tools/rebuild_master_index.py``) without importing pandas / the package
``__init__`` chain. ``batch.py`` imports these for full runs.
"""
from __future__ import annotations

import html
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any


@dataclass
class IndexRow:
    """One row in the master index, sourced from a compound ``manifest.json``."""

    control_id: str
    lib: str = ""
    name: str = ""
    formula: str = ""
    exact_mass: str = ""
    db_chemical_id: str = ""
    n_neighbors: int = 0
    status: str = ""
    report: str = ""
    msi_level: str = ""
    msi_label: str = ""


def _read_confidence_extra(report_dir: Path) -> dict:
    """Read msi_level/label + display_name from a report's confidence_inputs sidecar."""

    slug = report_dir.name.lower()
    p = report_dir / "data" / f"{slug}_confidence_inputs.json"
    if not p.exists():
        return {}
    try:
        d = json.loads(p.read_text(encoding="utf-8"))
    except (ValueError, OSError):
        return {}
    return {
        "msi_level": d.get("msi_level", ""),
        "msi_label": d.get("msi_label", ""),
        "display_name": d.get("display_name", ""),
    }


def render_master_index(rows: list[IndexRow], generated_at: str) -> str:
    """Render the self-contained master ``index.html`` (search + sortable table)."""

    def esc(value: Any) -> str:
        return html.escape(str(value))

    libs = sorted({r.lib for r in rows if r.lib})
    body_rows = []
    for idx, r in enumerate(rows, start=1):
        body_rows.append(
            "<tr>"
            f"<td class='num'>{idx}</td>"
            f"<td><a href='{esc(r.report)}'>{esc(r.control_id)}</a></td>"
            f"<td>{esc(r.name)}</td>"
            f"<td class='num' data-sort='{esc(r.msi_level)}' title='{esc(r.msi_label)}'>"
            f"{('L' + esc(r.msi_level)) if r.msi_level else ''}</td>"
            f"<td class='mono'>{esc(r.formula)}</td>"
            f"<td class='num'>{esc(r.exact_mass)}</td>"
            f"<td class='num'>{esc(r.lib)}</td>"
            f"<td class='mono'>{esc(r.db_chemical_id)}</td>"
            f"<td class='num'>{esc(r.n_neighbors)}</td>"
            f"<td class='status status-{esc(r.status)}'>{esc(r.status)}</td>"
            "</tr>"
        )
    lib_chips = " ".join(
        f"<span class='chip'>lib {esc(lib)}: {sum(1 for r in rows if r.lib == lib)}</span>" for lib in libs
    )
    rows_html = "\n".join(body_rows)
    total = len(rows)
    return f"""<!DOCTYPE html>
<html lang='en'>
<head>
<meta charset='utf-8'>
<meta name='viewport' content='width=device-width, initial-scale=1'>
<title>X-compound reports ({total})</title>
<style>
  :root {{ color-scheme: light dark; }}
  * {{ box-sizing: border-box; }}
  body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif;
         margin: 0; padding: 32px; line-height: 1.45; }}
  h1 {{ margin: 0 0 4px; font-size: 22px; }}
  .meta {{ color: #6a737d; font-size: 13px; margin-bottom: 16px; }}
  .chips {{ margin: 10px 0 18px; }}
  .chip {{ display: inline-block; padding: 3px 10px; margin: 2px 6px 2px 0; border-radius: 999px;
          background: rgba(127,127,127,0.14); font-size: 12px; }}
  #q {{ width: 100%; max-width: 420px; padding: 9px 12px; font-size: 14px; border-radius: 8px;
       border: 1px solid rgba(127,127,127,0.4); margin-bottom: 14px; }}
  table {{ border-collapse: collapse; width: 100%; font-size: 13px; }}
  th, td {{ text-align: left; padding: 7px 10px; border-bottom: 1px solid rgba(127,127,127,0.18); }}
  th {{ position: sticky; top: 0; cursor: pointer; user-select: none; background: rgba(127,127,127,0.10);
       backdrop-filter: blur(6px); }}
  th:hover {{ background: rgba(127,127,127,0.18); }}
  th[data-dir='asc']::after {{ content: ' \\25B2'; font-size: 10px; opacity: 0.7; }}
  th[data-dir='desc']::after {{ content: ' \\25BC'; font-size: 10px; opacity: 0.7; }}
  td.num {{ text-align: right; font-variant-numeric: tabular-nums; }}
  td.mono {{ font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace; font-size: 12px; }}
  tr:hover td {{ background: rgba(127,127,127,0.07); }}
  a {{ color: #0969da; text-decoration: none; }}
  a:hover {{ text-decoration: underline; }}
  .status {{ font-size: 12px; }}
  .status-failed {{ color: #cf222e; }}
  .status-skipped {{ color: #9a6700; }}
  .status-generated {{ color: #1a7f37; }}
</style>
</head>
<body>
  <h1>X-compound HTML reports</h1>
  <div class='meta'>{total} compounds &middot; generated {esc(generated_at)}</div>
  <div class='chips'>{lib_chips}</div>
  <input id='q' type='search' placeholder='Filter by id, name, formula, lib&hellip;' autocomplete='off'>
  <table id='t'>
    <thead><tr>
      <th data-k='n'>#</th><th>Control ID</th><th>Name</th>
      <th title='MSI identification level (1-5)'>MSI</th><th>Formula</th>
      <th>Exact mass</th><th>Lib</th><th>DB chem id</th><th>Neighbors</th><th>Status</th>
    </tr></thead>
    <tbody>
{rows_html}
    </tbody>
  </table>
<script>
  const q = document.getElementById('q');
  const tbody = document.querySelector('#t tbody');
  const rows = Array.from(tbody.rows);
  q.addEventListener('input', () => {{
    const needle = q.value.toLowerCase();
    for (const row of rows) {{
      row.style.display = row.textContent.toLowerCase().includes(needle) ? '' : 'none';
    }}
  }});
  const headers = document.querySelectorAll('#t thead th');
  headers.forEach((th, col) => {{
    let asc = true;
    th.addEventListener('click', () => {{
      headers.forEach(h => {{ if (h !== th) h.removeAttribute('data-dir'); }});
      th.setAttribute('data-dir', asc ? 'asc' : 'desc');
      const cellVal = (r) => {{
        const c = r.cells[col];
        return (c.dataset && c.dataset.sort !== undefined && c.dataset.sort !== '')
          ? c.dataset.sort : c.textContent.trim();
      }};
      const visible = rows.slice();
      visible.sort((a, b) => {{
        const x = cellVal(a), y = cellVal(b);
        const nx = parseFloat(x), ny = parseFloat(y);
        const cmp = (!isNaN(nx) && !isNaN(ny)) ? nx - ny : String(x).localeCompare(String(y));
        return asc ? cmp : -cmp;
      }});
      asc = !asc;
      visible.forEach(r => tbody.appendChild(r));
    }});
  }});
</script>
</body>
</html>
"""
