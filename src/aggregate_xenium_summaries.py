#!/usr/bin/env python3
"""
Aggregate Xenium analysis_summary.html files into one table.

Each analysis_summary.html embeds its report as a JSON object (`const data =
{...}`). The metric tables shown on the page are stored under keys ending in
`_table`, each with `rows` of [label, value]; alerts are under `alarms`. This
script pulls those rows from every summary and writes:

  - a CSV with one row per metric (section, metric) and one column per sample
  - the same table as a standalone HTML page

Values are kept as formatted in the report (e.g. "94.2%", "44,187").

Usage:
  python src/aggregate_xenium_summaries.py \
      --summaries data/**/output-*/analysis_summary.html \
      --csv output/xenium/analysis_summary_table.csv \
      --html output/xenium/analysis_summary_table.html
"""
import argparse
import html
import json
import os
import re

import pandas as pd

# Report tables in the order the page shows them; any other *_table found is
# appended after these.
TABLE_ORDER = [
    "run_info_table",
    "region_info_table",
    "panel_spec_table",
    "software_table",
    "decoding_yield_table",
    "segmentation_table",
    "stain_cell_segmentation_table",
    "neg_controls_table",
    "debug_info_table",
    "debug_info_metadata_table",
]


def load_report_data(path):
    """The JSON object embedded in an analysis_summary.html as `const data`."""
    with open(path) as handle:
        text = handle.read()
    m = re.search(r"const data = (\{.*?\})\s*;?\s*\n", text, re.S)
    if not m:
        raise ValueError(f"no embedded report data in {path}")
    return json.loads(m.group(1))


def cell_text(value):
    """Table cells are strings, or lists of strings for multi-part cells."""
    if isinstance(value, list):
        return "; ".join(cell_text(v) for v in value)
    return "" if value is None else str(value)


def table_title(key):
    return key.removesuffix("_table").replace("_", " ").capitalize()


def extract_metrics(data):
    """[(section, metric, value)] from every metric table and the alarms."""
    table_keys = [k for k in TABLE_ORDER if k in data]
    table_keys += sorted(k for k in data
                         if k.endswith("_table") and k not in TABLE_ORDER)

    records = []
    for key in table_keys:
        table = data[key]
        if not isinstance(table, dict):
            continue
        for row in table.get("rows") or []:
            if len(row) >= 2:
                records.append((table_title(key), cell_text(row[0]),
                                cell_text(row[1])))

    alarms = (data.get("alarms") or {}).get("alarms") or []
    raised = [a.get("title", a.get("parent", "")) for a in alarms if a.get("raised")]
    records.append(("Alerts", "Raised alerts", "; ".join(raised)))
    return records


def sample_label(data, path):
    """<slide>_<region>, from the region info table, else the output dir name."""
    region = dict((r[0], cell_text(r[1]))
                  for r in (data.get("region_info_table") or {}).get("rows") or [])
    if region.get("Region name"):
        slide = region.get("Slide ID")
        return f"{slide}_{region['Region name']}" if slide else region["Region name"]
    return os.path.basename(os.path.dirname(path))


def aggregate(paths):
    columns = {}
    order = []
    for path in sorted(paths):
        data = load_report_data(path)
        label = sample_label(data, path)
        if label in columns:
            raise ValueError(f"duplicate sample label {label} ({path})")
        series = {}
        for section, metric, value in extract_metrics(data):
            series[(section, metric)] = value
            if (section, metric) not in order:
                order.append((section, metric))
        columns[label] = series

    table = pd.DataFrame(columns)
    table = table.reindex(pd.MultiIndex.from_tuples(order, names=["section", "metric"]))
    return table[sorted(table.columns)]


def write_html(table, path):
    samples = list(table.columns)
    head = "".join(f"<th>{html.escape(s)}</th>" for s in samples)
    body = []
    current = None
    for (section, metric), row in table.iterrows():
        if section != current:
            body.append(f'<tr class="section"><th colspan="{len(samples) + 1}">'
                        f"{html.escape(section)}</th></tr>")
            current = section
        cells = "".join(f"<td>{html.escape('' if pd.isna(v) else str(v))}</td>"
                        for v in row)
        body.append(f'<tr><th class="metric">{html.escape(metric)}</th>{cells}</tr>')

    page = f"""<!DOCTYPE html>
<html><head><meta charset="utf-8"><title>Xenium analysis summaries</title>
<style>
body {{ font-family: sans-serif; font-size: 13px; margin: 16px; }}
table {{ border-collapse: collapse; }}
th, td {{ border: 1px solid #ccc; padding: 3px 6px; white-space: nowrap; }}
thead th {{ position: sticky; top: 0; background: #eee; }}
th.metric {{ position: sticky; left: 0; background: #fafafa; text-align: left;
            font-weight: normal; white-space: normal; min-width: 260px; }}
tr.section th {{ background: #dde6f0; text-align: left; }}
td {{ text-align: right; }}
</style></head><body>
<h1>Xenium analysis summaries ({len(samples)} samples)</h1>
<table><thead><tr><th>Metric</th>{head}</tr></thead>
<tbody>
{chr(10).join(body)}
</tbody></table></body></html>
"""
    with open(path, "w") as handle:
        handle.write(page)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--summaries", nargs="+", required=True,
                        help="analysis_summary.html files")
    parser.add_argument("--csv", required=True, help="output CSV")
    parser.add_argument("--html", help="output HTML table")
    args = parser.parse_args()

    table = aggregate(args.summaries)
    os.makedirs(os.path.dirname(os.path.abspath(args.csv)), exist_ok=True)
    table.to_csv(args.csv)
    if args.html:
        os.makedirs(os.path.dirname(os.path.abspath(args.html)), exist_ok=True)
        write_html(table, args.html)
    print(f"{table.shape[1]} samples x {table.shape[0]} metrics -> {args.csv}")


if __name__ == "__main__":
    main()
