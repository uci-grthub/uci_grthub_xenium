#!/usr/bin/env python3
"""
Generate a PDF project report for a Xenium run processed by this workflow.

Sections: project information, summary statistics, panel and run, pipeline
overview (with the parameters from config.yaml), per-sample Xenium QC and
filtering, merged objects (with each group's UMAP), and the output files.

Inputs:
  --metrics     output/xenium/metrics_summary.csv (rule xenium_metrics_summary)
  --stats       JSON from src/xenium_report_stats.py (cells, genes, clusters
                and preprocessing parameters of every AnnData object)
  --experiments every region's experiment.xenium (run, panel, software)
  --umaps       merged/plots/<group>/umap.png, one per merge group (optional)

Usage:
  python src/generate_xenium_report.py --metrics output/xenium/metrics_summary.csv \
      --stats output/xenium/report_stats.json \
      --experiments data/*/*/output-*/experiment.xenium \
      --umaps output/xenium/anndata/merged/plots/*/umap.png \
      --config config.yaml --out output/Xenium_Project_Report.pdf

Requires reportlab, pandas, pyyaml and Pillow (the snakemake_dfs env).
"""
import argparse
import json
import os
from datetime import datetime

import pandas as pd
import yaml
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.pagesizes import landscape, letter
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import inch
from reportlab.lib.utils import ImageReader
from reportlab.platypus import (BaseDocTemplate, Frame, Image, NextPageTemplate,
                                PageBreak, PageTemplate, Paragraph, Spacer, Table,
                                TableStyle)

BLUE = colors.HexColor("#1f4788")
MARGIN = 0.75 * inch

styles = getSampleStyleSheet()
TITLE = ParagraphStyle("Title", parent=styles["Heading1"], fontSize=24, textColor=BLUE,
                       spaceAfter=12, alignment=TA_CENTER, fontName="Helvetica-Bold")
HEADING = ParagraphStyle("Heading", parent=styles["Heading2"], fontSize=14, textColor=BLUE,
                         spaceAfter=12, spaceBefore=12, fontName="Helvetica-Bold")
SUBHEADING = styles["Heading3"]
BODY = ParagraphStyle("Body", parent=styles["BodyText"], fontSize=10, alignment=TA_LEFT,
                      spaceAfter=10, leading=14)
CELL = ParagraphStyle("Cell", parent=styles["BodyText"], fontSize=8, leading=10,
                      wordWrap="CJK")
PATH = ParagraphStyle("Path", parent=styles["BodyText"], fontSize=7, leading=9,
                      fontName="Courier", wordWrap="CJK")

KEY_VALUE_STYLE = [
    ("BACKGROUND", (0, 0), (0, -1), colors.HexColor("#E8EEF7")),
    ("FONTNAME", (0, 0), (0, -1), "Helvetica-Bold"),
    ("VALIGN", (0, 0), (-1, -1), "TOP"),
    ("FONTSIZE", (0, 0), (-1, -1), 10),
    ("BOTTOMPADDING", (0, 0), (-1, -1), 6),
    ("TOPPADDING", (0, 0), (-1, -1), 6),
    ("GRID", (0, 0), (-1, -1), 1, colors.grey),
]


def header_table_style(font_size=9, numeric_from=None):
    style = [
        ("BACKGROUND", (0, 0), (-1, 0), BLUE),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.whitesmoke),
        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("FONTSIZE", (0, 0), (-1, -1), font_size),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
        ("TOPPADDING", (0, 0), (-1, -1), 4),
        ("GRID", (0, 0), (-1, -1), 0.5, colors.grey),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#F5F5F5")]),
    ]
    if numeric_from is not None:
        style.append(("ALIGN", (numeric_from, 0), (-1, -1), "RIGHT"))
    return TableStyle(style)


def key_value_table(rows, widths=(2.3 * inch, 4.7 * inch)):
    table = Table(rows, colWidths=list(widths))
    table.setStyle(TableStyle(KEY_VALUE_STYLE))
    return table


def fmt_int(value):
    return "" if value is None or pd.isna(value) else f"{int(round(value)):,}"


def fmt_pct(value):
    return "" if value is None or pd.isna(value) else f"{100 * value:.1f}%"


def unique(values):
    return sorted({str(v) for v in values if v not in (None, "")})


def load_experiments(paths):
    """experiment.xenium per region, keyed by <slide_id>_<region_name>."""
    runs = {}
    for path in paths:
        with open(path) as fh:
            meta = json.load(fh)
        slide = meta.get("slide_id") or ""
        region = meta.get("region_name") or ""
        runs[f"{slide}_{region}" if slide else region] = meta
    return runs


class XeniumReport:
    def __init__(self, metrics, stats, experiments, umaps, config, author):
        self.metrics = metrics.set_index("sample_id") if "sample_id" in metrics else metrics
        self.samples = stats["samples"]
        self.merged = stats["merged"]
        self.runs = experiments
        self.umaps = umaps
        self.config = config
        self.xenium = config.get("xenium", {}) or {}
        self.author = author

    # -- sections ---------------------------------------------------------

    def project_information(self):
        runs = self.runs.values()
        rows = [
            ["Generation date", datetime.now().strftime("%Y-%m-%d %H:%M")],
            ["Author", self.author or ""],
            ["Run name(s)", ", ".join(unique(r.get("run_name") for r in runs))],
            ["Run start", ", ".join(unique((r.get("run_start_time") or "")[:10] for r in runs))],
            ["Slides", ", ".join(unique(r.get("slide_id") for r in runs))],
            ["Instrument(s)", ", ".join(unique(r.get("instrument_sn") for r in runs))],
        ]
        return [Paragraph("Project Information", HEADING), key_value_table(rows)]

    def summary_statistics(self):
        cells_raw = sum(s["preprocessing"].get("n_cells_raw", 0) for s in self.samples.values())
        cells_kept = sum(s["cells"] for s in self.samples.values())
        rows = [
            ["Metric", "Value"],
            ["Regions (samples)", str(len(self.samples))],
            ["Slides", str(len(unique(r.get("slide_id") for r in self.runs.values())))],
            ["Cells detected by Xenium", fmt_int(cells_raw)],
            ["Cells retained after QC", f"{fmt_int(cells_kept)} "
                                        f"({fmt_pct(cells_kept / cells_raw) if cells_raw else ''})"],
        ]
        if "median_transcripts_per_cell" in self.metrics:
            rows.append(["Median transcripts per cell (range)",
                         self._range(self.metrics["median_transcripts_per_cell"])])
        if "median_genes_per_cell" in self.metrics:
            rows.append(["Median genes per cell (range)",
                         self._range(self.metrics["median_genes_per_cell"])])
        if self.merged:
            rows.append(["Merged objects", ", ".join(self.merged)])
        table = Table(rows, colWidths=[3.0 * inch, 3.0 * inch])
        table.setStyle(header_table_style(10))
        return [Paragraph("Summary Statistics", HEADING), table]

    def panel_and_run(self):
        run = next(iter(self.runs.values()), {})
        n_targets = (run.get("panel_num_targets_predesigned") or 0) + \
            (run.get("panel_num_targets_custom") or 0)
        rows = [
            ["Item", "Value"],
            ["Panel", str(run.get("panel_name", ""))],
            ["Panel design id", str(run.get("panel_design_id", ""))],
            ["Organism / tissue", f"{run.get('panel_organism', '')} / {run.get('panel_tissue_type', '')}"],
            ["Targets", f"{n_targets:,} ({run.get('panel_num_targets_predesigned', 0):,} predesigned, "
                        f"{run.get('panel_num_targets_custom', 0):,} custom)"],
            ["Preservation", str(run.get("preservation_method", "")).upper()],
            ["Chemistry", str(run.get("chemistry_version", ""))],
            ["Instrument software", ", ".join(unique(r.get("instrument_sw_version") for r in self.runs.values()))],
            ["Onboard analysis", ", ".join(unique(r.get("analysis_sw_version") for r in self.runs.values()))],
            ["Pixel size", f"{run.get('pixel_size', '')} µm"],
        ]
        panels = unique(r.get("panel_name") for r in self.runs.values())
        table = Table(rows, colWidths=[1.9 * inch, 5.1 * inch], repeatRows=1)
        table.setStyle(header_table_style(9))
        out = [Paragraph("Panel and Run", HEADING), table]
        if len(panels) > 1:
            out.append(Paragraph(f"Note: regions use {len(panels)} panels: {', '.join(panels)}.", BODY))
        return out

    def pipeline_overview(self):
        x = self.xenium
        # describe what the merged objects were built with, not what config
        # says now (an object built before a config change keeps its setting)
        corrections = {g: m["preprocessing"].get("batch_correction", "none")
                       for g, m in self.merged.items()}
        harmony = sorted(g for g, c in corrections.items() if c == "harmony")
        uncorrected = sorted(g for g, c in corrections.items() if c != "harmony")
        n_pcs = next((s["preprocessing"].get("n_pcs") for s in self.samples.values()), 30)
        merge_text = ("Samples in each merge group were rebuilt from raw counts, outer-joined "
                      "on genes and re-normalized with the same cell-size constant before "
                      "embedding and clustering. ")
        if harmony:
            merge_text += (f"For {', '.join(harmony)}, Harmony was run on the principal "
                           "components over sample_id; the neighbor graph, UMAP and Leiden "
                           "clusters use the corrected components while expression values "
                           "remain uncorrected. ")
        if uncorrected:
            merge_text += (f"{', '.join(uncorrected)}: no batch correction; sample_id is kept "
                           "as the batch label.")
        onboard = ", ".join(unique(r.get("analysis_sw_version") for r in self.runs.values()))
        steps = [
            (f"Xenium Onboard Analysis ({onboard})",
             "Decoding, DAPI-based cell segmentation and transcript assignment were run on the "
             "instrument. Per-region QC is summarized from each region's metrics_summary.csv and "
             "analysis_summary.html."),
            ("Per-sample AnnData (scanpy)",
             "Each region's cell_feature_matrix.h5 was loaded with gene expression features only "
             "(negative control probes and codewords removed) and cell metadata from cells.csv.gz. "
             f"Cells with fewer than {x.get('min_genes', 10)} genes and genes detected in fewer than "
             f"{x.get('min_cells', 5)} cells were removed."),
            ("Normalization",
             "Counts were scaled by cell size: counts / cell_area x "
             f"{x.get('target_area', 100)} µm², then log1p. The same constant is used for every "
             "sample and merged object, so expression values are comparable across them."),
            ("Embedding and clustering",
             f"PCA ({n_pcs} components), k-nearest-neighbor graph, UMAP and Leiden clustering "
             f"(resolution {x.get('leiden_resolution', 0.3)}). Marker genes per cluster were ranked "
             "with a Wilcoxon rank-sum test."),
        ]
        if self.merged:
            steps.append(("Merged objects", merge_text))
        out = [Paragraph("Pipeline Overview", HEADING)]
        for title, text in steps:
            out += [Paragraph(f"<b>{title}</b>", SUBHEADING), Paragraph(text, BODY)]
        return out

    def sample_qc(self):
        header = ["Sample", "Cells\n(Xenium)", "Median\ntranscripts", "Median\ngenes",
                  "Q20\ndecoded", "Transcripts\nassigned", "Cells\nafter QC", "Retained",
                  "Genes", "Clusters"]
        rows = [header]
        for sample in sorted(self.samples):
            s = self.samples[sample]
            m = self.metrics.loc[sample] if sample in self.metrics.index else pd.Series(dtype=float)
            raw = s["preprocessing"].get("n_cells_raw")
            rows.append([
                Paragraph(sample, CELL),
                fmt_int(m.get("num_cells_detected", raw)),
                fmt_int(m.get("median_transcripts_per_cell")),
                fmt_int(m.get("median_genes_per_cell")),
                fmt_pct(m.get("fraction_transcripts_decoded_q20")),
                fmt_pct(m.get("fraction_transcripts_assigned")),
                fmt_int(s["cells"]),
                fmt_pct(s["cells"] / raw) if raw else "",
                fmt_int(s["genes"]),
                fmt_int(s["clusters"]),
            ])
        widths = [2.4 * inch] + [0.82 * inch] * 9
        table = Table(rows, colWidths=widths, repeatRows=1)
        table.setStyle(header_table_style(8, numeric_from=1))
        return [
            Paragraph("Per-Sample QC and Filtering", HEADING),
            table,
            Spacer(1, 0.1 * inch),
            Paragraph("Cells (Xenium), median transcripts/genes per cell, Q20 decoded and "
                      "transcripts assigned are from Xenium Onboard Analysis. Cells after QC, "
                      "genes and clusters are from the per-sample AnnData objects. Full metrics "
                      "are in output/xenium/analysis_summary_table.html.", BODY),
        ]

    def merged_objects(self, frame_width, frame_height):
        if not self.merged:
            return []
        rows = [["Group", "Samples", "Cells", "Genes", "Clusters", "Batch correction"]]
        for group, m in self.merged.items():
            rows.append([
                Paragraph(group, CELL),
                str(len(m.get("cells_per_sample", {}))),
                fmt_int(m["cells"]), fmt_int(m["genes"]), fmt_int(m["clusters"]),
                str(m["preprocessing"].get("batch_correction", "none")),
            ])
        table = Table(rows, colWidths=[2.4 * inch, 0.9 * inch, 1.2 * inch, 0.9 * inch,
                                       0.9 * inch, 1.6 * inch], repeatRows=1)
        table.setStyle(header_table_style(9, numeric_from=1))
        out = [Paragraph("Merged Objects", HEADING), table, Spacer(1, 0.1 * inch)]
        configured = self.xenium.get("batch_correction", "none")
        stale = [g for g, m in self.merged.items()
                 if m["preprocessing"].get("batch_correction", "none") != configured]
        if stale:
            out.append(Paragraph(
                f"Note: built with a different batch correction than the current config "
                f"(xenium.batch_correction: {configured}): {', '.join(stale)}. Rebuild with "
                "--forcerun xenium_merge to match.", BODY))
        for group in self.merged:
            umap = self.umaps.get(group)
            if not umap:
                continue
            out += [PageBreak(), Paragraph(f"Merged UMAP: {group}", HEADING),
                    self._image(umap, frame_width, frame_height - 1.0 * inch)]
        return out

    def output_files(self):
        out_dir = self.config.get("paths", {}).get("output", "output")
        rows = [
            ["Output", "Path"],
            ["Xenium QC (all metrics)", f"{out_dir}/xenium/analysis_summary_table.html, .csv"],
            ["Xenium metrics (numeric)", f"{out_dir}/xenium/metrics_summary.csv"],
            ["Per-region Xenium reports", f"{out_dir}/xenium/analysis_summaries/"],
            ["Per-sample AnnData", f"{out_dir}/xenium/anndata/<sample>.h5ad"],
            ["Per-sample plots", f"{out_dir}/xenium/anndata/plots/<sample>/"],
            ["Merged AnnData", f"{out_dir}/xenium/anndata/merged/<group>.h5ad"],
            ["Merged plots", f"{out_dir}/xenium/anndata/merged/plots/<group>/"],
        ]
        rows = [rows[0]] + [[r[0], Paragraph(r[1], PATH)] for r in rows[1:]]
        table = Table(rows, colWidths=[2.3 * inch, 4.7 * inch], repeatRows=1)
        table.setStyle(header_table_style(9))
        return [
            Paragraph("Output Files", HEADING), table, Spacer(1, 0.1 * inch),
            Paragraph("AnnData layout: raw counts in layers['counts'], normalized log values in "
                      "X, centroids in obsm['spatial'], Leiden clusters in obs['leiden'] "
                      "(merged objects also keep per-sample clusters as obs['leiden_sample']).",
                      BODY),
        ]

    # -- helpers ----------------------------------------------------------

    @staticmethod
    def _range(series):
        series = pd.to_numeric(series, errors="coerce").dropna()
        if series.empty:
            return ""
        return f"{series.median():,.0f} ({series.min():,.0f}-{series.max():,.0f})"

    @staticmethod
    def _image(path, max_width, max_height):
        width, height = ImageReader(path).getSize()
        scale = min(max_width / width, max_height / height)
        return Image(path, width=width * scale, height=height * scale)

    # -- document ---------------------------------------------------------

    def build(self, out_path):
        pw, ph = letter
        lw, lh = landscape(letter)
        portrait = Frame(MARGIN, MARGIN, pw - 2 * MARGIN, ph - 2 * MARGIN, id="portrait")
        wide = Frame(MARGIN, MARGIN, lw - 2 * MARGIN, lh - 2 * MARGIN, id="landscape")

        def page_number(canvas, doc):
            canvas.setFont("Helvetica", 8)
            canvas.drawRightString(canvas._pagesize[0] - MARGIN, 0.5 * inch, f"Page {doc.page}")

        def on_portrait(canvas, doc):
            canvas.setPageSize(letter)
            page_number(canvas, doc)

        def on_landscape(canvas, doc):
            canvas.setPageSize(landscape(letter))
            page_number(canvas, doc)

        doc = BaseDocTemplate(out_path, pagesize=letter, title="Xenium Project Report",
                              author=self.author or "",
                              pageTemplates=[PageTemplate("Portrait", [portrait], onPage=on_portrait),
                                             PageTemplate("Landscape", [wide], onPage=on_landscape)])
        elements = [Paragraph("Xenium Project Report", TITLE), Spacer(1, 0.2 * inch)]
        elements += self.project_information() + [Spacer(1, 0.2 * inch)]
        elements += self.summary_statistics()
        elements += [PageBreak()] + self.panel_and_run() + self.pipeline_overview()
        elements += [NextPageTemplate("Landscape"), PageBreak()] + self.sample_qc()
        if self.merged:
            elements += [PageBreak()] + self.merged_objects(wide._width, wide._height)
        elements += [NextPageTemplate("Portrait"), PageBreak()] + self.output_files()

        os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
        doc.build(elements)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--metrics", required=True, help="metrics_summary.csv")
    parser.add_argument("--stats", required=True, help="JSON from xenium_report_stats.py")
    parser.add_argument("--experiments", nargs="+", required=True, help="experiment.xenium files")
    parser.add_argument("--umaps", nargs="*", default=[],
                        help="merged/plots/<group>/umap.png files")
    parser.add_argument("--config", default="config.yaml")
    parser.add_argument("--author", default=None, help="default: config report.author")
    parser.add_argument("--out", required=True, help="output .pdf")
    args = parser.parse_args()

    with open(args.config) as fh:
        config = yaml.safe_load(fh) or {}
    with open(args.stats) as fh:
        stats = json.load(fh)
    author = args.author or (config.get("report") or {}).get("author", "")
    umaps = {os.path.basename(os.path.dirname(p)): p for p in args.umaps}

    report = XeniumReport(pd.read_csv(args.metrics), stats, load_experiments(args.experiments),
                          umaps, config, author)
    report.build(args.out)
    print(f"report -> {args.out}")


if __name__ == "__main__":
    main()
