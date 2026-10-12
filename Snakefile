# Snakemake workflow for Xenium runs (plus idle single-cell Cell Ranger rules)
# Xenium outputs are discovered under xenium.data_root; see README.md

import glob
import json
import os
import re

from snakemake.exceptions import WorkflowError

# Load configuration; config.yaml is per-project and gitignored
if not os.path.exists("config.yaml"):
    raise WorkflowError("config.yaml not found: cp config.example.yaml config.yaml and fill it in")
configfile: "config.yaml"

# start/done/error mail goes to config email; empty disables it
EMAIL = config.get("email") or ""
onstart:
    if EMAIL:
        shell("mail -s 'STARTED' {EMAIL} < {log}")
onsuccess:
    if EMAIL:
        shell("mail -s 'DONE' {EMAIL} < {log}")
onerror:
    if EMAIL:
        shell("mail -s 'ERROR' {EMAIL} < {log}")

SLURM_ACCOUNT = config.get("slurm_account") or ""
SLURM_GPU_ACCOUNT = config.get("slurm_gpu_account") or ""
MULTIQC_SIF = config.get("multiqc_sif", "multiqc-1.20.sif")

# Auto-detect samples from FASTQ directory
def get_samples_from_fastq_dir():
    """
    Detect sample names by inspecting the FASTQ directory structure.

    Supported layouts:
    - per-sample subfolders: <FASTQ_DIR>/<sample>/<sample>_S*_L*_R1_*.fastq.gz
    - flat files in FASTQ_DIR: <FASTQ_DIR>/<sample>_S*_L*_R1_*.fastq.gz

    Returns a sorted list of unique sample names, excluding 'Undetermined'.
    """
    # Prefer configured path; fall back to default
    fastq_path = config.get("paths", {}).get("fastqs", "data/FASTQ")
    if not os.path.exists(fastq_path):
        return []

    samples = set()

    # 1) Look inside one-level subdirectories (typical 10x demux output)
    subdirs = [d for d in glob.glob(os.path.join(fastq_path, "*")) if os.path.isdir(d)]
    for d in subdirs:
        sample_dir_name = os.path.basename(d)
        if sample_dir_name == "Undetermined":
            continue
        # R1 files inside the subdir; be lenient on prefix
        r1_in_dir = glob.glob(os.path.join(d, "*_R1_*.fastq.gz"))
        if r1_in_dir:
            samples.add(sample_dir_name)

    # 2) Also consider flat layout at root (no subdirs)
    r1_files_root = glob.glob(os.path.join(fastq_path, "*_R1_*.fastq.gz"))
    for file in r1_files_root:
        basename = os.path.basename(file)
        m = re.match(r'^(.+?)_S\d+_', basename)
        if m:
            name = m.group(1)
            if name != "Undetermined":
                samples.add(name)

    return sorted(samples)

# Extract sample list and configuration
# Use auto-detected samples if available, otherwise fall back to config
auto_samples = get_samples_from_fastq_dir()
SAMPLES = auto_samples if auto_samples else (config.get("samples") or [])
FASTQ_DIR = config.get("paths", {}).get("fastqs", "data/FASTQ")
TRANSCRIPTOME = config["references"]["transcriptome"]
OUTPUT_DIR = config["paths"]["output"]

# Xenium configuration
XENIUM = config.get("xenium", {})
# merged objects: group name -> glob over sample ids
XENIUM_MERGE_GROUPS = XENIUM.get("merge_groups") or {"merged": "*"}

# sample ids and merge groups are single path components, so
# anndata/{sample}.h5ad never matches anndata/merged/{group}.h5ad
wildcard_constraints:
    sample = "[^/]+",
    group = "[^/]+"

def get_xenium_samples():
    """
    Xenium output directories under xenium.data_root, keyed by sample id.

    Each region on a slide is delivered as one output-XETG* directory holding
    an experiment.xenium manifest; the instrument has already run the onboard
    analysis, so there is nothing to count. The sample id is
    <slide>_<region name>, which stays unique when a region name is reused on
    a second slide.
    """
    import json

    root = XENIUM.get("data_root", "data")
    samples = {}
    for manifest in sorted(glob.glob(os.path.join(root, "**", "experiment.xenium"),
                                     recursive=True)):
        outdir = os.path.dirname(manifest)
        with open(manifest) as handle:
            meta = json.load(handle)
        region = meta.get("region_name") or os.path.basename(outdir)
        slide = meta.get("slide_id") or ""
        sample_id = f"{slide}_{region}" if slide else region
        if sample_id in samples:
            raise ValueError(f"Duplicate Xenium sample id '{sample_id}': "
                             f"{samples[sample_id]['dir']} and {outdir}")
        samples[sample_id] = {"dir": outdir, "slide": slide, "region": region}
    return samples

XENIUM_SAMPLES = get_xenium_samples()

# Print detected samples for debugging
print(f"Xenium samples ({len(XENIUM_SAMPLES)}): {sorted(XENIUM_SAMPLES)}")
print(f"Auto-detected samples: {auto_samples}")
print(f"Using samples: {SAMPLES}")
print(f"FASTQ directory: {FASTQ_DIR}")
print(f"Number of samples to process: {len(SAMPLES)}")

def xenium_file(name):
    """Input function: <name> inside a Xenium sample's output directory."""
    return lambda wildcards: os.path.join(XENIUM_SAMPLES[wildcards.sample]["dir"], name)

# Rule all - defines final outputs for all samples
rule all:
    input:
        f"{OUTPUT_DIR}/xenium/analysis_summary_table.csv",
        f"{OUTPUT_DIR}/xenium/analysis_summary_table.html",
        f"{OUTPUT_DIR}/xenium/metrics_summary.csv",
        expand(f"{OUTPUT_DIR}/xenium/analysis_summaries/{{sample}}_analysis_summary.html",
               sample=sorted(XENIUM_SAMPLES))


# Target rule: the single-cell (Cell Ranger) outputs; idle unless FASTQs are
# present under paths.fastqs
rule cellranger_all:
    input:
        expand(f"{OUTPUT_DIR}/cellranger/{{sample}}/outs/web_summary.html", sample=SAMPLES),
        expand(f"{OUTPUT_DIR}/cellranger/{{sample}}/outs/filtered_feature_bc_matrix", sample=SAMPLES),
        f"{OUTPUT_DIR}/multi_sample_summary.csv",
        f"{OUTPUT_DIR}/scanpy/combined_harmony_integrated.h5ad",
        f"{OUTPUT_DIR}/scanpy/inspect_integrated_anndata_combined.ipynb",
        expand(f"{OUTPUT_DIR}/qc/filtering_timeline/{{sample}}.png", sample=SAMPLES),
        expand(f"{OUTPUT_DIR}/qc/violin/{{sample}}.png", sample=SAMPLES)


# Rule: one table of the metrics shown in every analysis_summary.html
rule xenium_summary_table:
    input:
        summaries = [os.path.join(XENIUM_SAMPLES[s]["dir"], "analysis_summary.html")
                     for s in sorted(XENIUM_SAMPLES)]
    output:
        csv = f"{OUTPUT_DIR}/xenium/analysis_summary_table.csv",
        html = f"{OUTPUT_DIR}/xenium/analysis_summary_table.html"
    params:
        script = "src/aggregate_xenium_summaries.py"
    shell:
        """
        python {params.script} --summaries {input.summaries:q} \
            --csv {output.csv} --html {output.html}
        """


# Rule: stack every metrics_summary.csv (numeric, one row per sample)
rule xenium_metrics_summary:
    input:
        metrics = [os.path.join(XENIUM_SAMPLES[s]["dir"], "metrics_summary.csv")
                   for s in sorted(XENIUM_SAMPLES)]
    output:
        csv = f"{OUTPUT_DIR}/xenium/metrics_summary.csv"
    run:
        import pandas as pd

        frames = []
        for sample_id, path in zip(sorted(XENIUM_SAMPLES), input.metrics):
            df = pd.read_csv(path)
            df.insert(0, "sample_id", sample_id)
            frames.append(df)
        pd.concat(frames, ignore_index=True).to_csv(output.csv, index=False)


# Rule: copy each analysis_summary.html next to the others, named by sample
rule collect_xenium_summaries:
    input:
        xenium_file("analysis_summary.html")
    output:
        f"{OUTPUT_DIR}/xenium/analysis_summaries/{{sample}}_analysis_summary.html"
    shell:
        "cp {input:q} {output:q}"

localrules: xenium_summary_table, xenium_metrics_summary, collect_xenium_summaries


# Target rule: one processed AnnData per Xenium sample
rule xenium_anndata_all:
    input:
        expand(f"{OUTPUT_DIR}/xenium/anndata/{{sample}}.h5ad", sample=sorted(XENIUM_SAMPLES)),
        expand(f"{OUTPUT_DIR}/xenium/anndata/plots/{{sample}}/markers.csv",
               sample=sorted(XENIUM_SAMPLES)),
        expand(f"{OUTPUT_DIR}/xenium/anndata/merged/{{group}}.h5ad",
               group=sorted(XENIUM_MERGE_GROUPS)),
        expand(f"{OUTPUT_DIR}/xenium/anndata/merged/plots/{{group}}/markers.csv",
               group=sorted(XENIUM_MERGE_GROUPS)),
        f"{OUTPUT_DIR}/xenium/anndata/plots/qc_violin_by_slide.png"


def merge_group_samples(group):
    """Samples whose id matches the xenium.merge_groups pattern for group."""
    import fnmatch

    samples = [s for s in sorted(XENIUM_SAMPLES)
               if fnmatch.fnmatch(s, XENIUM_MERGE_GROUPS[group])]
    if not samples:
        raise ValueError(f"merge group '{group}' matches no sample: "
                         f"{XENIUM_MERGE_GROUPS[group]}")
    return samples


# Rule: the samples of one merge group in one AnnData, re-normalized from raw
# counts with the same fixed xenium.target_area; xenium.batch_correction
# (none | harmony) sets whether the embedding is corrected over sample_id;
# see src/merge_xenium_anndata.py
rule xenium_merge:
    input:
        h5ads = lambda wildcards: expand(f"{OUTPUT_DIR}/xenium/anndata/{{sample}}.h5ad",
                                         sample=merge_group_samples(wildcards.group))
    output:
        h5ad = f"{OUTPUT_DIR}/xenium/anndata/merged/{{group}}.h5ad"
    log:
        "logs/xenium_merge/{group}.log"
    params:
        script = "src/merge_xenium_anndata.py",
        python = XENIUM.get("python", "python"),
        target_area = XENIUM.get("target_area", 100),
        min_cells = XENIUM.get("min_cells", 5),
        resolution = XENIUM.get("leiden_resolution", 0.3),
        batch_correction = XENIUM.get("batch_correction", "none")
    threads: 8
    resources:
        mem_mb = 96000
    shell:
        """
        PYTHONNOUSERSITE=1 {params.python} {params.script} \
            --h5ads {input.h5ads} --out {output.h5ad} \
            --min-cells {params.min_cells} --resolution {params.resolution} \
            --target-area {params.target_area} \
            --batch-correction {params.batch_correction} > {log} 2>&1
        """


XENIUM_MERGED_PLOTS = ["umap.png", "umap_by_sample.png", "umap_qc.png",
                       "cluster_composition.png", "markers_dotplot.png", "markers.csv"]

# Rule: plots for a merged object (UMAP by sample and cluster, per-sample
# UMAP panels to judge batch effects, cluster composition, marker dotplot)
rule xenium_merge_plots:
    input:
        h5ad = f"{OUTPUT_DIR}/xenium/anndata/merged/{{group}}.h5ad"
    output:
        [f"{OUTPUT_DIR}/xenium/anndata/merged/plots/{{group}}/{name}"
         for name in XENIUM_MERGED_PLOTS]
    log:
        "logs/xenium_merge_plots/{group}.log"
    params:
        script = "src/plot_xenium_anndata.py",
        outdir = f"{OUTPUT_DIR}/xenium/anndata/merged/plots/{{group}}",
        python = XENIUM.get("python", "python")
    threads: 2
    resources:
        mem_mb = 48000
    shell:
        """
        PYTHONNOUSERSITE=1 {params.python} {params.script} --merged \
            --h5ad {input.h5ad} --outdir {params.outdir} > {log} 2>&1
        """


XENIUM_PLOTS = ["qc_violin.png", "qc_spatial.png", "filtering_spatial.png",
                "pca_variance.png", "umap_qc.png",
                "leiden.png", "spatial_leiden.png", "markers_dotplot.png", "markers.csv"]

# Rule: standard per-sample plots (QC, UMAP, Leiden on tissue, marker dotplot);
# see src/plot_xenium_anndata.py. cells.csv.gz feeds filtering_spatial.png
# (cells removed by QC filtering highlighted on the tissue).
rule xenium_plots:
    input:
        h5ad = f"{OUTPUT_DIR}/xenium/anndata/{{sample}}.h5ad",
        cells = xenium_file("cells.csv.gz")
    output:
        [f"{OUTPUT_DIR}/xenium/anndata/plots/{{sample}}/{name}" for name in XENIUM_PLOTS]
    log:
        f"logs/xenium_plots/{{sample}}.log"
    params:
        script = "src/plot_xenium_anndata.py",
        outdir = f"{OUTPUT_DIR}/xenium/anndata/plots/{{sample}}",
        python = XENIUM.get("python", "python")
    threads: 2
    resources:
        mem_mb = 24000
    shell:
        """
        PYTHONNOUSERSITE=1 {params.python} {params.script} \
            --h5ad {input.h5ad} --outdir {params.outdir} \
            --cells-csv {input.cells:q} > {log} 2>&1
        """


# Rule: per-slide QC violin plots pooling every sample of a slide into one
# violin (genes, transcripts, cell area per cell); see
# src/plot_xenium_qc_by_slide.py
rule xenium_qc_by_slide:
    input:
        h5ads = expand(f"{OUTPUT_DIR}/xenium/anndata/{{sample}}.h5ad",
                       sample=sorted(XENIUM_SAMPLES))
    output:
        png = f"{OUTPUT_DIR}/xenium/anndata/plots/qc_violin_by_slide.png"
    log:
        "logs/xenium_qc_by_slide.log"
    params:
        script = "src/plot_xenium_qc_by_slide.py",
        python = XENIUM.get("python", "python"),
        slide_map = json.dumps({s: XENIUM_SAMPLES[s]["slide"]
                                for s in sorted(XENIUM_SAMPLES)})
    threads: 2
    resources:
        mem_mb = 24000
    shell:
        """
        PYTHONNOUSERSITE=1 {params.python} {params.script} \
            --h5ads {input.h5ads} --slide-map {params.slide_map:q} \
            --out {output.png} > {log} 2>&1
        """


# Rule: cells/genes/clusters and preprocessing parameters of every AnnData
# object, as JSON for the report (the report env has no anndata);
# see src/xenium_report_stats.py
rule xenium_report_stats:
    input:
        h5ads = expand(f"{OUTPUT_DIR}/xenium/anndata/{{sample}}.h5ad",
                       sample=sorted(XENIUM_SAMPLES)),
        merged = expand(f"{OUTPUT_DIR}/xenium/anndata/merged/{{group}}.h5ad",
                        group=sorted(XENIUM_MERGE_GROUPS))
    output:
        json = f"{OUTPUT_DIR}/xenium/report_stats.json"
    log:
        "logs/xenium_report_stats.log"
    params:
        script = "src/xenium_report_stats.py",
        python = XENIUM.get("python", "python")
    threads: 1
    resources:
        mem_mb = 16000
    shell:
        """
        PYTHONNOUSERSITE=1 {params.python} {params.script} \
            --h5ads {input.h5ads} --merged {input.merged} \
            --out {output.json} > {log} 2>&1
        """


# Rule: PDF project report (run, panel, pipeline parameters, per-sample QC
# and filtering, merged objects with their UMAPs, output files);
# see src/generate_xenium_report.py. report.python needs reportlab, pandas,
# pyyaml and Pillow.
rule xenium_report:
    input:
        metrics = f"{OUTPUT_DIR}/xenium/metrics_summary.csv",
        stats = f"{OUTPUT_DIR}/xenium/report_stats.json",
        experiments = [os.path.join(XENIUM_SAMPLES[s]["dir"], "experiment.xenium")
                       for s in sorted(XENIUM_SAMPLES)],
        umaps = expand(f"{OUTPUT_DIR}/xenium/anndata/merged/plots/{{group}}/umap.png",
                       group=sorted(XENIUM_MERGE_GROUPS)),
        config = "config.yaml"
    output:
        pdf = f"{OUTPUT_DIR}/Xenium_Project_Report.pdf"
    log:
        "logs/xenium_report.log"
    params:
        script = "src/generate_xenium_report.py",
        python = (config.get("report") or {}).get("python", "envs/report/.pixi/envs/default/bin/python")
    shell:
        """
        PYTHONNOUSERSITE=1 {params.python} {params.script} --metrics {input.metrics} --stats {input.stats} \
            --experiments {input.experiments:q} --umaps {input.umaps} \
            --config {input.config} --out {output.pdf} > {log} 2>&1
        """

localrules: xenium_report


# Rule: Xenium sample -> AnnData (QC filter, cell-area normalization,
# PCA/UMAP/Leiden); see src/xenium_to_anndata.py
rule xenium_anndata:
    input:
        matrix = xenium_file("cell_feature_matrix.h5"),
        cells = xenium_file("cells.csv.gz")
    output:
        h5ad = f"{OUTPUT_DIR}/xenium/anndata/{{sample}}.h5ad"
    log:
        f"logs/xenium_anndata/{{sample}}.log"
    params:
        script = "src/xenium_to_anndata.py",
        xenium_dir = lambda wildcards: XENIUM_SAMPLES[wildcards.sample]["dir"],
        python = XENIUM.get("python", "python"),
        min_genes = XENIUM.get("min_genes", 10),
        target_area = XENIUM.get("target_area", 100),
        min_cells = XENIUM.get("min_cells", 5),
        resolution = XENIUM.get("leiden_resolution", 0.3)
    threads: 4
    resources:
        mem_mb = 32000
    shell:
        """
        PYTHONNOUSERSITE=1 {params.python} {params.script} \
            --xenium-dir {params.xenium_dir:q} --sample-id {wildcards.sample} \
            --out {output.h5ad} --min-genes {params.min_genes} \
            --min-cells {params.min_cells} --resolution {params.resolution} \
            --target-area {params.target_area} > {log} 2>&1
        """


# Helper function to get FASTQ files for a sample
def get_fastq_files(wildcards, read):
    """Get R1 or R2 FASTQ file for a sample"""
    import glob
    # Support both flat and per-sample subfolder layouts
    pattern_flat = f"{FASTQ_DIR}/{wildcards.sample}_S*_L*_{read}_*.fastq.gz"
    pattern_dir = f"{FASTQ_DIR}/{wildcards.sample}/*_{read}_*.fastq.gz"
    files = glob.glob(pattern_flat) + glob.glob(pattern_dir)
    if not files:
        raise ValueError(f"No {read} FASTQ files found for sample {wildcards.sample}")
    return files[0]  # Return first match

# Rule: FastQC on raw FASTQ files
rule fastqc:
    input:
        r1 = lambda wildcards: get_fastq_files(wildcards, "R1"),
        r2 = lambda wildcards: get_fastq_files(wildcards, "R2")
    output:
        r1_html = f"fastqc/{{sample}}_r1_fastqc.html",
        r1_zip = f"fastqc/{{sample}}_r1_fastqc.zip",
        r2_html = f"fastqc/{{sample}}_r2_fastqc.html",
        r2_zip = f"fastqc/{{sample}}_r2_fastqc.zip"
    threads: 2
    resources:
        mem_mb = 4000,
        cpus = 2,
        partition = "standard",
        account = SLURM_ACCOUNT
    shell:
        """
        module load fastqc/0.11.9
        mkdir -p fastqc
        fastqc -o fastqc -t {threads} {input.r1} {input.r2}
        module unload fastqc/0.11.9
        """

# Rule: MultiQC report
rule multiqc:
    input:
        expand(f"fastqc/{{sample}}_r1_fastqc.html", sample=SAMPLES),
        expand(f"fastqc/{{sample}}_r2_fastqc.html", sample=SAMPLES)
    output:
        report = f"{OUTPUT_DIR}/multiqc_report.html"
    threads: 2
    resources:
        mem_mb = 4000,
        cpus = 2,
        partition = "standard",
        account = SLURM_ACCOUNT
    shell:
        """
        module load singularity/3.11.3
        singularity run {MULTIQC_SIF} multiqc {OUTPUT_DIR}/fastqc -o {OUTPUT_DIR}
        module unload singularity/3.11.3
        """

# Rule: CellRanger count for each sample
rule cellranger_count:
    input:
        transcriptome = TRANSCRIPTOME,
        fastq_dir = f"{FASTQ_DIR}/{{sample}}"
    output:
        web_summary = f"{OUTPUT_DIR}/cellranger/{{sample}}/outs/web_summary.html",
        filtered_matrix_h5 = f"{OUTPUT_DIR}/cellranger/{{sample}}/outs/filtered_feature_bc_matrix.h5",
        filtered_matrix_dir = directory(f"{OUTPUT_DIR}/cellranger/{{sample}}/outs/filtered_feature_bc_matrix"),
        bam_file = f"{OUTPUT_DIR}/cellranger/{{sample}}/outs/possorted_genome_bam.bam",
        molecule_info = f"{OUTPUT_DIR}/cellranger/{{sample}}/outs/molecule_info.h5"
    params:
        output_dir = OUTPUT_DIR,
        localcores = config["params"]["localcores"],
        localmem = config["params"]["localmem"],
        chemistry = config["params"].get("chemistry", "auto"),
        sample_id = "{sample}"
    threads: 8
    resources:
        mem_mb = 47663,  # 48GB in MB
        mem_mb_per_cpu=6000,
        cpus = 8,
        partition = "standard",
        account = SLURM_ACCOUNT
    shell:
        """
        module load cellranger/8.0.1
        
        rm  -rf {params.sample_id}
        rm  -rf {params.output_dir}/cellranger/{params.sample_id}

        cellranger count --id={params.sample_id} \
        --transcriptome={input.transcriptome} \
        --fastqs={input.fastq_dir} \
        --create-bam=true \
        --sample={params.sample_id} \
        --localcores={params.localcores} \
        --localmem={params.localmem} \
        --chemistry={params.chemistry}

        mv {params.sample_id} {params.output_dir}/cellranger/{params.sample_id}
        
        module unload cellranger/8.0.1
        """



# Rule: Generate multi-sample summary
rule multi_sample_summary:
    input:
        metrics_summaries = expand(f"{OUTPUT_DIR}/cellranger/{{sample}}/outs/metrics_summary.csv", sample=SAMPLES)
    output:
        summary = f"{OUTPUT_DIR}/multi_sample_summary.csv"
    params:
        script = "src/filter_cellranger_out.py",
        mt_thresh = config.get('mt_thresh', 5),
        min_genes = config.get('min_genes', 200),
        min_cells = config.get('min_cells', 5),
        organism = config.get('organism', 'auto')
    threads: 8
    resources:
        mem_mb = 47663,  # 128GB in MB
        mem_mb_per_cpu=6000,
        cpus = 8,
        partition = "standard",
        account = SLURM_ACCOUNT
    run:
        import pandas as pd
        import os
        import subprocess
        import sys

        dfs = []
        for path in input.metrics_summaries:
            # Ensure the file exists and is not empty
            if not os.path.exists(path):
                print(f"Warning: metrics file not found: {path}")
                continue
            try:
                df = pd.read_csv(path)
            except pd.errors.EmptyDataError:
                print(f"Warning: metrics file empty: {path}")
                continue

            # Extract sample name from the path: .../cellranger/{sample}/outs/metrics_summary.csv
            parts = path.split(os.sep)
            sample = None
            # Try to locate the 'cellranger' dir and pick the next component as sample
            if 'cellranger' in parts:
                idx = parts.index('cellranger')
                if idx + 1 < len(parts):
                    sample = parts[idx + 1]
            if sample is None:
                # Fallback: use parent directory name two levels up
                sample = os.path.basename(os.path.dirname(os.path.dirname(path)))

            df['sample_id'] = sample

            # Call the filter script and capture the number of cells remaining
            matrix_h5 = os.path.join(OUTPUT_DIR, 'cellranger', sample, 'outs', 'filtered_feature_bc_matrix.h5')
            cmd = [sys.executable, params.script,
                   '--adata_path', matrix_h5,
                   '--min_genes', str(params.min_genes),
                   '--min_cells', str(params.min_cells),
                   '--mt_thresh', str(params.mt_thresh),
                   '--organism', str(params.organism)]
            # Run the filter script and capture stdout/stderr separately so warnings/errors
            # printed to stderr are not mixed into the stdout result.
            proc = subprocess.run(cmd, capture_output=True, text=True)
            if proc.returncode != 0:
                # Log stderr for debugging, but do not include it in the stored output
                print(f"Warning: filter script failed for sample {sample}. stderr:\n{proc.stderr}")
                out = proc.stdout.strip()
            else:
                out = proc.stdout.strip()

            df['cells_after_filtering'] = out
            df['min_cells'] = params.min_cells
            df['min_genes'] = params.min_genes
            df['mt_thresh'] = params.mt_thresh
            dfs.append(df)

        if not dfs:
            # No data to write; create an empty file with a header
            os.makedirs(os.path.dirname(output.summary), exist_ok=True)
            with open(output.summary, 'w') as fh:
                fh.write('')
            print(f"No metrics summaries found; created empty {output.summary}")
            return

        combined = pd.concat(dfs, sort=False, ignore_index=True)
        # Pop the column
        popped_column = combined.pop('sample_id')
        # Insert the column at the front (index 0)
        combined.insert(0, 'sample_id', popped_column)
        os.makedirs(os.path.dirname(output.summary), exist_ok=True)
        combined.to_csv(output.summary, index=False)
        print(f"Wrote combined metrics to {output.summary}")


# Rule: collect web_summary.html files into a single directory for easy viewing
rule collect_web_summaries:
    input:
        web_summaries = expand(f"{OUTPUT_DIR}/cellranger/{{sample}}/outs/web_summary.html", sample=SAMPLES)
    output:
        directory(f"{OUTPUT_DIR}/web_summaries")
    params:
        outdir = f"{OUTPUT_DIR}/web_summaries"
    run:
        import os, shutil
        # Ensure output directory exists
        os.makedirs(params.outdir, exist_ok=True)
        for src in input.web_summaries:
            # Try to extract sample name from the path
            m = re.search(rf"{re.escape(OUTPUT_DIR)}/cellranger/([^/]+)/outs/web_summary.html", src)
            if m:
                sample = m.group(1)
            else:
                # Fallback: infer from parent directories
                sample = os.path.basename(os.path.dirname(os.path.dirname(src)))
            dest = os.path.join(params.outdir, f"{sample}_web_summary.html")
            shutil.copy(src, dest)
        # Touch a file to mark completion (optional; directory() is sufficient)
        open(os.path.join(params.outdir, ".done"), 'w').close()

# Rule: 10x scVI integration
rule tenx_scvi_integration:
    input:
       filtered_matrix_dirs = expand(f"{OUTPUT_DIR}/cellranger/{{sample}}/outs/filtered_feature_bc_matrix", sample=SAMPLES)
    output:
        combined_adata = f"{OUTPUT_DIR}/scanpy/combined.h5ad",
        integration_results = f"{OUTPUT_DIR}/scanpy/combined_integrated.h5ad"
    conda: "scvi-tools"
    params:
        script = "src/tenx_scvi_integration.py",
        input_dir = f"{OUTPUT_DIR}/cellranger",
        min_genes = config.get("min_genes", 300),
        min_cells = config.get("min_cells", 5),
        n_top_genes = config.get("n_top_genes", 2000),
        batch_key = config.get("batch_key", "batch"),
        output_prefix = f"{OUTPUT_DIR}/scanpy/combined"
    threads: 4
    resources:
        mem_mb = 48000,  # 32GB in MB
        cpus = 8,
        partition = "gpu",
        account = SLURM_GPU_ACCOUNT
    shell:
        """
        mkdir -p {OUTPUT_DIR}/scanpy
        python {params.script} \
            --filtered_matrix_dirs {input.filtered_matrix_dirs} \
            --output_prefix {params.output_prefix} \
            --min_genes {params.min_genes} \
            --min_cells {params.min_cells} \
            --n_top_genes {params.n_top_genes} \
            --batch_key {params.batch_key}
        """

# Rule: 10x harmony integration
rule tenx_harmony_integration:
    input:
       filtered_matrix_dirs = expand(f"{OUTPUT_DIR}/cellranger/{{sample}}/outs/filtered_feature_bc_matrix", sample=SAMPLES)
    output:
        combined_adata = f"{OUTPUT_DIR}/scanpy/combined.h5ad",
        integration_results = f"{OUTPUT_DIR}/scanpy/combined_harmony_integrated.h5ad"
    conda: "scvi-tools"
    params:
        script = "src/tenx_harmony_integration.py",
        input_dir = f"{OUTPUT_DIR}/cellranger",
        min_genes = config.get("min_genes", 300),
        min_cells = config.get("min_cells", 5),
        n_top_genes = config.get("n_top_genes", 2000),
        batch_key = config.get("batch_key", "batch"),
        output_prefix = f"{OUTPUT_DIR}/scanpy/combined",
        metadata = config.get("metadata", None)
    threads: 8
    resources:
        mem_mb = 48000,  # 32GB in MB
        cpus = 8,
        account = SLURM_ACCOUNT
    shell:
        """
        mkdir -p {OUTPUT_DIR}/scanpy
        python {params.script} \
            --filtered_matrix_dirs {input.filtered_matrix_dirs} \
            --output_prefix {params.output_prefix} \
            --min_genes {params.min_genes} \
            --min_cells {params.min_cells} \
            --n_top_genes {params.n_top_genes} \
            --batch_key {params.batch_key} \
            --metadata {params.metadata}
        """

# Rule: 10x harmony notebook
rule tenx_harmony_notebook:
    input:
       integration_results = f"{OUTPUT_DIR}/scanpy/combined_harmony_integrated.h5ad"
    output:
        integration_notebook = f"{OUTPUT_DIR}/scanpy/inspect_integrated_anndata_combined.ipynb"
    conda: "scvi-tools"
    params:
        script = "src/submit_harmony_integration.sh",
        input_dir = f"{OUTPUT_DIR}/cellranger",
        min_genes = config.get("min_genes", 300),
        min_cells = config.get("min_cells", 5),
        n_top_genes = config.get("n_top_genes", 2000),
        batch_key = config.get("batch_key", "batch"),
        output_prefix = f"{OUTPUT_DIR}/scanpy/combined"
    threads: 4
    resources:
        mem_mb = 32000,  # 32GB in MB
        cpus = 8,
        account = SLURM_ACCOUNT
    shell:
        """
        mkdir -p {OUTPUT_DIR}/scanpy
        echo {input.integration_results}
        {params.script} --no-integration \
        --min-genes {params.min_genes} \
        {params.output_prefix}
        """


# Rule: per-sample filtering timeline plot
rule plot_filtering_timeline:
    input:
        matrix = f"{OUTPUT_DIR}/cellranger/{{sample}}/outs/filtered_feature_bc_matrix.h5"
    output:
        png = f"{OUTPUT_DIR}/qc/filtering_timeline/{{sample}}.png"
    params:
        script = "src/plot_filtering_timeline.py",
        batch_key = config.get("batch_key", "batch"),
        min_genes = config.get("min_genes", 200),
        min_cells = config.get("min_cells", 3)
    threads: 4
    resources:
        mem_mb = 24000,  # 24GB in MB
        cpus = 4,
        partition = "standard",
        account = SLURM_ACCOUNT
    shell:
        """
        mkdir -p {OUTPUT_DIR}/qc
        python {params.script} \
            --input {input.matrix} \
            --batch-key {params.batch_key} \
            --batch-value {wildcards.sample} \
            --min-genes {params.min_genes} \
            --min-cells {params.min_cells} \
            --out {output.png}
        """

# Rule: per-sample QC violin plots (genes, transcripts, % mitochondrial
# counts per cell); see src/plot_cellranger_qc.py
rule plot_cellranger_qc:
    input:
        matrix = f"{OUTPUT_DIR}/cellranger/{{sample}}/outs/filtered_feature_bc_matrix.h5"
    output:
        png = f"{OUTPUT_DIR}/qc/violin/{{sample}}.png"
    params:
        script = "src/plot_cellranger_qc.py",
        sample_id = "{sample}"
    threads: 2
    resources:
        mem_mb = 8000,
        cpus = 2,
        partition = "standard",
        account = SLURM_ACCOUNT
    shell:
        """
        mkdir -p {OUTPUT_DIR}/qc/violin
        python {params.script} \
            --matrix {input.matrix} \
            --sample-id {params.sample_id} \
            --out {output.png}
        """


# Rule: Velocyto
rule velocyto:
    input:
        sample_folder = f"{OUTPUT_DIR}/{{sample}}"
    output:
        loom = f"{OUTPUT_DIR}/velocyto/{{sample}}.loom"
    threads: 4
    shell:
        """
        module load velocyto/0.17.17
        velocyto run10x -m config["references"]["repeat_mask"] {input.sample_folder} config["gtf"] > {output.loom}
        module unload velocyto/0.17.17
        """

# Rule: scVelo
rule scvelo:
    input:
        loom_file = f"{OUTPUT_DIR}/velocyto/{{sample}}.loom",
        anndata_file = f"{OUTPUT_DIR}/scanpy/{{sample}}.h5ad",
        seurat_embeddings = f"{OUTPUT_DIR}/seurat/{{sample}}_embeddings.csv"
    output:
        scvelo_h5ad = f"{OUTPUT_DIR}/scanpy/{{sample}}_scvelo.h5ad"
    threads: 2
    shell:
        '''python src/compute_velocity.py {input.anndata_file} {input.loom_file} {input.seurat_embeddings}'''

# Rule: pileup and phasing
rule pileup_and_phasing:
    input:
        bam = f"{OUTPUT_DIR}/{{sample}}/outs/possorted_genome_bam.bam",
        barcodes = f"{OUTPUT_DIR}/{{sample}}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz",
        script = "src/pileup_and_phase.R"
    output:
        allele_df = f"{OUTPUT_DIR}/numbat/{{sample}}_allele_counts.tsv.gz"
    threads: 4
    shell:
        '''Rscript {input.script} --label {wildcards.sample} --samples {wildcards.sample} --bams {input.bam} --barcodes {input.barcodes} --gmap config["gmap"] --snpvcf config["snpvcf"] --paneldir config["paneldir"] --outdir {OUTPUT_DIR}/numbat/{wildcards.sample} --ncores {threads}'''

# Rule: SCENIC analysis
include: "rules/scenic.smk"

# Rule: CellBender remove-background
rule cellbender:
    input:
        h5 = f"{OUTPUT_DIR}/cellranger/{{sample}}/outs/raw_feature_bc_matrix.h5"
    output:
        clean_h5 = f"{OUTPUT_DIR}/cellbender/{{sample}}_filtered.h5"
    conda:
        "envs/cellbender_env.yaml"
    threads: 8
    resources:
        mem_mb = 64000,
        cpus = 8,
        partition = "standard"
    params:
        expected_cells = config.get("cellbender_expected_cells", 3000),
        total_droplets = config.get("cellbender_total_droplets", 20000)
    shell:
        """
        cellbender remove-background \
            --input {input.h5} \
            --output {output.clean_h5} \
            --expected-cells {params.expected_cells} \
            --total-droplets-included {params.total_droplets} \
            --cpu-threads {threads}
        """

# Rule: Scrublet doublet detection
rule scrublet:
    input:
        matrix = f"{OUTPUT_DIR}/cellranger/{{sample}}/outs/filtered_feature_bc_matrix.h5"
    output:
        doublets = f"{OUTPUT_DIR}/scrublet/{{sample}}_doublets.csv"
    threads: 2
    resources:
        mem_mb = 8000,
        cpus = 2,
        partition = "standard"
    params:
        script = "src/run_scrublet.py"
    shell:
        """
        mkdir -p {OUTPUT_DIR}/scrublet
        python {params.script} --input {input.matrix} --output {output.doublets}
        """
