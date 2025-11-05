# Multi-sample CellRanger Snakemake workflow
# For processing multiple single-cell RNA-seq samples

EMAIL = "kstachel@uci.edu"
onstart:
   shell("mail -s 'STARTED' {EMAIL} < {log}")
onsuccess:
   shell("mail -s 'DONE' {EMAIL} < {log}")
onerror:
   shell("mail -s 'ERROR' {EMAIL} < {log}")

from datetime import datetime
import loompy
import glob
import os
import re

# Load configuration
configfile: "config.yaml"

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
SAMPLES = auto_samples if auto_samples else config.get("samples", [])
FASTQ_DIR = config.get("paths", {}).get("fastqs", "data/FASTQ")
TRANSCRIPTOME = config["references"]["transcriptome"]
OUTPUT_DIR = config["paths"]["output"]

# Print detected samples for debugging
print(f"Auto-detected samples: {auto_samples}")
print(f"Using samples: {SAMPLES}")
print(f"FASTQ directory: {FASTQ_DIR}")
print(f"Number of samples to process: {len(SAMPLES)}")

# Rule all - defines final outputs for all samples
rule all:
    input:
        expand(f"{OUTPUT_DIR}/cellranger/{{sample}}/outs/web_summary.html", sample=SAMPLES),
        expand(f"{OUTPUT_DIR}/cellranger/{{sample}}/outs/filtered_feature_bc_matrix", sample=SAMPLES),
        f"{OUTPUT_DIR}/multi_sample_summary.csv",
        f"{OUTPUT_DIR}/scanpy/combined_harmony_integrated.h5ad",
        f"{OUTPUT_DIR}/scanpy/inspect_integrated_anndata_combined.ipynb",
        # per-sample filtering timeline plots
        expand(f"{OUTPUT_DIR}/qc/filtering_timeline/{{sample}}.png", sample=SAMPLES)
        # loompy outputs
        # expand(f"{OUTPUT_DIR}/loom/{{sample}}.loom", sample=SAMPLES),
        # scenic outputs
        # expand(f"{OUTPUT_DIR}/scenic/{{sample}}.loom", sample=SAMPLES),
        # # CellBender outputs
        # expand(f"{OUTPUT_DIR}/cellbender/{{sample}}_filtered.h5", sample=SAMPLES),
        # Scrublet outputs
        # expand(f"{OUTPUT_DIR}/scrublet/{{sample}}_doublets.csv", sample=SAMPLES)
        # # Seurat outputs
        # expand(f"{OUTPUT_DIR}/seurat/{{sample}}_seu.rds", sample=SAMPLES),
        # expand(f"{OUTPUT_DIR}/seurat/{{sample}}_embeddings.csv", sample=SAMPLES),
        # # Velocyto outputs
        # expand(f"{OUTPUT_DIR}/velocyto/{{sample}}.loom", sample=SAMPLES)
        # # scVelo outputs
        # expand(f"{OUTPUT_DIR}/scanpy/{{sample}}_scvelo.h5ad", sample=SAMPLES),
        # # Pileup and phasing outputs
        # expand(f"{OUTPUT_DIR}/numbat/{{sample}}_allele_counts.tsv.gz", sample=SAMPLES),
        # # Numbat outputs
        # expand(f"{OUTPUT_DIR}/numbat/{{sample}}/done.txt", sample=SAMPLES)


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
        account = "sbsandme_lab"
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
        account = "sbsandme_lab"
    shell:
        """
        module load singularity/3.11.3
        singularity run /dfs9/ucightf-lab/kstachel/TOOLS/multiqc-1.20.sif multiqc {OUTPUT_DIR}/fastqc -o {OUTPUT_DIR}
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
        account = "sbsandme_lab"
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
        account = "sbsandme_lab"
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
        account = "sbsandme_lab_gpu"
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
        account = "sbsandme_lab"
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
        account = "sbsandme_lab"
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
        account = "sbsandme_lab"
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

# Rule: loompy
rule loompy:
    input:
        folder = f"{OUTPUT_DIR}/cellranger/{{sample}}"
    output:
        loom_file = f"{OUTPUT_DIR}/loom/{{sample}}.loom"
    resources:
        mem_mb = 8000,  # 8GB for loompy
        cpus = 4,
        partition = "standard"
    threads: 4
    run:
        loompy.create_from_cellranger(input.folder, f"{OUTPUT_DIR}/loom", genome="hg38")

# Rule: Seurat
rule seurat:
    input:
        matrix_dir = f"{OUTPUT_DIR}/{{sample}}/outs/filtered_feature_bc_matrix/",
        script = "src/process_seurat.R"
    output:
        seu_path = f"{OUTPUT_DIR}/seurat/{{sample}}_seu.rds"
    log:
        f"{OUTPUT_DIR}/logs/seurat_{{sample}}.log"
    threads: 4
    shell:
        '''Rscript {input.script} matrix_dir="{input.matrix_dir}" seu_path="{output.seu_path}" celltype_ref="config['celltype_ref']" > {log} 2>&1'''

# Rule: Save Seurat embeddings
rule seurat_embeddings:
    input:
        seu_path = f"{OUTPUT_DIR}/seurat/{{sample}}_seu.rds",
        nb_path = f"{OUTPUT_DIR}/numbat/{{sample}}_numbat.rds",
        script = "src/save_seurat_embeddings.R"
    output:
        seurat_embeddings = f"{OUTPUT_DIR}/seurat/{{sample}}_embeddings.csv"
    log:
        f"{OUTPUT_DIR}/logs/seurat_{{sample}}_embeddings.log"
    threads: 2
    shell:
        '''Rscript {input.script} seu_path="{input.seu_path}" nb_path="{input.nb_path}" > {log} 2>&1'''

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

# Rule: Numbat
rule numbat:
    input:
        allele_df = f"{OUTPUT_DIR}/numbat/{{sample}}_allele_counts.tsv.gz",
        matrix_file = f"{OUTPUT_DIR}/{{sample}}/outs/filtered_feature_bc_matrix/matrix.mtx.gz",
        seu_path = f"{OUTPUT_DIR}/seurat/{{sample}}_seu.rds",
        script = "src/run_numbat.R"
    output:
        done_file = f"{OUTPUT_DIR}/numbat/{{sample}}/done.txt"
    threads: 4
    shell:
        '''Rscript {input.script} seu_path="{input.seu_path}" ref_path=config["ref_path"] tau=config["tau"] read_prop=config["read_prop"] max_iter=config["max_iter"] min_LLR=config["min_LLR"] t=config["numbat_t"] cell_ceiling=config["cell_ceiling"] max_entropy=config["max_entropy"] allele_df="{input.allele_df}" matrix_file="{input.matrix_file}" out_dir="{OUTPUT_DIR}/numbat/{{wildcards.sample}}" ncores={threads} rprof_out={OUTPUT_DIR}/numbat/{{wildcards.sample}}/log.prof > {output.done_file} 2>&1'''

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
