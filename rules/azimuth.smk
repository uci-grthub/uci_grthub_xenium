# Azimuth cell type annotation of the per-sample Xenium AnnData objects.
#
# Azimuth ships no mouse skin reference, so one is built from a labeled
# CELLxGENE scRNA-seq atlas (azimuth.reference_url), restricted to the Xenium
# panel genes and to the cell types in azimuth.keep_labels. Each sample's raw
# counts are mapped onto it with RunAzimuth. R steps run in the pixi env in
# envs/azimuth (azimuth.rscript); Python steps in xenium.python. Counts pass
# from Python to R as Matrix Market directories (src/h5ad_to_mtx.py).
#
# Included from the Snakefile; uses XENIUM, XENIUM_SAMPLES and OUTPUT_DIR.

AZIMUTH = config.get("azimuth") or {}
AZIMUTH_LABELS = AZIMUTH.get("label_columns") or ["cell_type"]
AZIMUTH_DIR = f"{OUTPUT_DIR}/xenium/azimuth"
AZIMUTH_REFERENCE_H5AD = "data/reference/azimuth/reference.h5ad"
AZIMUTH_PLOTS = ["azimuth_spatial.png", "azimuth_umap.png", "azimuth_scores.png",
                 "azimuth_vs_leiden.png"]


def xenium_gene_panel():
    """gene_panel.json of the first sample; every sample shares one panel."""
    return os.path.join(XENIUM_SAMPLES[sorted(XENIUM_SAMPLES)[0]]["dir"], "gene_panel.json")


# Target rule: Azimuth labels for every Xenium sample
rule xenium_azimuth_all:
    input:
        expand(f"{OUTPUT_DIR}/xenium/anndata/annotated/{{sample}}.h5ad",
               sample=sorted(XENIUM_SAMPLES)),
        expand(f"{OUTPUT_DIR}/xenium/anndata/plots/{{sample}}/{{name}}",
               sample=sorted(XENIUM_SAMPLES), name=AZIMUTH_PLOTS),
        f"{AZIMUTH_DIR}/reference/reference_umap.png"


# Rule: the reference atlas as published on CELLxGENE
rule azimuth_reference_download:
    output:
        AZIMUTH_REFERENCE_H5AD
    params:
        url = AZIMUTH.get("reference_url", "")
    shell:
        """
        if [ -z "{params.url}" ]; then
            echo "azimuth.reference_url is not set in config.yaml" >&2; exit 1
        fi
        wget -q -O {output:q}.part {params.url:q} && mv {output:q}.part {output:q}
        """

localrules: azimuth_reference_download


# Rule: reference raw counts for the kept cell types (at most
# azimuth.max_cells_per_label each), panel genes only, renamed to the panel's
# symbols so they match the Xenium feature names
rule azimuth_reference_export:
    input:
        h5ad = AZIMUTH_REFERENCE_H5AD,
        panel = xenium_gene_panel()
    output:
        temp(directory(f"{AZIMUTH_DIR}/reference_mtx"))
    log:
        "logs/azimuth/reference_export.log"
    params:
        script = "src/h5ad_to_mtx.py",
        python = XENIUM.get("python", "python"),
        label_column = AZIMUTH_LABELS[0],
        keep_labels = json.dumps(AZIMUTH.get("keep_labels") or []),
        max_cells = AZIMUTH.get("max_cells_per_label", 2000),
        labels = AZIMUTH_LABELS
    threads: 2
    resources:
        mem_mb = 48000
    shell:
        """
        PYTHONNOUSERSITE=1 {params.python} {params.script} \
            --h5ad {input.h5ad} --layer raw --panel {input.panel:q} \
            --label-column {params.label_column:q} --keep-labels {params.keep_labels:q} \
            --max-cells-per-label {params.max_cells} --obs-columns {params.labels:q} \
            --outdir {output} > {log} 2>&1
        """


# Rule: SCTransform/PCA/UMAP on the reference and package it for RunAzimuth;
# see src/build_azimuth_reference.R
rule azimuth_reference_build:
    input:
        f"{AZIMUTH_DIR}/reference_mtx"
    output:
        ref = f"{AZIMUTH_DIR}/reference/ref.Rds",
        annoy = f"{AZIMUTH_DIR}/reference/idx.annoy",
        umap = f"{AZIMUTH_DIR}/reference/reference_umap.png"
    log:
        "logs/azimuth/reference_build.log"
    params:
        script = "src/build_azimuth_reference.R",
        rscript = AZIMUTH.get("rscript", "envs/azimuth/.pixi/envs/default/bin/Rscript"),
        outdir = f"{AZIMUTH_DIR}/reference",
        labels = ",".join(AZIMUTH_LABELS),
        dims = AZIMUTH.get("dims", 50)
    threads: 4
    resources:
        mem_mb = 64000
    shell:
        """
        {params.rscript} {params.script} --mtx-dir {input} --outdir {params.outdir} \
            --label-columns {params.labels:q} --dims {params.dims} > {log} 2>&1
        """


# Rule: a sample's raw counts (layers["counts"]) for R
rule xenium_azimuth_export:
    input:
        h5ad = f"{OUTPUT_DIR}/xenium/anndata/{{sample}}.h5ad"
    output:
        temp(directory(f"{AZIMUTH_DIR}/query_mtx/{{sample}}"))
    log:
        "logs/azimuth/export/{sample}.log"
    params:
        script = "src/h5ad_to_mtx.py",
        python = XENIUM.get("python", "python")
    threads: 1
    resources:
        mem_mb = 16000
    shell:
        """
        PYTHONNOUSERSITE=1 {params.python} {params.script} \
            --h5ad {input.h5ad} --layer counts --outdir {output} > {log} 2>&1
        """


# Rule: map one sample onto the reference; see src/run_azimuth.R
rule xenium_azimuth:
    input:
        mtx = f"{AZIMUTH_DIR}/query_mtx/{{sample}}",
        ref = f"{AZIMUTH_DIR}/reference/ref.Rds",
        annoy = f"{AZIMUTH_DIR}/reference/idx.annoy"
    output:
        csv = f"{AZIMUTH_DIR}/{{sample}}_predictions.csv"
    log:
        "logs/azimuth/{sample}.log"
    params:
        script = "src/run_azimuth.R",
        rscript = AZIMUTH.get("rscript", "envs/azimuth/.pixi/envs/default/bin/Rscript"),
        reference = f"{AZIMUTH_DIR}/reference"
    threads: 4
    resources:
        mem_mb = 64000
    shell:
        """
        {params.rscript} {params.script} --mtx-dir {input.mtx} \
            --reference {params.reference} --out {output.csv} > {log} 2>&1
        """


# Rule: predictions into obs of a copy of the sample AnnData (cells below
# azimuth.min_prediction_score labeled unassigned), plus plots;
# see src/add_azimuth_labels.py
rule xenium_azimuth_annotate:
    input:
        h5ad = f"{OUTPUT_DIR}/xenium/anndata/{{sample}}.h5ad",
        csv = f"{AZIMUTH_DIR}/{{sample}}_predictions.csv"
    output:
        h5ad = f"{OUTPUT_DIR}/xenium/anndata/annotated/{{sample}}.h5ad",
        plots = [f"{OUTPUT_DIR}/xenium/anndata/plots/{{sample}}/{name}"
                 for name in AZIMUTH_PLOTS]
    log:
        "logs/azimuth/annotate/{sample}.log"
    params:
        script = "src/add_azimuth_labels.py",
        python = XENIUM.get("python", "python"),
        outdir = f"{OUTPUT_DIR}/xenium/anndata/plots/{{sample}}",
        labels = AZIMUTH_LABELS,
        min_score = AZIMUTH.get("min_prediction_score", 0.5)
    threads: 2
    resources:
        mem_mb = 24000
    shell:
        """
        PYTHONNOUSERSITE=1 {params.python} {params.script} \
            --h5ad {input.h5ad} --predictions {input.csv} --labels {params.labels:q} \
            --min-score {params.min_score} --out {output.h5ad} --outdir {params.outdir} > {log} 2>&1
        """
