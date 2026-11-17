# Cell types of every Xenium section on the tissue, one panel per annotation
# method: Azimuth (per sample, rules/azimuth.smk) and the label transfers run
# per merge group outside Snakemake (src/submit_tangram.sh,
# src/submit_seurat_transfer.sh). The transfer CSVs are inputs only; no rule
# here makes them.
#
# Included from the Snakefile after rules/azimuth.smk; uses XENIUM,
# XENIUM_SAMPLES, XENIUM_MERGE_GROUPS, OUTPUT_DIR, AZIMUTH, AZIMUTH_DIR and
# AZIMUTH_LABELS.

CELLTYPE_METHODS = config.get("celltype_methods") or {}
# method -> predictions CSV of a merge group ({group} is filled in)
CELLTYPE_TRANSFER = CELLTYPE_METHODS.get("transfer") or {
    "tangram": f"{OUTPUT_DIR}/xenium/tangram/{{group}}_predictions.csv",
    "seurat_sct": f"{OUTPUT_DIR}/xenium/seurat_transfer/{{group}}_predictions_sct.csv",
}
CELLTYPE_PLOT = "celltype_methods_spatial.png"


def sample_merge_group(sample):
    """The merge group of a sample; the most specific pattern when several match."""
    import fnmatch

    groups = [g for g, pattern in XENIUM_MERGE_GROUPS.items()
              if fnmatch.fnmatch(sample, pattern)]
    if not groups:
        raise ValueError(f"sample {sample} is in no xenium.merge_groups group")
    return max(groups, key=lambda g: len(XENIUM_MERGE_GROUPS[g]))


def celltype_transfer_csvs(wildcards):
    group = sample_merge_group(wildcards.sample)
    return {method: path.format(group=group) for method, path in CELLTYPE_TRANSFER.items()}


# Target rule: the methods-comparison plot for every Xenium section
rule xenium_celltype_methods_all:
    input:
        expand(f"{OUTPUT_DIR}/xenium/anndata/plots/{{sample}}/{CELLTYPE_PLOT}",
               sample=sorted(XENIUM_SAMPLES))


# Rule: one section's cell types per method on the tissue (score below
# celltype_methods.min_score = unassigned); see
# src/plot_celltype_methods_spatial.py
rule xenium_celltype_methods_spatial:
    input:
        unpack(lambda wildcards: {
            "h5ad": f"{OUTPUT_DIR}/xenium/anndata/{wildcards.sample}.h5ad",
            "azimuth": f"{AZIMUTH_DIR}/{wildcards.sample}_predictions.csv",
            **celltype_transfer_csvs(wildcards)})
    output:
        f"{OUTPUT_DIR}/xenium/anndata/plots/{{sample}}/{CELLTYPE_PLOT}"
    log:
        "logs/celltype_methods/{sample}.log"
    params:
        script = "src/plot_celltype_methods_spatial.py",
        python = XENIUM.get("python", "python"),
        transfer = lambda wildcards: " ".join(
            f"{method}={path}" for method, path in celltype_transfer_csvs(wildcards).items()),
        azimuth_label = AZIMUTH_LABELS[0],
        label = CELLTYPE_METHODS.get("label", "cell_type_015"),
        min_score = CELLTYPE_METHODS.get("min_score", 0.5)
    threads: 1
    resources:
        mem_mb = 16000
    shell:
        """
        PYTHONNOUSERSITE=1 {params.python} {params.script} --h5ad {input.h5ad} \
            --azimuth {input.azimuth} --azimuth-label {params.azimuth_label:q} \
            --transfer {params.transfer} --label {params.label} \
            --min-score {params.min_score} --out {output} > {log} 2>&1
        """


# obs column written per transfer method (a key of CELLTYPE_TRANSFER)
TRANSFER_OBS_COLUMNS = {"tangram": "tangram_cell_type", "seurat_sct": "seurat_celltype"}

# merged-object directories under anndata/ that get transfer obs columns:
# merged/ (rule xenium_merge) and merged_nocorrection/ (same cells, made
# outside Snakemake)
TRANSFER_OBS_MERGED_DIRS = ["merged", "merged_nocorrection"]


def transfer_obs_done(method):
    """.done flags of every per-sample and merged object for one method."""
    return (expand(f"{OUTPUT_DIR}/xenium/anndata/{method}_obs/{{sample}}.done",
                   sample=sorted(XENIUM_SAMPLES))
            + expand(f"{OUTPUT_DIR}/xenium/anndata/{method}_obs/{{mdir}}/{{group}}.done",
                     mdir=TRANSFER_OBS_MERGED_DIRS, group=sorted(XENIUM_MERGE_GROUPS)))


# Target rule: Tangram cell types in obs of every per-sample and merged AnnData
rule xenium_tangram_obs_all:
    input:
        transfer_obs_done("tangram")


# Target rule: Seurat (SCT) cell types in obs of every per-sample and merged
# AnnData
rule xenium_seurat_obs_all:
    input:
        transfer_obs_done("seurat_sct")


# Rule: add <column> and <column>_score (TRANSFER_OBS_COLUMNS) to obs of
# <sample>.h5ad in place (mtime kept); the .done flag is the output. See
# src/add_transfer_obs.py
rule xenium_transfer_obs:
    input:
        h5ad = f"{OUTPUT_DIR}/xenium/anndata/{{sample}}.h5ad",
        predictions = lambda wildcards: celltype_transfer_csvs(wildcards)[wildcards.method]
    output:
        touch(f"{OUTPUT_DIR}/xenium/anndata/{{method}}_obs/{{sample}}.done")
    wildcard_constraints:
        method = "|".join(TRANSFER_OBS_COLUMNS),
        sample = "[^/]+"
    log:
        "logs/{method}_obs/{sample}.log"
    params:
        script = "src/add_transfer_obs.py",
        python = XENIUM.get("python", "python"),
        column = lambda wildcards: TRANSFER_OBS_COLUMNS[wildcards.method],
        label = CELLTYPE_METHODS.get("label", "cell_type_015"),
        min_score = CELLTYPE_METHODS.get("min_score", 0.5)
    threads: 1
    resources:
        mem_mb = 8000
    shell:
        """
        PYTHONNOUSERSITE=1 {params.python} {params.script} --h5ad {input.h5ad} \
            --method {wildcards.method} --predictions {input.predictions} \
            --column {params.column} --label {params.label} \
            --min-score {params.min_score} > {log} 2>&1
        """


# Rule: as xenium_transfer_obs, for anndata/<mdir>/<group>.h5ad
rule xenium_transfer_obs_merged:
    input:
        h5ad = f"{OUTPUT_DIR}/xenium/anndata/{{mdir}}/{{group}}.h5ad",
        predictions = lambda wildcards: CELLTYPE_TRANSFER[wildcards.method].format(group=wildcards.group)
    output:
        touch(f"{OUTPUT_DIR}/xenium/anndata/{{method}}_obs/{{mdir}}/{{group}}.done")
    wildcard_constraints:
        method = "|".join(TRANSFER_OBS_COLUMNS),
        mdir = "|".join(TRANSFER_OBS_MERGED_DIRS)
    log:
        "logs/{method}_obs/{mdir}/{group}.log"
    params:
        script = "src/add_transfer_obs.py",
        python = XENIUM.get("python", "python"),
        column = lambda wildcards: TRANSFER_OBS_COLUMNS[wildcards.method],
        label = CELLTYPE_METHODS.get("label", "cell_type_015"),
        min_score = CELLTYPE_METHODS.get("min_score", 0.5)
    threads: 1
    resources:
        mem_mb = 16000
    shell:
        """
        PYTHONNOUSERSITE=1 {params.python} {params.script} --h5ad {input.h5ad} \
            --method {wildcards.method} --predictions {input.predictions} \
            --column {params.column} --label {params.label} \
            --min-score {params.min_score} > {log} 2>&1
        """
