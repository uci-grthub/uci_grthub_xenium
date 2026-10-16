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


# Target rule: Tangram cell types in obs of every per-sample AnnData
rule xenium_tangram_obs_all:
    input:
        expand(f"{OUTPUT_DIR}/xenium/anndata/tangram_obs/{{sample}}.done",
               sample=sorted(XENIUM_SAMPLES))


# Rule: add tangram_cell_type and tangram_cell_type_score to obs of
# <sample>.h5ad in place (mtime kept); the .done flag is the output. See
# src/add_tangram_obs.py
rule xenium_tangram_obs:
    input:
        h5ad = f"{OUTPUT_DIR}/xenium/anndata/{{sample}}.h5ad",
        tangram = lambda wildcards: celltype_transfer_csvs(wildcards)["tangram"]
    output:
        touch(f"{OUTPUT_DIR}/xenium/anndata/tangram_obs/{{sample}}.done")
    log:
        "logs/tangram_obs/{sample}.log"
    params:
        script = "src/add_tangram_obs.py",
        python = XENIUM.get("python", "python"),
        label = CELLTYPE_METHODS.get("label", "cell_type_015"),
        min_score = CELLTYPE_METHODS.get("min_score", 0.5)
    threads: 1
    resources:
        mem_mb = 8000
    shell:
        """
        PYTHONNOUSERSITE=1 {params.python} {params.script} --h5ad {input.h5ad} \
            --tangram {input.tangram} --label {params.label} \
            --min-score {params.min_score} > {log} 2>&1
        """


# merged-object directories under anndata/ that get Tangram obs columns:
# merged/ (rule xenium_merge) and merged_nocorrection/ (same cells, made
# outside Snakemake)
TANGRAM_OBS_MERGED_DIRS = ["merged", "merged_nocorrection"]


# Target rule: Tangram cell types in obs of every merged object
rule xenium_tangram_obs_merged_all:
    input:
        expand(f"{OUTPUT_DIR}/xenium/anndata/tangram_obs/{{mdir}}/{{group}}.done",
               mdir=TANGRAM_OBS_MERGED_DIRS, group=sorted(XENIUM_MERGE_GROUPS))


# Rule: as xenium_tangram_obs, for anndata/<mdir>/<group>.h5ad
rule xenium_tangram_obs_merged:
    input:
        h5ad = f"{OUTPUT_DIR}/xenium/anndata/{{mdir}}/{{group}}.h5ad",
        tangram = lambda wildcards: CELLTYPE_TRANSFER["tangram"].format(group=wildcards.group)
    output:
        touch(f"{OUTPUT_DIR}/xenium/anndata/tangram_obs/{{mdir}}/{{group}}.done")
    wildcard_constraints:
        mdir = "|".join(TANGRAM_OBS_MERGED_DIRS)
    log:
        "logs/tangram_obs/{mdir}/{group}.log"
    params:
        script = "src/add_tangram_obs.py",
        python = XENIUM.get("python", "python"),
        label = CELLTYPE_METHODS.get("label", "cell_type_015"),
        min_score = CELLTYPE_METHODS.get("min_score", 0.5)
    threads: 1
    resources:
        mem_mb = 16000
    shell:
        """
        PYTHONNOUSERSITE=1 {params.python} {params.script} --h5ad {input.h5ad} \
            --tangram {input.tangram} --label {params.label} \
            --min-score {params.min_score} > {log} 2>&1
        """
