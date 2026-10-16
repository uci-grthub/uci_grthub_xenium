#!/bin/bash
#SBATCH --job-name=celltype_methods
#SBATCH -A sbsandme_lab
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --error=logs/celltype_methods_batch.%j.err
#SBATCH --output=logs/celltype_methods_batch.%j.out

# Cell types per annotation method (Azimuth, Tangram, Seurat SCT) on the
# tissue for every Xenium section (rule xenium_celltype_methods_all,
# rules/celltype_methods.smk). The Tangram and Seurat predictions must exist
# (src/submit_tangram.sh, src/submit_seurat_transfer.sh). Usage:
#   sbatch src/submit_celltype_methods.sh
# Extra snakemake options go in SNAKEMAKE_ARGS, e.g. to redraw every plot:
#   SNAKEMAKE_ARGS="--forcerun xenium_celltype_methods_spatial" sbatch src/submit_celltype_methods.sh
# Submit from the repository root.
cd "${SLURM_SUBMIT_DIR:-.}"

snakemake --cores 4 --resources mem_mb=64000 --rerun-triggers mtime --rerun-incomplete \
    ${SNAKEMAKE_ARGS:-} -- xenium_celltype_methods_all
