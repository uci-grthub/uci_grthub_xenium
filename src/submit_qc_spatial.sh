#!/bin/bash
#SBATCH --job-name=xenium_qc_spatial
#SBATCH -A sbsandme_lab
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --error=logs/qc_spatial_batch.%j.err
#SBATCH --output=logs/qc_spatial_batch.%j.out

# Regenerate the per-sample plot sets: filtering_spatial.png is new (missing
# everywhere, so it triggers every sample's xenium_plots rerun), which also
# rewrites qc_spatial.png in the reverted single-row style.
cd "${SLURM_SUBMIT_DIR:-.}"

TARGETS=""
for f in output/xenium/anndata/*.h5ad; do
    s=$(basename "$f" .h5ad)
    TARGETS="$TARGETS output/xenium/anndata/plots/$s/filtering_spatial.png"
done

echo "Targets:$TARGETS"
snakemake --cores 2 $TARGETS --rerun-triggers mtime --rerun-incomplete
