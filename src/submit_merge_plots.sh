#!/bin/bash
#SBATCH --job-name=xenium_merge_plots
#SBATCH -A sbsandme_lab
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --time=08:00:00
#SBATCH --error=logs/merge_plots_batch.%j.err
#SBATCH --output=logs/merge_plots_batch.%j.out

# Build merged AnnData + merged plots for one or more merge groups
# (config.yaml xenium.merge_groups). Usage:
#   sbatch src/submit_merge_plots.sh <group> [group ...]
# Extra snakemake options go in SNAKEMAKE_ARGS, e.g. to rebuild existing merged
# objects after a config change:
#   SNAKEMAKE_ARGS="--forcerun xenium_merge" sbatch src/submit_merge_plots.sh <group>
# Submit from the repository root.
cd "${SLURM_SUBMIT_DIR:-.}"

if [ $# -eq 0 ]; then
    echo "usage: sbatch $0 <group> [group ...]" >&2
    exit 1
fi

TARGETS=""
for g in "$@"; do
    TARGETS="$TARGETS output/xenium/anndata/merged/plots/$g/markers.csv"
done

echo "Targets:$TARGETS"
snakemake --cores 8 --resources mem_mb=100000 --rerun-triggers mtime --rerun-incomplete \
    ${SNAKEMAKE_ARGS:-} -- $TARGETS
