#!/bin/bash
#SBATCH --job-name=build_squidpy_sif
#SBATCH -A sbsandme_lab
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/build_squidpy_sif.%j.out
#SBATCH --error=logs/build_squidpy_sif.%j.err

# Build envs/squidpy/squidpy.sif from envs/squidpy/squidpy.def, which packs the
# pixi `squidpy-dfs` environment from pixi.lock. Run from the project root:
#   sbatch src/submit_build_squidpy_sif.sh   (or: pixi run build-sif)
set -euo pipefail

source /etc/profile.d/modules.sh
module load apptainer/1.4.5

# build scratch on node-local disk, not /dfs9
scratch="/tmp/${USER}/apptainer_${SLURM_JOB_ID:-$$}"
mkdir -p "$scratch/tmp" "$scratch/cache"
trap 'rm -rf "$scratch"' EXIT
export APPTAINER_TMPDIR="$scratch/tmp"
export APPTAINER_CACHEDIR="$scratch/cache"
unset SINGULARITY_CACHEDIR
# the site bind list (/dfs9,...) would be mounted into %post before those
# directories exist in the image; binds are only needed at run time
unset SINGULARITY_BINDPATH APPTAINER_BINDPATH SINGULARITY_BIND APPTAINER_BIND

# base images are public: pull anonymously so stale ghcr.io credentials in
# ~/.apptainer/docker-config.json cannot make the registry deny the pull
echo '{}' > "$scratch/authfile.json"

# build to a temp name so a failed build never replaces a working image
out=envs/squidpy/squidpy.sif
apptainer build --fakeroot --force --authfile "$scratch/authfile.json" \
    "$out.tmp" envs/squidpy/squidpy.def
mv -f "$out.tmp" "$out"

apptainer exec "$out" python -c "import scanpy, squidpy, anndata; print('scanpy', scanpy.__version__, 'squidpy', squidpy.__version__, 'anndata', anndata.__version__)"
