#!/bin/bash
#SBATCH --job-name=azimuth_env
#SBATCH -A sbsandme_lab
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/azimuth_env.%j.out
#SBATCH --error=logs/azimuth_env.%j.err

# Build the R env for rules/azimuth.smk (envs/azimuth): conda packages from
# pixi.lock, then Azimuth and its GitHub-only dependencies.
set -euo pipefail

export MAKEFLAGS="-j${SLURM_CPUS_PER_TASK}"
pixi install --manifest-path envs/azimuth/pixi.toml
pixi run --manifest-path envs/azimuth/pixi.toml install-azimuth
envs/azimuth/.pixi/envs/default/bin/Rscript -e 'library(Azimuth); packageVersion("Azimuth")'
