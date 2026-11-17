#!/bin/bash
#SBATCH --job-name=tangram_env
#SBATCH -A sbsandme_lab
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/tangram_env.%j.out
#SBATCH --error=logs/tangram_env.%j.err

# Build the Python env for src/run_tangram.py (envs/tangram).
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"

pixi install --manifest-path envs/tangram/pixi.toml
PYTHONNOUSERSITE=1 envs/tangram/.pixi/envs/default/bin/python -c \
    'import tangram, torch, scanpy; print("tangram", tangram.__version__, "torch", torch.__version__)'
