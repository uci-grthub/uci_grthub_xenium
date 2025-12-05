#!/bin/bash
#SBATCH --job-name=integrate_w_scvi
#SBATCH -A SBSANDME_LAB_GPU
#SBATCH -p gpu
#SBATCH --gres=gpu:1
#SBATCH --nodes=4
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=08:00:00
#SBATCH --error=logs/scvi_integration_master.%j.err
#SBATCH --output=logs/scvi_integration_master.%j.out

# Small RNA processing workflow submission script
python /scripts/tenx_scvi_integration.py --output_prefix output/scanpy/combined_filtered_tables
