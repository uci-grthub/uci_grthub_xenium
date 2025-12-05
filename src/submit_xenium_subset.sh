#!/bin/bash
#SBATCH --job-name=xenium_combine
#SBATCH -A SBSANDME_LAB
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=6G
#SBATCH --time=08:00:00
#SBATCH --error=logs/subset_xenium.%j.err
#SBATCH --output=logs/subset_xenium.%j.out

. ~/.mymambainit-24.3.0
mamba activate samui

python src/subset_xenium_by_coordinates_and_combine.py -x data/output-XETG00221__0069982__Aged_cohort__20251003__181018
python src/subset_xenium_by_coordinates_and_combine.py -x data/output-XETG00221__0069979__Adult_cohort__20251003__181017
