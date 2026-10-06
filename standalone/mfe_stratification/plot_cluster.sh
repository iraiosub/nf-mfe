#!/usr/bin/env bash
#SBATCH --job-name=mfe-violins
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
set -euo pipefail
ANALYSIS=/camp/lab/ulej/home/users/luscomben/users/iosubi/projects/structurome_blencowe/mfe_stratification_2026-10-05
PYTHON="${MFE_PYTHON:-$HOME/miniconda3/envs/comp-hiclip-dev/bin/python}"
"$PYTHON" "$ANALYSIS/scripts/plot_observed_mfe.py" --analysis-dir "$ANALYSIS/results"
