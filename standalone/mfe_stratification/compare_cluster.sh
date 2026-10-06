#!/usr/bin/env bash
#SBATCH --job-name=mfe-methods
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
set -euo pipefail
BASE=/camp/lab/ulej/home/users/luscomben/users/iosubi/projects/structurome_blencowe
ANALYSIS="$BASE/mfe_stratification_2026-10-05"
PYTHON="${MFE_PYTHON:-$HOME/miniconda3/envs/comp-hiclip-dev/bin/python}"
"$PYTHON" "$ANALYSIS/scripts/compare_methods.py" --analysis-dir "$ANALYSIS/results"
