#!/usr/bin/env bash
#SBATCH --job-name=mfe-all
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
set -euo pipefail
BASE=/camp/lab/ulej/home/users/luscomben/users/iosubi/projects/structurome_blencowe
ANALYSIS="$BASE/mfe_stratification_2026-10-05"
PYTHON="${MFE_PYTHON:-$HOME/miniconda3/envs/comp-hiclip-dev/bin/python}"
"$PYTHON" "$ANALYSIS/scripts/prepare_all_samples.py" \
 --results-dir "$BASE/results_human" --existing-analysis-dir "$ANALYSIS/results" \
 --gtf "$ANALYSIS/reference/gencode.v33.annotation.gtf.gz" --outdir "$ANALYSIS/results_all_samples"
"$PYTHON" "$ANALYSIS/scripts/plot_observed_mfe.py" --analysis-dir "$ANALYSIS/results_all_samples"
"$PYTHON" "$ANALYSIS/scripts/plot_pooled_proportions.py" --plot-dir "$ANALYSIS/results_all_samples/observed_mfe_violins"
"$PYTHON" "$ANALYSIS/scripts/plot_mfe_with_shuffled.py" --analysis-dir "$ANALYSIS/results_all_samples"
