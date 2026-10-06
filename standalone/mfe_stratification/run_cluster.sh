#!/usr/bin/env bash
#SBATCH --job-name=mfe-stratify
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
set -euo pipefail
BASE=/camp/lab/ulej/home/users/luscomben/users/iosubi/projects/structurome_blencowe
ANALYSIS="$BASE/mfe_stratification_2026-10-05"
PYTHON="${MFE_PYTHON:-$HOME/miniconda3/envs/comp-hiclip-dev/bin/python}"
"$PYTHON" "$ANALYSIS/scripts/stratify_mfe.py" \
  --samplesheet "$BASE/samplesheet_ira_data.csv" \
  --results-dir "$BASE/results_human" \
  --gtf "$ANALYSIS/reference/gencode.v33.annotation.gtf.gz" \
  --date 2026-10-05 \
  --outdir "$ANALYSIS/results"
