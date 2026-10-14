# Replicate-pooled method and condition plots

```bash
python standalone/method_condition_plots/plot_conditions.py \
  --pipeline-plots /path/to/results_human/plot_method_figures \
  --outdir /path/to/condition_plots
```

Uses the pipeline's existing `plot_source.tsv.gz` and count exports. No folding
or GC recalculation is performed, and old files elsewhere in `results_human`
are not added. Replicate hybrid rows are pooled, rather than averaging sample
densities or sample medians. Source populations and normalization definitions
are inherited from the pipeline exports.

SPLASH is separated into GM12878 cytoplasm, HEK293 cytoplasm, HeLa cytoplasm,
and K562 chromatin. Replicates are pooled within each condition. CAR-SPLASH,
KARR-seq, PARIS and RIC-seq each pool their matched-condition replicates;
STAU1 high and low remain separate conditions. Unknown sample naming patterns
fail instead of silently joining unrelated conditions.

Two figures are saved as PDF/PNG: `condition_overview` (32 inches tall for the
current ten groups) and `condition_span` (20 inches tall, one group per column).
The grey styling, solid observed/dashed shuffled curves, −70 kcal/mol raw MFE
limit, and MFE per observed base-pair normalization are retained. The span
figure pools any exon/intron categories. The overview inherits the pipeline's
stratification setting. No main title or bottom explanatory text is added.

`sample_condition_mapping.tsv`, `sample_counts.tsv`, `condition_counts.tsv`,
and `plot_group_counts.tsv` describe membership and full-data counts/statistics.
`manifest.json` records the pipeline source path, size and timestamp, figure
sizes and validation. The row-level source is the existing pipeline source
joined to the mapping on `sample`; it is not duplicated. Each sample's
source/GC/observed/matched counts must agree with the pipeline, and every
cohort's span counts must sum to its pooled energy population.
