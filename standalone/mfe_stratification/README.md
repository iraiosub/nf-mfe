# Standalone MFE stratification

## Preferred figures: observed MFE violin plots

The current requested plots include **all 22 completed samples**, regardless of
date. They are in `results_all_samples/observed_mfe_violins/` locally and in
the cluster analysis folder. The all-sample run completed successfully (job `59548382`), covering
**8,837,915 input pairs**; sample counts, MFE inclusion totals and stacked-bar
fractions were verified. The earlier `results_2026-10-05/` and cluster
`results/` folders contain the initial 13-sample run. They show
**observed MFE only**, with a box plot inside each violin:

- `mfe_by_method_exon_class.pdf/png`: one panel per exon/intron/mixed/unannotated
  class, comparing methods.
- `mfe_by_method_loop_size.pdf/png`: one panel per genomic loop bin, comparing
  methods on observed MFE.
- `mfe_by_method_loop_and_exon.pdf/png`: one panel per method, comparing exonic
  and intronic MFE within each loop bin.
- `mfe_by_method.pdf/png`: overall observed MFE comparison.
- `loop_sizes_by_method.pdf/png`: genomic loop-size violins across all methods
  on a log10(1 + nt) axis, with `method_loop_size_summary.tsv` in raw nt.
- `mfe_by_method_special_pairs.pdf/png`: special pair categories.
- `samples/SAMPLE/mfe_by_annotation_class.pdf/png` and
  `samples/SAMPLE/mfe_by_loop_size_bin.pdf/png`: per-sample stratifications.
- `exon_intron_proportions_by_sample.pdf/png`: stacked proportions for every sample.
- `exon_intron_proportions_by_method.pdf/png`: stacked proportions pooled by assay.
- Matching `_counts.tsv` and `_fractions.tsv`: exact counts and denominators.
- `observed_mfe_stratified_summary.tsv`: full-data numerical summaries.

Stacked bars include all input rows, independently of the MFE filter, and
retain mixed/partial and unannotated categories. The older eight samples map
to SPLASH using their source paths in `samplesheet_human.csv`; `paris_2016`
maps to PARIS.

Create the all-sample analysis from every top-level `*_mfe.tsv` final result:

```bash
sbatch --output=all-samples-%j.log scripts/run_all_samples_cluster.sh
```

`prepare_all_samples.py` inventories final tables, checks the annotation hash
and cached input metadata, reuses matching annotated tables through symlinks,
and annotates only missing samples. It writes a complete inventory and manifest.
The all-sample output folder must initially be empty. No RNA folding is run.

Regenerate these figures without annotation or folding:

```bash
python plot_observed_mfe.py --analysis-dir /path/to/completed/results
# Or, from the cluster analysis folder:
sbatch --output=violins-%j.log scripts/plot_cluster.sh
```

Violin widths are equal across groups. KDE uses a fixed seed and up to 10,000
pairs per group with bandwidth factor 0.25. Box medians, quartiles, whiskers and
tables use every retained pair. These figures supersede the initial multi-metric
boxplot figures described below; those initial files are retained as prior outputs.

Postprocess existing final `SAMPLE_mfe.tsv` files. No Nextflow changes, sequence
extraction, RNA folding, or shuffle generation. Dependencies: Python >=3.8 and
`requirements.txt` (pandas supplies python-dateutil).

## Current analysis

Date: **2026-10-05**, interpreted in Europe/London. Samples are the intersection
of `sample_id` entries in `samplesheet_ira_data.csv` and final result files whose
modification date matches that date. The `.csv` samplesheet is actually tab
separated; the script detects its delimiter. Missing listed results cause an
error, and older results are recorded as skipped. Modification date is a proxy
for generation date; copied or touched files may have misleading dates.

Cluster folder:
`/camp/lab/ulej/home/users/luscomben/users/iosubi/projects/structurome_blencowe/mfe_stratification_2026-10-05/`

It contains `scripts/`, `reference/`, and a separate `results/` folder. Annotation
is the full `gencode.v33.annotation.gtf.gz` copied from Dropbox's
`Ira Iosub/comp_hiclip/ref/`. Use GRCh38 coordinates; the script does not lift over.

Submit from the cluster analysis folder:

```bash
sbatch --output=analysis-%j.log scripts/run_cluster.sh
```

For another location, run directly:

```bash
python stratify_mfe.py \
  --samplesheet /path/to/samplesheet_ira_data.csv \
  --results-dir /path/to/results_human \
  --gtf /path/to/gencode.v33.annotation.gtf.gz \
  --date 2026-10-05 \
  --outdir /path/to/new/empty/results
```

The output directory must be empty. The script reads only final merged tables,
never intermediate chunks, and processes input in chunks of 100,000 rows.

## Definitions

Input arms are BED-style zero-based, half-open `[ll, lr)` and `[rl, rr)` intervals,
as used by this repo's PREPARE_BED process. GTF coordinates are converted from
one-based inclusive. `chr` prefixes and mitochondrial M/MT are normalized.
Annotation is strand-specific and includes all GENCODE gene/transcript types.

- **Exonic:** each entire arm fits within an individual exon. Arms can occupy
  different exons, transcripts, or genes. Adjacent/overlapping exons are not
  merged to manufacture containment. This is an annotation-wide definition;
  it does not require both arms to belong to one common transcript.
- **Intronic:** each entire arm fits inside a gene span on its strand and has
  no overlap with any exon on that strand across all annotated isoforms.
  This conservative definition excludes alternatively exonic segments.
- **Mixed or partial:** partial exon overlap, exon/intron pairs, or pairs with
  one annotated and one unannotated arm.
- **Unannotated:** neither arm has gene containment or exon overlap on its
  strand. This includes intergenic regions and unknown annotation contigs.

Loop size is the genomic inter-arm gap:
`max(ll, rl) - min(lr, rr)` for same-chromosome, same-strand, nonoverlapping
arms. It is independent of arm ordering. This is not a spliced transcript
loop length or the duplexfold paired-site gap. Bins are `[0,50)`, `[50,100)`,
`[100,500)`, `[500,1000)`, `[1000,5000)`, `[5000,10000)`, `[10000,50000)`,
and `50000+` nt. Overlapping, opposite-strand and trans-chromosomal pairs
have separate categories and no numerical loop size.

Plots and summaries exclude nonfinite observed MFE and positive observed MFE,
matching the existing plotting script's positive-MFE filter. Zero MFE is kept.
All input rows remain in annotated tables with an `included_in_plots` flag.
No MFE values are recalculated or normalized by arm length; comparisons are
descriptive and may reflect differences in arm length or assay composition.

## Outputs

Each sample folder contains:

- `annotated_pairs.tsv.gz`: coordinates, identifiers, existing MFE metrics,
  arm annotations, pair class, loop size/bin and inclusion flag. Sequence and
  structure columns are omitted to keep this analysis compact.
- `loop_size_bin.png/pdf`: separate loop-size boxplots.
- `annotation_class.png/pdf`: separate exonic/intronic/mixed/unannotated boxplots.
- `summary.tsv`: counts, means, medians and quartiles by loop bin, annotation
  class, and their joint stratification.

Available observed, shuffled-mean, delta, z-score and flipped-arm MFE columns
receive individual panels. Boxes show median/IQR and 1.5-IQR whiskers;
visual outliers are hidden but remain in summaries. Absent categories are
omitted from plots. No significance tests or pooled cross-sample tests are run.

Top-level `all_samples_summary.tsv`, `sample_qc.tsv` and `manifest.json` record
sample-level results, exclusions, exact selected/skipped files, arguments,
completion timestamp and annotation SHA256. Keep the manifest with results.

Validation:

```bash
python -m unittest discover -s standalone/mfe_stratification -v
```

## Method comparison

After annotation finishes:

```bash
python compare_methods.py --analysis-dir /path/to/completed/results
# Cluster alternative, from the analysis folder:
sbatch --output=methods-%j.log scripts/compare_cluster.sh
```

`results/method_comparison/` contains separate `mfe_by_method.png/pdf` and
`loop_sizes_by_method.png/pdf`, sample and method statistics, loop-bin fractions,
and the explicit sample-to-method mapping. Methods are CAR-SPLASH, KARR-seq,
PARIS, RIC-seq and STAU1 hiCLIP. Both pooled pair distributions and individual
sample medians are shown. Method tables also report the mean/median of sample
medians, so sample contributions can be compared without pair-count weighting.

Loop comparisons use **all input pairs**, independently of MFE filtering;
numerical loop distributions use same-chromosome, same-strand nonoverlapping
pairs only. The stacked-bin panel retains all special categories. The ECDF and
sample-median panel use `log10(1 + loop size)` to retain zero gaps. Genomic gap
calculation is tested on minus-strand genes and both arm orders.

RIC-seq is HeLa while the other assays are HEK293/HEK293T; cell type, arm lengths
and assay composition confound method comparisons. STAU1 high/low are subsets,
not assumed independent biological replicates. Comparisons are descriptive.

## Completed run: 5 October 2026

Both cluster jobs completed successfully (annotation `59546870`, comparison
`59547855`). Analysed 13 samples, **7,802,735 input pairs**. Annotation boundary
and minus-strand tests passed; output QC totals and method-bin fractions were
checked. Plots and summary tables are also available locally in
`standalone/mfe_stratification/results_2026-10-05/` (git-ignored). Full annotated
pair tables remain in the cluster results folder.

Pooled descriptive medians:

| Assay | Observed MFE (kcal/mol) | Genomic loop size (nt) |
| --- | ---: | ---: |
| CAR-SPLASH | -14.4 | 3 |
| KARR-seq | -28.4 | 227 |
| PARIS | -17.0 | 92 |
| RIC-seq | -27.9 | 600 |
| STAU1 hiCLIP | -25.4 | 139 |

CAR-SPLASH is dominated by short gaps; RIC-seq has the largest median gap.
KARR-seq and RIC-seq have the most negative pooled observed MFE medians.
These are descriptive results, not estimates of intrinsic assay performance:
raw MFE depends on arm length/composition, cell types differ, and pooled
statistics weight datasets by their number of pairs. The sample-median panels
provide a complementary view. MFE and loop distributions use the different
inclusion rules described above.

## Replicates pooled within sample type

`observed_mfe_violins/exon_intron_proportions_pooled_replicates.pdf/png` pools
raw counts across replicate entries, then calculates proportions. The 22 samples
form 11 groups. Method, cell type and compartment remain distinct: GM12878,
HEK293 and HeLa cytoplasmic SPLASH each have their own group; K562 chromatin
is separate. STAU1 high/low subsets and the older PARIS dataset also remain
separate. Mapping, counts and fractions are saved alongside the plot.

```bash
python plot_pooled_proportions.py --plot-dir /path/to/results_all_samples/observed_mfe_violins
```

## Observed MFE with shuffled means beside it

`results_all_samples/mfe_with_shuffled_violins/` contains updated overall,
exon-class, loop-bin, joint exon/loop and per-sample MFE figures. Each method
has two side-by-side violins in the same method colour: observed on the left
(darker), per-pair `mean_shuffled_mfe` on the right (lighter), each with a box
plot inside. No flipped-arm values are plotted.

Both distributions use exactly the same pairs: finite observed MFE <=0 and
finite precomputed shuffled mean. `paired_plot_qc.tsv` records any missing
shuffled means. `mfe_stratified_summary.tsv` reports separate observed and
shuffled-mean statistics. This shows the distribution of per-pair means across
pairs, not the distribution of individual shuffles. It requires no folding.

```bash
python plot_mfe_with_shuffled.py --analysis-dir /path/to/results_all_samples
# On the cluster, from the analysis folder:
sbatch --output=shuffled-%j.log scripts/plot_shuffled_cluster.sh
```

Paired-violin run completed successfully (job `59549551`): 22 samples,
8,798,960 matched pairs, no missing shuffled means among retained observed
pairs. Observed/shuffled counts match within every plotted stratum.
