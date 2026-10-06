#!/usr/bin/env python3
"""Annotate and plot existing duplex MFE tables; never compute MFE."""
import argparse
import gzip
import hashlib
import json
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from dateutil.tz import gettz

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

METRICS = ['mfe', 'mean_shuffled_mfe', 'delta_mfe', 'zscore_mfe',
           'mfe_lseq_flipped', 'mfe_rseq_flipped']
CLASSES = ['exonic', 'intronic', 'mixed_or_partial', 'unannotated']
LOOPS = ['0–49', '50–99', '100–499', '500–999', '1000–4999',
         '5000–9999', '10000–49999', '50000+', 'overlapping',
         'opposite_strand', 'trans_chromosomal']


def chrom(value):
    value = str(value)[3:] if str(value).startswith('chr') else str(value)
    return 'MT' if value in ('M', 'MT') else value


class IntervalIndex:
    """Individual intervals: prefix maxima allow containment/overlap lookup."""
    def __init__(self, intervals):
        intervals = sorted(set(intervals))
        self.starts = np.array([s for s, e in intervals], dtype=np.int64)
        self.ends = np.maximum.accumulate([e for s, e in intervals]) if intervals else np.array([], dtype=np.int64)

    def query(self, starts, ends, containment):
        i = np.searchsorted(self.starts, starts if containment else ends,
                            side='right' if containment else 'left') - 1
        out = np.zeros(len(starts), dtype=bool)
        valid = i >= 0
        out[valid] = (self.ends[i[valid]] >= ends[valid] if containment
                      else self.ends[i[valid]] > starts[valid])
        return out


def load_gtf(path):
    raw = defaultdict(list)
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 9 or fields[2] not in ('exon', 'gene'):
                continue
            raw[(fields[2], chrom(fields[0]), fields[6])].append((int(fields[3]) - 1, int(fields[4])))
    if not any(k[0] == 'exon' for k in raw) or not any(k[0] == 'gene' for k in raw):
        raise ValueError('GTF must contain exon and gene features (use full GENCODE annotation).')
    return {k: IntervalIndex(v) for k, v in raw.items()}


def annotate(df, index):
    df = df.copy()
    for col in ['ll', 'lr', 'rl', 'rr']:
        x = pd.to_numeric(df[col], errors='raise')
        if x.isna().any() or (x != np.floor(x)).any():
            raise ValueError(f'Invalid coordinates in {col}')
        df[col] = x.astype('int64')
    for arm, start, end in [('l', 'll', 'lr'), ('r', 'rl', 'rr')]:
        if ((df[start] < 0) | (df[end] <= df[start])).any():
            raise ValueError(f'Invalid {arm} arm BED intervals')
        df[arm + 'chr'] = df[arm + 'chr'].map(chrom)
        if not df[arm + 'strand'].isin(['+', '-']).all():
            raise ValueError('Annotation requires + or - arm strands')
        exon, overlap, gene = [np.zeros(len(df), dtype=bool) for _ in range(3)]
        for (c, s), positions in df.groupby([arm + 'chr', arm + 'strand']).indices.items():
            starts, ends = df[start].to_numpy()[positions], df[end].to_numpy()[positions]
            for feature, dest, contained in [('exon', exon, True), ('exon', overlap, False), ('gene', gene, True)]:
                idx = index.get((feature, c, s))
                if idx is not None:
                    dest[positions] = idx.query(starts, ends, contained)
        df[arm + '_annotation'] = np.select([exon, gene & ~overlap, overlap | gene],
                                            ['exonic', 'intronic', 'mixed_or_partial'], default='unannotated')
    l, r = df.l_annotation, df.r_annotation
    df['annotation_class'] = np.select([(l == 'exonic') & (r == 'exonic'),
                                         (l == 'intronic') & (r == 'intronic'),
                                         (l == 'unannotated') & (r == 'unannotated')],
                                        ['exonic', 'intronic', 'unannotated'], default='mixed_or_partial')
    same_chr = df.lchr == df.rchr
    same_strand = df.lstrand == df.rstrand
    gap = np.maximum(df.ll, df.rl) - np.minimum(df.lr, df.rr)
    valid = same_chr & same_strand & (gap >= 0)
    df['loop_size_nt'] = gap.where(valid)
    bins = pd.cut(df.loop_size_nt, [0, 50, 100, 500, 1000, 5000, 10000, 50000, np.inf],
                  labels=LOOPS[:8], right=False).astype('object')
    df['loop_size_bin'] = bins.where(valid, np.select([~same_chr, ~same_strand],
                                  ['trans_chromosomal', 'opposite_strand'], default='overlapping'))
    return df


def summarize(df, sample, group_cols, metrics):
    rows = []
    for key, part in df.groupby(group_cols, observed=True, sort=False):
        key = key if isinstance(key, tuple) else (key,)
        for metric in metrics:
            vals = part[metric].to_numpy(dtype=float)
            vals = vals[np.isfinite(vals)]
            rows.append(dict(sample=sample, **dict(zip(group_cols, key)), metric=metric,
                             n_rows=len(part), n_valid=len(vals),
                             mean=np.mean(vals) if len(vals) else np.nan,
                             median=np.median(vals) if len(vals) else np.nan,
                             q25=np.quantile(vals, .25) if len(vals) else np.nan,
                             q75=np.quantile(vals, .75) if len(vals) else np.nan))
    return rows


def plot(df, sample, column, order, outdir, metrics):
    fig, axes = plt.subplots(len(metrics), 1, figsize=(13, 3.5 * len(metrics)), squeeze=False)
    for ax, metric in zip(axes[:, 0], metrics):
        data, labels = [], []
        for label in order:
            vals = df.loc[df[column] == label, metric].to_numpy(dtype=float)
            vals = vals[np.isfinite(vals)]
            if len(vals):
                data.append(vals)
                labels.append(f'{label}\nn={len(vals):,}')
        if data:
            ax.boxplot(data, labels=labels, showfliers=False)
        else:
            ax.text(.5, .5, 'No finite values', ha='center', transform=ax.transAxes)
        ax.set_ylabel(metric + (' (kcal/mol)' if metric != 'zscore_mfe' else ''))
        ax.tick_params(axis='x', labelsize=8)
        ax.grid(axis='y', alpha=.2)
    fig.suptitle(f'{sample}: {column}\nBoxes: median and IQR; whiskers: 1.5×IQR; outliers hidden')
    fig.tight_layout(rect=(0, 0, 1, .96))
    for ext in ('png', 'pdf'):
        fig.savefig(outdir / f'{column}.{ext}', dpi=160, bbox_inches='tight')
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--samplesheet', required=True)
    p.add_argument('--results-dir', required=True)
    p.add_argument('--gtf', required=True)
    p.add_argument('--outdir', required=True)
    p.add_argument('--date', default='2026-10-05', help='Final result modification date in Europe/London, or all')
    p.add_argument('--chunksize', type=int, default=100000)
    args = p.parse_args()
    if args.date != 'all':
        datetime.strptime(args.date, '%Y-%m-%d')
    out = Path(args.outdir)
    if out.exists() and any(out.iterdir()):
        raise ValueError('Output directory must be empty to avoid mixing runs')
    sheet = pd.read_csv(args.samplesheet, sep=None, engine='python', dtype=str)
    if 'sample_id' not in sheet or sheet.sample_id.duplicated().any():
        raise ValueError('Samplesheet needs unique sample_id values')
    selected, skipped = [], []
    for sample in sheet.sample_id:
        if not sample or Path(sample).name != sample or sample in ('.', '..'):
            raise ValueError('Unsafe sample_id')
        path = Path(args.results_dir) / f'{sample}_mfe.tsv'
        if not path.exists():
            raise FileNotFoundError(path)
        date = datetime.fromtimestamp(path.stat().st_mtime, gettz('Europe/London')).date().isoformat()
        entry = dict(sample=sample, path=str(path.resolve()), modification_date=date, bytes=path.stat().st_size)
        (selected if args.date == 'all' or date == args.date else skipped).append(entry)
    if not selected:
        raise ValueError('No samplesheet-listed final results match the requested date')
    out.mkdir(parents=True, exist_ok=True)
    index = load_gtf(args.gtf)
    all_stats, qc = [], []
    for entry in selected:
        sample, path = entry['sample'], entry['path']
        print(f'Annotating {sample}', flush=True)
        dest = out / sample
        dest.mkdir()
        columns = pd.read_csv(path, sep='\t', nrows=0).columns
        required = {'lchr', 'll', 'lr', 'lstrand', 'rchr', 'rl', 'rr', 'rstrand', 'mfe'}
        if required - set(columns):
            raise ValueError(f'{sample}: missing {sorted(required - set(columns))}')
        metrics = [m for m in METRICS if m in columns]
        usecols = list(required | set(metrics) | (set(columns) & {'name', 'row_id', 'empirical_p_lower'}))
        compact = []
        first = True
        total = positive = missing = 0
        for chunk in pd.read_csv(path, sep='\t', usecols=usecols, chunksize=args.chunksize,
                                  dtype={'lchr': str, 'rchr': str}):
            annotated = annotate(chunk, index)
            for metric in metrics:
                annotated[metric] = pd.to_numeric(annotated[metric], errors='raise')
            finite = np.isfinite(annotated.mfe)
            positive += int((finite & (annotated.mfe > 0)).sum())
            missing += int((~finite).sum())
            total += len(annotated)
            annotated['included_in_plots'] = finite & (annotated.mfe <= 0)
            annotated.to_csv(dest / 'annotated_pairs.tsv.gz', sep='\t', index=False,
                             compression='gzip', mode='w' if first else 'a', header=first)
            first = False
            compact.append(annotated.loc[annotated.included_in_plots,
                           ['loop_size_bin', 'annotation_class'] + metrics])
        df = pd.concat(compact, ignore_index=True)
        stats = []
        for groups in [['loop_size_bin'], ['annotation_class'], ['loop_size_bin', 'annotation_class']]:
            rows = summarize(df, sample, groups, metrics)
            for row in rows:
                row['stratification'] = '+'.join(groups)
            stats.extend(rows)
        pd.DataFrame(stats).to_csv(dest / 'summary.tsv', sep='\t', index=False)
        all_stats.extend(stats)
        for column, order in [('loop_size_bin', LOOPS), ('annotation_class', CLASSES)]:
            plot(df, sample, column, order, dest, metrics)
        counts = df.annotation_class.value_counts().to_dict()
        qc.append(dict(sample=sample, total_rows=total, plotted_rows=len(df),
                       positive_mfe_excluded=positive, nonfinite_mfe_excluded=missing, **counts))
    pd.DataFrame(all_stats).to_csv(out / 'all_samples_summary.tsv', sep='\t', index=False)
    pd.DataFrame(qc).fillna(0).to_csv(out / 'sample_qc.tsv', sep='\t', index=False)
    manifest = dict(arguments=vars(args), selected=selected, skipped=skipped,
                    gtf_sha256=hashlib.sha256(Path(args.gtf).read_bytes()).hexdigest(),
                    completed=datetime.now(gettz('Europe/London')).isoformat())
    (out / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(f'Completed {len(selected)} samples: {out}', flush=True)


if __name__ == '__main__':
    main()
