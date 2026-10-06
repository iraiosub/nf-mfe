#!/usr/bin/env python3
"""Compare assay methods using completed standalone annotations, without folding."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from stratify_mfe import METRICS, LOOPS

METHODS = ['CAR-SPLASH', 'KARR-seq', 'PARIS', 'RIC-seq', 'STAU1_hiCLIP', 'SPLASH']


def method_for(sample):
    if sample == 'paris_2016':
        return 'PARIS'
    if sample.startswith(('GM12878_cyt_', 'HEK293_cyt_', 'HeLa_cyt_', 'k562_chr_')):
        return 'SPLASH'
    for method in METHODS:
        if sample.startswith(method + '_'):
            return method
    raise ValueError(f'Unknown assay method: {sample}; extend METHODS explicitly')


def descriptive(values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    return dict(n_valid=len(values), mean=values.mean() if len(values) else np.nan,
                median=np.median(values) if len(values) else np.nan,
                q25=np.quantile(values, .25) if len(values) else np.nan,
                q75=np.quantile(values, .75) if len(values) else np.nan)


def save(fig, dest, name):
    fig.tight_layout(rect=(0, 0, 1, .94))
    for ext in ('pdf', 'png'):
        fig.savefig(dest / f'{name}.{ext}', dpi=160, bbox_inches='tight')
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--analysis-dir', required=True, help='Completed stratify_mfe.py results folder')
    args = parser.parse_args()
    root = Path(args.analysis_dir)
    manifest = json.loads((root / 'manifest.json').read_text())
    dest = root / 'method_comparison'
    if dest.exists() and any(dest.iterdir()):
        raise ValueError('method_comparison must be empty')
    dest.mkdir(exist_ok=True)
    pooled = {method: {} for method in METHODS}
    samples, bins, mapping = [], [], []
    for entry in manifest['selected']:
        sample = entry['sample']
        method = method_for(sample)
        mapping.append(dict(sample=sample, method=method))
        print(f'Comparing {sample}', flush=True)
        path = root / sample / 'annotated_pairs.tsv.gz'
        columns = pd.read_csv(path, sep='\t', nrows=0).columns
        metrics = [m for m in METRICS if m in columns]
        cols = metrics + ['loop_size_nt', 'loop_size_bin', 'included_in_plots']
        arrays = {m: [] for m in metrics + ['loop_size_nt']}
        counts = {label: 0 for label in LOOPS}
        n_rows = 0
        for chunk in pd.read_csv(path, sep='\t', usecols=cols, chunksize=100000):
            n_rows += len(chunk)
            for label, n in chunk.loop_size_bin.value_counts().items():
                counts[label] += int(n)
            # Loop distributions describe all coordinate-valid input pairs,
            # independently of observed MFE inclusion.
            arrays['loop_size_nt'].append(chunk.loop_size_nt.dropna().to_numpy(dtype=float))
            for metric in metrics:
                arrays[metric].append(chunk.loc[chunk.included_in_plots, metric].to_numpy(dtype=float))
        for metric, pieces in arrays.items():
            vals = np.concatenate(pieces)
            pooled[method].setdefault(metric, []).append(vals)
            samples.append(dict(sample=sample, method=method, metric=metric, **descriptive(vals)))
        for label, n in counts.items():
            bins.append(dict(sample=sample, method=method, loop_size_bin=label,
                             n_rows=n, fraction=n / n_rows if n_rows else np.nan))
    sample_stats = pd.DataFrame(samples)
    sample_stats.to_csv(dest / 'sample_statistics.tsv', sep='\t', index=False)
    pd.DataFrame(mapping).to_csv(dest / 'sample_method_mapping.tsv', sep='\t', index=False)
    bin_df = pd.DataFrame(bins)
    bin_df.to_csv(dest / 'sample_loop_bin_fractions.tsv', sep='\t', index=False)
    method_stats = []
    combined = {}
    for method, metrics in pooled.items():
        combined[method] = {}
        for metric, pieces in metrics.items():
            vals = np.concatenate(pieces)
            vals = vals[np.isfinite(vals)]
            combined[method][metric] = vals
            medians = sample_stats.loc[(sample_stats.method == method) & (sample_stats.metric == metric), 'median']
            method_stats.append(dict(method=method, metric=metric, n_samples=len(pieces),
                                     mean_sample_median=medians.mean(),
                                     median_sample_median=medians.median(), **descriptive(vals)))
    pd.DataFrame(method_stats).to_csv(dest / 'method_statistics.tsv', sep='\t', index=False)
    # Pair-weighted distributions and one point per sample median are separate.
    metrics = [m for m in METRICS if any(m in values for values in combined.values())]
    fig, axes = plt.subplots(len(metrics), 2, figsize=(16, 3.6 * len(metrics)), squeeze=False)
    for row, metric in enumerate(metrics):
        methods = [m for m in METHODS if len(combined[m].get(metric, []))]
        vals = [combined[m][metric] for m in methods]
        axes[row, 0].boxplot(vals, labels=[f'{m}\nn={len(v):,}' for m, v in zip(methods, vals)], showfliers=False)
        axes[row, 0].set_title('Pooled pairs; each pair has equal weight')
        axes[row, 0].set_ylabel(metric + (' (kcal/mol)' if metric != 'zscore_mfe' else ''))
        ax = axes[row, 1]
        for pos, method in enumerate(methods, 1):
            part = sample_stats.loc[(sample_stats.method == method) & (sample_stats.metric == metric)]
            jitter = np.linspace(-.12, .12, len(part)) if len(part) > 1 else [0]
            ax.scatter(pos + np.asarray(jitter), part['median'], s=35)
            for x, (_, record) in zip(pos + np.asarray(jitter), part.iterrows()):
                # Mapping table retains the full sample identity.
                ax.annotate(record['sample'].removeprefix(method + '_') if hasattr(str, 'removeprefix')
                            else record['sample'][len(method) + 1:], (x, record['median']),
                            xytext=(3, 3), textcoords='offset points', fontsize=6, rotation=15)
        ax.set_xticks(range(1, len(methods) + 1))
        ax.set_xticklabels(methods)
        ax.set_ylabel('Sample median ' + metric)
        ax.set_title('One point per samplesheet entry')
        for ax in axes[row]:
            ax.grid(axis='y', alpha=.2)
            ax.tick_params(axis='x', labelsize=8)
    fig.suptitle('MFE comparison across assay methods (observed MFE ≤ 0; finite)')
    save(fig, dest, 'mfe_by_method')

    fig, axes = plt.subplots(1, 3, figsize=(20, 6))
    for pos, method in enumerate(METHODS, 1):
        vals = combined[method].get('loop_size_nt', np.array([]))
        if not len(vals):
            continue
        x = np.sort(vals)
        # Log1p preserves zero-sized gaps and makes the full range visible.
        axes[0].plot(np.log10(1 + x), np.arange(1, len(x) + 1) / len(x), label=f'{method} (n={len(x):,})')
        part = sample_stats.loc[(sample_stats.method == method) & (sample_stats.metric == 'loop_size_nt')]
        jitter = np.linspace(-.12, .12, len(part)) if len(part) > 1 else [0]
        axes[1].scatter(pos + np.asarray(jitter), np.log10(1 + part['median']), s=35)
    axes[0].set_xlabel('log10(1 + genomic loop size in nt)')
    axes[0].set_ylabel('Cumulative fraction of pairs')
    axes[0].legend(fontsize=8)
    axes[0].set_title('Pooled genomic gap distributions')
    axes[1].set_xticks(range(1, len(METHODS) + 1))
    axes[1].set_xticklabels(METHODS, rotation=25)
    axes[1].set_ylabel('log10(1 + sample median loop size in nt)')
    axes[1].set_title('One point per samplesheet entry')
    grouped = bin_df.groupby(['method', 'loop_size_bin']).n_rows.sum().unstack(fill_value=0).reindex(index=METHODS, columns=LOOPS, fill_value=0)
    fractions = grouped.div(grouped.sum(axis=1), axis=0)
    fractions.to_csv(dest / 'method_loop_bin_fractions.tsv', sep='\t')
    bottom = np.zeros(len(METHODS))
    colors = plt.cm.tab20(np.linspace(0, 1, len(LOOPS)))
    for label, color in zip(LOOPS, colors):
        axes[2].bar(METHODS, fractions[label], bottom=bottom, label=label, color=color)
        bottom += fractions[label].to_numpy()
    axes[2].tick_params(axis='x', rotation=25)
    axes[2].set_ylabel('Fraction of all pairs')
    axes[2].set_title('Loop bins and ineligible pair categories')
    axes[2].legend(fontsize=7, loc='upper left', bbox_to_anchor=(1.02, 1))
    fig.suptitle('Loop-size comparison across assay methods; all input pairs')
    save(fig, dest, 'loop_sizes_by_method')
    (dest / 'comparison_notes.txt').write_text(
        'Methods are assays, mapped explicitly from sample prefixes. Pooled summaries weight pairs equally; '
        'sample-median panels weight samples equally. STAU1 high/low entries are subsets, not assumed independent replicates. '
        'HeLa RIC-seq versus HEK293/HEK293T methods is confounded by cell type and assay/arm-length differences. '
        'MFE panels require finite observed MFE <=0; loops use all coordinate-valid pairs. '
        'Genomic gaps are strand invariant, not spliced transcript distances. No significance tests.\n')
    print(f'Completed method comparison: {dest}', flush=True)


if __name__ == '__main__':
    main()
