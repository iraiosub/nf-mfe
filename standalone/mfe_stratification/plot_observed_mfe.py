#!/usr/bin/env python3
"""Observed MFE violin/box plots by assay, exon class and genomic loop size."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from compare_methods import METHODS, method_for
from stratify_mfe import CLASSES, LOOPS

COLORS = dict(zip(METHODS, ['#3478A5', '#E69F00', '#279E73', '#CC6677', '#8B6BB1', '#56B4B9']))
CLASS_COLORS = {'exonic': '#3478A5', 'intronic': '#E69F00'}


def violin_box(ax, values, position, color, width=.72, alpha=.65):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if not len(values):
        return
    # KDE is reproducibly sampled for speed; box statistics use every pair.
    density_values = values
    if len(values) > 10000:
        density_values = np.random.default_rng(33).choice(values, 10000, replace=False)
    if len(np.unique(density_values)) > 1:
        violin = ax.violinplot([density_values], positions=[position], widths=width,
                              showextrema=False, points=150, bw_method=.25)
        for body in violin['bodies']:
            body.set_facecolor(color)
            body.set_edgecolor(color)
            body.set_alpha(alpha)
    q1, med, q3 = np.quantile(values, [.25, .5, .75])
    inside = values[(values >= q1 - 1.5 * (q3-q1)) & (values <= q3 + 1.5 * (q3-q1))]
    stats = dict(q1=q1, med=med, q3=q3, whislo=inside.min(), whishi=inside.max(), fliers=[])
    ax.bxp([stats], positions=[position], widths=width*.22, showfliers=False,
           patch_artist=True, manage_ticks=False,
           boxprops=dict(facecolor='white', edgecolor='#333333', linewidth=1),
           medianprops=dict(color='#222222', linewidth=1.5),
           whiskerprops=dict(color='#333333'), capprops=dict(color='#333333'))


def style(ax):
    ax.set_ylabel('Observed MFE (kcal/mol)')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', alpha=.16)
    ax.set_axisbelow(True)


def save(fig, dest, name, title, footer=None):
    fig.suptitle(title, fontsize=16, fontweight='bold')
    fig.text(.5, .01, footer or 'Observed MFE ≤ 0 • box: median/IQR, whiskers: 1.5×IQR • violin density: up to 10,000 pairs/group',
             ha='center', fontsize=9, color='#555555')
    fig.tight_layout(rect=(0, .035, 1, .95))
    for ext in ['png', 'pdf']:
        fig.savefig(dest / f'{name}.{ext}', dpi=180, bbox_inches='tight')
    plt.close(fig)


def facets(df, column, categories, dest, name, title):
    ncols = 2
    fig, axes = plt.subplots((len(categories)+1)//2, ncols,
                              figsize=(14, 4*((len(categories)+1)//2)), squeeze=False, sharey=True)
    for ax, category in zip(axes.flat, categories):
        part = df.loc[df[column] == category]
        labels = []
        for pos, method in enumerate(METHODS, 1):
            vals = part.loc[part.method == method, 'mfe'].to_numpy()
            violin_box(ax, vals, pos, COLORS[method])
            labels.append(f'{method.replace("STAU1_hiCLIP", "STAU1 hiCLIP")}\nn={len(vals):,}')
        ax.set_xticks(range(1, len(METHODS)+1))
        ax.set_xticklabels(labels, fontsize=8)
        ax.set_xlim(.4, len(METHODS)+.6)
        ax.set_title(category.replace('_', ' ') + (' nt' if column == 'loop_size_bin' and category in LOOPS[:8] else ''), fontweight='bold')
        style(ax)
    for ax in list(axes.flat)[len(categories):]:
        ax.set_visible(False)
    save(fig, dest, name, title)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--analysis-dir', required=True)
    args = p.parse_args()
    root = Path(args.analysis_dir)
    dest = root / 'observed_mfe_violins'
    dest.mkdir(exist_ok=True)
    manifest = json.loads((root / 'manifest.json').read_text())
    parts = []
    class_counts = []
    loop_parts = []
    for entry in manifest['selected']:
        sample = entry['sample']
        method = method_for(sample)
        print(f'Plotting {sample}', flush=True)
        frames = []
        for chunk in pd.read_csv(root / sample / 'annotated_pairs.tsv.gz', sep='\t',
                                 usecols=['mfe', 'annotation_class', 'loop_size_bin', 'loop_size_nt', 'included_in_plots'], chunksize=100000):
            valid_loops=chunk.loc[chunk.loop_size_nt.notna(), ['loop_size_nt']].copy()
            valid_loops['method']=method
            valid_loops['sample']=sample
            loop_parts.append(valid_loops)
            for category, count in chunk.annotation_class.value_counts().items():
                class_counts.append(dict(sample=sample, method=method, annotation_class=category, n=int(count)))
            chunk = chunk.loc[chunk.included_in_plots & np.isfinite(chunk.mfe)].copy()
            frames.append(chunk.drop(columns='included_in_plots'))
        df = pd.concat(frames, ignore_index=True)
        sample_dest = dest / 'samples' / sample
        sample_dest.mkdir(parents=True, exist_ok=True)
        for column, categories, title in [('annotation_class', CLASSES, 'Exon/intron stratification'),
                                          ('loop_size_bin', LOOPS, 'Loop-size stratification')]:
            fig, ax = plt.subplots(figsize=(max(9, len(categories)*1.1), 5))
            labels = []
            for pos, category in enumerate(categories, 1):
                vals = df.loc[df[column] == category, 'mfe'].to_numpy()
                violin_box(ax, vals, pos, COLORS[method])
                labels.append(f'{category.replace("_", " ")}\nn={len(vals):,}')
            ax.set_xticks(range(1, len(categories)+1))
            ax.set_xticklabels(labels, fontsize=8, rotation=15 if len(categories)>5 else 0)
            ax.set_xlabel('Genomic inter-arm gap (nt)' if column == 'loop_size_bin' else 'GENCODE v33 annotation class')
            style(ax)
            save(fig, sample_dest, 'mfe_by_'+column, sample + ': ' + title)
        df['method'] = method
        df['sample'] = sample
        parts.append(df)
    df = pd.concat(parts, ignore_index=True)
    for column in ['method', 'sample', 'annotation_class', 'loop_size_bin']:
        df[column] = df[column].astype('category')
    facets(df, 'annotation_class', CLASSES, dest, 'mfe_by_method_exon_class', 'Observed MFE by method — exon/intron stratification')
    facets(df, 'loop_size_bin', LOOPS[:8], dest, 'mfe_by_method_loop_size', 'Observed MFE by method — genomic loop-size stratification')
    facets(df, 'loop_size_bin', LOOPS[8:], dest, 'mfe_by_method_special_pairs', 'Observed MFE — overlapping, opposite-strand and trans pairs')
    fig, ax = plt.subplots(figsize=(10, 6))
    labels=[]
    for pos, method in enumerate(METHODS, 1):
        vals=df.loc[df.method==method, 'mfe'].to_numpy()
        violin_box(ax, vals, pos, COLORS[method])
        labels.append(f'{method}\nn={len(vals):,}')
    ax.set_xticks(range(1,len(METHODS)+1)); ax.set_xticklabels(labels)
    style(ax)
    save(fig,dest,'mfe_by_method','Observed MFE by assay method')
    # Joint stratification: compare exon/intron MFE within each loop bin and method.
    fig, axes = plt.subplots(3, 2, figsize=(16, 14), sharey=True)
    for ax, method in zip(axes.flat, METHODS):
        for pos, loop in enumerate(LOOPS[:8], 1):
            for category, offset in [('exonic', -.19), ('intronic', .19)]:
                vals=df.loc[(df.method==method)&(df.loop_size_bin==loop)&(df.annotation_class==category),'mfe'].to_numpy()
                violin_box(ax, vals, pos+offset, CLASS_COLORS[category], width=.34)
        ax.set_xticks(range(1,9)); ax.set_xticklabels(LOOPS[:8], fontsize=8, rotation=20)
        ax.set_xlabel('Genomic inter-arm gap (nt)'); ax.set_title(method, fontweight='bold'); style(ax)
    for ax in list(axes.flat)[len(METHODS):]:
        ax.set_visible(False)
    fig.legend(handles=[Patch(facecolor=c, label=k.title()) for k,c in CLASS_COLORS.items()],
               loc='lower right', bbox_to_anchor=(.85,.18), frameon=False)
    save(fig,dest,'mfe_by_method_loop_and_exon','Observed MFE — exon/intron comparison within genomic loop bins')
    summaries=[]
    for groups in [['method'], ['method','annotation_class'], ['method','loop_size_bin'],
                   ['method','loop_size_bin','annotation_class'], ['sample','method','annotation_class'], ['sample','method','loop_size_bin']]:
        for key, part in df.groupby(groups, observed=True):
            key=key if isinstance(key,tuple) else (key,)
            vals=part.mfe.to_numpy()
            q1,med,q3=np.quantile(vals,[.25,.5,.75])
            summaries.append(dict(zip(groups,key),stratification='+'.join(groups),n=len(vals),mean=vals.mean(),median=med,q25=q1,q75=q3))
    pd.DataFrame(summaries).to_csv(dest/'observed_mfe_stratified_summary.tsv',sep='\t',index=False)
    loop_df=pd.concat(loop_parts,ignore_index=True)
    fig,ax=plt.subplots(figsize=(11,6))
    loop_stats=[]
    for pos,method in enumerate(METHODS,1):
        values=loop_df.loc[loop_df.method==method,'loop_size_nt'].to_numpy()
        violin_box(ax,np.log10(1+values),pos,COLORS[method])
        loop_stats.append(dict(method=method,n=len(values),mean=values.mean(),median=np.median(values),q25=np.quantile(values,.25),q75=np.quantile(values,.75)))
    ax.set_xticks(range(1,len(METHODS)+1)); ax.set_xticklabels(METHODS)
    style(ax); ax.set_ylabel('log10(1 + genomic loop size in nt)')
    save(fig,dest,'loop_sizes_by_method','Genomic loop sizes by assay method',footer='All same-chromosome, same-strand nonoverlapping pairs • box within violin • zero gaps retained')
    pd.DataFrame(loop_stats).to_csv(dest/'method_loop_size_summary.tsv',sep='\t',index=False)
    counts = pd.DataFrame(class_counts)
    for group, filename in [('sample', 'exon_intron_proportions_by_sample'), ('method', 'exon_intron_proportions_by_method')]:
        table = counts.groupby([group, 'annotation_class']).n.sum().unstack(fill_value=0).reindex(columns=CLASSES, fill_value=0)
        fractions = table.div(table.sum(axis=1), axis=0)
        table.to_csv(dest/(filename+'_counts.tsv'),sep='\t')
        fractions.to_csv(dest/(filename+'_fractions.tsv'),sep='\t')
        fig, ax = plt.subplots(figsize=(13, max(5, len(table)*.38)))
        left=np.zeros(len(table))
        for category,color in zip(CLASSES,['#3478A5','#E69F00','#9B9B9B','#E4E4E4']):
            ax.barh(range(len(table)), fractions[category], left=left, label=category.replace('_',' ').title(), color=color)
            left += fractions[category].to_numpy()
        ax.set_yticks(range(len(table)))
        ax.set_yticklabels([f'{label}  (n={int(n):,})' for label,n in zip(table.index,table.sum(axis=1))],fontsize=8)
        ax.invert_yaxis(); ax.set_xlim(0,1)
        ax.set_xlabel('Fraction of all input pairs'); ax.legend(loc='upper center',bbox_to_anchor=(.5,-.08),ncol=4,frameon=False)
        ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
        save(fig,dest,filename,'GENCODE v33 exon/intron proportions by '+group,
             footer='All input pairs • both complete arms must be exon-contained for exonic classification')
    (dest/'plot_notes.txt').write_text('Observed MFE only; finite observed MFE <=0. Full data used for box statistics and summaries. '
        'Violin KDE uses fixed seed 33 and at most 10,000 pairs per group; bandwidth factor 0.25. '
        'Widths are equal across groups; counts shown in method facets and per-sample plots. '
        'Pooled method comparisons weight pairs equally. No shuffled or flipped values are plotted. '
        'Loop size is the strand-invariant genomic gap, not spliced transcript length.\n')
    print(f'Completed observed MFE plots: {dest}',flush=True)


if __name__ == '__main__':
    main()
