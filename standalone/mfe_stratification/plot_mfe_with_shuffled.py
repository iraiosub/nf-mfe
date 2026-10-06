#!/usr/bin/env python3
"""Side-by-side observed and per-pair mean shuffled MFE violins."""
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
from plot_observed_mfe import COLORS, violin_box, style, save

FOOTER = 'Left/darker: observed • right/lighter: mean shuffled • same finite pairs • box statistics use all pairs'


def paired(ax, part, position, color):
    violin_box(ax, part.mfe.to_numpy(), position-.19, color, width=.34, alpha=.8)
    violin_box(ax, part.mean_shuffled_mfe.to_numpy(), position+.19, color, width=.34, alpha=.3)


def finish(fig, dest, name, title):
    fig.legend(handles=[Patch(facecolor='#555555',alpha=.8,label='Observed MFE (left)'),
                        Patch(facecolor='#555555',alpha=.3,label='Mean shuffled MFE (right)')],
               loc='upper center',bbox_to_anchor=(.5,.95),ncol=2,frameon=False)
    save(fig,dest,name,title,footer=FOOTER)


def facets(df,column,categories,dest,name,title):
    fig,axes=plt.subplots((len(categories)+1)//2,2,figsize=(16,4.3*((len(categories)+1)//2)),squeeze=False,sharey=True)
    for ax,category in zip(axes.flat,categories):
        labels=[]
        part=df[df[column]==category]
        for pos,method in enumerate(METHODS,1):
            group=part[part.method==method]
            paired(ax,group,pos,COLORS[method])
            labels.append(f'{method.replace("STAU1_hiCLIP","STAU1 hiCLIP")}\nn={len(group):,}')
        ax.set_xticks(range(1,len(METHODS)+1)); ax.set_xticklabels(labels,fontsize=8)
        ax.set_xlim(.4,len(METHODS)+.6)
        ax.set_title(category.replace('_',' ')+(' nt' if column=='loop_size_bin' and category in LOOPS[:8] else ''),fontweight='bold')
        style(ax); ax.set_ylabel('MFE (kcal/mol)')
    for ax in list(axes.flat)[len(categories):]:
        ax.set_visible(False)
    finish(fig,dest,name,title)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--analysis-dir',required=True)
    args=p.parse_args()
    root=Path(args.analysis_dir)
    dest=root/'mfe_with_shuffled_violins'; dest.mkdir(exist_ok=True)
    manifest=json.loads((root/'manifest.json').read_text())
    parts=[]; qc=[]
    for entry in manifest['selected']:
        sample=entry['sample']; method=method_for(sample)
        print('Plotting paired MFE: '+sample,flush=True)
        frames=[]; eligible=missing=0
        path=root/sample/'annotated_pairs.tsv.gz'
        header=pd.read_csv(path,sep='\t',nrows=0).columns
        if 'mean_shuffled_mfe' not in header:
            raise ValueError(sample+': no precomputed mean_shuffled_mfe')
        for chunk in pd.read_csv(path,sep='\t',usecols=['mfe','mean_shuffled_mfe','annotation_class','loop_size_bin','included_in_plots'],chunksize=100000):
            chunk=chunk.loc[chunk.included_in_plots & np.isfinite(chunk.mfe)]
            eligible+=len(chunk)
            missing+=int((~np.isfinite(chunk.mean_shuffled_mfe)).sum())
            frames.append(chunk.loc[np.isfinite(chunk.mean_shuffled_mfe)].drop(columns='included_in_plots'))
        df=pd.concat(frames,ignore_index=True)
        qc.append(dict(sample=sample,method=method,observed_eligible=eligible,missing_shuffled_mean=missing,paired_n=len(df)))
        sample_dest=dest/'samples'/sample; sample_dest.mkdir(parents=True,exist_ok=True)
        for column,categories,title in [('annotation_class',CLASSES,'Exon/intron stratification'),('loop_size_bin',LOOPS,'Genomic loop-size stratification')]:
            fig,ax=plt.subplots(figsize=(max(10,len(categories)*1.15),5.5))
            labels=[]
            for pos,category in enumerate(categories,1):
                group=df[df[column]==category]
                paired(ax,group,pos,COLORS[method]); labels.append(f'{category.replace("_"," ")}\nn={len(group):,}')
            ax.set_xticks(range(1,len(categories)+1)); ax.set_xticklabels(labels,fontsize=8,rotation=15 if len(categories)>5 else 0)
            style(ax); ax.set_ylabel('MFE (kcal/mol)')
            finish(fig,sample_dest,'mfe_by_'+column,sample+': '+title)
        df['sample']=sample; df['method']=method; parts.append(df)
    df=pd.concat(parts,ignore_index=True)
    for column in ['method','sample','annotation_class','loop_size_bin']:
        df[column]=df[column].astype('category')
    facets(df,'annotation_class',CLASSES,dest,'mfe_by_method_exon_class','Observed and mean shuffled MFE — exon/intron stratification')
    facets(df,'loop_size_bin',LOOPS[:8],dest,'mfe_by_method_loop_size','Observed and mean shuffled MFE — genomic loop-size stratification')
    facets(df,'loop_size_bin',LOOPS[8:],dest,'mfe_by_method_special_pairs','Observed and mean shuffled MFE — special pair categories')
    fig,ax=plt.subplots(figsize=(12,6))
    labels=[]
    for pos,method in enumerate(METHODS,1):
        group=df[df.method==method]; paired(ax,group,pos,COLORS[method]); labels.append(f'{method}\nn={len(group):,}')
    ax.set_xticks(range(1,len(METHODS)+1)); ax.set_xticklabels(labels,fontsize=9)
    style(ax); ax.set_ylabel('MFE (kcal/mol)')
    finish(fig,dest,'mfe_by_method','Observed and mean shuffled MFE by assay method')
    fig,axes=plt.subplots(len(METHODS),2,figsize=(16,3.4*len(METHODS)),sharey=True)
    for row,method in enumerate(METHODS):
        for col,category in enumerate(['exonic','intronic']):
            ax=axes[row,col]
            part=df[(df.method==method)&(df.annotation_class==category)]
            for pos,loop in enumerate(LOOPS[:8],1):
                paired(ax,part[part.loop_size_bin==loop],pos,COLORS[method])
            ax.set_xticks(range(1,9)); ax.set_xticklabels(LOOPS[:8],fontsize=8,rotation=15)
            ax.set_title(method+' | '+category,fontweight='bold'); ax.set_xlabel('Genomic inter-arm gap (nt)')
            style(ax); ax.set_ylabel('MFE (kcal/mol)')
    finish(fig,dest,'mfe_by_method_loop_and_exon','Observed and mean shuffled MFE — loop bins within exon/intron classes')
    rows=[]
    for groups in [['method'],['method','annotation_class'],['method','loop_size_bin'],['method','loop_size_bin','annotation_class'],['sample','method','annotation_class'],['sample','method','loop_size_bin']]:
        for key,part in df.groupby(groups,observed=True):
            key=key if isinstance(key,tuple) else (key,)
            for metric in ['mfe','mean_shuffled_mfe']:
                values=part[metric].to_numpy(); q1,median,q3=np.quantile(values,[.25,.5,.75])
                rows.append(dict(zip(groups,key),stratification='+'.join(groups),metric=metric,n=len(values),mean=values.mean(),median=median,q25=q1,q75=q3))
    pd.DataFrame(rows).to_csv(dest/'mfe_stratified_summary.tsv',sep='\t',index=False)
    pd.DataFrame(qc).to_csv(dest/'paired_plot_qc.tsv',sep='\t',index=False)
    (dest/'plot_notes.txt').write_text('Left/darker observed; right/lighter per-pair mean_shuffled_mfe. Both use the method colour. '
        'Only pairs with finite observed MFE <=0 and finite mean_shuffled_mfe enter both distributions. '
        'Shuffled means are precomputed per-pair means, not individual shuffle values. '
        'Boxes/summaries use all pairs; KDE uses seed 33, up to 10,000 pairs per group and bandwidth 0.25. '
        'All 22 completed samples are included. No flipped-arm values or new folding.\n')
    print('Completed paired MFE plots: '+str(dest),flush=True)


if __name__=='__main__':
    main()
