#!/usr/bin/env python3
"""Replot pipeline source data after pooling replicates within method/condition."""
import argparse
from collections import Counter
import json
from pathlib import Path
import re
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'bin'))
from plot_method_figures import SPAN_BINS, describe, draw_bands, legend, save, style, violin_box

METRICS=[('percent_shortest_arm_paired','Shortest arm paired (%)'),
         ('paired_gc_percent','GC among paired nucleotides (%)'),
         ('mfe','MFE (kcal/mol)'),('mfe_per_bp','MFE / base-pair count (kcal/mol/bp)')]
ORDER=['CAR-SPLASH | HEK293','KARR-seq | HEK293','PARIS | HEK293T','RIC-seq | HeLa',
       'SPLASH | GM12878 | cytoplasm','SPLASH | HEK293 | cytoplasm',
       'SPLASH | HeLa | cytoplasm','SPLASH | K562 | chromatin',
       'STAU1 hiCLIP | HEK293 | high','STAU1 hiCLIP | HEK293 | low']


def sample_metadata(sample):
    match=re.fullmatch(r'SPLASH_(GM12878|HEK293|HeLa|k562)_(cyt|chr)_rep(\d+)',sample)
    if match:
        cell,compartment,replicate=match.groups()
        cell='K562' if cell=='k562' else cell
        compartment={'cyt':'cytoplasm','chr':'chromatin'}[compartment]
        method='SPLASH';condition=cell+' | '+compartment
    elif sample.startswith('STAU1_hiCLIP_HEK293_'):
        match=re.fullmatch(r'STAU1_hiCLIP_HEK293_(high|low)_(.+)',sample)
        if not match:raise ValueError('Unknown STAU1 condition: '+sample)
        method='STAU1 hiCLIP';condition='HEK293 | '+match.group(1);replicate=''
    else:
        match=re.fullmatch(r'(CAR-SPLASH|KARR-seq|PARIS|RIC-seq)_([^_]+)_(\d+)_(.+)',sample)
        if not match:raise ValueError('Unknown sample grouping: '+sample)
        method,condition,replicate,accession=match.groups()
    cohort=method+' | '+condition
    if cohort not in ORDER:raise ValueError('Unexpected cohort: '+cohort)
    return dict(sample=sample,method=method,condition=condition,cohort=cohort,replicate=replicate)


def display(cohort):
    parts=cohort.split(' | ')
    return parts[0]+'\n'+' · '.join(parts[1:])


def render(data,groups,dest,categories):
    rows=[];shuffled=bool(data.matched_eligible.any())
    # Materialize each cohort once, rather than repeatedly copying the full source.
    cohorts={label:part for label,part in data.groupby('cohort',observed=True,sort=False)}
    fig,axes=plt.subplots(4,len(categories),figsize=(16 if len(categories)>1 else 12,max(26,3.2*len(groups))),
                          squeeze=False,gridspec_kw={'hspace':.32,'wspace':.25})
    for row,(metric,xlabel) in enumerate(METRICS):
        bands={}
        for col,category in enumerate(categories):
            ax=axes[row,col]
            for pos,label in enumerate(groups,1):
                part=cohorts[label]
                if category!='all':part=part[part.annotation_class==category]
                if row>=2:part=part[part.matched_eligible if shuffled else part.observed_eligible]
                sources=[('observed',metric)]
                if row>=2 and shuffled:sources.append(('shuffled','mean_shuffled_'+metric))
                for source,column in sources:
                    values=part[column].to_numpy()
                    rows.append(dict(figure='overview',cohort=label,annotation_class=category,
                                     span_bin='all',metric=metric,source=source,**describe(values)))
                    if row<2:violin_box(ax,values,pos)
                    else:bands[col,pos,source]=values
            style(ax,[display(g) for g in groups],xlabel)
            if row<2:ax.set_xlim(0,100)
            if row==0 and category!='all':ax.set_title(category.replace('_',' ').title(),fontsize=14,fontweight='bold')
        if row>=2:draw_bands(axes[row],bands,metric,len(groups))
    legend(fig,shuffled)
    fig.subplots_adjust(top=.97,bottom=.04,left=.25 if len(categories)==1 else .18,right=.98)
    save(fig,dest,'condition_overview')
    fig,axes=plt.subplots(2,len(groups),figsize=(4.3*len(groups),20),squeeze=False,
                          gridspec_kw={'hspace':.22,'wspace':.12})
    for row,(metric,xlabel) in enumerate(METRICS[2:]):
        bands={}
        for col,label in enumerate(groups):
            part=cohorts[label]
            part=part[part.matched_eligible if shuffled else part.observed_eligible]
            for pos,span in enumerate(SPAN_BINS,1):
                band=part[part.span_bin==span]
                sources=[('observed',metric)]
                if shuffled:sources.append(('shuffled','mean_shuffled_'+metric))
                for source,column in sources:
                    values=band[column].to_numpy();bands[col,pos,source]=values
                    rows.append(dict(figure='span',cohort=label,annotation_class='all',span_bin=span,
                                     metric=metric,source=source,**describe(values)))
            style(axes[row,col],SPAN_BINS if col==0 else ['']*len(SPAN_BINS),xlabel)
            if row==0:axes[row,col].set_title(display(label),fontsize=12,fontweight='bold',pad=18)
            if col==0:axes[row,col].set_ylabel('Hybrid span (nt)',fontsize=12)
        draw_bands(axes[row],bands,metric,len(SPAN_BINS))
    legend(fig,shuffled);fig.subplots_adjust(top=.93,bottom=.05,left=.04,right=.995)
    save(fig,dest,'condition_span')
    pd.DataFrame(rows).to_csv(dest/'plot_group_counts.tsv',sep='\t',index=False)
    return rows


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--pipeline-plots',required=True);p.add_argument('--outdir',required=True)
    args=p.parse_args();root=Path(args.pipeline_plots);dest=Path(args.outdir);dest.mkdir(parents=True,exist_ok=True)
    source=root/'plot_source.tsv.gz';before=source.stat()
    notes=json.loads((root/'plot_notes.json').read_text())
    expected=pd.read_csv(root/'sample_group_counts.tsv',sep='\t')
    metadata=pd.DataFrame([sample_metadata(s) for s in sorted(expected['sample'].unique())])
    metadata.to_csv(dest/'sample_condition_mapping.tsv',sep='\t',index=False)
    mapping=metadata.set_index('sample').cohort.to_dict()
    groups=[label for label in ORDER if label in metadata.cohort.unique()]
    numeric=['percent_shortest_arm_paired','paired_gc_percent','mfe','mean_shuffled_mfe','mfe_per_bp','mean_shuffled_mfe_per_bp']
    columns=['sample','annotation_class','span_bin','observed_eligible','matched_eligible']+numeric
    frames=[];counts=Counter();observed=Counter();matched=Counter();gc=Counter()
    for i,chunk in enumerate(pd.read_csv(source,sep='\t',usecols=columns,chunksize=500000,
                               dtype={**{k:float for k in numeric},'sample':str,'observed_eligible':bool,'matched_eligible':bool}),1):
        unknown=set(chunk['sample'].unique())-set(mapping)
        if unknown:raise ValueError('Unmapped source samples: '+str(unknown))
        counts.update(chunk['sample'].value_counts().to_dict())
        observed.update(chunk.loc[chunk.observed_eligible,'sample'].value_counts().to_dict())
        matched.update(chunk.loc[chunk.matched_eligible,'sample'].value_counts().to_dict())
        gc.update(chunk.loc[np.isfinite(chunk.paired_gc_percent),'sample'].value_counts().to_dict())
        chunk['cohort']=pd.Categorical(chunk['sample'].map(mapping),categories=groups,ordered=True)
        chunk['annotation_class']=pd.Categorical(chunk.annotation_class,categories=['all','pure_exon','intron_overlapping'])
        chunk['span_bin']=pd.Categorical(chunk.span_bin,categories=SPAN_BINS,ordered=True)
        if chunk.annotation_class.isna().any() or chunk.span_bin.isna().any():raise ValueError('Invalid source facet')
        frames.append(chunk.drop(columns='sample'))
        if i%10==0:print('Loaded '+str(sum(counts.values()))+' pipeline source rows',flush=True)
    after=source.stat()
    if (before.st_size,before.st_mtime_ns)!=(after.st_size,after.st_mtime_ns):raise ValueError('Source changed while reading')
    assert sum(counts.values())==notes['source_rows']
    expected_totals=expected.groupby('sample')[['source_rows','gc_rows','observed_mfe_rows','matched_mfe_rows']].sum()
    for sample,row in expected_totals.iterrows():
        assert counts[sample]==row.source_rows and gc[sample]==row.gc_rows
        assert observed[sample]==row.observed_mfe_rows and matched[sample]==row.matched_mfe_rows
    data=pd.concat(frames,ignore_index=True);frames.clear()
    print('Rendering '+str(len(groups))+' pooled method/condition groups',flush=True)
    categories=['pure_exon','intron_overlapping'] if notes['stratified'] else ['all']
    rows=render(data,groups,dest,categories)
    sample_counts=metadata.copy()
    for column,counter in [('source_rows',counts),('gc_rows',gc),('observed_mfe_rows',observed),('matched_mfe_rows',matched)]:
        sample_counts[column]=sample_counts['sample'].map(counter).fillna(0).astype('int64')
    sample_counts.to_csv(dest/'sample_counts.tsv',sep='\t',index=False)
    totals=sample_counts.groupby(['cohort','method','condition'],sort=False).agg(
        n_samples=('sample','size'),source_rows=('source_rows','sum'),gc_rows=('gc_rows','sum'),
        observed_mfe_rows=('observed_mfe_rows','sum'),matched_mfe_rows=('matched_mfe_rows','sum')).reset_index()
    totals.to_csv(dest/'condition_counts.tsv',sep='\t',index=False)
    summary=pd.DataFrame(rows)
    for (cohort,metric,who),part in summary[summary.figure=='span'].groupby(['cohort','metric','source']):
        assert part.n.sum()==int(totals.loc[totals.cohort==cohort,'matched_mfe_rows' if notes['matched_controls'] else 'observed_mfe_rows'].iloc[0])
    manifest=dict(pipeline_source=str(source.resolve()),source_bytes=before.st_size,source_mtime_ns=before.st_mtime_ns,
                  samples=len(metadata),cohorts=len(groups),source_rows=len(data),matched_mfe_rows=int(data.matched_eligible.sum()),
                  grouping='Pool hybrid rows across replicates within method and condition; STAU1 high/low remain separate',
                  source_join='Join pipeline plot_source.tsv.gz to sample_condition_mapping.tsv on sample',
                  pipeline_plot_notes=notes,overview_inches=[16 if len(categories)>1 else 12,max(26,3.2*len(groups))],
                  span_inches=[4.3*len(groups),20],validation='All per-sample source/GC/observed/matched counts equal pipeline exports; span counts reconcile per cohort')
    (dest/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('Completed: '+str(dest),flush=True)

if __name__=='__main__':main()
