#!/usr/bin/env python3
"""Render only the grey method overview and MFE-by-span figures, plus source tables."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

SPAN_BINS = ['0–49','50–99','100–499','500–999','1,000–4,999','5,000–9,999','10,000–49,999','≥50,000']
COLOR = '#444444'


def finite(values):
    values = np.asarray(values,dtype=float)
    return values[np.isfinite(values)]


def describe(values):
    values = finite(values)
    return dict(n=len(values),mean=values.mean() if len(values) else np.nan,
                median=np.median(values) if len(values) else np.nan,
                q25=np.quantile(values,.25) if len(values) else np.nan,
                q75=np.quantile(values,.75) if len(values) else np.nan)


def kde(values):
    values = finite(values)
    if len(values)<2 or values.std()==0: return None
    if len(values)>10000: values=np.random.default_rng(33).choice(values,10000,replace=False)
    bw=.25*values.std()
    if bw==0: return None
    grid=np.linspace(values.min()-3*bw,values.max()+3*bw,350); density=np.zeros(len(grid))
    for block in np.array_split(values,max(1,int(np.ceil(len(values)/500)))):
        density+=np.exp(-.5*((grid[:,None]-block)/bw)**2).sum(axis=1)
    return grid,density/(len(values)*bw*np.sqrt(2*np.pi))


def violin_box(ax, values, pos):
    values=finite(values)
    if not len(values):return
    sample=np.random.default_rng(33).choice(values,min(10000,len(values)),replace=False)
    if sample.std()>0:
        violin=ax.violinplot([sample],positions=[pos],vert=False,widths=.78,
                             showextrema=False,bw_method=.25,points=150)
        for body in violin['bodies']:
            body.set_facecolor('#bdbdbd');body.set_edgecolor('#777777');body.set_alpha(1)
    q1,med,q3=np.quantile(values,[.25,.5,.75])
    inside=values[(values>=q1-1.5*(q3-q1))&(values<=q3+1.5*(q3-q1))]
    ax.bxp([dict(q1=q1,med=med,q3=q3,whislo=inside.min(),whishi=inside.max(),fliers=[])],
           positions=[pos],vert=False,widths=.19,showfliers=False,patch_artist=True,
           boxprops=dict(facecolor='#eeeeee',edgecolor='#333333'),
           medianprops=dict(color='#222222',linewidth=1.5),
           whiskerprops=dict(color='#333333'),capprops=dict(color='#333333'))


def style(ax, labels, xlabel):
    ax.set_yticks(range(1,len(labels)+1));ax.set_yticklabels(labels,fontsize=10)
    ax.set_ylim(len(labels)+.55,.25);ax.set_xlabel(xlabel,fontsize=10)
    ax.grid(axis='x',color='#eeeeee');ax.set_axisbelow(True)
    ax.spines['top'].set_visible(False);ax.spines['right'].set_visible(False)


def legend(fig, shuffled):
    handles=[Line2D([],[],color=COLOR,label='Observed')]
    if shuffled:handles.append(Line2D([],[],color=COLOR,linestyle='--',label='Mean shuffled'))
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,.99),ncol=len(handles),frameon=False)


def draw_bands(axes, values, metric, nlabels):
    curves={key:kde(value) for key,value in values.items()}
    actual=[curve for curve in curves.values() if curve is not None]
    if actual:
        peak=max(d.max() for g,d in actual)
        lo=min(g.min() for g,d in actual);hi=max(g.max() for g,d in actual)
    else:
        nonempty=[finite(v) for v in values.values() if len(finite(v))]
        lo=min(v.min() for v in nonempty)-1 if nonempty else -1
        hi=max(v.max() for v in nonempty)+1 if nonempty else 1
        peak=1
    if metric=='mfe':lo=-70
    hi=max(hi,0)
    for (col,pos,source),curve in curves.items():
        ax=axes[col]
        if curve is not None:
            grid,density=curve
            ax.plot(grid,pos-.68*density/peak,color=COLOR,lw=1.3,linestyle='-' if source=='observed' else '--')
        else:
            value=finite(values[col,pos,source])
            if len(value):ax.vlines(value[0],pos-.5,pos,color=COLOR,linestyle='-' if source=='observed' else '--')
    for ax in axes:
        ax.set_xlim(lo,hi);ax.set_ylim(nlabels+.25,.15)
        for pos in range(1,nlabels+1):ax.axhline(pos,color='#dddddd',lw=.5,zorder=0)


def save(fig, outdir, stem):
    for ext in ['pdf','png']:fig.savefig(outdir/(stem+'.'+ext),dpi=150)
    plt.close(fig)


def render(data, outdir, stratify):
    methods=sorted(data.method.unique());categories=['pure_exon','intron_overlapping'] if stratify else ['all']
    metrics=[('percent_shortest_arm_paired','Shortest arm paired (%)'),
             ('paired_gc_percent','GC among paired nucleotides (%)'),
             ('mfe','MFE (kcal/mol)'),('mfe_per_bp','MFE / base-pair count (kcal/mol/bp)')]
    shuffled=bool(data.matched_eligible.any())
    energy=data[data.matched_eligible if shuffled else data.observed_eligible]
    summaries=[]
    fig,axes=plt.subplots(4,len(categories),figsize=(13 if stratify else 8,max(18,len(methods)*2.8)),
                           squeeze=False,gridspec_kw={'hspace':.40,'wspace':.25})
    for row,(metric,xlabel) in enumerate(metrics):
        bands={}
        for col,category in enumerate(categories):
            part=data[data.annotation_class==category] if stratify else data
            ep=energy[energy.annotation_class==category] if stratify else energy
            for pos,method in enumerate(methods,1):
                sub=(part if row<2 else ep);sub=sub[sub.method==method]
                sources=[('observed',metric)]
                if row>=2 and shuffled:sources.append(('shuffled','mean_shuffled_'+metric))
                for source,column in sources:
                    values=sub[column].to_numpy()
                    summaries.append(dict(figure='overview',method=method,group=category,span_bin='all',
                                          metric=metric,source=source,**describe(values)))
                    if row<2:violin_box(axes[row,col],values,pos)
                    else:bands[col,pos,source]=values
            style(axes[row,col],methods,xlabel)
            if row<2:axes[row,col].set_xlim(0,100)
            if row==0:axes[row,col].set_title(category.replace('_',' ').title(),fontsize=13,fontweight='bold',pad=12)
        if row>=2:draw_bands(axes[row],bands,metric,len(methods))
        axes[row,0].text(-.14,1.08,chr(65+row),transform=axes[row,0].transAxes,fontweight='bold',fontsize=12)
    legend(fig,shuffled);fig.subplots_adjust(top=.95,bottom=.06,left=.16,right=.98)
    save(fig,outdir,'method_overview')
    fig,axes=plt.subplots(2,len(methods),figsize=(max(8,4.4*len(methods)),13),squeeze=False,
                          gridspec_kw={'hspace':.28,'wspace':.13})
    for row,(metric,xlabel) in enumerate(metrics[2:]):
        bands={}
        for col,method in enumerate(methods):
            part=energy[energy.method==method]
            for pos,span in enumerate(SPAN_BINS,1):
                sub=part[part.span_bin==span]
                sources=[('observed',metric)]
                if shuffled:sources.append(('shuffled','mean_shuffled_'+metric))
                for source,column in sources:
                    values=sub[column].to_numpy();bands[col,pos,source]=values
                    summaries.append(dict(figure='span',method=method,group='all',span_bin=span,
                                          metric=metric,source=source,**describe(values)))
            style(axes[row,col],SPAN_BINS if col==0 else ['']*len(SPAN_BINS),xlabel)
            if row==0:axes[row,col].set_title(method,fontsize=13,fontweight='bold',pad=15)
            if col==0:axes[row,col].set_ylabel('Hybrid span (nt)',fontsize=12)
        draw_bands(axes[row],bands,metric,len(SPAN_BINS))
    legend(fig,shuffled);fig.subplots_adjust(top=.93,bottom=.07,left=.10,right=.99)
    save(fig,outdir,'method_span')
    pd.DataFrame(summaries).to_csv(outdir/'plot_group_counts.tsv',sep='\t',index=False)
    return shuffled


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--inputs',nargs='+',required=True);p.add_argument('--outdir',default='.')
    p.add_argument('--stratify',action='store_true');args=p.parse_args()
    outdir=Path(args.outdir);outdir.mkdir(parents=True,exist_ok=True)
    data=pd.concat([pd.read_csv(path,sep='\t',dtype={'sample':str,'method':str,'name':str,'row_id':str,
                   'lchr':str,'rchr':str}) for path in sorted(args.inputs)],ignore_index=True)
    if data.empty:raise ValueError('No valid hybrids remain after filtering')
    if data.duplicated(['sample','source_row_index']).any():raise ValueError('Duplicate sample/source-row identities')
    data['span_bin']=pd.cut(data.span_nt,[0,50,100,500,1000,5000,10000,50000,np.inf],labels=SPAN_BINS,right=False)
    if data.span_bin.isna().any():raise ValueError('Invalid span')
    data.to_csv(outdir/'plot_source.tsv.gz',sep='\t',index=False,compression='gzip')
    shuffled=render(data,outdir,args.stratify)
    counts=data.groupby(['sample','method','annotation_class','span_bin'],observed=True).agg(
        source_rows=('row_id','size'),gc_rows=('paired_gc_percent','count'),
        observed_mfe_rows=('observed_eligible','sum'),matched_mfe_rows=('matched_eligible','sum')).reset_index()
    counts.to_csv(outdir/'sample_group_counts.tsv',sep='\t',index=False)
    (outdir/'plot_notes.json').write_text(json.dumps(dict(
        stratified=args.stratify,source_rows=len(data),matched_controls=shuffled,
        span='max(lr,rr)-min(ll,rl)',normalization='MFE / observed base_pair_count, same denominator for shuffled mean',
        pairing='100 * base_pair_count / shortest genomic arm length',
        gc='100 * paired G+C / paired nucleotide count across both arms, excluding unpaired positions',
        raw_mfe_xmin=-70,density_sample_cap=10000,density_seed=33,density_bandwidth_sd=.25,
        composition_population='All retained valid structures; zero-pair GC undefined',
        mfe_population='Finite observed MFE <=0 and positive paired count; finite shuffled mean also required when controls exist',
        density_scaling='Common x limits and density heights across methods and groups within each metric'),indent=2)+'\n')

if __name__=='__main__':main()
