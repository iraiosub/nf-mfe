#!/usr/bin/env python3
"""Pool annotation counts across replicates of the same assay/cell/compartment."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from stratify_mfe import CLASSES


def group_for(sample):
    if sample.startswith('CAR-SPLASH_HEK293_'):
        return 'CAR-SPLASH | HEK293'
    if sample.startswith('KARR-seq_HEK293_'):
        return 'KARR-seq | HEK293'
    if sample.startswith('PARIS_HEK293T_'):
        return 'PARIS | HEK293T'
    if sample.startswith('RIC-seq_HeLa_'):
        return 'RIC-seq | HeLa'
    if sample.startswith('STAU1_hiCLIP_HEK293_high_'):
        return 'STAU1 hiCLIP | HEK293 | high'
    if sample.startswith('STAU1_hiCLIP_HEK293_low_'):
        return 'STAU1 hiCLIP | HEK293 | low'
    if sample == 'paris_2016':
        return 'PARIS | 2016 dataset'
    for cell in ['GM12878', 'HEK293', 'HeLa']:
        if sample.startswith(cell+'_cyt_rep'):
            return 'SPLASH | '+cell+' | cytoplasm'
    if sample.startswith('k562_chr_rep'):
        return 'SPLASH | K562 | chromatin'
    raise ValueError('Unmapped sample: '+sample)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plot-dir',required=True)
    args=p.parse_args()
    dest=Path(args.plot_dir)
    counts=pd.read_csv(dest/'exon_intron_proportions_by_sample_counts.tsv',sep='\t',index_col=0).reindex(columns=CLASSES,fill_value=0)
    mapping=pd.DataFrame({'sample':counts.index,'pooled_group':[group_for(s) for s in counts.index]})
    mapping.to_csv(dest/'pooled_replicate_group_mapping.tsv',sep='\t',index=False)
    counts['pooled_group']=mapping.pooled_group.to_numpy()
    pooled=counts.groupby('pooled_group',sort=False)[CLASSES].sum()
    fractions=pooled.div(pooled.sum(axis=1),axis=0)
    sizes=mapping.groupby('pooled_group').size()
    stem='exon_intron_proportions_pooled_replicates'
    pooled.to_csv(dest/(stem+'_counts.tsv'),sep='\t')
    fractions.to_csv(dest/(stem+'_fractions.tsv'),sep='\t')
    fig,ax=plt.subplots(figsize=(13,7))
    left=np.zeros(len(pooled))
    for category,color in zip(CLASSES,['#3478A5','#E69F00','#9B9B9B','#E4E4E4']):
        ax.barh(range(len(pooled)),fractions[category],left=left,label=category.replace('_',' ').title(),color=color)
        left+=fractions[category].to_numpy()
    ax.set_yticks(range(len(pooled)))
    ax.set_yticklabels([f'{label}  (samples={sizes[label]}, n={int(n):,})' for label,n in zip(pooled.index,pooled.sum(axis=1))],fontsize=9)
    ax.invert_yaxis(); ax.set_xlim(0,1)
    ax.set_xlabel('Fraction of all input pairs (replicate counts pooled)')
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    ax.legend(loc='upper center',bbox_to_anchor=(.5,-.12),ncol=4,frameon=False)
    fig.suptitle('Exon/intron proportions — replicates pooled within assay and sample type',fontsize=15,fontweight='bold')
    fig.text(.5,.015,'Cell types, cytoplasm/chromatin, STAU1 subsets and PARIS datasets kept separate',ha='center',fontsize=9,color='#555555')
    fig.tight_layout(rect=(0,.065,1,.94))
    for ext in ['pdf','png']:
        fig.savefig(dest/(stem+'.'+ext),dpi=180,bbox_inches='tight')
    plt.close(fig)
    assert len(mapping)==22 and len(pooled)==11
    assert pooled.to_numpy().sum()==counts[CLASSES].to_numpy().sum()
    assert np.allclose(fractions.sum(axis=1),1)
    print('Verified 22 samples pooled into 11 groups; counts conserved; proportions sum to 1.')


if __name__=='__main__':
    main()
