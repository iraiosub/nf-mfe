#!/usr/bin/env python3
"""Prepare row-identified method plot features from a final nf-mfe TSV."""
import argparse
import gzip
import json
import re
from collections import defaultdict, Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import numpy as np
import pandas as pd
from method_plot_gc import paired_gc_percent


def chromosome(value):
    value = str(value)
    value = value[3:] if value.lower().startswith('chr') else value
    return 'MT' if value.upper() in ('M', 'MT') else value


def comma_list(value):
    return [item.strip() for item in value.split(',') if item.strip()]


class Intervals:
    def __init__(self, values):
        values = sorted(set(values))
        self.starts = np.array([s for s,e in values], dtype=np.int64)
        self.ends = np.maximum.accumulate([e for s,e in values])

    def contains(self, starts, ends):
        ix = np.searchsorted(self.starts, starts, side='right')-1
        out = np.zeros(len(starts), dtype=bool); valid = ix >= 0
        out[valid] = self.ends[ix[valid]] >= ends[valid]
        return out


def load_gtf(path, gene_types):
    raw = defaultdict(list); found = set(); selected = 0
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as handle:
        for line in handle:
            if line.startswith('#'): continue
            f = line.rstrip('\n').split('\t')
            if len(f) != 9 or f[2] not in ('gene', 'exon'): continue
            attrs = dict(re.findall(r'(\w+) "([^"]*)"', f[8]))
            biotype = attrs.get('gene_type', attrs.get('gene_biotype', ''))
            if f[2] == 'gene': found.add(biotype)
            if gene_types and biotype not in gene_types: continue
            if f[2] == 'gene': selected += 1
            raw[f[2], chromosome(f[0]), f[6]].append((int(f[3])-1,int(f[4])))
    unknown = set(gene_types)-found
    if unknown: raise ValueError('Gene types absent from GTF: '+', '.join(sorted(unknown)))
    if not selected: raise ValueError('GTF has no selected gene features')
    if not any(k[0] == 'exon' for k in raw): raise ValueError('GTF must contain selected exon features')
    return {k:Intervals(v) for k,v in raw.items()}


def annotate(data, index):
    if index is None:
        return np.ones(len(data), dtype=bool), np.full(len(data), 'all', dtype=object)
    same_gene = np.zeros(len(data), dtype=bool)
    exonic = {arm:np.zeros(len(data), dtype=bool) for arm in ['l','r']}
    for (c,strand), pos in data.groupby(['lchr','lstrand']).indices.items():
        idx = index.get(('gene',c,strand))
        if idx is not None:
            start = np.minimum(data.ll.to_numpy()[pos], data.rl.to_numpy()[pos])
            end = np.maximum(data.lr.to_numpy()[pos], data.rr.to_numpy()[pos])
            same_gene[pos] = idx.contains(start,end)
    for arm,start,end in [('l','ll','lr'),('r','rl','rr')]:
        for (c,strand),pos in data.groupby([arm+'chr',arm+'strand']).indices.items():
            idx = index.get(('exon',c,strand))
            if idx is not None:
                exonic[arm][pos] = idx.contains(data[start].to_numpy()[pos],data[end].to_numpy()[pos])
    return same_gene, np.where(exonic['l'] & exonic['r'], 'pure_exon', 'intron_overlapping')


def prepare(args):
    if args.gene_types and not args.gtf: raise ValueError('--gene-types requires --gtf')
    index = load_gtf(args.gtf, comma_list(args.gene_types)) if args.gtf else None
    excluded = {chromosome(c) for c in comma_list(args.exclude_chromosomes)}
    required = ['name','lchr','ll','lr','lstrand','rchr','rl','rr','rstrand',
                'lseq','rseq','dot_bracket','mfe']
    columns = pd.read_csv(args.input, sep='\t', nrows=0).columns.tolist()
    missing = set(required)-set(columns)
    if missing: raise ValueError('Missing required columns: '+', '.join(sorted(missing)))
    readcols = required + (['row_id'] if 'row_id' in columns else []) + (['mean_shuffled_mfe'] if 'mean_shuffled_mfe' in columns else [])
    qc = Counter(); statuses = Counter(); first = True
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with ProcessPoolExecutor(max_workers=args.processes) as pool:
        for data in pd.read_csv(args.input, sep='\t', usecols=readcols, chunksize=100000,
                                dtype={k:str for k in ['name','row_id','lchr','rchr','lseq','rseq','dot_bracket']}):
            data = data.reset_index(drop=True)
            data['source_row_index'] = np.arange(qc['input_rows'],qc['input_rows']+len(data))
            if 'row_id' not in data: data['row_id'] = data.source_row_index.astype(str)
            qc['input_rows'] += len(data)
            for arm,start,end in [('l','ll','lr'),('r','rl','rr')]:
                data[arm+'chr'] = data[arm+'chr'].map(chromosome)
                for col in [start,end]:
                    values = pd.to_numeric(data[col],errors='raise')
                    if values.isna().any() or (values != np.floor(values)).any():
                        raise ValueError('Noninteger coordinates in '+col)
                    data[col] = values.astype('int64')
                if ((data[start]<0)|(data[end]<=data[start])).any(): raise ValueError('Invalid arm interval')
                if not data[arm+'strand'].isin(['+','-']).all(): raise ValueError('Invalid strand')
            excluded_chr = data.lchr.isin(excluded) | data.rchr.isin(excluded)
            cis = (data.lchr == data.rchr) & (data.lstrand == data.rstrand)
            qc['excluded_chromosome_rows'] += int(excluded_chr.sum())
            qc['excluded_trans_or_opposite_strand_rows'] += int((~excluded_chr & ~cis).sum())
            data = data.loc[~excluded_chr & cis].reset_index(drop=True)
            gene_keep, classes = annotate(data, index)
            qc['excluded_gene_rows'] += int((~gene_keep).sum())
            data['annotation_class'] = classes
            data = data.loc[gene_keep].reset_index(drop=True)
            structures = data.dot_bracket.fillna('')
            opens = structures.str.count(r'\('); closes = structures.str.count(r'\)')
            valid = structures.str.fullmatch(r'[.(]*&[.)]*').fillna(False) & (opens == closes)
            qc['invalid_structure_rows'] += int((~valid).sum())
            data['base_pair_count'] = opens
            data = data.loc[valid].copy()
            data['paired_nucleotide_count'] = 2*data.base_pair_count
            data['shortest_arm_length_nt'] = np.minimum(data.lr-data.ll,data.rr-data.rl)
            if (data.base_pair_count > data.shortest_arm_length_nt).any(): raise ValueError('Paired count exceeds shortest arm')
            data['span_nt'] = np.maximum(data.lr,data.rr)-np.minimum(data.ll,data.rl)
            data['percent_shortest_arm_paired'] = 100*data.base_pair_count/data.shortest_arm_length_nt
            data['mfe'] = pd.to_numeric(data.mfe,errors='raise')
            if 'mean_shuffled_mfe' not in data: data['mean_shuffled_mfe'] = np.nan
            data['mean_shuffled_mfe'] = pd.to_numeric(data.mean_shuffled_mfe,errors='raise')
            tasks = list(data[['lseq','rseq','dot_bracket','mfe']].itertuples(index=False,name=None))
            results = list(pool.map(paired_gc_percent,tasks,chunksize=250))
            data['paired_gc_percent'] = [r[0] for r in results]
            data['gc_mapping_status'] = [r[1] for r in results]
            statuses.update(data.gc_mapping_status)
            if not np.array_equal([r[2] for r in results],data.paired_nucleotide_count): raise ValueError('Paired counts differ')
            denominator = data.base_pair_count.where(data.base_pair_count > 0)
            for source in ['mfe','mean_shuffled_mfe']:
                data[source+'_per_bp'] = data[source]/denominator
            data['observed_eligible'] = np.isfinite(data.mfe) & (data.mfe<=0) & (data.base_pair_count>0)
            data['matched_eligible'] = data.observed_eligible & np.isfinite(data.mean_shuffled_mfe)
            data['sample'] = args.sample; data['method'] = args.method
            qc['source_rows'] += len(data); qc['observed_mfe_rows'] += int(data.observed_eligible.sum())
            qc['matched_mfe_rows'] += int(data.matched_eligible.sum())
            qc['gc_rows'] += int(np.isfinite(data.paired_gc_percent).sum())
            data.drop(columns=['lseq','rseq']).to_csv(args.output,sep='\t',index=False,
                compression='gzip', mode='w' if first else 'a',header=first)
            first = False
    (Path(args.output).parent/(args.sample+'.plot_qc.json')).write_text(json.dumps(dict(
        sample=args.sample,method=args.method,gtf=args.gtf,gene_types=comma_list(args.gene_types),
        exclude_chromosomes=sorted(excluded),counts=dict(qc),gc_mapping=dict(statuses)),indent=2)+'\n')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ['input','output','sample','method']: p.add_argument('--'+key,required=True)
    p.add_argument('--gtf'); p.add_argument('--gene-types',default='')
    p.add_argument('--exclude-chromosomes',default=''); p.add_argument('--processes',type=int,default=1)
    args = p.parse_args()
    if args.processes<1: p.error('--processes must be positive')
    prepare(args)

if __name__ == '__main__': main()
