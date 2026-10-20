#!/usr/bin/env python3
"""Plot the multi-species clustered duplex TSV on the assumed common human rDNA.

Requires numpy and matplotlib. Inputs use 1-based inclusive coordinates.
Coverage is weighted by `count` by default, or one per cluster with --weight rows.
RPM is deliberately not used: this rRNA-only table is not a complete library.
"""
import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
import numpy as np

TRACKS = [('intra', 'Within RNA', '#7752a3'),
          ('other_rrna', 'Other rRNA', '#138c81'),
          ('exon', 'Exon partner', '#ca791b'),
          ('intron', 'Intron partner', '#3579b1'),
          ('other', 'Other partner', '#888888'),
          ('clippy_direct', 'rRNA arm with\nClippy annotation', '#b33c86'),
          ('clippy_partner', 'rRNA footprint:\npartner has Clippy', '#986832')]
GROUPS = ['human', 'mouse_lifted', 'rat_lifted', 'merged', 'merged_both_arms_available']
ALIASES = {'rDNA', 'rRNA_45S'}


def arm(row, side):
    """Never substitute rodent source coordinates for failed mappings."""
    original = row['species'] == 'human'
    mapped = row[side+'_mapped'] == 'TRUE'
    prefix = side + ('_genomic_' if original else '_LO_target_')
    available = original or mapped
    return dict(seq=row[prefix+'seqnames'] if available else '',
                start=int(row[prefix+'start'])-1 if available else None,
                end=int(row[prefix+'end']) if available else None,
                region=row[side+'_region'] if original or not mapped else row[side+'_LO_target_region'],
                original_region=row[side+'_region'],
                available=available, peak=row[side+'_clippy_has_peak']=='TRUE')


def owner(a, bounds):
    if a['seq'] not in ALIASES:
        return None
    return next((r for r, (lo, hi) in bounds.items()
                 if lo <= a['start'] < a['end'] <= hi), None)


def add_interval(delta, intervals, weight):
    """Union the two arms so overlapping same-RNA arms count once per base."""
    merged = []
    for lo, hi in sorted(intervals):
        if merged and lo <= merged[-1][1]:
            merged[-1][1] = max(hi, merged[-1][1])
        else:
            merged.append([lo, hi])
    for lo, hi in merged:
        delta[lo] += weight
        delta[hi] -= weight


def binned(v, size):
    return np.array([v[i:i+size].mean() for i in range(0, len(v), size)])


def collect(path, bounds, size, weighting):
    data = {g: {r: dict(delta={k: np.zeros(hi-lo+1) for k, _, _ in TRACKS},
                       contacts=np.zeros(((hi-lo+size-1)//size,)*2),
                       stats=Counter(), regions=Counter())
                for r, (lo, hi) in bounds.items()} for g in GROUPS}
    audit = Counter()
    peaks = {}
    with path.open() as f:
        for row in csv.DictReader(f, delimiter='\t'):
            sp = row['species']
            if sp not in {'human', 'mouse', 'rat'}:
                raise ValueError(f'Unexpected species: {sp}')
            audit[sp+'_rows'] += 1
            w = int(row['count']) if weighting == 'count' else 1
            if w < 0:
                raise ValueError('Negative count')
            arms = [arm(row, s) for s in ['L', 'R']]
            owners = [owner(a, bounds) for a in arms]
            groups = ['human' if sp == 'human' else sp+'_lifted', 'merged']
            if all(a['available'] for a in arms):
                groups += ['merged_both_arms_available']
                audit[sp+'_both_arms_available_rows'] += 1
            for s, a in zip(['L', 'R'], arms):
                if row[s+'_region']=='rRNA':
                    audit[sp+'_rrna_arms'] += 1
                    if not a['available']:
                        audit[sp+'_rrna_arms_without_mapping'] += 1
                if a['peak']:
                    key = (sp, row[s+'_clippy_peak_ids'])
                    peaks[key] = dict(species=sp, peak_id=key[1], region=row[s+'_clippy_peak_regions'],
                                      cDNA_score=row[s+'_clippy_peak_cDNA_scores'])
                if a['seq'] in ALIASES and owners[['L','R'].index(s)] is None:
                    audit[sp+'_rdna_arms_outside_or_crossing_mature_bounds'] += 1
            for rna, (lo, hi) in bounds.items():
                if rna not in owners:
                    continue
                indices = [i for i, o in enumerate(owners) if o == rna]
                intervals = [(arms[i]['start']-lo, arms[i]['end']-lo) for i in indices]
                region = None
                if owners == [rna, rna]:
                    kind = 'intra'
                else:
                    partner = arms[1-indices[0]]
                    # Original region supports classification even if partner failed LO.
                    region = partner['region'] or partner['original_region'] or 'unknown'
                    if partner['seq'] in ALIASES:
                        region = 'rRNA'
                    kind = ('other_rrna' if region == 'rRNA' else 'exon' if region in
                            {'CDS', 'UTR5', 'UTR3', 'exon'} else 'intron' if region == 'intron' else 'other')
                for g in groups:
                    d = data[g][rna]
                    d['stats']['rows'] += 1
                    d['stats'][kind+'_rows'] += 1
                    d['stats'][kind+'_weight'] += w
                    add_interval(d['delta'][kind], intervals, w)
                    if region is not None:
                        d['regions'][region] += w
                        if not arms[1-indices[0]]['available']:
                            d['stats']['partner_unmapped_rows'] += 1
                    if kind == 'intra':
                        i, j = sorted((min(d['contacts'].shape[0]-1,
                                         int(((a['start']+a['end']-1)/2-lo)//size)) for a in arms))
                        d['contacts'][i,j] += w
                    direct = [(arms[i]['start']-lo, arms[i]['end']-lo)
                              for i in indices if arms[i]['peak']]
                    linked = [(arms[i]['start']-lo, arms[i]['end']-lo)
                              for i in indices if arms[1-i]['peak']]
                    for key, iv in [('clippy_direct',direct), ('clippy_partner',linked)]:
                        if iv:
                            add_interval(d['delta'][key], iv, w)
                            d['stats'][key+'_rows'] += 1
                            d['stats'][key+'_weight'] += w
    for group in data.values():
        for d in group.values():
            d['coverage'] = {}
            for k, delta in d.pop('delta').items():
                assert delta.sum() == 0
                d['coverage'][k] = np.cumsum(delta)[:-1]
                assert np.all(d['coverage'][k] >= 0)
            assert d['contacts'].sum() == d['stats']['intra_weight']
    for rna in bounds:
        for k, _, _ in TRACKS:
            assert np.array_equal(data['merged'][rna]['coverage'][k], sum(
                data[g][rna]['coverage'][k] for g in GROUPS[:3]))
        assert np.array_equal(data['merged'][rna]['contacts'], sum(data[g][rna]['contacts'] for g in GROUPS[:3]))
    return data, dict(audit), list(peaks.values())


def expansion(ax, segments, rna, length):
    selected = [s for s in segments if s['rna']==rna]
    centres = np.array([(int(s['start_nt'])+int(s['end_nt']))/2 for s in selected])
    labels = centres.copy()
    for i in range(1,len(labels)):
        labels[i] = max(labels[i],labels[i-1]+length*.04)
    if len(labels):
        labels[-1] = min(labels[-1],length*.98)
        for i in range(len(labels)-2,-1,-1):
            labels[i] = min(labels[i],labels[i+1]-length*.04)
    for s, c, y in zip(selected,centres,labels):
        ax.axhspan(int(s['start_nt']),int(s['end_nt']),xmin=.72,xmax=.98,color='#8585a4')
        ax.plot([.55,.72],[y,c],color='#8585a4',lw=.5)
        ax.text(.52,y,s['segment'],ha='right',va='center',fontsize=6,color='#555578')
    ax.set(xlim=(0,1),xticks=[],title='Expansion\nsegments')
    ax.tick_params(left=False,labelleft=False)
    for spine in ax.spines.values(): spine.set_visible(False)


def plot_group(group, data, assets, segments, out, size, weighting, arc_min):
    fig = plt.figure(figsize=(29,10))
    grid = fig.add_gridspec(2,12,width_ratios=[.9,3,3,1.1,.75,.75,.75,.75,.75,.75,.75,2.6],
                          left=.025,right=.985,bottom=.13,top=.80,hspace=.47,wspace=.65)
    maximum = max(d['contacts'].max() for d in data.values())
    norm = Normalize(0,max(1,np.log2(1+maximum)))
    row_axes = []
    for ri,(rna,d) in enumerate(data.items()):
        a = assets[rna]; length = a['length']
        edges = np.r_[np.arange(0,length,size)+.5,length+.5]
        centres = (edges[:-1]+edges[1:])/2
        axes = [fig.add_subplot(grid[ri,j]) for j in range(12)]
        row_axes.append(axes)
        for ax in axes[:11]: ax.set_ylim(length+.5,.5)
        expansion(axes[0],segments,rna,length)
        physical = axes[1].imshow(a['distance'],extent=(.5,length+.5,length+.5,.5),
                                  cmap='coolwarm_r',vmin=0,vmax=120,aspect='equal',interpolation='nearest')
        axes[1].set(title=f'{rna} · cryo-EM · 4V6X/{a["chain"]}\nPhysical distance',xlabel='RNA position (nt)',ylabel='RNA position (nt)')
        matrix = d['contacts'] + d['contacts'].T - np.diag(np.diag(d['contacts']))
        contact = axes[2].pcolormesh(edges,edges,np.log2(1+matrix),cmap='Reds',norm=norm,rasterized=True)
        axes[1].set_aspect('equal', adjustable='box')
        axes[2].set_aspect('equal')
        axes[2].set_aspect('equal', adjustable='box')
        axes[2].set(title=f'{rna} · cluster contacts\n{d["stats"]["intra_rows"]:,} rows',xlabel='RNA position (nt)')
        if not matrix.any(): axes[2].text(.5,.5,'No within-RNA records',transform=axes[2].transAxes,ha='center')
        i,j = np.nonzero(np.triu(d['contacts'],1) >= arc_min)
        theta = np.linspace(-np.pi/2,np.pi/2,50)
        radius = (centres[j]-centres[i])/2
        if len(i):
            seg = np.stack([radius[:,None]*np.cos(theta),
                            (centres[i]+centres[j])[:,None]/2+radius[:,None]*np.sin(theta)],axis=2)
            lines = LineCollection(seg,cmap='Reds',norm=norm,lw=.6)
            lines.set_array(np.log2(1+d['contacts'][i,j])); axes[3].add_collection(lines)
        axes[3].set(xlim=(0,length/2),xticks=[],title=f'Arcs ≥{arc_min:g}\nweight/bin pair')
        for ax in axes[2:11]: ax.tick_params(labelleft=False)
        for ci,(kind,label,color) in enumerate(TRACKS,4):
            ax = axes[ci]; values = binned(d['coverage'][kind],size)
            ax.fill_betweenx(centres,0,values,color=color,alpha=.7)
            ax.plot(values,centres,color=color,lw=.6)
            ax.set(xlim=(0,values.max()*1.05 or 1),
                   title=f'{label}\nn={d["stats"][kind+"_rows"]:,}',xlabel='Mean coverage')
            ax.title.set_color(color)
            ax.title.set_fontsize(8)
            ax.locator_params(axis='x',nbins=2)
        regions = sorted(d['regions'],key=lambda x: ({'rRNA':0,'CDS':1,'UTR5':2,'UTR3':3,'intron':4}.get(x,5),x))
        counts = [d['regions'][r] for r in regions]; total = sum(counts)
        axes[11].barh(regions,counts,color=[{'rRNA':'#138c81','CDS':'#ca791b','UTR5':'#ca791b','UTR3':'#ca791b','intron':'#3579b1'}.get(r,'#888888') for r in regions])
        axes[11].invert_yaxis()
        axes[11].set(title=f'{rna} partner regions\nTotal weight={total:,}',xlabel='Support count' if weighting=='count' else 'Cluster rows',xlim=(0,max(counts,default=1)*1.65))
        for n,value in enumerate(counts):
            axes[11].text(value+max(counts)*.02,n,f'{value:,} ({value/total:.1%})',va='center',fontsize=7)
    fig.suptitle(f'STAU clustered duplexes · {group.replace("_"," ")}',fontsize=19,weight='bold',y=.97)
    fig.text(.025,.91,f'{size}-nt bins · weight: {weighting} · raw coverage (rRNA-only input; no library RPM) · independent coverage scales',fontsize=11)
    fig.text(.025,.87,'Human rDNA = rRNA_45S assumed; rodents use successful target mappings only. Clippy tracks mark annotated arms / partners, not independent peak coverage.',fontsize=10)
    fig.colorbar(physical,cax=fig.add_axes([.08,.045,.25,.015]),orientation='horizontal',label='Physical distance (Å)')
    fig.colorbar(contact,cax=fig.add_axes([.48,.045,.35,.015]),orientation='horizontal',label='log₂(1 + contact weight)')
    # Match all vertical profiles to the square heatmap's actual drawing area.
    fig.canvas.draw()
    for axes in row_axes:
        position = axes[1].get_position()
        for ax in [axes[0]] + axes[3:]:
            box = ax.get_position()
            ax.set_position([box.x0, position.y0, box.width, position.height])
        for ax in axes[1:3]:
            box = ax.get_window_extent()
            assert abs(box.width-box.height) < .01, 'Heatmap must be square'
    for ext in ['png','pdf']:
        fig.savefig(out/f'{group}.profiles.{ext}',dpi=170)
    plt.close(fig)


def plot_facets(data, assets, segments, out, size, weighting):
    """One source per row; each coverage facet has its own range."""
    from matplotlib.backends.backend_pdf import PdfPages
    names = {'human': 'Human', 'mouse_lifted': 'Mouse lifted',
             'rat_lifted': 'Rat lifted', 'merged': 'All merged',
             'merged_both_arms_available': 'Merged: both arms available'}
    with PdfPages(out/'sources_and_merged.pdf') as pdf:
        for rna, asset in assets.items():
            length = asset['length']
            edges = np.r_[np.arange(0,length,size)+.5,length+.5]
            y = (edges[:-1]+edges[1:])/2
            fig = plt.figure(figsize=(22,18))
            grid = fig.add_gridspec(len(GROUPS), len(TRACKS)+1,
                width_ratios=[.8]+[1]*len(TRACKS), left=.07,right=.99,
                bottom=.05,top=.90,hspace=.70,wspace=.45)
            for ri, group in enumerate(GROUPS):
                es = fig.add_subplot(grid[ri,0])
                es.set_ylim(length+.5,.5)
                expansion(es,segments,rna,length)
                es.set_title('Expansion segments',fontsize=8)
                es.text(-.45,.5,names[group],transform=es.transAxes,
                        rotation=90,ha='center',va='center',fontsize=11,weight='bold')
                for ci,(kind,label,color) in enumerate(TRACKS,1):
                    ax = fig.add_subplot(grid[ri,ci])
                    d = data[group][rna]
                    values = binned(d['coverage'][kind],size)
                    ax.fill_betweenx(y,0,values,color=color,alpha=.65)
                    ax.plot(values,y,color=color,lw=.6)
                    ax.set(ylim=(length+.5,.5),xlim=(0,values.max()*1.05 or 1),
                           title=f'{label}\nn={d["stats"][kind+"_rows"]:,} rows',
                           xlabel='Mean coverage')
                    ax.title.set_fontsize(9)
                    ax.locator_params(axis='x',nbins=3)
                    if ci==1: ax.set_ylabel('RNA position (nt)')
                    else: ax.tick_params(labelleft=False)
                    if not values.any():
                        ax.text(.5,.5,'No positive annotation' if kind=='clippy_direct'
                                else 'No records',transform=ax.transAxes,
                                ha='center',va='center',fontsize=8,color='#777777')
            fig.suptitle(f'{rna} · human, lifted and merged source facets',fontsize=19,weight='bold',y=.98)
            fig.text(.07,.945,f'{size}-nt bins · raw {weighting} weights · independent coverage range per facet · expansion segments repeated per facet',fontsize=11)
            fig.text(.07,.923,'Clippy partner column = rRNA interaction footprints with a peak-annotated partner; it is not a peak track on rRNA.',fontsize=10)
            fig.savefig(out/f'{rna}.sources_and_merged.png',dpi=170)
            fig.savefig(out/f'{rna}.sources_and_merged.pdf')
            pdf.savefig(fig)
            plt.close(fig)
    # Remove the superseded overlay, which is replaced by the source facets.
    (out/'sources_and_merged.png').unlink(missing_ok=True)


def plot_aligned_species(data, assets, segments, out, size, weighting):
    """Narrow adjacent species tracks, with only one row per mature RNA."""
    sources = [('mouse_lifted','Mouse'), ('rat_lifted','Rat'), ('human','Human')]
    fig = plt.figure(figsize=(24,10))
    grid = fig.add_gridspec(2,len(TRACKS)+1,width_ratios=[1]+[2]*len(TRACKS),
        left=.035,right=.99,bottom=.10,top=.82,hspace=.62,wspace=.40)
    for ri,(rna,asset) in enumerate(assets.items()):
        length = asset['length']
        edges = np.r_[np.arange(0,length,size)+.5,length+.5]
        y = (edges[:-1]+edges[1:])/2
        es = fig.add_subplot(grid[ri,0])
        es.set_ylim(length+.5,.5)
        expansion(es,segments,rna,length)
        es.set_title(f'{rna}\nExpansion segments',fontsize=8)
        for ci,(kind,label,color) in enumerate(TRACKS,1):
            subgrid = grid[ri,ci].subgridspec(1,3,wspace=.13)
            for si,(group,name) in enumerate(sources):
                ax = fig.add_subplot(subgrid[0,si])
                d = data[group][rna]
                values = binned(d['coverage'][kind],size)
                ax.fill_betweenx(y,0,values,color=color,alpha=.65)
                ax.plot(values,y,color=color,lw=.55)
                ax.set(ylim=(length+.5,.5),xlim=(0,values.max()*1.05 or 1),
                       title=f'{name}\nn={d["stats"][kind+"_rows"]:,}')
                ax.title.set_fontsize(6)
                ax.tick_params(axis='x',labelsize=6)
                ax.locator_params(axis='x',nbins=2)
                if ci==1 and si==0:
                    ax.set_ylabel(f'{rna} position (nt)',fontsize=8)
                    ax.tick_params(axis='y',labelsize=7)
                else:
                    ax.tick_params(axis='y',left=False,labelleft=False)
                if si==1:
                    ax.set_xlabel('Mean coverage',fontsize=7)
                if not values.any():
                    ax.text(.5,.5,'No annotation' if kind=='clippy_direct' else 'No records',
                        transform=ax.transAxes,rotation=90,ha='center',va='center',
                        fontsize=6,color='#777777')
            box = grid[ri,ci].get_position(fig)
            fig.text((box.x0+box.x1)/2,box.y1+.065,label,ha='center',va='bottom',
                     fontsize=9,color=color)
    fig.suptitle('rRNA coverage · mouse, rat and human aligned',fontsize=19,weight='bold',y=.98)
    fig.text(.035,.935,f'{size}-nt bins · raw {weighting} weights · independent coverage ranges · mouse and rat use lifted human coordinates',fontsize=11)
    fig.text(.035,.900,'Each category: Mouse | Rat | Human. Clippy partner signal marks rRNA footprints linked to a peak-annotated partner.',fontsize=10)
    for ext in ['png','pdf']:
        fig.savefig(out/f'species_aligned.{ext}',dpi=170)
    plt.close(fig)


def plot_interaction_types(data, out, weighting):
    """Partner annotation bars, faceted by species and mature rRNA."""
    sources = [('human','Human'), ('mouse_lifted','Mouse lifted'),
               ('rat_lifted','Rat lifted')]
    preferred = ['rRNA','CDS','UTR5','UTR3','exon','intron','lncRNA','ncRNA','snoRNA','tRNA','intergenic','unknown']
    regions = set().union(*(data[g][rna]['regions'].keys()
                            for g,_ in sources for rna in ['18S','28S']))
    regions = sorted(regions,key=lambda r:(preferred.index(r) if r in preferred else len(preferred),r))
    colors = {'rRNA':'#138c81','CDS':'#ca791b','UTR5':'#ca791b',
              'UTR3':'#ca791b','exon':'#ca791b','intron':'#3579b1'}
    fig,axes = plt.subplots(2,3,figsize=(17,10))
    fig.subplots_adjust(left=.10,right=.98,bottom=.08,top=.80,hspace=.45,wspace=.62)
    for ri,rna in enumerate(['18S','28S']):
        for ci,(group,name) in enumerate(sources):
            ax = axes[ri,ci]
            counts = [data[group][rna]['regions'][r] for r in regions]
            total = sum(counts)
            ax.barh(range(len(regions)),counts,color=[colors.get(r,'#888888') for r in regions])
            ax.set(yticks=range(len(regions)),yticklabels=['Other rRNA' if r=='rRNA' else r for r in regions],
                   title=f'{rna} · {name}\nTotal weight = {total:,}',
                   xlabel='Support count' if weighting=='count' else 'Cluster rows',
                   xlim=(0,max(counts,default=0)*1.65 or 1))
            ax.invert_yaxis()
            for i,value in enumerate(counts):
                ax.text(value+max(counts,default=0)*.02,i,
                        f'{value:,} ({value/total:.1%})' if total else '0 (0.0%)',
                        fontsize=8,va='center')
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.tick_params(labelsize=9)
            ax.locator_params(axis='x',nbins=3)
    fig.suptitle('rRNA partner interaction types · separate species facets',fontsize=18,weight='bold',y=.97)
    fig.text(.10,.91,f'Raw {weighting} weights · independent count ranges · percentages within each rRNA/species facet',fontsize=11)
    fig.text(.10,.875,'Mouse and rat use lifted rRNA coordinates; all external rDNA partners are included in Other rRNA.',fontsize=10)
    for ext in ['png','pdf']:
        fig.savefig(out/f'interaction_types.species.{ext}',dpi=170)
    plt.close(fig)


def write_tsv(path, rows, fields):
    with path.open('w') as f:
        w = csv.DictWriter(f,fieldnames=fields,delimiter='\t'); w.writeheader(); w.writerows(rows)


def main():
    p = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--input',type=Path,default=Path.home()/'Downloads/data_duplex_events_all_species_with_selected_liftover_with_clippy_peak_info_rRNA_only.tsv')
    p.add_argument('--reference-repo',type=Path,default=Path.home()/'github/pdb-rna-benchmarking',help='Existing PARIS/KARR structural assets and expansion annotations')
    p.add_argument('--out',type=Path,default=Path(__file__).resolve().parent/'results')
    p.add_argument('--weight',choices=['count','rows'],default='count')
    p.add_argument('--bin-size',type=int,default=10)
    p.add_argument('--arc-min',type=float,default=2,help='Minimum raw weight per off-diagonal contact bin pair')
    args = p.parse_args()
    if args.bin_size < 1 or args.arc_min <= 0: p.error('Bin size and arc minimum must be positive')
    args.out.mkdir(parents=True,exist_ok=True)
    cache = args.reference_repo/'rrna_contacts/outputs/paper_gamma_rpm5'
    assets = {}; bounds = {}
    for rna in ['18S','28S']:
        folder = cache/rna
        m = json.loads((folder/'mapping_audit.json').read_text())
        bounds[rna] = (m['reference_start']-1,m['reference_end'])
        distance = np.load(folder/'centroid_distance.npy').copy()
        modelled = np.load(folder/'modelled_fraction.npy')
        distance[modelled<.5,:] = np.nan; distance[:,modelled<.5] = np.nan
        assets[rna] = dict(length=m['reference_end']-m['reference_start']+1,chain=m['chain'],distance=distance)
    with (args.reference_repo/'stau1_hiclip_profiles/expansion_segments.tsv').open() as f:
        segments = list(csv.DictReader(f,delimiter='\t'))
    data,audit,peaks = collect(args.input,bounds,args.bin_size,args.weight)
    records=[]; regions=[]; summary=[]
    for group in GROUPS:
        for rna,d in data[group].items():
            summary.append(dict(group=group,rna=rna,**dict(d['stats'])))
            np.savez_compressed(args.out/f'{group}.{rna}.npz',contacts_upper=d['contacts'],**d['coverage'])
            for kind,_,_ in TRACKS:
                for i,v in enumerate(binned(d['coverage'][kind],args.bin_size)):
                    records.append(dict(group=group,rna=rna,track=kind,start_nt=i*args.bin_size+1,
                                        end_nt=min((i+1)*args.bin_size,assets[rna]['length']),mean_coverage=v))
            regions.extend(dict(group=group,rna=rna,region=r,weight=w) for r,w in d['regions'].items())
        plot_group(group,data[group],assets,segments,args.out,args.bin_size,args.weight,args.arc_min)
        print('Wrote',group,flush=True)
    plot_facets(data, assets, segments, args.out, args.bin_size, args.weight)
    plot_aligned_species(data, assets, segments, args.out, args.bin_size, args.weight)
    plot_interaction_types(data, args.out, args.weight)
    write_tsv(args.out/'profiles.tsv',records,['group','rna','track','start_nt','end_nt','mean_coverage'])
    write_tsv(args.out/'partner_regions.tsv',regions,['group','rna','region','weight'])
    write_tsv(args.out/'clippy_peaks.tsv',peaks,['species','peak_id','region','cDNA_score'])
    fields=['group','rna']+sorted(set().union(*(r.keys() for r in summary))-{'group','rna'})
    write_tsv(args.out/'summary.tsv',summary,fields)
    (args.out/'audit.json').write_text(json.dumps(dict(input=str(args.input.resolve()),
        sha256=hashlib.sha256(args.input.read_bytes()).hexdigest(),counts=audit,bounds_0based_halfopen=bounds,
        weight=args.weight,bin_size=args.bin_size,arc_min=args.arc_min,
        assumptions=['Human rDNA and rRNA_45S share the U13369 coordinate system.',
                     'Input intervals are 1-based inclusive.',
                     'Count is support weight, not independently verified deduplicated reads.'],
        groups={'human':'native human coordinates','mouse_lifted':'mouse rRNA arms with successful mappings',
                'rat_lifted':'rat rRNA arms with successful mappings','merged':'sum of those three sources',
                'merged_both_arms_available':'human plus rodent rows where BOTH arms mapped'},
        partner_annotation='Mapped annotation when available; otherwise original annotation. Unmapped partners retained in main coverage, excluded from strict merged group.',
        clippy='Direct: rRNA arm has peak annotation. Partner: opposite arm has peak annotation. Coverage is interaction footprint, not peak width or cDNA score. Unique peak IDs exported separately, per source species.',
        contact_map='Within-mature-RNA cluster arm midpoint pairs; upper triangle stored once; displayed symmetrically. Arcs use raw bin-pair weight.',
        normalization='Raw support/row weights; no full-library denominator available.'),indent=2)+'\n')
    print('Results:',args.out.resolve())


if __name__=='__main__': main()
