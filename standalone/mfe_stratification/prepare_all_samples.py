#!/usr/bin/env python3
"""Inventory all final tables, reuse matching annotations and annotate missing samples."""
import argparse
import hashlib
import json
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from dateutil.tz import gettz
import pandas as pd


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--results-dir',required=True)
    p.add_argument('--existing-analysis-dir',required=True)
    p.add_argument('--gtf',required=True)
    p.add_argument('--outdir',required=True)
    args=p.parse_args()
    out=Path(args.outdir)
    if out.exists() and any(out.iterdir()):
        raise ValueError('All-sample output folder must be empty')
    out.mkdir(parents=True,exist_ok=True)
    old=Path(args.existing_analysis_dir)
    manifest=json.loads((old/'manifest.json').read_text())
    sha=hashlib.sha256(Path(args.gtf).read_bytes()).hexdigest()
    if manifest['gtf_sha256'] != sha:
        raise ValueError('Cached annotation GTF does not match')
    cached={e['sample']: e for e in manifest['selected']}
    selected=[]; missing=[]; sources={}
    for path in sorted(Path(args.results_dir).glob('*_mfe.tsv')):
        sample=path.name[:-len('_mfe.tsv')]
        entry=dict(sample=sample,path=str(path.resolve()),bytes=path.stat().st_size,
                   modification_date=datetime.fromtimestamp(path.stat().st_mtime,gettz('Europe/London')).date().isoformat())
        selected.append(entry)
        cache=cached.get(sample)
        if cache and all(cache[k]==entry[k] for k in ['path','bytes','modification_date']) and (old/sample/'annotated_pairs.tsv.gz').exists():
            sources[sample]=old/sample
        else:
            missing.append(entry)
    if not selected:
        raise ValueError('No final MFE tables found')
    pd.DataFrame(selected).to_csv(out/'all_processed_samples.tsv',sep='\t',index=False)
    if missing:
        sheet=out/'additional_samples.tsv'
        pd.DataFrame([dict(sample_id=e['sample'],file_path=e['path']) for e in missing]).to_csv(sheet,sep='\t',index=False)
        extra=out/'additional_annotation'
        subprocess.run([sys.executable,str(Path(__file__).with_name('stratify_mfe.py')),
                        '--samplesheet',str(sheet),'--results-dir',args.results_dir,
                        '--gtf',args.gtf,'--date','all','--outdir',str(extra)],check=True)
        sources.update({e['sample']:extra/e['sample'] for e in missing})
    for sample,source in sources.items():
        (out/sample).symlink_to(source.resolve(),target_is_directory=True)
    qc=[pd.read_csv(old/'sample_qc.tsv',sep='\t')]
    if missing:
        qc.append(pd.read_csv(out/'additional_annotation/sample_qc.tsv',sep='\t'))
    q=pd.concat(qc,ignore_index=True)
    q=q[q['sample'].isin(sources)].drop_duplicates('sample',keep='last')
    q.to_csv(out/'sample_qc.tsv',sep='\t',index=False)
    summaries=[pd.read_csv(sources[e['sample']]/'summary.tsv',sep='\t') for e in selected]
    pd.concat(summaries,ignore_index=True).to_csv(out/'all_samples_summary.tsv',sep='\t',index=False)
    result=dict(arguments=vars(args),selected=selected,skipped=[],gtf_sha256=sha,
                annotation_sources={k:str(v.resolve()) for k,v in sources.items()},
                completed=datetime.now(gettz('Europe/London')).isoformat())
    (out/'manifest.json').write_text(json.dumps(result,indent=2)+'\n')
    print(f'Prepared {len(selected)} samples ({len(missing)} newly annotated): {out}',flush=True)


if __name__=='__main__':
    main()
