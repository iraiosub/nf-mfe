import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'bin'))
from prepare_method_plot_data import annotate,load_gtf


def fixture(directory, controls=True):
    rows=[]
    for i,(c,ll,rl) in enumerate([('chr1',110,210),('1',160,210),('chr1',510,710),
                                  ('chr1',1200,1250),('rDNA',10,20),('chr1',20,30)]):
        rows.append(dict(name='hybrid'+str(i),row_id=str(i),lchr=c,ll=ll,lr=ll+2,lstrand='+',
                         rchr='chr2' if i==5 else c,rl=rl,rr=rl+2,rstrand='+',
                         lseq='GG',rseq='CC',dot_bracket='((&))',mfe=-2,
                         **({'mean_shuffled_mfe':-1} if controls else {})))
    table=directory/'input.tsv'
    frame=pd.DataFrame(rows)
    if not controls: frame=frame.drop(columns='row_id')
    frame.to_csv(table,sep='\t',index=False)
    gtf=directory/'genes.gtf'
    gtf.write_text(''.join('chr1\ttest\t%s\t%s\t%s\t.\t+\t.\tgene_id "%s"; gene_type "%s";\n'%r for r in [
        ('gene',101,400,'a','protein_coding'),('exon',101,150,'a','protein_coding'),
        ('exon',201,250,'a','protein_coding'),('gene',501,800,'b','lncRNA'),
        ('exon',501,550,'b','lncRNA'),('exon',701,750,'b','lncRNA')]))
    return table,gtf


class MethodFiguresTests(unittest.TestCase):
    def run_prepare(self,directory,extra=(),controls=True):
        table,gtf=fixture(directory,controls)
        output=directory/'sample.plot_features.tsv.gz'
        subprocess.run([sys.executable,str(ROOT/'bin/prepare_method_plot_data.py'),'--input',str(table),
                        '--output',str(output),'--sample','sample','--method','Test method',*extra],check=True,
                       stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        return pd.read_csv(output,sep='\t'),output,gtf

    def test_no_gtf_uses_all_cis_and_excludes_either_arm_chromosome(self):
        with tempfile.TemporaryDirectory() as tmp:
            d=Path(tmp);data,path,gtf=self.run_prepare(d,['--exclude-chromosomes','rDNA'])
            self.assertEqual(len(data),4)
            self.assertEqual(set(data.annotation_class),{'all'})
            self.assertIn('hybrid3',set(data.name))  # Intergenic is unknowable without GTF.
            self.assertTrue((data.base_pair_count==2).all())
            self.assertTrue((data.mfe_per_bp==-1).all())
            self.assertTrue((data.mean_shuffled_mfe_per_bp==-.5).all())
            self.assertTrue((data.paired_gc_percent==100).all())
            self.assertEqual(data.span_nt.iloc[0],102)

    def test_gtf_selects_same_gene_and_exact_biotypes(self):
        with tempfile.TemporaryDirectory() as tmp:
            d=Path(tmp);table,gtf=fixture(d)
            data,path,gtf=self.run_prepare(d,['--gtf',str(gtf),'--gene-types','protein_coding'])
            self.assertEqual(data.name.tolist(),['hybrid0','hybrid1'])
            self.assertEqual(data.annotation_class.tolist(),['pure_exon','intron_overlapping'])
            self.assertEqual(json.loads((d/'sample.plot_qc.json').read_text())['counts']['excluded_gene_rows'],3)
            with self.assertRaisesRegex(ValueError,'absent'):load_gtf(gtf,['protein_codin'])

    def test_gene_filter_requires_gtf(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(subprocess.CalledProcessError):
                self.run_prepare(Path(tmp),['--gene-types','protein_coding'])

    def test_observed_only_and_exported_plot_counts(self):
        with tempfile.TemporaryDirectory() as tmp:
            d=Path(tmp);data,path,gtf=self.run_prepare(d,controls=False)
            self.assertEqual(int(data.matched_eligible.sum()),0)
            self.assertEqual(data.row_id.astype(int).tolist(),data.source_row_index.tolist())
            out=d/'plots'
            subprocess.run([sys.executable,str(ROOT/'bin/plot_method_figures.py'),'--inputs',str(path),
                            '--outdir',str(out)],check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            self.assertEqual(len(list(out.glob('*.pdf'))),2)
            self.assertEqual(len(list(out.glob('*.png'))),2)
            source=pd.read_csv(out/'plot_source.tsv.gz',sep='\t');self.assertEqual(len(source),5)
            counts=pd.read_csv(out/'plot_group_counts.tsv',sep='\t')
            self.assertEqual(set(counts.source),{'observed'})
            for metric in ['mfe','mfe_per_bp']:
                self.assertEqual(counts[(counts.figure=='span')&(counts.metric==metric)].n.sum(),5)
            self.assertEqual(pd.read_csv(out/'sample_group_counts.tsv',sep='\t').source_rows.sum(),5)

    def test_overlapping_genes_cannot_fabricate_same_gene(self):
        from prepare_method_plot_data import Intervals
        index={('gene','1','+'):Intervals([(0,20),(10,30)])}
        data=pd.DataFrame(dict(lchr=['1'],rchr=['1'],lstrand=['+'],rstrand=['+'],ll=[1],lr=[3],rl=[25],rr=[28]))
        keep,classes=annotate(data,index)
        self.assertFalse(keep[0])

if __name__=='__main__':unittest.main()
