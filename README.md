# Chimeric Interaction Sequence Extraction Pipeline

## This Nextflow pipeline processes chimeric interaction tables by:
1. Splitting large tables into 10k row chunks (preserving headers)
2. Extracting FASTA sequences for left and right coordinates using bedtools
3. Adding extracted sequences as new columns (lseq, rseq)
4. Concatenating all chunks back into a single table with one header

## Requirements

- Nextflow >= 23.04.0
- bedtools >= 2.31.0
- Python 3

## Quick Start

### Option 1: Using a samplesheet

```bash
nextflow run main.nf \
  --input samplesheet.tsv \
  --fasta /path/to/genome.fa \
  --outdir results
```

**Samplesheet format** (TSV with header):
```
sample_id	file_path
GM12878_rep1	/path/to/GM12878_rep1.chim.txt
GM12878_rep2	/path/to/GM12878_rep2.chim.txt
```

### Option 2: Using a glob pattern

```bash
nextflow run main.nf \
  --input_pattern 'data/*.chim.txt' \
  --fasta /path/to/genome.fa \
  --outdir results
```

## Parameters

### Required
- `--fasta`: Path to reference genome FASTA file

### Input (choose one)
- `--input`: Path to samplesheet (TSV format)
- `--input_pattern`: Glob pattern for input files (e.g., "data/*.txt")

### Optional
- `--outdir`: Output directory (default: './results')
- `--chunk_size`: Number of rows per chunk (default: 10000)

## Input File Format

Input files should be tab-delimited with the following columns:
```
lchr	ll	lr	lstrand	rchr	rl	rr	rstrand	name	mapq
```

## Output

For each sample, the pipeline produces:
- `{sample_id}_with_sequences.txt`: Final table with added sequence columns

Output columns include all original columns plus:
- `lseq`: Extracted sequence for left coordinate (strand-aware)
- `rseq`: Extracted sequence for right coordinate (strand-aware)

## Execution Profiles

### Docker
```bash
nextflow run main.nf -profile docker --input ... --fasta ...
```

### Singularity
```bash
nextflow run main.nf -profile singularity --input ... --fasta ...
```

### Conda
```bash
nextflow run main.nf -profile conda --input ... --fasta ...
```

## Pipeline Overview

```
Input Tables
    ↓
SPLIT_TABLE (10k rows per chunk)
    ↓
PREPARE_BED (convert to BED format)
    ↓
EXTRACT_SEQUENCES (bedtools getfasta -s)
    ↓
ADD_SEQUENCES (add lseq, rseq columns)
    ↓
CONCATENATE_TABLES (merge chunks)
    ↓
Final Output
```

## Example

```bash
# Using samplesheet
nextflow run main.nf \
  --input samplesheet.tsv \
  --fasta hg38.fa \
  --chunk_size 10000 \
  --outdir results \
  -profile docker
```
### Settings for local testing
`conda activate nfcore_tools_34`

`nextflow run main.nf --fasta /Volumes/lab-ulej/home/users/luscomben/users/iosubi/projects/structurome_blencowe/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa --input samplesheet.example.tsv --outdir results -profile docker`

#### Data flow with chunk matching
```
SPLIT_TABLE → [meta, chunk_0001.txt]
              [meta, chunk_0002.txt]
              [meta, chunk_0003.txt]
                     ↓
PREPARE_BED → [meta, chunk_0001.txt, chunk_0001_left.bed, chunk_0001_right.bed]
              [meta, chunk_0002.txt, chunk_0002_left.bed, chunk_0002_right.bed]
              [meta, chunk_0003.txt, chunk_0003_left.bed, chunk_0003_right.bed]
                     ↓
EXTRACT_SEQUENCES → [meta, chunk_0001.txt, left.fa, right.fa]
                    [meta, chunk_0002.txt, left.fa, right.fa]
                    [meta, chunk_0003.txt, left.fa, right.fa]
                     ↓
ADD_SEQUENCES → [meta, chunk_0001_with_sequences.txt]
                [meta, chunk_0002_with_sequences.txt]
                [meta, chunk_0003_with_sequences.txt]
```
