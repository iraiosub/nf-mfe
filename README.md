# Chimeric Interaction Sequence Extraction and MFE Computation Pipeline

## This Nextflow pipeline processes chimeric interaction tables by:
1. Splitting large tables into 10k row chunks (preserving headers)
2. Extracting FASTA sequences for left and right coordinates using bedtools
3. Adding extracted sequences as new columns (lseq, rseq)
4. Performing MFE calculation on the extracted sequences ([RNAduplex](https://www.tbi.univie.ac.at/RNA/ViennaRNA/doc/html/man/RNAduplex.html))
5. Optionally generating shuffled control sequences for MFE comparisons ([uShuffle](https://link.springer.com/article/10.1186/1471-2105-9-192))
6. Concatenating all chunks back into a single table with one header

## Requirements

- Nextflow >= 25.10.2
- Docker or Singularity

## Quick Start

> Important: we assume BED-style 0-based half-open coordinates are given in the input. 

### Option 1: Using a samplesheet

```bash
nextflow run main.nf \
  --input samplesheet.tsv \
  --fasta /path/to/genome.fa \
  --fai /path/to/genome.fa.fai \
  --chunk_size <n> \
  --shuffled_mfe \
  --n_shuffles <n> \
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
  --outdir results \
  <...>
```

## Parameters

### Required
- `--fasta`: Path to reference genome FASTA file
- `--fai`: Path to reference genome FASTA index

### Input (choose one)
- `--input`: Path to samplesheet (TSV format)
- `--input_pattern`: Glob pattern for input files (e.g., "data/*.txt")

### Optional
- `--outdir`: Output directory (default: './results')
- `--chunk_size`: Number of rows per chunk (default: 10000)
- `--shuffled_mfe`: Enable MFE calculations for shuffled control sequences (default: false)
- `--n_shuffles`: Number of times the sequence is shuffled (default: 100)
- `--klet_shuffles`: klet for sequence shuffling (default: 2)

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
- `mfe` : MFE (kcal/mol)
- `dot_bracket`: Dot-bracket notation corresponding to the duplex
- `mean_shuffled_mfe`: Mean MFE across shuffles
- `sd_shuffled_mfe`
- `delta_mfe`: Difference between Observed MFE and the mean MFE across shuffles
- `zscore_mfe`
- `empirical_p_lower`
- `n_shuffles_ok`: Number of shuffles for which MFE was succesfully computed


## Execution Profiles

### Docker
```bash
nextflow run main.nf -profile docker --input ... --fasta ...
```

### Singularity
```bash
nextflow run main.nf -profile singularity --input ... --fasta ...
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
CALCULATE_SHUFFLED_MFE or CALCULATE_SHUFFLED_MFE (MFE-related columns)
    ↓
CONCATENATE_TABLES (merge chunks)
    ↓
PLOT_MFE_SUMMARY (per sample)
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

## Notes to self

### Settings for local testing
`conda activate nfcore_tools_34`

`nextflow run main.nf --fasta /Volumes/lab-ulej/home/users/luscomben/users/iosubi/projects/structurome_blencowe/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa --input samplesheet.example.tsv --outdir results -profile docker --chunk_size 10 -resume --fai /Volumes/lab-ulej/home/users/luscomben/users/iosubi/projects/structurome_blencowe/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai --shuffled_mfe --n_shuffles 4`

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
ADD_SEQUENCES → [meta, chunk_0001_sequences.tsv]
                [meta, chunk_0002_sequences.tsv]
                [meta, chunk_0003_sequences.tsv]

etc.
```
