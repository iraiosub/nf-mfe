# Pipeline Summary

## Overview
This pipeline extracts FASTA sequences for chimeric interaction coordinates and adds them as new columns to your input tables.

## Pipeline Structure

```
chim-seq-pipeline/
├── main.nf                              # Main workflow script
├── nextflow.config                      # Configuration and parameters
├── README.md                            # Usage documentation
├── samplesheet.example.tsv              # Example samplesheet
└── modules/
    └── local/
        ├── split_table.nf               # Split tables into 10k row chunks
        ├── prepare_bed.nf               # Convert chunks to BED format
        ├── extract_sequences.nf         # Extract sequences with bedtools
        ├── add_sequences.nf             # Add sequences as new columns
        └── concatenate_tables.nf        # Merge chunks with single header
```

## Modules Description

### 1. SPLIT_TABLE
- **Input**: Table file with chimeric interactions
- **Output**: Multiple chunk files (10k rows each, with headers)
- **Function**: Splits large input tables into manageable chunks

### 2. PREPARE_BED
- **Input**: Table chunk
- **Output**: Two BED files (left and right coordinates)
- **Function**: Converts table coordinates to BED format for bedtools

### 3. EXTRACT_SEQUENCES
- **Input**: BED files + Reference FASTA
- **Output**: FASTA files with extracted sequences
- **Function**: Uses bedtools getfasta with strand awareness (-s flag)
- **Features**: 
  - Respects strand orientation
  - Extracts sequences for both left and right coordinates

### 4. ADD_SEQUENCES
- **Input**: Table chunk + FASTA sequences
- **Output**: Table chunk with added lseq and rseq columns
- **Function**: Parses FASTA and merges sequences back to table

### 5. CONCATENATE_TABLES
- **Input**: All processed chunks for a sample
- **Output**: Single table with all rows and one header
- **Function**: Merges all chunks, keeping only the first header

## Data Flow

```
Input: chimeric_interactions.txt
├── lchr, ll, lr, lstrand  (left coordinate)
├── rchr, rl, rr, rstrand  (right coordinate)
├── name, mapq
↓
Split into chunks (10k rows each)
↓
For each chunk:
  ├── Create left.bed:  lchr  ll  lr  name  mapq  lstrand
  ├── Create right.bed: rchr  rl  rr  name  mapq  rstrand
  ↓
  ├── bedtools getfasta -s → left.fa
  ├── bedtools getfasta -s → right.fa
  ↓
  Add sequences to table:
  ├── lseq (sequence from left coordinate)
  ├── rseq (sequence from right coordinate)
↓
Concatenate all chunks
↓
Output: sample_with_sequences.txt
├── All original columns
├── lseq (new column)
└── rseq (new column)
```

## Usage Examples

### Basic usage with samplesheet
```bash
nextflow run main.nf \
  --input samplesheet.tsv \
  --fasta genome.fa \
  --outdir results
```

### With Docker
```bash
nextflow run main.nf \
  -profile docker \
  --input samplesheet.tsv \
  --fasta genome.fa \
  --outdir results
```

### Using glob pattern
```bash
nextflow run main.nf \
  --input_pattern 'data/*.chim.txt' \
  --fasta genome.fa \
  --chunk_size 5000 \
  --outdir results
```

### With custom chunk size
```bash
nextflow run main.nf \
  --input samplesheet.tsv \
  --fasta genome.fa \
  --chunk_size 20000 \
  --outdir results
```

## Parameters

| Parameter | Required | Default | Description |
|-----------|----------|---------|-------------|
| `--input` | Either this or input_pattern | null | Path to samplesheet (TSV) |
| `--input_pattern` | Either this or input | null | Glob pattern for input files |
| `--fasta` | Yes | null | Reference genome FASTA file |
| `--outdir` | No | ./results | Output directory |
| `--chunk_size` | No | 10000 | Rows per chunk |

## Output Files

```
results/
├── {sample_id}_with_sequences.txt    # Final output with sequences
└── pipeline_info/
    ├── timeline.html
    ├── report.html
    ├── trace.txt
    └── dag.html
```

## Performance Notes

- **Parallelization**: Each chunk is processed independently
- **Memory**: Most processes use 2-4GB
- **Chunk size**: Default 10k rows balances parallelization and overhead
- **Scaling**: Processing time scales linearly with number of chunks

## Key Features

✅ Preserves headers in all chunks  
✅ Strand-aware sequence extraction  
✅ Automatic parallelization per chunk  
✅ Support for multiple samples via samplesheet  
✅ Configurable chunk size  
✅ Container support (Docker/Singularity/Conda)  
✅ Comprehensive pipeline reports
