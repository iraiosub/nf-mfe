process ADD_SEQUENCES {
    tag "${meta.id}_${chunk.baseName}"
    label 'process_low'

    container 'community.wave.seqera.io/library/bedtools_python:88ef333ca5f1123a'

    input:
    tuple val(meta), path(chunk), path(left_fa), path(right_fa)

    output:
    tuple val(meta), path(chunk), path("*_sequences.tsv"), emit: sequence_table

    script:
    def prefix = "${chunk.baseName}"
    """
    #!/usr/bin/env python3

    import csv
    import sys
    from collections import OrderedDict

    # Parse FASTA files into dictionaries
    def parse_fasta(fasta_file):
        sequences = {}
        current_id = None
        current_seq = []

        with open(fasta_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_id:
                        if current_id in sequences:
                            raise ValueError(f"Duplicate sequence name found in FASTA: '{current_id}'. Names must be unique within each chunk for sequence merging.")
                        sequences[current_id] = ''.join(current_seq)
                    # Extract ID from header (remove '>' and strand info)
                    # bedtools -nameOnly adds strand as suffix: readID(+) or readID(-)
                    raw_id = line[1:].split()[0]
                    # Strip strand suffix if present: remove (+) or (-)
                    if raw_id.endswith('(+)') or raw_id.endswith('(-)'):
                        current_id = raw_id[:-3]
                    else:
                        current_id = raw_id
                    current_seq = []
                else:
                    current_seq.append(line)
            
            if current_id:
                if current_id in sequences:
                    raise ValueError(f"Duplicate sequence name found in FASTA: '{current_id}'. Names must be unique within each chunk for sequence merging.")
                sequences[current_id] = ''.join(current_seq)
        
        return sequences
    
    # Parse both FASTA files
    left_seqs = parse_fasta('${left_fa}')
    right_seqs = parse_fasta('${right_fa}')
    
    # DEBUG: Print what we found
    print(f"DEBUG: Parsed {len(left_seqs)} left sequences", file=sys.stderr)
    print(f"DEBUG: Parsed {len(right_seqs)} right sequences", file=sys.stderr)
    if left_seqs:
        first_left_id = list(left_seqs.keys())[0]
        print(f"DEBUG: First left ID: '{first_left_id}'", file=sys.stderr)
        print(f"DEBUG: First left seq: '{left_seqs[first_left_id][:50]}...'", file=sys.stderr)
    
    # Read input table and add sequence columns
    matched_count = 0
    total_count = 0
    seen_ids = set()
    
    with open('${chunk}', 'r') as infile, \\
         open('${prefix}_sequences.tsv', 'w') as outfile:
        
        reader = csv.DictReader(infile, delimiter='\\t')
        
        # Add new columns to fieldnames
        fieldnames = reader.fieldnames + ['lseq', 'rseq']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter='\\t')
        writer.writeheader()
        
        for row in reader:
            read_id = row['name']
            total_count += 1

            if read_id in seen_ids:
                raise ValueError(f"Duplicate name found in chunk '${chunk}': '{read_id}'. Names must be unique within each chunk for sequence merging.")
            seen_ids.add(read_id)
            
            # DEBUG: Print first few IDs we're looking for
            if total_count <= 3:
                print(f"DEBUG: Looking for read ID: '{read_id}'", file=sys.stderr)
            
            # Add sequences to the row
            row['lseq'] = left_seqs.get(read_id, '')
            row['rseq'] = right_seqs.get(read_id, '')
            
            if row['lseq'] and row['rseq']:
                matched_count += 1
            
            writer.writerow(row)
    
    print(f"DEBUG: Matched {matched_count}/{total_count} rows with sequences", file=sys.stderr)
    """
}
