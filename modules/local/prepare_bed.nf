process PREPARE_BED {
    tag "${meta.id}_${chunk.baseName}"
    label 'process_low'

    container 'community.wave.seqera.io/library/bedtools_python:88ef333ca5f1123a'

    input:
    tuple val(meta), path(chunk)
    
    output:
    tuple val(meta), path(chunk), path("*_left.bed"), path("*_right.bed"), emit: beds
    
    script:
    def prefix = "${chunk.baseName}"
    """
    #!/usr/bin/env python3
    
    import csv
    
    input_file = "${chunk}"
    prefix = "${prefix}"
    
    left_bed = open(f"{prefix}_left.bed", 'w')
    right_bed = open(f"{prefix}_right.bed", 'w')
    
    with open(input_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\\t')
        has_mapq = reader.fieldnames and 'mapq' in reader.fieldnames
        
        for row in reader:
            score = row['mapq'] if has_mapq and row.get('mapq', '') else '0'

            # Left BED entry
            left_bed.write(f"{row['lchr']}\\t{row['ll']}\\t{row['lr']}\\t{row['name']}\\t{score}\\t{row['lstrand']}\\n")
            
            # Right BED entry
            right_bed.write(f"{row['rchr']}\\t{row['rl']}\\t{row['rr']}\\t{row['name']}\\t{score}\\t{row['rstrand']}\\n")
    
    left_bed.close()
    right_bed.close()
    """
}
