process PREPARE_BED {
    tag "${meta.id}_${chunk.baseName}"
    
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
        
        for row in reader:
            # Left BED entry
            left_bed.write(f"{row['lchr']}\\t{row['ll']}\\t{row['lr']}\\t{row['name']}\\t{row['mapq']}\\t{row['lstrand']}\\n")
            
            # Right BED entry
            right_bed.write(f"{row['rchr']}\\t{row['rl']}\\t{row['rr']}\\t{row['name']}\\t{row['mapq']}\\t{row['rstrand']}\\n")
    
    left_bed.close()
    right_bed.close()
    """
}
