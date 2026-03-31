process CONCATENATE_TABLES {
    tag "${meta.id}"
    label 'process_high'

    container 'community.wave.seqera.io/library/bedtools_python:88ef333ca5f1123a'

    publishDir params.outdir, mode: 'copy'
    
    input:
    tuple val(meta), path(chunks)
    
    output:
    tuple val(meta), path("${meta.id}_mfe.tsv"), emit: final_table
    
    script:
    """
    #!/usr/bin/env python3
    
    import glob
    import os
    
    sample_id = "${meta.id}"
    
    # Get all chunk files and sort them
    chunk_files = sorted(glob.glob("*_mfe.tsv"))
    
    if not chunk_files:
        print("Error: No chunk files found!")
        exit(1)
    
    with open(f"{sample_id}_mfe.tsv", 'w') as outfile:
        first_file = True
        
        for chunk_file in chunk_files:
            with open(chunk_file, 'r') as infile:
                if first_file:
                    # Write header from first file
                    outfile.write(infile.readline())
                    first_file = False
                else:
                    # Skip header for subsequent files
                    infile.readline()
                
                # Write all data lines
                for line in infile:
                    outfile.write(line)
    
    print(f"Concatenated {len(chunk_files)} chunks into {sample_id}_mfe.tsv")
    """
}
