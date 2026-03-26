process SPLIT_TABLE {
    tag "${meta.id}"
    
    input:
    tuple val(meta), path(table)
    val chunk_size
    
    output:
    tuple val(meta), path("${meta.id}_chunk_*.txt"), emit: chunks
    
    script:
    """
    #!/usr/bin/env python3
    
    import sys
    
    chunk_size = ${chunk_size}
    input_file = "${table}"
    sample_id = "${meta.id}"
    
    with open(input_file, 'r') as f:
        header = f.readline()
        
        chunk_num = 1
        line_count = 0
        chunk_file = None
        
        for line in f:
            if line_count % chunk_size == 0:
                if chunk_file:
                    chunk_file.close()
                chunk_file = open(f"{sample_id}_chunk_{chunk_num:04d}.txt", 'w')
                chunk_file.write(header)
                chunk_num += 1
            
            chunk_file.write(line)
            line_count += 1
        
        if chunk_file:
            chunk_file.close()
    
    print(f"Split {line_count} rows into {chunk_num - 1} chunks", file=sys.stderr)
    """
}
