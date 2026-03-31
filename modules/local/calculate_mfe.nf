process CALCULATE_MFE {
    tag "${meta.id}_${chunk.baseName}"
    label 'process_medium'

    input:
    tuple val(meta), path(chunk), path(sequence_table)

    output:
    tuple val(meta), path(chunk), path("${meta.id}_${chunk.baseName}_mfe.tsv"), emit: mfe

    script:
    def prefix = "${meta.id}_${chunk.baseName}"
    """
    python mfe_chunk.py \\
        --input ${sequence_table} \\
        --output ${prefix}_mfe.tsv \\
        --processes ${task.cpus}
    """
}