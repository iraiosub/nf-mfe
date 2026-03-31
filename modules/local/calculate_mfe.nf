process CALCULATE_MFE {
    tag "${meta.id}_${chunk.baseName}"
    label 'process_medium'

    container 'community.wave.seqera.io/library/pandas_viennarna:d01fcc9dd6f10930'

    input:
    tuple val(meta), path(chunk), path(sequence_table)

    output:
    tuple val(meta), path("*_mfe.tsv"), emit: mfe

    script:
    def prefix = "${chunk.baseName}"
    """
    mfe_chunk.py \\
        --input ${sequence_table} \\
        --output ${prefix}_mfe.tsv \\
        --processes ${task.cpus}
    """
}