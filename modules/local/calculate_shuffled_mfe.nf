process CALCULATE_SHUFFLED_MFE {
    tag "${meta.id}_${chunk.baseName}"
    label 'process_low'

    container 'community.wave.seqera.io/library/pandas_ushuffle_viennarna:aa51db9ac2370318'

    input:
    tuple val(meta), path(chunk), path(mfe_table)

    output:
    tuple val(meta), path("*shuffled_mfe.tsv"), emit: shuffled_mfe
    tuple val(meta), path("*_shuffled_mfe_all.tsv"), emit: shuffled_mfe_all

    script:
    def prefix = "${chunk.baseName}"
    """
    mfe_shuffle_chunk.py \
        --input ${mfe_table} \
        --output-summary ${prefix}_shuffled_mfe.tsv \
        --output-detail ${prefix}_shuffled_mfe_all.tsv \
        --processes ${task.cpus} \
        --n-shuffles ${params.n_shuffles} \
        --klet ${params.klet_shuffles}
    """
}