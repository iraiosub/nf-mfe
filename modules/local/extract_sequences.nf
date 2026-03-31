process EXTRACT_SEQUENCES {
    tag "${meta.id}_${chunk.baseName}"
    label 'process_medium'

    container 'community.wave.seqera.io/library/bedtools_python:88ef333ca5f1123a'

    input:
    tuple val(meta), path(chunk), path(left_bed), path(right_bed)
    path fasta
    path fai

    output:
    tuple val(meta), path(chunk), path("*_left.fa"), path("*_right.fa"), emit: sequences

    script:
    def prefix = "${chunk.baseName}"
    """
    # Extract left sequences with strand awareness
    bedtools getfasta -s -nameOnly \\
        -fi ${fasta} \\
        -bed ${left_bed} \\
        -fo ${prefix}_left.fa

    # Extract right sequences with strand awareness
    bedtools getfasta -s -nameOnly \\
        -fi ${fasta} \\
        -bed ${right_bed} \\
        -fo ${prefix}_right.fa
    """
}
