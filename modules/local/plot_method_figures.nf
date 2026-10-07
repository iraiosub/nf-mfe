process PREPARE_METHOD_PLOT_DATA {
    tag "${meta.id}"
    publishDir "${params.outdir}/prepare_method_plot_data", mode: 'copy', pattern: '*.plot_qc.json'
    label 'process_low'
    container 'community.wave.seqera.io/library/pandas_ushuffle_viennarna:aa51db9ac2370318'

    input:
    tuple val(meta), path(final_table)
    path gtf
    val gene_types
    val exclude_chromosomes

    output:
    tuple val(meta), path("${meta.id}.plot_features.tsv.gz"), emit: features
    tuple val(meta), path("${meta.id}.plot_qc.json"), emit: qc

    script:
    def quote = { value -> "'" + value.toString().replace("'", "'\"'\"'") + "'" }
    def gtf_arg = gtf ? "--gtf ${quote(gtf)}" : ''
    """
    prepare_method_plot_data.py \\
        --input ${quote(final_table)} \\
        --output ${quote(meta.id + '.plot_features.tsv.gz')} \\
        --sample ${quote(meta.id)} \\
        --method ${quote(meta.method ?: meta.id)} \\
        --processes ${task.cpus} \\
        --gene-types ${quote(gene_types ?: '')} \\
        --exclude-chromosomes ${quote(exclude_chromosomes ?: '')} \\
        ${gtf_arg}
    """
}

process PLOT_METHOD_FIGURES {
    label 'process_medium'
    container 'community.wave.seqera.io/library/matplotlib_numpy_pandas_scipy:3a411aa680dcde7e'

    input:
    path feature_tables, stageAs: 'features/*'
    val stratify

    output:
    path 'method_overview.png', emit: overview_png
    path 'method_overview.pdf', emit: overview_pdf
    path 'method_span.png', emit: span_png
    path 'method_span.pdf', emit: span_pdf
    path 'plot_source.tsv.gz', emit: source
    path '*group_counts.tsv', emit: counts
    path 'plot_notes.json', emit: notes

    script:
    def stratify_arg = stratify ? '--stratify' : ''
    """
    plot_method_figures.py --inputs features/*.tsv.gz --outdir . ${stratify_arg}
    """
}
