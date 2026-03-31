process PLOT_MFE_SUMMARY {
    tag "${meta.id}"
    label 'process_low'

    container 'community.wave.seqera.io/library/matplotlib_numpy_pandas_scipy:3a411aa680dcde7e'

    input:
    tuple val(meta), path(final_table)

    output:
    tuple val(meta), path("${meta.id}.mfe_summary_plots.png"), emit: png
    tuple val(meta), path("${meta.id}.mfe_summary_plots.pdf"), emit: pdf
    tuple val(meta), path("${meta.id}.mfe_summary_stats.tsv"), emit: stats

    script:
    """
    plot_mfe_summary.py \\
        --input ${final_table} \\
        --outdir . \\
        --sample-name ${meta.id}
    """
}