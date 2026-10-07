include { PREPARE_METHOD_PLOT_DATA; PLOT_METHOD_FIGURES } from '../modules/local/plot_method_figures'

workflow METHOD_PLOTS {
    take:
    final_tables

    main:
    if (params.gene_types && !params.gtf) {
        error '--gene_types requires --gtf'
    }
    def gtf = params.gtf ? file(params.gtf, checkIfExists: true) : []
    PREPARE_METHOD_PLOT_DATA(final_tables, gtf, params.gene_types, params.exclude_chromosomes)
    def features = PREPARE_METHOD_PLOT_DATA.out.features.map { meta, table -> table }.collect()
    PLOT_METHOD_FIGURES(features, params.gtf != null && params.gtf.toString() != '')

    emit:
    overview_pdf = PLOT_METHOD_FIGURES.out.overview_pdf
    overview_png = PLOT_METHOD_FIGURES.out.overview_png
    span_pdf = PLOT_METHOD_FIGURES.out.span_pdf
    span_png = PLOT_METHOD_FIGURES.out.span_png
    source = PLOT_METHOD_FIGURES.out.source
    counts = PLOT_METHOD_FIGURES.out.counts
    qc = PREPARE_METHOD_PLOT_DATA.out.qc
}
