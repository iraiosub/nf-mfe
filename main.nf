#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
========================================================================================
    Chimeric Interaction Sequence Extraction Pipeline
========================================================================================
    Splits chimeric interaction tables, extracts sequences, and concatenates results
----------------------------------------------------------------------------------------
*/

// Include modules
include { SPLIT_TABLE } from './modules/local/split_table'
include { PREPARE_BED } from './modules/local/prepare_bed'
include { EXTRACT_SEQUENCES } from './modules/local/extract_sequences'
include { ADD_SEQUENCES } from './modules/local/add_sequences'
include { CALCULATE_MFE } from './modules/local/calculate_mfe'
include { CALCULATE_MFE_CONTROLS } from './modules/local/calculate_mfe_controls'
include { CONCATENATE_TABLES } from './modules/local/concatenate_tables'
include { PLOT_MFE_SUMMARY } from './modules/local/plot_mfe_summary'
include { METHOD_PLOTS } from './workflows/method_figures'

/*
========================================================================================
    Main Workflow
========================================================================================
*/

workflow {
    main:

    // Create input channel from samplesheet or glob pattern
    def ch_input
    if (params.input) {
        // Read samplesheet (TSV format: sample_id, file_path)
        ch_input = Channel.fromPath(params.input)
            .splitCsv(header: true, sep: '\t')
            .map { row ->
                def meta = [id: row.sample_id, method: row.method ?: row.sample_id]
                [meta, file(row.file_path, checkIfExists: true)]
            }
    } else if (params.input_pattern) {
        // Use glob pattern
        ch_input = Channel.fromPath(params.input_pattern, checkIfExists: true)
            .map { file_item ->
                def meta = [id: file_item.baseName.replaceAll(/\..*$/, '')]
                [meta, file_item]
            }
    } else {
        error "Please provide either --input (samplesheet) or --input_pattern (glob pattern)"
    }

    // Load reference genome FASTA
    def ch_fasta = file(params.fasta, checkIfExists: true)
    def ch_fai = file(params.fai, checkIfExists: true)

    // Step 1: Split tables into chunks
    SPLIT_TABLE(
        ch_input,
        params.chunk_size
    )

    // Flatten chunks to process each one separately
    def ch_chunks = SPLIT_TABLE.out.chunks
        .transpose()

    // Step 2: Prepare BED files for each chunk
    PREPARE_BED(ch_chunks)

    // Step 3: Extract sequences using bedtools getfasta
    EXTRACT_SEQUENCES(
        PREPARE_BED.out.beds,
        ch_fasta,
        ch_fai
    )

    // Step 4: Add sequences as new columns
    ADD_SEQUENCES(EXTRACT_SEQUENCES.out.sequences)

    // Step 5: Calculate MFE
    def run_controls = params.shuffled_mfe || params.flipped_arm_mfe

    def ch_mfe
    if (!run_controls) {

        CALCULATE_MFE(ADD_SEQUENCES.out.sequence_table)
        ch_mfe = CALCULATE_MFE.out.mfe

    } else {
        CALCULATE_MFE_CONTROLS(ADD_SEQUENCES.out.sequence_table)
        ch_mfe = CALCULATE_MFE_CONTROLS.out.mfe

    }

    // Step 6: Group chunks by sample and concatenate
    def ch_grouped = ch_mfe.groupTuple(by: 0)

    CONCATENATE_TABLES(ch_grouped)

    // Step 7: Plot concatenated shuffled-MFE tables
    if (run_controls && !params.plot_method_figures) {
        PLOT_MFE_SUMMARY(CONCATENATE_TABLES.out.final_table)
    }

    if (params.plot_method_figures) {
        METHOD_PLOTS(CONCATENATE_TABLES.out.final_table)
    }

    // Emit final outputs
    emit:
    final_tables = CONCATENATE_TABLES.out.final_table
}

/*
========================================================================================
    Workflow Event Handlers
========================================================================================
*/

workflow.onComplete {
    println "Pipeline completed at: ${workflow.complete}"
    println "Execution status: ${workflow.success ? 'OK' : 'failed'}"
    println "Duration: ${workflow.duration}"
}

workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}

// Plot existing final MFE tables without extracting sequences or recalculating MFE.
workflow METHOD_FIGURES {
    main:
    if (!params.input) {
        error 'METHOD_FIGURES requires --input TSV with sample_id, file_path and method columns'
    }
    def tables = Channel.fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            if (!row.sample_id || !row.file_path || !row.method) {
                error 'METHOD_FIGURES samplesheet requires nonempty sample_id, file_path and method'
            }
            if (!(row.sample_id ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/)) {
                error "Invalid sample_id: ${row.sample_id}; use letters, numbers, underscores, dots or hyphens"
            }
            tuple([id: row.sample_id, method: row.method], file(row.file_path, checkIfExists: true))
        }
    METHOD_PLOTS(tables)
}
