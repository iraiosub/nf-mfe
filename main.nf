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
include { CONCATENATE_TABLES } from './modules/local/concatenate_tables'

/*
========================================================================================
    Main Workflow
========================================================================================
*/

workflow {

    // Create input channel from samplesheet or glob pattern
    def ch_input
    if (params.input) {
        // Read samplesheet (TSV format: sample_id, file_path)
        ch_input = channel
            .fromPath(params.input)
            .splitCsv(header: true, sep: '\t')
            .map { row ->
                def meta = [id: row.sample_id]
                [meta, file(row.file_path, checkIfExists: true)]
            }
    } else if (params.input_pattern) {
        // Use glob pattern
        ch_input = channel
            .fromPath(params.input_pattern, checkIfExists: true)
            .map { file_item ->
                def meta = [id: file_item.baseName.replaceAll(/\..*$/, '')]
                [meta, file_item]
            }
    } else {
        error "Please provide either --input (samplesheet) or --input_pattern (glob pattern)"
    }

    // Load reference genome FASTA
    def ch_fasta = channel.fromPath(params.fasta, checkIfExists: true)
    ch_fasta = ch_fasta.first()

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
        ch_fasta
    )

    // Step 4: Add sequences as new columns
    ADD_SEQUENCES(EXTRACT_SEQUENCES.out.sequences)

    // Step 5: Group chunks by sample and concatenate
    def ch_grouped = ADD_SEQUENCES.out.table
        .groupTuple(by: 0)

    CONCATENATE_TABLES(ch_grouped)

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
