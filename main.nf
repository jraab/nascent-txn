#!/usr/bin/env nextflow

/*
========================================================================================
    TT-seq Processing Pipeline
========================================================================================
    A Nextflow DSL2 implementation of the Raab Lab TT-seq processing pipeline
    
    Authors: Jesse Raab Lab, UNC Chapel Hill
    Based on: crickbabs/DRB_TT-seq workflow
    Structure: raab-lab/cut-n-run template
========================================================================================
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/

WorkflowMain.initialise(workflow, params, log)

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

include { TTSEQ } from './subworkflows/ttseq'

//
// WORKFLOW: Run main TT-seq analysis workflow
//
workflow {
    
    if (params.help) {
        WorkflowMain.helpMessage(workflow, params, log)
        exit 0
    }

    //
    // Create channel from input file
    //
    Channel
        .fromPath(params.sample_sheet, checkIfExists: true)
        .set { ch_samplesheet }

    //
    // WORKFLOW: Run pipeline
    //
    TTSEQ (
        ch_samplesheet
    )
}

/*
========================================================================================
    FUNCTIONS
========================================================================================
*/

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id           = row.sample
    meta.single_end   = row.containsKey("single_end") ? row.single_end.toBoolean() : false
    meta.strandedness = row.containsKey("strandedness") ? row.strandedness : 'unstranded'
    
    // add path(s) of the fastq file(s) to the meta map
    def fastq_meta = []
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }
    if (meta.single_end) {
        fastq_meta = [ meta, [ file(row.fastq_1) ] ]
    } else {
        if (!file(row.fastq_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
        }
        fastq_meta = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
    }
    return fastq_meta
}

/*
========================================================================================
    THE END
========================================================================================
*/