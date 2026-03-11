#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { GERMLINE_VARIANT_CALLING } from './workflows/main.nf'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow {

    //
    // Create input channel from samplesheet or input parameters
    //
    // Read from samplesheet CSV
    channel.fromPath(params.input)
        .splitCsv(header: true)
        .map { row ->
            def meta = [:]
            meta.id = row.sample
            meta.sample = row.sample
            // Support both new format (with lane) and old format (without lane)
            meta.lane = row.lane ?: "L001"
            meta.read_group = "${row.sample}_${meta.lane}"
            def reads = []
            reads.add(file(row.fastq_1))
            if (row.fastq_2) {
                reads.add(file(row.fastq_2))
            }
            return [meta, reads]
        }
        .set { ch_input }


    GERMLINE_VARIANT_CALLING(
        ch_input
    )
}
