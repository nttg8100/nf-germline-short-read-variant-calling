#!/usr/bin/env nextflow

/*
========================================================================================
    PREPROCESSING SUBWORKFLOW
========================================================================================
    Deduplication, BQSR, Quality Metrics
========================================================================================
*/

include { GATKSPARK_MARKDUPLICATES } from '../../../modules/local/gatkspark/markduplicates/main'
include { GATKSPARK_BASERECALIBRATOR } from '../../../modules/local/gatkspark/baserecalibrator/main'
include { GATKSPARK_APPLYBQSR } from '../../../modules/local/gatkspark/applybqsr/main'
include { GATK_COLLECTMETRICS } from '../../../modules/local/gatk/collectmetrics/main'
include { SAMBAMBA_MARKDUP } from '../../../modules/local/sambamba/markdup/main'

workflow PREPROCESSING {
    take:
    preprocessor // value: gatk or sambamba
    aligned_bams // channel: [ val(meta), path(bam) ]
    ref_fasta // path: reference FASTA
    ref_fai // path: reference FAI
    ref_dict // path: reference dict
    dbsnp_vcf // path: dbSNP VCF
    dbsnp_tbi // path: dbSNP TBI
    known_indels_vcf // path: known indels VCF
    known_indels_tbi // path: known indels TBI

    main:
    ch_versions = channel.empty()

    if (preprocessor=="gatk"){
        //
        // STEP 5: Mark Duplicates with GATK Spark
        //
        GATKSPARK_MARKDUPLICATES(
        aligned_bams
        )
        ch_versions = ch_versions.mix(GATKSPARK_MARKDUPLICATES.out.versions)

        //
        // STEP 6: Base Quality Score Recalibration - Generate table
        //
        GATKSPARK_BASERECALIBRATOR(
            GATKSPARK_MARKDUPLICATES.out.bam.join(GATKSPARK_MARKDUPLICATES.out.bai),
            ref_fasta,
            ref_fai,
            ref_dict,
            dbsnp_vcf,
            dbsnp_tbi,
            known_indels_vcf,
            known_indels_tbi,
        )
        ch_versions = ch_versions.mix(GATKSPARK_BASERECALIBRATOR.out.versions)

        //
        // STEP 7: Apply BQSR
        //
        GATKSPARK_APPLYBQSR(
            GATKSPARK_MARKDUPLICATES.out.bam.join(GATKSPARK_MARKDUPLICATES.out.bai).join(GATKSPARK_BASERECALIBRATOR.out.table),
            ref_fasta,
            ref_fai,
            ref_dict,
        )
        ch_versions = ch_versions.mix(GATKSPARK_APPLYBQSR.out.versions)

        //
        // STEP 8: Collect Alignment Metrics
        //
        GATK_COLLECTMETRICS(
            GATKSPARK_APPLYBQSR.out.bam.join(GATKSPARK_APPLYBQSR.out.bai),
            ref_fasta,
            ref_fai,
            ref_dict,
        )
        ch_versions = ch_versions.mix(GATK_COLLECTMETRICS.out.versions)
        ch_processed_bam = GATKSPARK_APPLYBQSR.out.bam
        ch_processed_bai = GATKSPARK_APPLYBQSR.out.bai
        ch_alignment_summary = GATK_COLLECTMETRICS.out.alignment_summary
        ch_alignment_insert_metrics = GATK_COLLECTMETRICS.out.insert_metrics

    } else if (preprocessor=="sambamba"){
        SAMBAMBA_MARKDUP(
            aligned_bams
        )
        ch_versions = ch_versions.mix(SAMBAMBA_MARKDUP.out.versions)
        ch_processed_bam = SAMBAMBA_MARKDUP.out.bam
        ch_processed_bai = SAMBAMBA_MARKDUP.out.bai
        ch_alignment_summary = channel.empty() // no alignment summary metrics with sambamba
        ch_alignment_insert_metrics = channel.empty() // no insert metrics with sambamba
        
    }
    else {
        error "Unsupported preprocessor: ${preprocessor}"
    }
   

    emit:
    bam = ch_processed_bam // channel: [ val(meta), path(bam) ]
    bai = ch_processed_bai // channel: [ val(meta), path(bai) ]
    alignment_summary = ch_alignment_summary // channel: [ val(meta), path(metrics) ]
    insert_metrics = ch_alignment_insert_metrics // channel: [ val(meta), path(metrics) ]
    versions = ch_versions // channel: path(versions.yml)
}
