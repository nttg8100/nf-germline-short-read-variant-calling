#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
========================================================================================
    ANNOTATION SUBWORKFLOW
========================================================================================
    Steps 14-16: SnpEff Annotation, Statistics, and Visualization
========================================================================================
*/

// Include modules
include { SNPEFF } from '../../../modules/local/snpeff/annotate/main'
include { BCFTOOLS_STATS } from '../../../modules/local/bcftools/stats/main'
include { BCFTOOLS_QUERY } from '../../../modules/local/bcftools/query/main'
include { BEDTOOLS_GENOMECOV } from '../../../modules/local/bedtools/genomecov/main'

workflow ANNOTATION {
    
    take:
    vcf             // channel: [ val(meta), path(vcf) ]
    tbi             // channel: [ val(meta), path(tbi) ]
    bam             // channel: [ val(meta), path(bam) ]
    bai             // channel: [ val(meta), path(bai) ]
    snpeff_genome   // value: String (e.g., 'GRCh38.mane.1.0.refseq')
    
    main:
    ch_versions = Channel.empty()
    
    //
    // STEP 14: Functional Annotation with SnpEff
    //
    SNPEFF (
        vcf.join(tbi),
        snpeff_genome
    )
    ch_versions = ch_versions.mix(SNPEFF.out.versions)
    
    //
    // STEP 15: Variant Statistics with bcftools
    //
    BCFTOOLS_STATS (
        vcf.join(tbi)
    )
    ch_versions = ch_versions.mix(BCFTOOLS_STATS.out.versions)
    
    //
    // STEP 16: Create Visualization Files
    //
    // 16a: Create BED file from VCF with bcftools
    BCFTOOLS_QUERY (
        vcf.join(tbi)
    )
    ch_versions = ch_versions.mix(BCFTOOLS_QUERY.out.versions)
    
    // 16b: Generate coverage track with bedtools
    BEDTOOLS_GENOMECOV (
        bam.join(bai)
    )
    ch_versions = ch_versions.mix(BEDTOOLS_GENOMECOV.out.versions)
    
    emit:
    annotated_vcf = SNPEFF.out.vcf
    snpeff_log    = SNPEFF.out.log
    stats         = BCFTOOLS_STATS.out.stats
    bed           = BCFTOOLS_QUERY.out.bed
    bedgraph      = BEDTOOLS_GENOMECOV.out.bedgraph
    versions      = ch_versions
}
