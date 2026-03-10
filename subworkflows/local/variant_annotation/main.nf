#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
========================================================================================
    ANNOTATION SUBWORKFLOW
========================================================================================
*/

// Include modules
include { BCFTOOLS_NORM } from '../../../modules/local/bcftools/norm/main'
include { SNPEFF } from '../../../modules/local/snpeff/main'
include { VEP } from '../../../modules/local/vep/main'


workflow VARIANT_ANNOTATION {
    take:
    vcf // channel: [ val(meta), path(vcf) ]
    vcf_tbi // channel: [ val(meta), path(tbi) ]
    snpeff_cache // value 
    vep_cache // value
    vep_cache_version // value
    vep_genome // value
    vep_species // value
    fasta // path to reference fasta


    main:
    ch_versions = channel.empty()
    ch_snpeff_cache_db = channel.fromPath(snpeff_cache, checkIfExists: true).collect()
    ch_vep_cache_db = channel.fromPath(vep_cache, checkIfExists: true).collect()

    BCFTOOLS_NORM(
        vcf.join(vcf_tbi),
        fasta
    )
    ch_versions = ch_versions.mix(BCFTOOLS_NORM.out.versions)

    SNPEFF(
        BCFTOOLS_NORM.out.vcf.join(BCFTOOLS_NORM.out.tbi),
        ch_snpeff_cache_db
    )
    ch_versions = ch_versions.mix(SNPEFF.out.versions)

    VEP(
        SNPEFF.out.vcf,
        ch_vep_cache_db,
        vep_cache_version,
        vep_genome,
        vep_species,
        fasta
    )
    ch_versions = ch_versions.mix(VEP.out.versions)

    emit:
    annotated_vcf = VEP.out.vcf // channel: [ val(meta), path(annotated_vcf), path(annotated_tbi) ]
    versions = ch_versions
}
