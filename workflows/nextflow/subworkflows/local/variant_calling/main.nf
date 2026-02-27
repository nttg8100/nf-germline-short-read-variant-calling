#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
========================================================================================
    VARIANT CALLING SUBWORKFLOW
========================================================================================
    Steps 9-13: HaplotypeCaller, GenotypeGVCFs, Variant Filtering, Merging
========================================================================================
*/

// Include modules
include { GATK_HAPLOTYPECALLER } from '../../../modules/local/gatk/haplotypecaller/main'
include { GATK_GENOTYPEGVCFS } from '../../../modules/local/gatk/genotypegvcfs/main'
include { GATK_SELECTVARIANTS_SNP } from '../../../modules/local/gatk/selectvariants_snp/main'
include { GATK_VARIANTFILTRATION_SNP } from '../../../modules/local/gatk/variantfiltration_snp/main'
include { GATK_SELECTVARIANTS_INDEL } from '../../../modules/local/gatk/selectvariants_indel/main'
include { GATK_VARIANTFILTRATION_INDEL } from '../../../modules/local/gatk/variantfiltration_indel/main'
include { GATK_MERGEVCFS } from '../../../modules/local/gatk/mergevcfs/main'

workflow VARIANT_CALLING {
    
    take:
    bam             // channel: [ val(meta), path(bam) ]
    bai             // channel: [ val(meta), path(bai) ]
    ref_fasta       // value: path(fasta)
    ref_fai         // value: path(fai)
    ref_dict        // value: path(dict)
    dbsnp_vcf       // value: path(vcf)
    dbsnp_tbi       // value: path(tbi)
    
    main:
    ch_versions = Channel.empty()
    
    //
    // STEP 9: Variant Calling with HaplotypeCaller (GVCF mode)
    //
    GATK_HAPLOTYPECALLER (
        bam.join(bai),
        ref_fasta,
        ref_fai,
        ref_dict,
        dbsnp_vcf,
        dbsnp_tbi
    )
    ch_versions = ch_versions.mix(GATK_HAPLOTYPECALLER.out.versions)
    
    //
    // STEP 10: Genotype GVCFs
    //
    GATK_GENOTYPEGVCFS (
        GATK_HAPLOTYPECALLER.out.gvcf.join(GATK_HAPLOTYPECALLER.out.tbi),
        ref_fasta,
        ref_fai,
        ref_dict
    )
    ch_versions = ch_versions.mix(GATK_GENOTYPEGVCFS.out.versions)
    
    //
    // STEP 11: Select and Filter SNPs
    //
    GATK_SELECTVARIANTS_SNP (
        GATK_GENOTYPEGVCFS.out.vcf.join(GATK_GENOTYPEGVCFS.out.tbi),
        ref_fasta,
        ref_fai,
        ref_dict
    )
    ch_versions = ch_versions.mix(GATK_SELECTVARIANTS_SNP.out.versions)
    
    GATK_VARIANTFILTRATION_SNP (
        GATK_SELECTVARIANTS_SNP.out.vcf.join(GATK_SELECTVARIANTS_SNP.out.tbi),
        ref_fasta,
        ref_fai,
        ref_dict
    )
    ch_versions = ch_versions.mix(GATK_VARIANTFILTRATION_SNP.out.versions)
    
    //
    // STEP 12: Select and Filter Indels
    //
    GATK_SELECTVARIANTS_INDEL (
        GATK_GENOTYPEGVCFS.out.vcf.join(GATK_GENOTYPEGVCFS.out.tbi),
        ref_fasta,
        ref_fai,
        ref_dict
    )
    ch_versions = ch_versions.mix(GATK_SELECTVARIANTS_INDEL.out.versions)
    
    GATK_VARIANTFILTRATION_INDEL (
        GATK_SELECTVARIANTS_INDEL.out.vcf.join(GATK_SELECTVARIANTS_INDEL.out.tbi),
        ref_fasta,
        ref_fai,
        ref_dict
    )
    ch_versions = ch_versions.mix(GATK_VARIANTFILTRATION_INDEL.out.versions)
    
    //
    // STEP 13: Merge filtered SNPs and Indels
    //
    GATK_MERGEVCFS (
        GATK_VARIANTFILTRATION_SNP.out.vcf
            .join(GATK_VARIANTFILTRATION_SNP.out.tbi)
            .join(GATK_VARIANTFILTRATION_INDEL.out.vcf)
            .join(GATK_VARIANTFILTRATION_INDEL.out.tbi),
        ref_fasta,
        ref_fai,
        ref_dict
    )
    ch_versions = ch_versions.mix(GATK_MERGEVCFS.out.versions)
    
    emit:
    gvcf     = GATK_HAPLOTYPECALLER.out.gvcf
    vcf      = GATK_MERGEVCFS.out.vcf
    tbi      = GATK_MERGEVCFS.out.tbi
    versions = ch_versions
}
