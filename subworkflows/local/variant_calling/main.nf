#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
========================================================================================
    VARIANT CALLING SUBWORKFLOW
========================================================================================
    Steps 9-13: HaplotypeCaller, GenotypeGVCFs, Variant Filtering, Merging
========================================================================================
*/

// Include modules
// GATK-based variant calling
include { GATK_HAPLOTYPECALLER } from '../../../modules/local/gatk/haplotypecaller/main'
include { GATK_GENOTYPEGVCFS } from '../../../modules/local/gatk/genotypegvcfs/main'
include { GATK_SELECTVARIANTS_SNP } from '../../../modules/local/gatk/selectvariants_snp/main'
include { GATK_VARIANTFILTRATION_SNP } from '../../../modules/local/gatk/variantfiltration_snp/main'
include { GATK_SELECTVARIANTS_INDEL } from '../../../modules/local/gatk/selectvariants_indel/main'
include { GATK_VARIANTFILTRATION_INDEL } from '../../../modules/local/gatk/variantfiltration_indel/main'
include { GATK_MERGEVCFS } from '../../../modules/local/gatk/mergevcfs/main'

// FreeBayes-based variant calling
include { FREEBAYES } from '../../../modules/local/freebayes/main'

// DeepVariant-based variant calling
include { DEEPVARIANT } from '../../../modules/local/deepvariant/main'  

workflow VARIANT_CALLING {
    take:
    variant_caller // value: variant calling variant_caller (e.g. "GATK")
    bam // channel: [ val(meta), path(bam) ]
    bai // channel: [ val(meta), path(bai) ]
    ref_fasta // value: path(fasta)
    ref_fai // value: path(fai)
    ref_dict // value: path(dict)
    dbsnp_vcf // value: path(vcf)
    dbsnp_tbi // value: path(tbi)

    main:
    ch_versions = channel.empty()

    if (variant_caller == "gatk"){
         //
        // STEP 9: Variant Calling with HaplotypeCaller (GVCF mode)
        //
        GATK_HAPLOTYPECALLER(
            bam.join(bai),
            ref_fasta,
            ref_fai,
            ref_dict,
            dbsnp_vcf,
            dbsnp_tbi,
        )
        ch_versions = ch_versions.mix(GATK_HAPLOTYPECALLER.out.versions)

        //
        // STEP 10: Genotype GVCFs
        //
        GATK_GENOTYPEGVCFS(
            GATK_HAPLOTYPECALLER.out.gvcf.join(GATK_HAPLOTYPECALLER.out.gvcf_tbi),
            ref_fasta,
            ref_fai,
            ref_dict,
        )
        ch_versions = ch_versions.mix(GATK_GENOTYPEGVCFS.out.versions)

        //
        // STEP 11: Select and Filter SNPs
        //
        GATK_SELECTVARIANTS_SNP(
            GATK_GENOTYPEGVCFS.out.vcf.join(GATK_GENOTYPEGVCFS.out.vcf_tbi),
            ref_fasta,
            ref_fai,
            ref_dict,
        )
        ch_versions = ch_versions.mix(GATK_SELECTVARIANTS_SNP.out.versions)

        GATK_VARIANTFILTRATION_SNP(
            GATK_SELECTVARIANTS_SNP.out.vcf.join(GATK_SELECTVARIANTS_SNP.out.vcf_tbi),
            ref_fasta,
            ref_fai,
            ref_dict,
        )
        ch_versions = ch_versions.mix(GATK_VARIANTFILTRATION_SNP.out.versions)

        //
        // STEP 12: Select and Filter Indels
        //
        GATK_SELECTVARIANTS_INDEL(
            GATK_GENOTYPEGVCFS.out.vcf.join(GATK_GENOTYPEGVCFS.out.vcf_tbi),
            ref_fasta,
            ref_fai,
            ref_dict,
        )
        ch_versions = ch_versions.mix(GATK_SELECTVARIANTS_INDEL.out.versions)

        GATK_VARIANTFILTRATION_INDEL(
            GATK_SELECTVARIANTS_INDEL.out.vcf.join(GATK_SELECTVARIANTS_INDEL.out.vcf_tbi),
            ref_fasta,
            ref_fai,
            ref_dict,
        )
        ch_versions = ch_versions.mix(GATK_VARIANTFILTRATION_INDEL.out.versions)

        //
        // STEP 13: Merge filtered SNPs and Indels
        //
        GATK_MERGEVCFS(
            GATK_VARIANTFILTRATION_SNP.out.vcf
                .join(GATK_VARIANTFILTRATION_SNP.out.vcf_tbi)
                .join(GATK_VARIANTFILTRATION_INDEL.out.vcf)
                .join(GATK_VARIANTFILTRATION_INDEL.out.vcf_tbi),
            ref_fasta,
            ref_fai,
            ref_dict,
        )
        ch_versions = ch_versions.mix(GATK_MERGEVCFS.out.versions)
        ch_out_vcf = GATK_MERGEVCFS.out.vcf
        ch_out_vcf_tbi = GATK_MERGEVCFS.out.vcf_tbi
        ch_out_gvcf = GATK_HAPLOTYPECALLER.out.gvcf
        ch_out_gvcf_tbi = GATK_HAPLOTYPECALLER.out.gvcf_tbi
    }
    else if (variant_caller == "freebayes") {
        // FREEBAYES is included as proof of concept for variant calling with multipler tools, less accurate than GATK and deepvariant
        // Its output will not be processed as gvcf for cohort joint genotyping, but it can be used for single-sample variant calling
        FREEBAYES(
            bam.join(bai),
            ref_fasta,
            ref_fai,
        )
        ch_versions = ch_versions.mix(FREEBAYES.out.versions)
        ch_out_vcf = FREEBAYES.out.vcf
        ch_out_vcf_tbi = FREEBAYES.out.vcf_tbi
        ch_out_gvcf = channel.empty() // FreeBayes does not produce gVCF output
        ch_out_gvcf_tbi = channel.empty()
    }
    else if (variant_caller == "deepvariant") {
        DEEPVARIANT(
            bam.join(bai),
            ref_fasta,
            ref_fai,
        )
        ch_versions = ch_versions.mix(DEEPVARIANT.out.versions)
        ch_out_vcf = DEEPVARIANT.out.vcf
        ch_out_vcf_tbi = DEEPVARIANT.out.vcf_tbi
        ch_out_gvcf = DEEPVARIANT.out.gvcf
        ch_out_gvcf_tbi = DEEPVARIANT.out.gvcf_tbi
    }
    else {
        error "Unsupported variant calling variant_caller: ${variant_caller}"
    }

    emit:
    // VCF Channel - independent outputs
    vcf = ch_out_vcf          // channel: [ val(meta), path(vcf.gz) ]
    vcf_tbi = ch_out_vcf_tbi  // channel: [ val(meta), path(vcf.gz.tbi) ]
    
    // GVCF Channel - independent outputs
    gvcf = ch_out_gvcf            // channel: [ val(meta), path(gvcf.gz) ] (null if not applicable)
    gvcf_tbi = ch_out_gvcf_tbi    // channel: [ val(meta), path(gvcf.gz.tbi) ] (null if not applicable)
    versions = ch_versions
}
