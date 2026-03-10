// Include subworkflows
include { ALIGNMENT } from '../subworkflows/local/alignment/main'
include { PREPROCESSING } from '../subworkflows/local/alignment_preprocessing/main'
include { VARIANT_CALLING } from '../subworkflows/local/variant_calling/main'
include { VARIANT_ANNOTATION } from '../subworkflows/local/variant_annotation/main'
include { VARIANT_ALIGNMENT_QUALITY_CONTROL } from '../subworkflows/local/variant_alignment_quality_control/main'

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

workflow GERMLINE_VARIANT_CALLING {
    take:
    reads_ch // channel: [ val(meta), [ path(read1), path(read2) ] ]

    main:
    ch_versions = channel.empty()

    // Automatically set reference resources from iGenomes if --igenomes_base and --genome are set
    if (params.genome && params.genomes.containsKey(params.genome)) {
        def igenome_ref = params.genomes[params.genome]
        if (igenome_ref) {
            if (igenome_ref.fasta) {
                params.reference = igenome_ref.fasta
            }
            if (igenome_ref.fasta_fai) {
                params.reference_index = igenome_ref.fasta_fai
            }
            if (igenome_ref.dict) {
                params.reference_dict = igenome_ref.dict
            }
            if (igenome_ref.index_bwa2_reference) {
                params.index_bwa2_reference = igenome_ref.index_bwa2_reference
            }
            if (igenome_ref.bwa2_index) {
                params.bwa2_index = igenome_ref.bwa2_index
            }
            if (igenome_ref.dbsnp) {
                params.dbsnp = igenome_ref.dbsnp
            }
            if (igenome_ref.dbsnp_tbi) {
                params.dbsnp_tbi = igenome_ref.dbsnp_tbi
            }
            if (igenome_ref.known_indels) {
                params.known_indels = igenome_ref.known_indels
            }
            if (igenome_ref.known_indels_tbi) {
                params.known_indels_tbi = igenome_ref.known_indels_tbi
            }
            if (igenome_ref.snpeff_cache) {
                params.snpeff_cache = igenome_ref.snpeff_cache
            }
            if (igenome_ref.vep_cache_version) {
                params.vep_cache_version = igenome_ref.vep_cache_version
            }
            if (igenome_ref.vep_genome) {
                params.vep_genome = igenome_ref.vep_genome
            }
            if (igenome_ref.vep_species) {
                params.vep_species = igenome_ref.vep_species
            }
            if (igenome_ref.vep_cache) {
                params.vep_cache = igenome_ref.vep_cache
            }
        }
    }

    log.info(
        """
    ==============================================================================================================================
    nf-germline-short-read-variant-calling:
     - Nextflow Version
     - Workflow: GATK_VARIANT_CALLING
     - Subworkflows: PREPROCESSING, VARIANT_CALLING, ANNOTATION
     - Loaded genomes set: ${params.genome ? params.genome : 'None'}
     - Reference Genome: ${params.reference}
     - dbSNP VCF: ${params.dbsnp}
     - Known Indels VCF: ${params.known_indels}
     - Input Samplesheet: ${params.input}
     - Output Directory: ${params.outdir}
    ==============================================================================================================================
    """.stripIndent()
    )

    //
    // Prepare reference genome channels
    // Values from nextflow.config params block, override via CLI as needed
    ref_fasta = channel.fromPath(params.reference, checkIfExists: true).collect()
    ref_fai = channel.fromPath(params.reference_index, checkIfExists: true).collect()
    ref_dict = channel.fromPath(params.reference_dict, checkIfExists: true).collect()

    // Prepare known sites channels
    // Values from nextflow.config params block, override via CLI as needed
    dbsnp_vcf = channel.fromPath(params.dbsnp, checkIfExists: true).collect()
    dbsnp_tbi = channel.fromPath(params.dbsnp_tbi, checkIfExists: true).collect()
    known_indels_vcf = channel.fromPath(params.known_indels, checkIfExists: true).collect()
    known_indels_tbi = channel.fromPath(params.known_indels_tbi, checkIfExists: true).collect()

    //
    // SUBWORKFLOW: PREPROCESSING (Steps 1-3)
    // Includes: FASTP, BWA-MEM2, Sorting, Merging
    //
    ALIGNMENT(
        reads_ch,
        ref_fasta,
        ref_fai,
        ref_dict,
        params.bwa2_index,
        params.index_bwa2_reference,
    )
    ch_versions = ch_versions.mix(ALIGNMENT.out.versions)
    //
    // SUBWORKFLOW: PREPROCESSING (Steps 4-8)
    // Includes: FASTP, BWA-MEM2, Sorting, Merging, MarkDuplicates, BQSR, Metrics
    //
    if(!params.skip_preprocessing){
        PREPROCESSING(
            params.preprocessor,
            ALIGNMENT.out.bam,
            ref_fasta,
            ref_fai,
            ref_dict,
            dbsnp_vcf,
            dbsnp_tbi,
            known_indels_vcf,
            known_indels_tbi,
        )
        
        ch_versions = ch_versions.mix(PREPROCESSING.out.versions)
        ch_final_bam = PREPROCESSING.out.bam
        ch_final_bai = PREPROCESSING.out.bai
    }

   
    else {
        ch_final_bam = ALIGNMENT.out.bam
        ch_final_bai = ALIGNMENT.out.bai
    }
   
    VARIANT_CALLING(
        params.variant_caller,
        ch_final_bam,
        ch_final_bai,
        ref_fasta,
        ref_fai,
        ref_dict,
        dbsnp_vcf,
        dbsnp_tbi,
    )
    ch_versions = ch_versions.mix(VARIANT_CALLING.out.versions)

    VARIANT_ANNOTATION(
        VARIANT_CALLING.out.vcf,
        VARIANT_CALLING.out.vcf_tbi,
        params.snpeff_cache,
        params.vep_cache,
        params.vep_cache_version,
        params.vep_genome,
        params.vep_species,
        ref_fasta
    )
    // ch_versions = ch_versions.mix(VARIANT_ANNOTATION.out.versions)


    VARIANT_ALIGNMENT_QUALITY_CONTROL(
        VARIANT_CALLING.out.vcf,
        VARIANT_CALLING.out.vcf_tbi,
        ch_final_bam,
        ch_final_bai,
    )
}
