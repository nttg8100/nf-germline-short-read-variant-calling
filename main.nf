#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
========================================================================================
    GATK Variant Calling Pipeline - Complete 16-Step Nextflow Version
========================================================================================
    Based on the bash workflow - Full migration from Part 1
========================================================================================
*/

// Automatically set reference resources from iGenomes if --igenomes_base and --genome are set
if (params.genome && params.genomes.containsKey(params.genome)) {
    def igenome_ref = params.genomes[params.genome]
    if (igenome_ref) {
        if (igenome_ref.fasta   ) params.reference        = igenome_ref.fasta
        if (igenome_ref.fasta_fai) params.reference_index   = igenome_ref.fasta_fai
        if (igenome_ref.dict    ) params.reference_dict    = igenome_ref.dict
        if (igenome_ref.index_bwa2_reference) params.index_bwa2_reference = igenome_ref.index_bwa2_reference
        if (igenome_ref.bwa2_index) params.bwa2_index      = igenome_ref.bwa2_index
        if (igenome_ref.dbsnp   ) params.dbsnp            = igenome_ref.dbsnp
        if (igenome_ref.dbsnp_tbi) params.dbsnp_tbi       = igenome_ref.dbsnp_tbi
        if (igenome_ref.known_indels) params.known_indels  = igenome_ref.known_indels
        if (igenome_ref.known_indels_tbi) params.known_indels_tbi = igenome_ref.known_indels_tbi
        // Add/override more as needed if your config provides more assets
    }
}


// Include subworkflows
include { PREPROCESSING } from './subworkflows/local/preprocessing/main'
include { VARIANT_CALLING } from './subworkflows/local/variant_calling/main'
include { ANNOTATION } from './subworkflows/local/annotation/main'

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

workflow GATK_VARIANT_CALLING {
    
    take:
    reads_ch        // channel: [ val(meta), [ path(read1), path(read2) ] ]
    ref_fasta       // channel: path(fasta)
    ref_fai         // channel: path(fai)
    ref_dict        // channel: path(dict)
    bwa2_index_ch 
    index_bwa2_reference // channel: Optional [ path(amb), path(ann), ...) ] - if empty, will be generated
    dbsnp_vcf      // channel: path(dbsnp vcf)
    dbsnp_tbi      // channel: path(dbsnp tbi)
    known_indels_vcf // channel: path(known indel vcf)
    known_indels_tbi // channel: path(known indel tbi)
    
    main:
    ch_versions = Channel.empty()
    //
    // SUBWORKFLOW: PREPROCESSING (Steps 1-8)
    // Includes: FASTP, BWA-MEM2, Sorting, Merging, MarkDuplicates, BQSR, Metrics
    //
    PREPROCESSING (
        reads_ch,
        ref_fasta,
        ref_fai,
        ref_dict,
        bwa2_index_ch,
        index_bwa2_reference,
        dbsnp_vcf,
        dbsnp_tbi,
        known_indels_vcf,
        known_indels_tbi
    )
    ch_versions = ch_versions.mix(PREPROCESSING.out.versions)
    
    //
    // SUBWORKFLOW: VARIANT_CALLING (Steps 9-13)
    // Includes: HaplotypeCaller, GenotypeGVCFs, Variant Filtering, Merging
    //
    VARIANT_CALLING (
        PREPROCESSING.out.bam,
        PREPROCESSING.out.bai,
        ref_fasta,
        ref_fai,
        ref_dict,
        dbsnp_vcf,
        dbsnp_tbi
    )
    ch_versions = ch_versions.mix(VARIANT_CALLING.out.versions)
    
    //
    // SUBWORKFLOW: ANNOTATION (Steps 14-16)
    // Includes: SnpEff, bcftools stats, Visualization
    //
    ANNOTATION (
        VARIANT_CALLING.out.vcf,
        VARIANT_CALLING.out.tbi,
        PREPROCESSING.out.bam,
        PREPROCESSING.out.bai,
        params.snpeff_genome
    )
    ch_versions = ch_versions.mix(ANNOTATION.out.versions)
    
    emit:
    alignment_summary = PREPROCESSING.out.alignment_summary
    insert_metrics = PREPROCESSING.out.insert_metrics
    final_bam = PREPROCESSING.out.bam
    final_vcf = VARIANT_CALLING.out.vcf
    annotated_vcf = ANNOTATION.out.annotated_vcf
    versions = ch_versions
}

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
    Channel
        .fromPath(params.input)
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
            return [ meta, reads ]
        }
        .set { ch_input }
    
    log.info """
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
    
    //
    // Prepare reference genome channels
    // Values from nextflow.config params block, override via CLI as needed
    ref_fasta_ch = Channel.fromPath(params.reference, checkIfExists: true).collect()
    ref_fai_ch   = Channel.fromPath(params.reference_index, checkIfExists: true).collect()
    ref_dict_ch  = Channel.fromPath(params.reference_dict, checkIfExists: true).collect()

    // Prepare known sites channels
    // Values from nextflow.config params block, override via CLI as needed
    dbsnp_vcf_ch       = Channel.fromPath(params.dbsnp, checkIfExists: true).collect()
    dbsnp_tbi_ch       = Channel.fromPath(params.dbsnp_tbi, checkIfExists: true).collect()
    known_indels_vcf_ch = Channel.fromPath(params.known_indels, checkIfExists: true).collect()
    known_indels_tbi_ch = Channel.fromPath(params.known_indels_tbi, checkIfExists: true).collect()

    // Prepare BWA-MEM2 index files channel
    // Try to find existing index files, if not found, channel will be empty and index will be generated
    //

    //
    // RUN WORKFLOW
    //
    GATK_VARIANT_CALLING (
        ch_input,
        ref_fasta_ch,
        ref_fai_ch,
        ref_dict_ch,
        params.bwa2_index,
        params.index_bwa2_reference,
        dbsnp_vcf_ch,
        dbsnp_tbi_ch,
        known_indels_vcf_ch,
        known_indels_tbi_ch
    )
}

/*
========================================================================================
    COMPLETION SUMMARY
========================================================================================
*/

workflow.onComplete {
    log.info ( workflow.success ? """
        ==========================================
        Pipeline Completed Successfully!
        ==========================================
        Completed at : ${workflow.complete}
        Duration     : ${workflow.duration}
        Success      : ${workflow.success}
        WorkDir      : ${workflow.workDir}
        Exit status  : ${workflow.exitStatus}
        Results      : ${params.outdir}
        ==========================================
        """ : """
        ==========================================
        Pipeline Failed
        ==========================================
        Failed: ${workflow.errorReport}
        Exit status : ${workflow.exitStatus}
        ==========================================
        """
    )
}
