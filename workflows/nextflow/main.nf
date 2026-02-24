#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
========================================================================================
    GATK Variant Calling Pipeline - Complete 16-Step Nextflow Version
========================================================================================
    Based on the bash workflow - Full migration from Part 1
========================================================================================
*/

// Include modules
include { SAVE_REFERENCE } from './modules/save/reference/main'

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
    reference_ch    // channel: [ path(fasta), path(fai), path(dict) ]
    bwa_index_ch    // channel: Optional [ path(amb), path(ann), ...) ] - if empty, will be generated
    dbsnp_ch        // channel: [ path(vcf), path(tbi) ]
    known_indels_ch // channel: [ path(vcf), path(tbi) ]
    
    main:
    ch_versions = Channel.empty()
    
    // Convert input channels to value channels (collect once, reuse everywhere)
    ref_fasta = reference_ch.map { it[0] }.first()
    ref_fai   = reference_ch.map { it[1] }.first()
    ref_dict  = reference_ch.map { it[2] }.first()
    
    dbsnp_vcf = dbsnp_ch.map { it[0] }.first()
    dbsnp_tbi = dbsnp_ch.map { it[1] }.first()
    
    known_indels_vcf = known_indels_ch.map { it[0] }.first()
    known_indels_tbi = known_indels_ch.map { it[1] }.first()
    
    //
    // SUBWORKFLOW: PREPROCESSING (Steps 1-8)
    // Includes: FASTP, BWA-MEM2, Sorting, Merging, MarkDuplicates, BQSR, Metrics
    //
    PREPROCESSING (
        reads_ch,
        ref_fasta,
        ref_fai,
        ref_dict,
        bwa_index_ch,
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
        params.snpeff_genome ?: 'GRCh38.mane.1.0.refseq'
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
    if (params.input) {
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
    } else {
        // Direct parameters
        def meta = [:]
        meta.id = params.sample ?: 'sample1'
        
        def reads = []
        reads.add(file(params.fastq_r1))
        if (params.fastq_r2) {
            reads.add(file(params.fastq_r2))
        }
        
        ch_input = Channel.of([ meta, reads ])
    }
    
    //
    // Prepare reference files - either use iGenomes S3 paths or local files
    // Nextflow automatically handles S3 file fetching
    //
    if (params.use_igenomes) {
        // Get genome-specific paths from config
        def genome_config = params.genomes[params.genome]
        
        log.info """
        ==========================================
        Using iGenomes Reference Files
        ==========================================
        Genome        : ${params.genome}
        iGenomes Base : ${params.igenomes_base}
        Save Reference: ${params.save_reference}
        Variant Caller : ${params.variant_caller}
        ==========================================
        FASTA         : ${genome_config.fasta}
        BWA Index     : ${genome_config.bwa_index}
        dbSNP         : ${genome_config.dbsnp}
        Known Indels  : ${genome_config.known_indels}
        ==========================================
        """.stripIndent()
        
        //
        // Prepare reference genome channel from S3
        // Nextflow will automatically fetch these files when needed
        //
        reference_ch = Channel.of([
            file(genome_config.fasta),
            file(genome_config.fasta_fai),
            file(genome_config.dict)
        ])
        
        //
        // Prepare BWA-MEM2 index files channel from S3
        // BWA-MEM2 uses different index format: .0123, .bwt.2bit.64, .amb, .ann, .pac, .alt
        // Try to find existing index files, if not found, channel will be empty and index will be generated
        //
        bwa_index_ch = Channel.fromPath("${genome_config.bwa_index}.{0123,amb,ann,pac,bwt.2bit.64,bwt,sa,alt}", checkIfExists: false)
            .collect()
            .ifEmpty { Channel.empty() }
        
        //
        // Prepare known sites channels from S3
        //
        dbsnp_ch = Channel.of([
            file(genome_config.dbsnp),
            file(genome_config.dbsnp_tbi)
        ])
        
        known_indels_ch = Channel.of([
            file(genome_config.known_indels),
            file(genome_config.known_indels_tbi)
        ])
        
        //
        // Optionally save reference files to output directory for reuse
        //
        if (params.save_reference) {
            SAVE_REFERENCE(
                reference_ch,
                bwa_index_ch,
                dbsnp_ch,
                known_indels_ch
            )
        }
        
    } else {
        log.info """
        ==========================================
        Using Local Reference Files
        ==========================================
        Reference     : ${params.reference}
        dbSNP         : ${params.dbsnp}
        Known Indels  : ${params.known_indels}
        ==========================================
        """.stripIndent()
        
        //
        // Prepare reference genome channel
        //
        reference_ch = Channel.of([
            file(params.reference),
            file("${params.reference}.fai"),
            file(params.reference.toString().replace('.fasta', '.dict').replace('.fa', '.dict'))
        ])
        
        //
        // Prepare BWA-MEM2 index files channel
        // Try to find existing index files, if not found, channel will be empty and index will be generated
        //
        bwa_index_ch = Channel.fromPath("${params.reference}.{0123,amb,ann,pac,bwt.2bit.64,bwt,sa,alt}", checkIfExists: false)
            .collect()
            .ifEmpty { Channel.empty() }
        
        //
        // Prepare known sites channels
        //
        dbsnp_ch = Channel.of([
            file(params.dbsnp),
            file("${params.dbsnp}.tbi")
        ])
        
        known_indels_ch = Channel.of([
            file(params.known_indels),
            file("${params.known_indels}.tbi")
        ])
    }
    
    //
    // RUN WORKFLOW
    //
    GATK_VARIANT_CALLING (
        ch_input,
        reference_ch,
        bwa_index_ch,
        dbsnp_ch,
        known_indels_ch
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
