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
include { SAVE_REFERENCE } from './modules/save_reference'
include { FASTP } from './modules/fastp'
include { BWA_MEM2 } from './modules/bwa_mem2'
include { SAMTOOLS_SORT } from './modules/samtools_sort'
include { SAMTOOLS_MERGE } from './modules/samtools_merge'
include { GATKSPARK_MARKDUPLICATES } from './modules/gatkspark_markduplicates'
include { GATKSPARK_BASERECALIBRATOR } from './modules/gatkspark_baserecalibrator'
include { GATKSPARK_APPLYBQSR } from './modules/gatkspark_applybqsr'
include { GATK_COLLECTMETRICS } from './modules/gatk_collectmetrics'
include { GATK_HAPLOTYPECALLER } from './modules/gatk_haplotypecaller'
include { GATK_GENOTYPEGVCFS } from './modules/gatk_genotypegvcfs'
include { GATK_SELECTVARIANTS_SNP } from './modules/gatk_selectvariants_snp'
include { GATK_VARIANTFILTRATION_SNP } from './modules/gatk_variantfiltration_snp'
include { GATK_SELECTVARIANTS_INDEL } from './modules/gatk_selectvariants_indel'
include { GATK_VARIANTFILTRATION_INDEL } from './modules/gatk_variantfiltration_indel'
include { GATK_MERGEVCFS } from './modules/gatk_mergevcfs'
include { SNPEFF } from './modules/snpeff'
include { BCFTOOLS_STATS } from './modules/bcftools_stats'
include { BCFTOOLS_QUERY } from './modules/bcftools_query'
include { BEDTOOLS_GENOMECOV } from './modules/bedtools_genomecov'

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

workflow GATK_VARIANT_CALLING {
    
    take:
    reads_ch        // channel: [ val(meta), [ path(read1), path(read2) ] ]
    reference_ch    // channel: [ path(fasta), path(fai), path(dict) ]
    bwa_index_ch    // channel: [ path(amb), path(ann), path(bwt), path(pac), path(sa) ]
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
    // STEP 1: Adapter Trimming, Quality Filtering, and QC with fastp
    //
    FASTP (
        reads_ch
    )
    ch_versions = ch_versions.mix(FASTP.out.versions)
    
    //
    // STEP 2: Read Alignment with BWA-MEM
    //
    BWA_MEM2 (
        FASTP.out.reads,
        ref_fasta,
        ref_fai,
        ref_dict,
        bwa_index_ch
    )
    ch_versions = ch_versions.mix(BWA_MEM2.out.versions)
    
    //
    // STEP 3: Sort BAM file
    //
    SAMTOOLS_SORT (
        BWA_MEM2.out.bam
    )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)
    
    //
    // STEP 4: Merge lanes per sample (if multiple lanes exist)
    // Group by sample ID and collect all BAMs for each sample
    //
    SAMTOOLS_SORT.out.bam
        .map { meta, bam ->
            // Create new meta with only sample ID for merging
            def new_meta = [:]
            new_meta.id = meta.sample
            new_meta.sample = meta.sample
            return [ new_meta.id, new_meta, bam ]
        }
        .groupTuple()
        .map { sample_id, metas, bams ->
            // Take the first meta (they should all have same sample info)
            return [ metas[0], bams ]
        }
        .set { ch_bams_to_merge }
    
    SAMTOOLS_MERGE (
        ch_bams_to_merge
    )
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions)
    
    //
    // STEP 5: Mark Duplicates
    //
    GATKSPARK_MARKDUPLICATES (
        SAMTOOLS_MERGE.out.bam
    )
    ch_versions = ch_versions.mix(GATKSPARK_MARKDUPLICATES.out.versions)
    
    //
    // STEP 6: Base Quality Score Recalibration - Generate table
    //
    GATKSPARK_BASERECALIBRATOR (
        GATKSPARK_MARKDUPLICATES.out.bam.join(GATKSPARK_MARKDUPLICATES.out.bai),
        ref_fasta,
        ref_fai,
        ref_dict,
        dbsnp_vcf,
        dbsnp_tbi,
        known_indels_vcf,
        known_indels_tbi
    )
    ch_versions = ch_versions.mix(GATKSPARK_BASERECALIBRATOR.out.versions)
    
    //
    // STEP 7: Apply BQSR
    //
    GATKSPARK_APPLYBQSR (
        GATKSPARK_MARKDUPLICATES.out.bam
            .join(GATKSPARK_MARKDUPLICATES.out.bai)
            .join(GATKSPARK_BASERECALIBRATOR.out.table),
        ref_fasta,
        ref_fai,
        ref_dict
    )
    ch_versions = ch_versions.mix(GATKSPARK_APPLYBQSR.out.versions)
    
    //
    // STEP 8: Alignment Quality Assessment
    //
    GATK_COLLECTMETRICS (
        GATKSPARK_APPLYBQSR.out.bam.join(GATKSPARK_APPLYBQSR.out.bai),
        ref_fasta,
        ref_fai,
        ref_dict
    )
    ch_versions = ch_versions.mix(GATK_COLLECTMETRICS.out.versions)
    
    //
    // STEP 9: Variant Calling with HaplotypeCaller (GVCF mode) or FreeBayes
    //
    GATK_HAPLOTYPECALLER (
        GATKSPARK_APPLYBQSR.out.bam.join(GATKSPARK_APPLYBQSR.out.bai),
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
    
    // Set final VCF channel for downstream processes
    ch_final_vcf = GATK_MERGEVCFS.out.vcf
    ch_final_tbi = GATK_MERGEVCFS.out.tbi
    
    //
    // STEP 14: Functional Annotation with SnpEff
    //
    SNPEFF (
        ch_final_vcf.join(ch_final_tbi),
        params.snpeff_genome ?: 'GRCh38.mane.1.0.refseq'
    )
    ch_versions = ch_versions.mix(SNPEFF.out.versions)
    
    //
    // STEP 15: Variant Statistics with bcftools
    //
    BCFTOOLS_STATS (
        ch_final_vcf.join(ch_final_tbi)
    )
    ch_versions = ch_versions.mix(BCFTOOLS_STATS.out.versions)
    
    //
    // STEP 16: Create Visualization Files
    //
    // 16a: Create BED file from VCF with bcftools
    BCFTOOLS_QUERY (
        ch_final_vcf.join(ch_final_tbi)
    )
    ch_versions = ch_versions.mix(BCFTOOLS_QUERY.out.versions)
    
    // 16b: Generate coverage track with bedtools
    BEDTOOLS_GENOMECOV (
        GATKSPARK_APPLYBQSR.out.bam.join(GATKSPARK_APPLYBQSR.out.bai)
    )
    ch_versions = ch_versions.mix(BEDTOOLS_GENOMECOV.out.versions)
    
    emit:
    fastp_html = FASTP.out.html
    fastp_json = FASTP.out.json
    trimmed_reads = FASTP.out.reads
    final_bam = GATKSPARK_APPLYBQSR.out.bam
    final_vcf = ch_final_vcf
    annotated_vcf = SNPEFF.out.vcf
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
        //
        bwa_index_ch = Channel.fromPath("${genome_config.bwa_index}.{0123,amb,ann,pac,bwt.2bit.64,bwt,sa,alt}").collect()
        
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
        bwa_index_ch = Channel.fromPath("${params.reference}.{0123,amb,ann,pac,bwt.2bit.64,bwt,sa,alt}")
            .collect()
        
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
