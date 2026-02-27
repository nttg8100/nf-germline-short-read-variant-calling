#!/usr/bin/env nextflow

/*
========================================================================================
    PREPROCESSING SUBWORKFLOW
========================================================================================
    Steps 1-8: QC, Trimming, Alignment, Deduplication, BQSR, Quality Metrics
========================================================================================
*/

include { FASTP } from '../../../modules/local/fastp/trim/main'
include { BWAMEM2_INDEX } from '../../../modules/local/bwa/index/main'
include { BWA_MEM2 } from '../../../modules/local/bwa/mem2/main'
include { SAMTOOLS_SORT } from '../../../modules/local/samtools/sort/main'
include { SAMTOOLS_MERGE } from '../../../modules/local/samtools/merge/main'
include { GATKSPARK_MARKDUPLICATES } from '../../../modules/local/gatkspark/markduplicates/main'
include { GATKSPARK_BASERECALIBRATOR } from '../../../modules/local/gatkspark/baserecalibrator/main'
include { GATKSPARK_APPLYBQSR } from '../../../modules/local/gatkspark/applybqsr/main'
include { GATK_COLLECTMETRICS } from '../../../modules/local/gatk/collectmetrics/main'

workflow PREPROCESSING {
    
    take:
    reads_ch        // channel: [ val(meta), [ path(read1), path(read2) ] ]
    ref_fasta       // path: reference FASTA
    ref_fai         // path: reference FAI
    ref_dict        // path: reference dict
    bwa2_index   // channel: Optional BWA index files
    index_bwa2_reference // channel: Optional BWA index files
    dbsnp_vcf       // path: dbSNP VCF
    dbsnp_tbi       // path: dbSNP TBI
    known_indels_vcf // path: known indels VCF
    known_indels_tbi // path: known indels TBI
    
    main:
    ch_versions = Channel.empty()
    
    // Check if BWA index needs to be generated
    // If empty, generate it; otherwise use provided index
    if (index_bwa2_reference){
        BWAMEM2_INDEX(ref_fasta)
        ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions)
        bwa2_index_ch = BWAMEM2_INDEX.out.index
    }
    else{
        bwa2_index_ch = Channel.fromPath(bwa2_index)
    }
    //
    // STEP 1: Adapter Trimming, Quality Filtering, and QC with fastp
    //
    FASTP (
        reads_ch
    )
    ch_versions = ch_versions.mix(FASTP.out.versions)
    
    //
    // STEP 2: Read Alignment with BWA-MEM2
    //
    BWA_MEM2 (
        FASTP.out.reads,
        ref_fasta,
        ref_fai,
        ref_dict,
        bwa2_index_ch
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
    //
    SAMTOOLS_SORT.out.bam
        .map { meta, bam ->
            def new_meta = [:]
            new_meta.id = meta.sample
            new_meta.sample = meta.sample
            return [ new_meta.id, new_meta, bam ]
        }
        .groupTuple()
        .map { sample_id, metas, bams ->
            return [ metas[0], bams ]
        }
        .set { ch_bams_to_merge }
    
    SAMTOOLS_MERGE (
        ch_bams_to_merge
    )
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions)
    
    //
    // STEP 5: Mark Duplicates with GATK Spark
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
    // STEP 8: Collect Alignment Metrics
    //
    GATK_COLLECTMETRICS (
        GATKSPARK_APPLYBQSR.out.bam.join(GATKSPARK_APPLYBQSR.out.bai),
        ref_fasta,
        ref_fai,
        ref_dict
    )
    ch_versions = ch_versions.mix(GATK_COLLECTMETRICS.out.versions)
    
    emit:
    bam                = GATKSPARK_APPLYBQSR.out.bam                  // channel: [ val(meta), path(bam) ]
    bai                = GATKSPARK_APPLYBQSR.out.bai                  // channel: [ val(meta), path(bai) ]
    alignment_summary  = GATK_COLLECTMETRICS.out.alignment_summary    // channel: [ val(meta), path(metrics) ]
    insert_metrics     = GATK_COLLECTMETRICS.out.insert_metrics       // channel: [ val(meta), path(metrics) ]
    versions           = ch_versions                                   // channel: path(versions.yml)
}
