/*
========================================================================================
    ALIGNMENT SUBWORKFLOW
========================================================================================
    QC, Trimming, Alignment
========================================================================================
*/


include { FASTP } from '../../../modules/local/fastp/trim/main'
include { BWAMEM2_INDEX } from '../../../modules/local/bwa/index/main'
include { BWA_MEM2 } from '../../../modules/local/bwa/mem2/main'
include { SAMTOOLS_SORT } from '../../../modules/local/samtools/sort/main'
include { SAMTOOLS_MERGE } from '../../../modules/local/samtools/merge/main'


workflow ALIGNMENT {
    take:
    reads_ch // channel: [ val(meta), [ path(read1), path(read2) ] ]
    ref_fasta // path: reference FASTA
    ref_fai // path: reference FAI
    ref_dict // path: reference dict
    bwa2_index // channel: Optional BWA index files
    index_bwa2_reference // channel: Optional BWA index files

    main:
    ch_versions = channel.empty()

    // Check if BWA index needs to be generated
    // If empty, generate it; otherwise use provided index
    if (index_bwa2_reference) {
        BWAMEM2_INDEX(ref_fasta)
        ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions)
        bwa2_index_ch = BWAMEM2_INDEX.out.index
    }
    else {
        bwa2_index_ch = channel.fromPath(bwa2_index)
    }
    //
    // STEP 1: Adapter Trimming, Quality Filtering, and QC with fastp
    //
    FASTP(
        reads_ch
    )
    ch_versions = ch_versions.mix(FASTP.out.versions)

    //
    // STEP 2: Read Alignment with BWA-MEM2
    //
    BWA_MEM2(
        FASTP.out.reads,
        ref_fasta,
        ref_fai,
        ref_dict,
        bwa2_index_ch,
    )
    ch_versions = ch_versions.mix(BWA_MEM2.out.versions)

    //
    // STEP 3: Sort BAM file
    //
    SAMTOOLS_SORT(
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
            return [new_meta.id, new_meta, bam]
        }
        .groupTuple()
        .map { _sample_id, metas, bams ->
            return [metas[0], bams]
        }
        .set { ch_bams_to_merge }

    SAMTOOLS_MERGE(
        ch_bams_to_merge
    )
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions)


    emit:
    bam = SAMTOOLS_MERGE.out.bam // channel: [ val(meta), path(bam) ]
    bai = SAMTOOLS_MERGE.out.bai // channel: [ val(meta), path(bai) ]
    versions = ch_versions // channel: path(versions.yml)
}