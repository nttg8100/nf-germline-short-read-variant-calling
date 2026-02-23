process GATK_COLLECTMETRICS {
    tag "$meta.id"
    label 'process_medium'
    container 'broadinstitute/gatk:4.4.0.0'
    
    
    input:
    tuple val(meta), path(bam), path(bai)
    path(reference)
    path(fai)
    path(dict)
    
    output:
    tuple val(meta), path("*_alignment_summary.txt"), emit: alignment_summary
    tuple val(meta), path("*_insert_size_metrics.txt"), emit: insert_metrics, optional: true
    tuple val(meta), path("*_insert_size_histogram.pdf"), emit: insert_histogram, optional: true
    path("versions.yml"), emit: versions
    
    script:
    def prefix = "${meta.id}"
    """
    # CollectAlignmentSummaryMetrics
    gatk CollectAlignmentSummaryMetrics \\
        -R ${reference} \\
        -I ${bam} \\
        -O ${prefix}_alignment_summary.txt
    
    # CollectInsertSizeMetrics with R for histogram generation
    gatk CollectInsertSizeMetrics \\
        -I ${bam} \\
        -O ${prefix}_insert_size_metrics.txt \\
        -H ${prefix}_insert_size_histogram.pdf
    
    # Generate versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | grep -oP 'GATK [0-9.]+' | sed 's/GATK //' || echo "4.4.0.0")
    END_VERSIONS
    """
}
