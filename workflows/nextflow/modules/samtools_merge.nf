process SAMTOOLS_MERGE {
    tag "$meta.id"
    label 'process_low'
    
    
    container 'quay.io/biocontainers/samtools:1.18--h50ea8bc_1'
    
    input:
    tuple val(meta), path(bams)
    
    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.bai"), emit: bai
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    // If only one BAM, just symlink and index it; otherwise merge
    if (bams instanceof List && bams.size() > 1) {
        """
        # Merge multiple BAM files
        samtools merge \\
            -@ $task.cpus \\
            $args \\
            ${prefix}_merged.bam \\
            ${bams.join(' ')}
        
        # Index the merged BAM
        samtools index \\
            -@ $task.cpus \\
            ${prefix}_merged.bam
        
        # Create versions file
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(samtools --version | head -1 | sed 's/samtools //')
        END_VERSIONS
        """
    } else {
        """
        # Only one BAM file, create symlink and index
        ln -s ${bams[0]} ${prefix}_merged.bam
        
        samtools index \\
            -@ $task.cpus \\
            ${prefix}_merged.bam
        
        # Create versions file
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(samtools --version | head -1 | sed 's/samtools //')
        END_VERSIONS
        """
    }
}
