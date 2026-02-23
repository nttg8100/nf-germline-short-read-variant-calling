process GATK_MERGEVCFS {
    tag "$meta.id"
    label 'process_low'
    
    
    container 'broadinstitute/gatk:4.6.1.0'
    
    input:
    tuple val(meta), path(snp_vcf), path(snp_tbi), path(indel_vcf), path(indel_tbi)
    path reference
    path reference_fai
    path reference_dict
    
    output:
    tuple val(meta), path("*_filtered.vcf.gz"), emit: vcf
    tuple val(meta), path("*_filtered.vcf.gz.tbi"), emit: tbi
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    # Merge filtered SNPs and Indels
    gatk MergeVcfs \\
        -I $snp_vcf \\
        -I $indel_vcf \\
        $args \\
        -O ${prefix}_filtered.vcf.gz
    
    # Create versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | grep -oP 'The Genome Analysis Toolkit \\(GATK\\) v\\K[0-9.]+' || echo "4.4.0.0")
    END_VERSIONS
    """
}
