process GATK_VARIANTFILTRATION_INDEL {
    tag "$meta.id"
    label 'process_low'
    
    
    container 'broadinstitute/gatk:4.6.1.0'
    
    input:
    tuple val(meta), path(vcf), path(tbi)
    path reference
    path reference_fai
    path reference_dict
    
    output:
    tuple val(meta), path("*_filtered_indels.vcf.gz"), emit: vcf
    tuple val(meta), path("*_filtered_indels.vcf.gz.tbi"), emit: tbi
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Relaxed filters for test data (~1-2x coverage)
    def filters = task.ext.filters ?: '--filter-expression "QD < 1.0" --filter-name "QD1" --filter-expression "QUAL < 10.0" --filter-name "QUAL10" --filter-expression "FS > 300.0" --filter-name "FS300"'
    
    """
    # Apply hard filters for indels
    gatk VariantFiltration \\
        -R $reference \\
        -V $vcf \\
        $filters \\
        $args \\
        -O ${prefix}_filtered_indels.vcf.gz
    
    # Create versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | grep -oP 'The Genome Analysis Toolkit \\(GATK\\) v\\K[0-9.]+' || echo "4.4.0.0")
    END_VERSIONS
    """
}
