process GATK_SELECTVARIANTS_INDEL {
    tag "$meta.id"
    label 'process_low'
    container 'broadinstitute/gatk:4.6.1.0'
    
    input:
    tuple val(meta), path(vcf), path(tbi)
    path reference
    path reference_fai
    path reference_dict
    
    output:
    tuple val(meta), path("*_raw_indels.vcf.gz"), emit: vcf
    tuple val(meta), path("*_raw_indels.vcf.gz.tbi"), emit: tbi
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    # Select Indels
    gatk SelectVariants \\
        -R $reference \\
        -V $vcf \\
        --select-type-to-include INDEL \\
        $args \\
        -O ${prefix}_raw_indels.vcf.gz
    
    # Create versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | grep -oP 'The Genome Analysis Toolkit \\(GATK\\) v\\K[0-9.]+' || echo "4.4.0.0")
    END_VERSIONS
    """
}
