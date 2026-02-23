process GATK_HAPLOTYPECALLER {
    tag "$meta.id"
    label 'process_high'
    
    
    container 'broadinstitute/gatk:4.6.1.0'
    
    input:
    tuple val(meta), path(bam), path(bai)
    path reference
    path reference_fai
    path reference_dict
    path dbsnp
    path dbsnp_tbi
    
    output:
    tuple val(meta), path("*.g.vcf.gz"), emit: gvcf
    tuple val(meta), path("*.g.vcf.gz.tbi"), emit: tbi
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    # Variant Calling with HaplotypeCaller (GVCF mode)
    gatk HaplotypeCaller \\
        -R $reference \\
        -I $bam \\
        -O ${prefix}.g.vcf.gz \\
        -ERC GVCF \\
        --native-pair-hmm-threads ${task.cpus} \\
        --dbsnp $dbsnp \\
        $args
    
    # Create versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | grep -oP 'The Genome Analysis Toolkit \\(GATK\\) v\\K[0-9.]+' || echo "4.4.0.0")
    END_VERSIONS
    """
}
