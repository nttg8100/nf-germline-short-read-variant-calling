process BCFTOOLS_STATS {
    tag "$meta.id"
    label 'process_low'
    container 'quay.io/biocontainers/bcftools:1.17--haef29d1_0'
    
    input:
    tuple val(meta), path(vcf), path(tbi)
    
    output:
    tuple val(meta), path("*_variant_stats.txt"), emit: stats
    tuple val(meta), path("*_snp_count.txt"), emit: snp_count
    tuple val(meta), path("*_indel_count.txt"), emit: indel_count
    path("versions.yml"), emit: versions
    
    script:
    def prefix = "${meta.id}"
    """
    # Generate bcftools stats
    bcftools stats ${vcf} > ${prefix}_variant_stats.txt
    
    # Count variants by type
    bcftools view -v snps ${vcf} | bcftools query -f '.\\n' | wc -l > ${prefix}_snp_count.txt
    bcftools view -v indels ${vcf} | bcftools query -f '.\\n' | wc -l > ${prefix}_indel_count.txt
    
    # Generate versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/bcftools //')
    END_VERSIONS
    """
}
