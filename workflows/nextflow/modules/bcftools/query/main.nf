process BCFTOOLS_QUERY {
    tag "$meta.id"
    label 'process_low'
    container 'quay.io/biocontainers/bcftools:1.17--haef29d1_0'
    
    input:
    tuple val(meta), path(vcf), path(tbi)
    
    output:
    tuple val(meta), path("*_variants.bed"), emit: bed
    path("versions.yml"), emit: versions
    
    script:
    def prefix = "${meta.id}"
    """
    # Create BED file from VCF for genome browsers
    bcftools query -f '%CHROM\\t%POS0\\t%END\\t%ID\\n' ${vcf} > ${prefix}_variants.bed
    
    # Generate versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/bcftools //')
    END_VERSIONS
    """
}
