process BEDTOOLS_GENOMECOV {
    tag "$meta.id"
    label 'process_low'
    container 'quay.io/biocontainers/bedtools:2.31.0--hf5e1c6e_2'
    
    
    input:
    tuple val(meta), path(bam), path(bai)
    
    output:
    tuple val(meta), path("*_coverage.bedgraph"), emit: bedgraph
    path("versions.yml"), emit: versions
    
    script:
    def prefix = "${meta.id}"
    """
    # Generate coverage track (bedGraph format)
    bedtools genomecov -ibam ${bam} -bg > ${prefix}_coverage.bedgraph
    
    # Generate versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed 's/bedtools v//')
    END_VERSIONS
    """
}
