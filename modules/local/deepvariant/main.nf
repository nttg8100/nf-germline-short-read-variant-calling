process DEEPVARIANT {
    tag "$meta.id"
    label 'process_high'
    container "docker.io/google/deepvariant:1.9.0"

    input:
    tuple val(meta), path(bam), path(bai)
    path(fasta)
    path(fai)

    output:
    tuple val(meta), path("${meta.id}.vcf.gz")             , emit: vcf
    tuple val(meta), path("${meta.id}.vcf.gz.{tbi,csi}")   , emit: vcf_index
    tuple val(meta), path("${meta.id}.g.vcf.gz")           , emit: gvcf
    tuple val(meta), path("${meta.id}.g.vcf.gz.{tbi,csi}") , emit: gvcf_index
    path "versions.yml"                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    /opt/deepvariant/bin/run_deepvariant \\
        --sample_name ${meta.id} \\
        --model_type=WGS \\
        --ref=${fasta} \\
        --reads=${bam} \\
        --output_vcf=${meta.id}.vcf.gz \\
        --output_gvcf=${meta.id}.g.vcf.gz \\
        --intermediate_results_dir=tmp \\
        --num_shards=${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepvariant: \$(echo \$(/opt/deepvariant/bin/run_deepvariant --version) | sed 's/^.*version //; s/ .*\$//' )
    END_VERSIONS
    """
}