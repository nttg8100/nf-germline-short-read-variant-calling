process GATKSPARK_MARKDUPLICATES {
    tag "${meta.id}"
    label 'process_high'
    container 'quay.io/biocontainers/gatk4-spark:4.6.1.0--hdfd78af_0'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}.bam"),     emit: bam
    tuple val(meta), path("${meta.id}.bam.bai"), emit: bai, optional: true
    tuple val(meta), path("*.metrics"),     emit: metrics,   optional: true
    path "versions.yml",                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def avail_mem = (task.memory.mega * 0.8).intValue()
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        MarkDuplicatesSpark \\
        --input ${bam} \\
        --output ${meta.id}.bam \\
        --spark-master local[${task.cpus}] \\
        --tmp-dir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}