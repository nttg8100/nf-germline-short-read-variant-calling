process GATKSPARK_BASERECALIBRATOR {
    tag "$meta.id"
    label 'process_high'
    container 'quay.io/biocontainers/gatk4-spark:4.6.1.0--hdfd78af_0'
    
    input:
    tuple val(meta), path(bam), path(bai)
    path reference
    path reference_fai
    path reference_dict
    path dbsnp
    path dbsnp_tbi
    path known_indels
    path known_indels_tbi
    
    output:
    tuple val(meta), path("*.table"), emit: table
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def avail_mem = (task.memory.mega * 0.8).intValue()
    """
    # Base Quality Score Recalibration - Generate table
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        BaseRecalibratorSpark \\
        -I $bam \\
        -R $reference \\
        --known-sites $dbsnp \\
        --known-sites $known_indels \\
        -O ${meta.id}_recal.table \\
        --spark-master local[${task.cpus}] \\
        --tmp-dir .

    # Create versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
