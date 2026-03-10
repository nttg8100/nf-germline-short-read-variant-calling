process SNPEFF {
    tag "${meta.id}"
    label 'process_medium'
    container 'quay.io/biocontainers/snpeff:5.4.0c--hdfd78af_0'

    input:
    tuple val(meta),  path(vcf), path(vcf_tbi)
    path(cache)

    output:
    tuple val(meta), path("*.ann.vcf"),     emit: vcf
    tuple val(meta), path("*.csv"),         emit: report
    tuple val(meta), path("*.html"),        emit: summary_html
    tuple val(meta), path("*.genes.txt"),   emit: genes_txt
    path "versions.yml",                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def avail_mem = 6144
    if (!task.memory) {
        log.info('[snpEff] Available memory not known - defaulting to 6GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    snpEff \\
        -Xmx${avail_mem}M \\
        ${cache} \\
        -nodownload -canon -v \\
        -csvStats ${prefix}.csv \\
        -dataDir \${PWD}/${cache} \\
        ${vcf} \\
        > ${prefix}.ann.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpeff: \$(snpEff -version 2>&1 | grep -oP 'version [0-9.]+' | sed 's/version //' || echo "5.1")
    END_VERSIONS
    """
}