process VEP {
    tag "${meta.id}"
    label 'process_medium'
    container 'community.wave.seqera.io/library/ensembl-vep_perl-math-cdf:1e13f65f931a6954'

    input:
    tuple val(meta), path(vcf)
    path cache
    val cache_version
    val genome
    val species
    path fasta

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf, optional: true
    tuple val(meta), path("*.vcf.gz.tbi"), emit: tbi, optional: true
    tuple val(meta), path("*.tab.gz"), emit: tab, optional: true
    tuple val(meta), path("*.json.gz"), emit: json, optional: true
    tuple val(meta), val("${task.process}"), val('ensemblvep'), path("*.html"), topic: multiqc_files, emit: report, optional: true
    tuple val("${task.process}"), val('ensemblvep'), eval("vep --help | sed -n '/ensembl-vep/s/.*: //p'"), topic: versions, emit: versions_ensemblvep
    path "versions.yml", emit: versions
   
    script:
    def temp_vep_cache = "temp_vep_cache/${species}/${cache_version}_${genome}"
    """
    mkdir -p ${temp_vep_cache}
    cache_dirname=\$(find -L $cache -name all_vars.gz -print -quit | xargs dirname) 
    mv \$cache_dirname ${temp_vep_cache}

    vep \\
        -i ${vcf} \\
        -o ${meta.id}.vcf.gz \\
        --stats_file ${meta.id}.html --vcf --everything --filter_common --per_gene --total_length --offline --format vcf \\
        --compress_output bgzip \\
        --assembly ${genome} \\
        --species ${species} \\
        --cache \\
        --fasta ${fasta} \\
        --cache_version ${cache_version} \\
        --dir_cache \$PWD/temp_vep_cache \\
        --fork ${task.cpus}

    tabix -p vcf  ${meta.id}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ensembl-vep: \$(vep --version 2>&1 | grep -oP 'ensembl-vep [0-9.]+' | sed 's/ensembl-vep //' || echo "104.3")
    END_VERSIONS
    """
    
}