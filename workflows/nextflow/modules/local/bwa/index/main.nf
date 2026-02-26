process BWAMEM2_INDEX {
    tag "$fasta"
    label 'process_high'
    container "community.wave.seqera.io/library/bwa-mem2_htslib_samtools:e1f420694f8e42bd"
    
    input:
    path fasta
    
    output:
    path "${fasta}.*", emit: index
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    """
    # Index the reference genome with BWA-MEM2
    bwa-mem2 index ${fasta}
    
    # Create versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwamem2: \$(echo \$(bwa-mem2 version 2>&1) | sed 's/.* //')
    END_VERSIONS
    """
    
    stub:
    """
    touch ${fasta}.0123
    touch ${fasta}.amb
    touch ${fasta}.ann
    touch ${fasta}.bwt.2bit.64
    touch ${fasta}.pac
    touch ${fasta}.alt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwamem2: \$(echo \$(bwa-mem2 version 2>&1) | sed 's/.* //')
    END_VERSIONS
    """
}
