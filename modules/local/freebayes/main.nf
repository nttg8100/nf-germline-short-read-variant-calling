process FREEBAYES{
    tag "$meta.id"
    label "process_high"
    container "quay.io/biocontainers/freebayes:1.3.10--hbefcdb2_0"
    
    input:
    tuple val(meta), path(bam), path(bai)
    path ref_fasta
    path ref_fai
    
    
    output:
    tuple val(meta), path("${meta.id}.vcf.gz"), path("${meta.id}.vcf.gz.tbi")

    script:
    """
    # call variant
    freebayes-parallel <(fasta_generate_regions.py $ref_fai 100000) $task.cpus -f $ref_fasta $bam > ${meta.id}.vcf
    # index
    bgzip -c ${meta.id}.vcf > ${meta.id}.vcf.gz
    tabix -p vcf ${meta.id}.vcf.gz
    """
}