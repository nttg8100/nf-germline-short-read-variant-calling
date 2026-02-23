process SAVE_REFERENCE {
    tag "Saving reference files for reuse"
    label 'process_low'
    
    
    input:
    tuple path(fasta), path(fai), path(dict)
    path bwa_index
    tuple path(dbsnp), path(dbsnp_tbi)
    tuple path(known_indels), path(known_indels_tbi)
    
    output:
    path fasta, emit: fasta
    path fai, emit: fai
    path dict, emit: dict
    path bwa_index, emit: bwa_index
    path dbsnp, emit: dbsnp
    path dbsnp_tbi, emit: dbsnp_tbi
    path known_indels, emit: known_indels
    path known_indels_tbi, emit: known_indels_tbi
    
    script:
    """
    # Reference files will be published to output directory
    # Files are already staged by Nextflow and will be published via publishDir
    echo "Reference files ready for publishing:"
    echo "  FASTA: ${fasta}"
    echo "  FAI: ${fai}"
    echo "  DICT: ${dict}"
    echo "  dbSNP: ${dbsnp}"
    echo "  Known Indels: ${known_indels}"
    echo ""
    echo "BWA Index files:"
    for f in ${bwa_index}; do
        echo "  \$f"
    done
    echo ""
    echo "All reference files will be copied to ${params.outdir}/reference/"
    """
}
