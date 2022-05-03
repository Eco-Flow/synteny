process GFFREAD {
    label 'gffread'
    tag 'gffread'
    container = 'dceoy/gffread'
             
    input:

    tuple val(sample_id), path(fasta), path(gff)

    output:
        
	tuple val(sample_id), path( "${sample_id}.prot.fa" ), path(gff), emit: proteins

    script:
    """
	gffread -w ${sample_id}.prot.fa -g ${fasta} ${gff}
    """
}

