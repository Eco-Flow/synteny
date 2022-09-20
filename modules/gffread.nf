process GFFREAD {
    label 'gffread'
    tag "$sample_id"
    container = 'chriswyatt/gffread_python3'
             
    input:

        tuple val(sample_id), path(fasta), path(gff)

    output:

        tuple val(sample_id), path( "${sample_id}.prot.fa" ), path( "${sample_id}.gff_for_jvci.gff3" ), emit: proteins

    script:
    """
	gffread_unzip.pl ${sample_id} ${fasta} ${gff}	
    """
}


    
