process JCVI {
    label 'jcvi'
    tag 'jcvi'
    //publishDir "$params.outdir/Jcvi_results/"
    container = 'chriswyatt/jcvi'
             
    input:

    tuple val(sample_id), path(fasta), path(gff)

    output:
        
	tuple val(sample_id), path( "${sample_id}.cds" ), path( "${sample_id}.bed" ) , emit: jcvi_result

    script:
    """
	python -m jcvi.formats.gff bed --type=mRNA --key=ID ${gff} -o ${sample_id}.bed
    python -m jcvi.formats.fasta format ${fasta} ${sample_id}.cds
    """
}

