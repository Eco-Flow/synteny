process JCVI {

    label 'jcvi'
    tag "$sample_id"
    container = 'chriswyatt/jcvi'

    input:
    tuple val(sample_id), path(fasta), path(gff)

    output:
    tuple val(sample_id), path( "${sample_id}.cds" ), path( "${sample_id}.bed" ) , emit: new_format
    path( "${sample_id}.bed" ) , emit: beds

    script:
    """
    #Run the basic transformation of gff to bed and fasta to cds conversions.

    python -m jcvi.formats.gff bed --type=mRNA --key=ID ${gff} -o ${sample_id}.bed
    python -m jcvi.formats.fasta format ${fasta} ${sample_id}.cds
    """
}
