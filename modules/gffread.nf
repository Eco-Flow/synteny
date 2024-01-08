process GFFREAD {

    label 'gffread'
    tag "$sample_id"
    container = 'chriswyatt/gffread_python3'
    publishDir "$params.outdir/Gffread_results" , mode: "copy"

    input:
    tuple val(sample_id), path(fasta), path(gff)

    output:
    tuple val(sample_id), path( "${sample_id}.nucl.fa" ), path( "${sample_id}.gff_for_jvci.gff3" ), emit: proteins
    path( "${sample_id}.gff_for_jvci.gff3" ), emit: gff

    script:
    """
    gffread_unzip.pl ${sample_id} ${fasta} ${gff}
    """
}
