process LONGEST {

    label 'longest'
    tag "$fasta"
    container = 'chriswyatt/bioseqio'
    publishDir "$params.outdir/" , mode: "copy"
    debug true

    input:
    tuple val(sample_id), path(fasta), path(gff)

    output:
    tuple val (sample_id), path( "${sample_id}.largestIsoform.fa" ), path( "${sample_id}.gff_for_jvci.gff3" ), emit: longest_proteins

    script:
    """
    ncbi_gffread_to_gene_universal.pl ${fasta} ${sample_id}
    """
}
