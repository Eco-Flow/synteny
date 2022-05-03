process SYNTENY {
    label 'syn'
    tag 'syn'
    publishDir "$params.outdir/Jcvi_results/"
    container = 'chriswyatt/jcvi'
             
    input:

    tuple val(sample_id), path(fasta), path(gff), val(sample_id2), path(fasta2), path(gff2)

    output:
        
	path("${sample_id}.${sample_id2}.anchors"), emit: anchors
    path("*.pdf"), emit: pdf

    script:
    """
        python -m jcvi.compara.catalog ortholog ${sample_id} ${sample_id2} --no_strip_names
        python -m jcvi.compara.synteny depth --histogram ${sample_id}.${sample_id2}.anchors
    """
}

