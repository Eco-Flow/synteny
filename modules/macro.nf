process MACRO {
    label 'mac'
    tag 'mac'
    publishDir "$params.outdir/Jcvi_results/"
    container = 'chriswyatt/jcvi'
             
    input:
    path(seqids) 
    path(layout)
    path(anchors)
    tuple val(sample_id), path(fasta), path(gff), val(sample_id2), path(fasta2), path(gff2)
    
    output:
    path("*.pdf"), emit: pdf_macro

    script:
    """
        python -m jcvi.compara.synteny screen --minspan=30 --simple ${anchors} ${anchors}.new 
        python -m jcvi.graphics.karyotype ${seqids} ${layout}
    """
}

