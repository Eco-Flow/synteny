process MACRO {
    label 'mac'
    tag "$sample_id\_$sample_id2"
    publishDir "$params.outdir/Macro_results" , mode: "copy"
    container = 'chriswyatt/jcvi'
             
    input:
    path(seqids) 
    path(layout)
    path(anchors)
    tuple val(sample_id), path(fasta), path(gff), val(sample_id2), path(fasta2), path(gff2)
    
    output:
    path("${sample_id}.${sample_id2}.macro.pdf"), emit: pdf_macro

    script:
        """
        python -m jcvi.compara.synteny screen --minspan=30 --simple ${sample_id}.${sample_id2}.anchors ${anchors}.new 

        python -m jcvi.graphics.karyotype $seqids $layout 

        mv karyotype.pdf ${sample_id}.${sample_id2}.macro.pdf
        """
}
