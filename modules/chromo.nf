process CHROMOPAINT {
    label 'chromo'
    tag "${sample_id}_VS_${sample_id2}"
    publishDir "$params.outdir/Jcvi_results" , mode: "copy"
    container = 'chriswyatt/jcvi'
             
    input:

        path (hex)
        each (anchors)
        path ('*')
    
    output:
        
        path("*.pdf"), emit: pdf

    script:
    """
        echo '${anchors}' | rev | cut -d'/' -f 1 | rev > Name

        A="\$(cut -d'.' -f1 Name)"
        B="\$(cut -d'.' -f2 Name)"
        anchor.pl \$A.bed \$B.bed ${anchors} 
        python -m jcvi.graphics.chromosome Chromopaint.txt colour.idmap
        mv Chromopaint.pdf "\$A.\$B.chromo.pdf"
    """
}
