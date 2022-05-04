process MACRO {
    label 'mac'
    tag 'mac'
    publishDir "$params.outdir/Jcvi_results" , mode: "copy"
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
        python -m jcvi.compara.synteny screen --minspan=30 --simple ${anchors} ${anchors}.new 

        #Make a default layout file
        echo "# y, xstart, xend, rotation, color, label, va,  bed\n.6,     .1,    .8,       0,      , ${sample_id}, top, ${sample_id}.bed\n.4,     .1,    .8,       0,      , ${sample_id2}, top, ${sample_id2}.bed\n# edges\ne, 0, 1, ${sample_id}.${sample_id2}.anchors.simple" > layout_default

        #Make a default seqids file
        grep '>' $fasta  | head -n 10 | cut -c 2- | tr '\n' ',' | sed 's/.\$//' >> seqids_default
        grep '>' $fasta2 | head -n 10 | cut -c 2- | tr '\n' ',' | sed 's/.\$//' >> seqids_default        

        #Check which files to run and execute karyotype analysis
        if [ $params.seqids = 'data/default1' ]; then 
            if [ $params.layout = 'data/default2' ]; then
                python -m jcvi.graphics.karyotype seqids_default layout_default
            else 
                python -m jcvi.graphics.karyotype seqids_default $layout 
            fi
        else 
            if [ $params.layout = 'data/default2' ]; then
                python -m jcvi.graphics.karyotype $seqids layout_default
            else 
                python -m jcvi.graphics.karyotype $seqids $layout 
            fi
        fi

        mv karyotype.pdf ${sample_id}.${sample_id2}.macro.pdf
        """
}
