process CONFIG {

    label 'config'
    tag 'config'
    container = 'dceoy/gffread'

    input:
    path(seqids)
    path(layout)
    tuple val(sample_id), path(fasta1), path(gff), val(sample_id2), path(fasta2), path(gff2)

    output:
    path("seqids_default"), emit: seqids_out
    path("layout_default"), emit: layout_out

    script:
    """
    #Make a default seqids file

    if [ $params.seqids = 'data/default1' ]; then
        grep '>' $fasta1  | head -n 10 | cut -c 2- | awk '{print \$1}' | tr '\n' ',' | sed 's/.\$//' >> seqids_default
        printf "\n" >> seqids_default
        grep '>' $fasta2  | head -n 10 | cut -c 2- | awk '{print \$1}' | tr '\n' ',' | sed 's/.\$//' >> seqids_default
    else
        cat $seqids > seqids_default
    fi

    #Make a default layout file

    if [ $params.layout = 'data/default2' ]; then
        echo "# y, xstart, xend, rotation, color, label, va,  bed\n.6,     .1,    .8,       0,      , ${sample_id}, top, ${sample_id}.bed\n.4,     .1,    .8,       0,      , ${sample_id2}, top, ${sample_id2}.bed\n# edges\ne, 0, 1, ${sample_id}.${sample_id2}.anchors.simple" > layout_default
    else
        cat $layout > layout_default
    fi
    """
}
