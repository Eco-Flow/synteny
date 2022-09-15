process SCORE {
    label 'score'
    tag "$sample_id"
    container = 'chriswyatt/jcvi'
    publishDir "$params.outdir/Jcvi_results" , mode: "copy"

    input:
    path(anchors)
    path(simularity)

    output:
    path("My_scores.tsv"), emit: score_combine
    path("My_sim_cores.tsv"), emit: simil_combine

    script:
    """
    #Run Score for each gene on how close it is to the edge of the syntenic block


    #Run score for genome X in terms of size of syntenic blacks to species Y.

    summarise_anchors.pl 

    summarise_simularity.pl

    syntenicVSsimilarity.pl

    """
}
