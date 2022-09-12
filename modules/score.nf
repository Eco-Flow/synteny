process SCORE {
    label 'score'
    tag "$sample_id"
    container = 'chriswyatt/jcvi'

    input:
    path(anchors)

    output:
    path("My_scores.tsv"), emit: score_combine

    script:
    """
    #Run Score for each gene on how close it is to the edge of the syntenic block


    #Run score for genome X in terms of size of syntenic blacks to species Y.

    summarise_anchors.pl 

    """
}
