process SCORE {
    label 'score'
    tag "$sample_id"
    container = 'chriswyatt/perl_r_e1071'
    publishDir "$params.outdir/Summary" , mode: "copy", pattern:"My_scores.tsv"
    publishDir "$params.outdir/Summary" , mode: "copy", pattern:"My_sim_cores.tsv"
    publishDir "$params.outdir/Summary" , mode: "copy", pattern:"My_comp_synteny_similarity.tsv"
    publishDir "$params.outdir/Summary" , mode: "copy", pattern:"My_pair_synteny_identity.pdf"
    publishDir "$params.outdir/Summary" , mode: "copy", pattern:"Synteny_matrix.tsv"
    publishDir "$params.outdir/Synt_gene_scores" , mode: "copy", pattern:"*geneScore.tsv"

    input:
    path(anchors)
    path(simularity)

    output:
    path("My_scores.tsv"), emit: score_combine
    path("My_sim_cores.tsv"), emit: simil_combine
    path("My_pair_synteny_identity.pdf"), emit: pairwiseplot
    path("My_comp_synteny_similarity.tsv"), emit: pairdata
    path("Synteny_matrix.tsv"), emit:synmat
    path("*geneScore.tsv"), emit: pairedgenescores
    path("*SpeciesScoreSummary.txt"), emit:speciesSummary 

    script:
    """
    #Run Score for each gene on how close it is to the edge of the syntenic block


    #Run score for genome X in terms of size of syntenic blacks to species Y.

    summarise_anchors.pl 

    summarise_similarity.pl

    syntenicVSsimilarity.pl

    Synteny_gene_score.pl

    SyntenyScoreSummary.pl 
    """
}
