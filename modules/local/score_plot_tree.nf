process SCORE_PLOT_TREE {

    label 'process_single'
    tag "All genes"
    container = 'quay.io/ecoflowucl/chopgo:r-4.3.2_python-3_perl-5.38'
    publishDir "$params.outdir/figures/synteny_comparisons" , mode: "${params.publish_dir_mode}", pattern:"*-all_treesort.pdf"
    publishDir "$params.outdir/figures/synteny_comparisons/" , mode: "${params.publish_dir_mode}", pattern:"Chart_of_break_types.pdf"

    input:
    path(trans_inver_summary), stageAs: "Summary_of_pairwise_comparisons.tsv"
    path(filec)
    path(species_order), stageAs: "species_order"

    output:
    path("*-all_treesort.pdf")
    path("Chart_of_break_types.pdf"), emit: pie

    """
    #First combine all the results into a single table (filec)
    combine_consensus_scores.pl 
    
    #Then merge all the pair of results (averages- from filec- save to filed).
    #As we do not want to plot pairs twice.
    merge_consensus_scores.pl 

    #Then plot the percent protein identity versus syntenic scores (filed):
    plot_scores_tree.R

    #Make a pie chart of predicted values 
    Rscript ${projectDir}/bin/break_types_to_pie.R
    """
}
