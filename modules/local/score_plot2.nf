process SCORE_PLOTS_2 {

    label 'process_single'
    tag "All genes"
    container = 'quay.io/ecoflowucl/chopgo:r-4.3.2_python-3_perl-5.38'
    publishDir "$params.outdir/figures/synteny_comparisons_2/" , mode: "${params.publish_dir_mode}", pattern:"*-all.pdf"
    publishDir "$params.outdir/figures/synteny_comparisons_2/" , mode: "${params.publish_dir_mode}", pattern:"Chart_of_break_types.pdf"

    input:
    path(trans_inver_summary)
    path(filec)

    output:
    path("*-all.pdf")
    path("Chart_of_break_types.pdf"), emit: pie

    script:
    """
    #First combine all the results into a single table (filec)
    combine_consensus_scores.pl 

    echo test-all.pdf > test-all.pdf
    echo test > Chart_of_break_types.pdf

    #Then merge all the pair of results (averages- from filec- save to filed).
    #As we do not want to plot pairs twice.
    #merge_consensus_scores.pl 

    #Then plot the percent protein identity versus syntenic scores (filed):
    #plot_scores.R

    #Make a pie chart of predicted values 
    #Rscript ${projectDir}/bin/break_types_to_pie.R
    """
}
