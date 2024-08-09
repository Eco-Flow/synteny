process SCORE_PLOTS {

    label 'process_single'
    tag "All genes"
    container = 'ecoflowucl/chopgo:r-4.3.2_python-3.10_perl-5.38'
    publishDir "$params.outdir/figures/synteny_comparisons" , mode: "${params.publish_dir_mode}", pattern:"*.pdf"

    input:
    path(filec), stageAs: "filec"

    output:
    path("*-all.pdf")

    '''
    #First merge all the pair of results (averages- from filec).
    #As we do not want to plot pairs twice.
    merge_consensus_scores.pl 

    #Then plot the percent protein identity versus syntenic scores (filed):
    plot_scores.R
    '''
}
