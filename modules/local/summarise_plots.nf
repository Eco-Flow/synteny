process SUMMARISE_PLOTS {

   label 'process_single'
   tag "$sample_id"
   container = 'ecoflowucl/chopgo:r-4.3.2_python-3.10_perl-5.38'
   publishDir "$params.outdir/go_results" , mode: "${params.publish_dir_mode}"

   input:
   tuple val(cutoff), path("*")

   output:
   path( "*.pdf" ), emit: go_summary_pdf

   """

   # Run R code:

   Summarise_go_plots.R $cutoff
   
   """

}
