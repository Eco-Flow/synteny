process SUMMARISE_PLOTS_OTHER {
   
   label 'process_single'
   tag "$sample_id"
   container = 'ecoflowucl/chopgo:r-4.3.2_python-3.10_perl-5.38'
   publishDir "$params.outdir/figures/go_results/summarise/other" , mode: "${params.publish_dir_mode}", pattern:"*.pdf"

   input:
   tuple val(cutoff), path("*")

   output:
   path( "*.pdf" ), emit: go_summary_pdf
   path "versions.yml", emit: versions
   
   """
   # Run R code:
   Summarise_go_plots_junction.R $cutoff
   
   cat <<-END_VERSIONS > versions.yml
   "${task.process}":
      R version: \$(R --version | grep "R version" | sed 's/[(].*//' | sed 's/ //g' | sed 's/[^0-9]*//')
   END_VERSIONS
   """

}
