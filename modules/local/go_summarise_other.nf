process GO_SUMMARISE_OTHER {

    label 'process_single'
    tag "summary of go s"
    container = 'ecoflowucl/chopgo:r-4.3.2_python-3.10_perl-5.38'
    publishDir "$params.outdir/output_data/go_results/summarise/other" , mode: "${params.publish_dir_mode}", pattern:"*.tsv"

    input:
    tuple val(cutoff), val(filePaths)

    output:
    tuple val(cutoff), path( "Go_summary*.tsv" ), emit: go_summary_table
    path "versions.yml", emit: versions

    script:
    """
    # ${task.process}
    
    # Read in names of files:
    unixfilePaths=\$(echo "$filePaths" )
    echo \$unixfilePaths > files_in

    # Run main summary script:
    perl ${projectDir}/bin/Summarise_go_junctions.pl

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R version: \$(R --version | grep "R version" | sed 's/[(].*//' | sed 's/ //g' | sed 's/[^0-9]*//')
        GO version: \$(Rscript -e "as.data.frame(installed.packages())[ ,c(1,3)]" | grep topGO | sed 's/[^0-9]*//')
        Perl version: \$(perl --version | grep "version" | sed 's/.*(//g' | sed 's/[)].*//')
    END_VERSIONS

    """
    
}
