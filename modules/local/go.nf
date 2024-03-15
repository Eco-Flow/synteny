process GO {

    label 'process_single'
    tag "$speciessummaries"
    container = 'ecoflowucl/chopgo:r-4.3.2_python-3.10_perl-5.38'
    publishDir "$params.outdir/output_data/go_results" , mode: "${params.publish_dir_mode}"

    input:
    tuple path(go, stageAs: 'Go'), path(speciessummaries)
    path(beds)

    output:
    path( "*.pdf" ), emit: go_pdf
    path( "*.tab" ), emit: go_table
    path "versions.yml", emit: versions

    script:
    """
    perl ${projectDir}/bin/Synteny_go.pl
    for tab_file in *.tab; do
      md5sum \$tab_file > \$tab_file.md5
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R version: \$(R --version | grep "R version" | sed 's/[(].*//' | sed 's/ //g' | sed 's/[^0-9]*//')
        GO version: \$(Rscript -e "as.data.frame(installed.packages())[ ,c(1,3)]" | grep topGO | sed 's/[^0-9]*//')
        Perl version: \$(perl --version | grep "version" | sed 's/.*(//g' | sed 's/[)].*//')
    END_VERSIONS
    """
}
