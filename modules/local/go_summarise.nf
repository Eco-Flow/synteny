process GO_SUMMARISE {

    label 'process_single'
    tag "summary of go s"
    container = 'ecoflowucl/chopgo:r-4.3.2_python-3.10'
    publishDir "$params.outdir/GO_results" , mode: "copy"

    input:
    path(go)

    output:
    path( "Go_summary*.tsv" ), emit: go_summary_table
    path( "*.pdf" ), emit: go_summary_pdf
    path "versions.yml", emit: versions

    script:
    """
    Summarise_go.pl
    Rscript "${projectDir}/bin/plotting-synteny_go.R"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R version: \$(R --version | grep "R version" | sed 's/[(].*//' | sed 's/ //g' | sed 's/[^0-9]*//')
        GO version: \$(Rscript -e "as.data.frame(installed.packages())[ ,c(1,3)]" | grep topGO | sed 's/[^0-9]*//')
        Perl version: \$(perl --version | grep "version" | sed 's/.*(//g' | sed 's/[)].*//')
    END_VERSIONS
    """
}
