process GO_SUMMARISE {

    label 'go_summarise'
    tag "summary of go s"
    container = 'ecoflowucl/chopgo:r-4.3.2_python-3.10'
    publishDir "$params.outdir/GO_results" , mode: "copy"

    input:
    path(go)

    output:
    path( "Go_summary*.tsv" ), emit: go_summary_table
    path( "*.pdf" ), emit: go_summary_pdf

    script:
    """
    Summarise_go.pl
    Rscript "${projectDir}/bin/plotting-synteny_go.R"
    """
}
