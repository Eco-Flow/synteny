process GO_SUMMARISE {

    label 'go_summarise'
    tag "summary of go s"
    container = 'chriswyatt/chopgo'
    publishDir "$params.outdir/GO_results" , mode: "copy"
    errorStrategy = 'ignore'

    input:
    path(go)

    output:
    path( "Go_summary*.tsv" ), emit: go_summary_table
    path( "*.pdf" ), emit: go_summary_pdf

    script:
    """
    Summarise_go.pl
    plotting-synteny_go.R
    """
}
