process GO {

    label 'process_single'
    tag "$speciessummaries"
    container = 'ecoflowucl/chopgo:r-4.3.2_python-3.10'
    publishDir "$params.outdir/GO_results" , mode: "copy"

    input:
    path(go), stageAs: 'Go'
    path(speciessummaries)
    path(beds)

    output:
    path( "*.pdf" ), emit: go_pdf
    path( "*.tab" ), emit: go_table

    script:
    """
    Synteny_go.pl
    """
}
