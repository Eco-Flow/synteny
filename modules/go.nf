process GO {
    label 'go'
    tag "$speciessummaries"
    container = 'chriswyatt/chopgo'
    publishDir "$params.outdir/GO_results" , mode: "copy"
    errorStrategy = 'ignore'
             
    input:

	path(go)
	path(speciessummaries)
	path(beds)

    output:

        path( "*.pdf" ), emit: go_pdf
	path( "*.tab" ), emit: go_table

    script:
    """
	#head $speciessummaries > ${speciessummaries}.gohits.pdf
	#head $go > ${speciessummaries}.gotable.tsv
	Synteny_go.pl	
    """
}

