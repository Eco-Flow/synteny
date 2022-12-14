process GO {
    label 'go'
    tag "$sample_id"
    container = 'chriswyatt/chopgo'
    publishDir "$params.outdir/GO_results" , mode: "copy"
             
    input:

	path(go)
	path(speciessummaries)

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

