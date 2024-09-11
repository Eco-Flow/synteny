process GO {

    label 'process_single'
    tag "$speciessummaries ($cutoff percent)"
    container = 'ecoflowucl/chopgo:r-4.3.2_python-3.10_perl-5.38'
    publishDir "$params.outdir/output_data/go_results/individual/all" , mode: "${params.publish_dir_mode}", pattern:"*.tab"
    publishDir "$params.outdir/figures/go_results/individual/all" , mode: "${params.publish_dir_mode}", pattern:"*.pdf"

    input:
    tuple path(go, stageAs: 'Go'), path(speciessummaries), val (cutoff)
    path(beds)

    output:
    path( "*.pdf" ), emit: go_pdf
    path( "*ALL.tab" ), emit: go_pvals
    tuple val(cutoff), path( "*results_ALL.tab" ), emit: go_table
    path "versions.yml", emit: versions

    script:
    """
    #Run GO on a variety of cutoffs, plus a user parameter (default 10, call params.cutoff)
    perl ${projectDir}/bin/Synteny_go.pl ${cutoff}

    #Calculate md5 sums for output
    for tab_file in *ALL.tab; do
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
