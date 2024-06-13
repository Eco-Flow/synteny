process GO {

    label 'process_single'
    tag "$speciessummaries"
    container = 'ecoflowucl/chopgo:r-4.3.2_python-3.10_perl-5.38'
    publishDir "$params.outdir/output_data/go_results" , mode: "${params.publish_dir_mode}"

    input:
    tuple path(go, stageAs: 'Go'), path(speciessummaries)
    path(beds)
    val (cutoff)

    output:
    path( "*.pdf" ), emit: go_pdf
    path( "*.tab" ), emit: go_table
    path "versions.yml", emit: versions

    script:
    """
    #Run GO on a variety of cutoffs, plus a user parameter (default 10, call params.cutoff)
 
    perl ${projectDir}/bin/Synteny_go.pl 5
    mkdir Res5
    mv *.pdf Res5
    mv *.tab Res5
    perl ${projectDir}/bin/Synteny_go.pl 15
    mkdir Res15
    mv *.pdf Res15
    mv *.tab Res15
    perl ${projectDir}/bin/Synteny_go.pl 20
    mkdir Res20
    mv *.pdf Res20
    mv *.tab Res20
    perl ${projectDir}/bin/Synteny_go.pl ${cutoff}


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
