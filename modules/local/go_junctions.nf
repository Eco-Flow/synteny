process GO_JUNCTIONS {

    label 'process_single'
    tag "Run Junction GO ($cutoff percent)"
    container = 'ecoflowucl/chopgo:r-4.3.2_python-3.10_perl-5.38'
    publishDir "$params.outdir/output_data/go_junction_results" , mode: "${params.publish_dir_mode}"

    input:
    tuple path(go, stageAs: 'Go'), path(inversion_distances), val (cutoff)
    path(beds)
    //path(translocation_distances)

    output:
    path( "*.pdf" ), emit: go_pdf
    path( "*ALL.tab" ), emit: go_pvals
    tuple val(cutoff), path( "*results_ALL.tab" ), emit: go_table
    path "versions.yml", emit: versions

    script:
    """
    for file in *gene_scores.txt; do
    # Check if any files match the pattern
    if [ -e "\$file" ]; then
        perl "${projectDir}/bin/Junction_go.pl" ${cutoff} "\$file"
    else
        echo "No files found matching the pattern *gene_scores.txt"
        break
    fi
    done

    #Run GO on junction lists:
    #perl ${projectDir}/bin/Junction_go.pl ${cutoff} ${inversion_distances}

    #Calculate md5 sums for output
    #for tab_file in *ALL.tab; do
    #   md5sum \$tab_file > \$tab_file.md5
    #done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R version: \$(R --version | grep "R version" | sed 's/[(].*//' | sed 's/ //g' | sed 's/[^0-9]*//')
        GO version: \$(Rscript -e "as.data.frame(installed.packages())[ ,c(1,3)]" | grep topGO | sed 's/[^0-9]*//')
        Perl version: \$(perl --version | grep "version" | sed 's/.*(//g' | sed 's/[)].*//')
    END_VERSIONS
    """
}
