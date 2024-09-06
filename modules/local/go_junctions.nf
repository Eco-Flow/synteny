process GO_JUNCTIONS {

    label 'process_single'
    tag "Run $junction_score : ($cutoff percent)"
    container = 'ecoflowucl/chopgo:r-4.3.2_python-3.10_perl-5.38'
    publishDir "$params.outdir/figures/go_junction_results" , mode: "${params.publish_dir_mode}"

    input:
    tuple path(go, stageAs: 'Go'), path(junction_score), val (cutoff)
    path(beds)

    output:
    path( "*.pdf" ), emit: go_pdf, optional:true
    path( "*ALL.tab" ), emit: go_pvals, optional:true
    tuple val(cutoff), path( "*results_ALL.tab" ), emit: go_table, optional:true
    path "versions.yml", emit: versions

    script:
    """
    # Check if any files match the pattern
    for file in *gene_scores.txt; do
    if [ -e "\$file" ]; then
        perl "${projectDir}/bin/Junction_go.pl" ${cutoff} "\$file"
    else
        echo "No files found matching the pattern *gene_scores.txt"
        break
    fi
    done

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
