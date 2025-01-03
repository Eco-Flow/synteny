process SCORE3 {

    label 'process_single'
    tag "${anchors}"
    container = 'quay.io/ecoflowucl/chopgo:r-4.3.2_python-3.10_perl-5.38_pandas' 
    publishDir "$params.outdir/tables/junction_details" , mode: "${params.publish_dir_mode}", pattern:"*junction_details.tsv"
    publishDir "$params.outdir/tables/junction_summary" , mode: "${params.publish_dir_mode}", pattern:"*junction_summary.tsv"

    input:
    path(anchors)
    path(simularity)
    path(gffs)
    path(beds)
    path(last)
    path(unfilteredlast)

    output:
    path("*junction_details.tsv"), emit: junction_types
    path("*_junction_summary.tsv"), emit: junction_summary
    path "versions.yml", emit: versions

    script:
    """
    
    #Refined junction scores:
    grep -Ev 'L\$' ${anchors} > ${anchors}_lifted_removed
    perl ${projectDir}/bin/Junction_focal_classifier.pl ${anchors}_lifted_removed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R version: \$(R --version | grep "R version" | sed 's/[(].*//' | sed 's/ //g' | sed 's/[^0-9]*//')
        Perl version: \$(perl --version | grep "version" | sed 's/.*(//g' | sed 's/[)].*//')
    END_VERSIONS

    """
}
