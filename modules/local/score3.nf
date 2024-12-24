process SCORE3 {

    label 'process_single'
    tag "${anchors}"
    container = 'quay.io/ecoflowucl/chopgo:r-4.3.2_python-3.10_perl-5.38_pandas' 
    publishDir "$params.outdir/tables/junction_details" , mode: "${params.publish_dir_mode}", pattern:"*junction_details.tsv"

    input:
    path(anchors)
    path(simularity)
    path(gffs)
    path(beds)
    path(last)
    path(unfilteredlast)

    output:
    path("*junction_details.tsv"), emit: score3_junction_types
    path "versions.yml", emit: versions

    script:
    """
    #Refined junction scores:
    grep -Ev 'L\$' ${anchors} > ${anchors}_lifted_removed
    perl ${projectDir}/bin/Junction_focal_classifier.pl ${anchors}_lifted_removed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R version: \$(R --version | grep "R version" | sed 's/[(].*//' | sed 's/ //g' | sed 's/[^0-9]*//')
        ggplot2 version: \$(Rscript -e "as.data.frame(installed.packages())[ ,c(1,3)]" | grep ggplot2 | sed 's/  */ /g' | cut -f 3 -d " ")
        ggstar version: \$(Rscript -e "as.data.frame(installed.packages())[ ,c(1,3)]" | grep ggstar | sed 's/[^0-9]*//')
        pheatmap version: \$(Rscript -e "as.data.frame(installed.packages())[ ,c(1,3)]" | grep pheatmap | sed 's/[^0-9]*//')
        Perl version: \$(perl --version | grep "version" | sed 's/.*(//g' | sed 's/[)].*//')
    END_VERSIONS


    """
}
