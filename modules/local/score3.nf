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
    path("*junction_summary.tsv"), emit: junction_summary
    path("versions.yml"), emit: versions
    path("*inver.gene_scores.txt"), optional: true, emit: inver_distancescores
    path("*inter.gene_scores.txt"), optional: true, emit: inter_distancescores
    path("*indel_large.gene_scores.txt"), optional: true, emit: large_indel_distancescores
    path("*indel_tiny.gene_scores.txt"), optional: true, emit: tiny_indel_distancescores
    path("*indel_small.gene_scores.txt"), optional: true, emit: small_indel_distancescores

    script:
    """
    # A script to calculate numbers of different types of break, and the distance from each gene to such breaks.

    #Refined junction scores:
    grep -Ev 'L\$' ${anchors} > ${anchors}_lifted_removed
    perl ${projectDir}/bin/sort_anchors.pl ${anchors}_lifted_removed
    perl ${projectDir}/bin/Junction_focal_classifier.pl ${anchors}_lifted_removed_anchors_sorted.tsv

    #List all genes (in all types of break), on boundary of break:
    perl ${projectDir}/bin/strip_genes_near_breaks.pl *_junction_details.tsv

    #Calculate gene scores for inversion and translocation junction distance
    perl ${projectDir}/bin/Calculate_distance_to_inver.score3.pl
    perl ${projectDir}/bin/Calculate_distance_to_trans.score3.pl
    perl ${projectDir}/bin/Calculate_distance_to_indel_large.score3.pl
    perl ${projectDir}/bin/Calculate_distance_to_indel_small.score3.pl
    perl ${projectDir}/bin/Calculate_distance_to_indel_tiny.score3.pl

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R version: \$(R --version | grep "R version" | sed 's/[(].*//' | sed 's/ //g' | sed 's/[^0-9]*//')
        Perl version: \$(perl --version | grep "version" | sed 's/.*(//g' | sed 's/[)].*//')
    END_VERSIONS

    """
}
