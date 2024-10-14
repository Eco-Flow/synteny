process SCORE2 {

    label 'process_single'
    tag "${anchors}"
    container = 'ecoflowucl/chopgo:r-4.3.2_python-3.10_perl-5.38' 
    publishDir "$params.outdir/tables" , mode: "${params.publish_dir_mode}", pattern:"Trans_Inversion_junction_merged.txt"
    publishDir "$params.outdir/tables/paired_anchor_change_junction_prediction" , mode: "${params.publish_dir_mode}", pattern:"*Classification_summary.tsv"
    publishDir "$params.outdir/tables/junctionscores/" , mode: "${params.publish_dir_mode}", pattern:"*_gene_scores.txt"

    input:
    path(anchors)
    path(simularity)
    path(gffs)
    path(beds)
    path(last)
    path(unfilteredlast)

    output:
    path("${anchors}.classified"), emit: filec
    path("*Classification_summary.tsv"), emit:classifications
    tuple env(myValue), path("*.translocation_gene_scores.txt"), emit: genetransdistancescores
    tuple env(myValue), path("*.inversion_gene_scores.txt"), emit: geneinverdistancescores
    tuple env(myValue), path("*.other_gene_scores.txt"), emit: geneotherdistancescores
    path "versions.yml", emit: versions

    script:
    """

    #Refined junction scores:
    perl ${projectDir}/bin/Best_synteny_classifier_v6.pl
    perl ${projectDir}/bin/Best_synteny_classifier_v6.classify.pl

    #Calculate gene scores for inversion and translocation junction distance
    perl ${projectDir}/bin/Calculate_distance_to_inver.pl
    perl ${projectDir}/bin/Calculate_distance_to_trans.pl
    perl ${projectDir}/bin/Calculate_distance_to_other.pl

    #Merge the two outputs
    cat Trans_Inversion_junction_count.txt > ${anchors}.classified

    # Extract the species name
    myValue=\$(cat species_tested*)

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
