process SCORE_TREE {

    label 'process_low'
    tag "All genes"
    container = 'quay.io/ecoflowucl/chopgo:r-4.3.2_python-3.10_perl-5.38'
    publishDir "$params.outdir/tables" , mode: "${params.publish_dir_mode}", pattern:"My_scores.tsv"
    publishDir "$params.outdir/tables" , mode: "${params.publish_dir_mode}", pattern:"My_sim_cores.tsv"
    publishDir "$params.outdir/tables" , mode: "${params.publish_dir_mode}", pattern:"My_comp_synteny_similarity.tsv"
    publishDir "$params.outdir/figures/synteny_comparisons" , mode: "${params.publish_dir_mode}", pattern:"My_pair_synteny_identity.pdf"
    publishDir "$params.outdir/tables" , mode: "${params.publish_dir_mode}", pattern:"Synteny_matrix.tsv"
    publishDir "$params.outdir/tables/synt_gene_scores" , mode: "${params.publish_dir_mode}", pattern:"*geneScore.tsv"
    publishDir "$params.outdir/tables/synt_species_scores_summaries" , mode: "${params.publish_dir_mode}", pattern:"*SpeciesScoreSummary.txt"
    publishDir "$params.outdir/tables" , mode: "${params.publish_dir_mode}", pattern:"Summary_of_pairwise_comparisons.withDistance.tsv"
    publishDir "$params.outdir/tables" , mode: "${params.publish_dir_mode}", pattern:"tree_distance_matrix.tsv"
    publishDir "$params.outdir/figures/synteny_comparisons" , mode: "${params.publish_dir_mode}", pattern:"*-all_treesort.pdf"

    input:
    path(anchors)
    path(simularity)
    path(gffs)
    path(beds)
    path(last)
    path(unfilteredlast)
    path(tree)


    output:
    path("My_scores.tsv"), emit: score_combine
    path("My_sim_cores.tsv"), emit: simil_combine
    path("My_pair_synteny_identity.pdf"), emit: pairwiseplot
    path("My_comp_synteny_similarity.tsv"), emit: pairdata
    path("Synteny_matrix.tsv"), emit:synmat
    path("*geneScore.tsv"), emit: pairedgenescores
    path("*SpeciesScoreSummary.txt"), emit:speciesSummary
    path("Summary_of_pairwise_comparisons.withDistance.tsv"), emit:trans_inver_summary
    path("species_order"), emit: species_order
    path("tree_distance_matrix.tsv"), emit: tree_distance_matrix
    path "versions.yml", emit: versions

    script:
    """
    # If gff files are compressed, decompress them (useful in testing)
    if [ "\$(ls -A | grep -i \\.*.gff3.gz\$)" ]; then
       for gff in *.gff3.gz; do zcat \$gff > "\${gff%.gz}"; done
    fi

    # Calculate the pairwise, syntenic gene #, max syntenic and average block length.
    perl ${projectDir}/bin/summarise_anchors.pl

    # Calculate average percentage protein identity (pairwise).
    perl ${projectDir}/bin/summarise_similarity.pl

    # Compare syntenic lengths versus protein similarity.
    perl ${projectDir}/bin/syntenicVSsimilarity.pl

    # Calculates gene scores for presence in syntenic block and distance to syntenic break (Genus_species	gene_name_if syntenic	#Distance in genes to break#):
    perl ${projectDir}/bin/Synteny_gene_score.pl
    
    # Summarises gene counts of multiple species to calculate average distance to break:
    perl ${projectDir}/bin/SyntenyScoreSummary.pl

    # Now we run newick tools to get the correct species order, based on the phylogenetic tree in nw format.

    # This is the file called tree from Goatee, used to sort the species order
    nw_labels $tree | grep -v 'N[0-9]' > species_order

    # Calculate tree distances from matrix
    nw_distance -n -m m $tree > tree_distance_matrix.tsv

    # Attempts to summarise all the steps run so far to produce a table (Summary_of_pairwise_comparisons.tsv)
    perl ${projectDir}/bin/Trans_location_Inversion_score_treeSort.pl

    # Adds on the tree distance score:
    perl ${projectDir}/bin/add_tree_distance_measure_to_table.pl Summary_of_pairwise_comparisons.tsv tree_distance_matrix.tsv  > Summary_of_pairwise_comparisons.withDistance.tsv

    # md5sums for nf test
    md5sum My_scores.tsv > My_scores.tsv.md5
    md5sum My_sim_cores.tsv > My_sim_cores.tsv.md5
    md5sum My_comp_synteny_similarity.tsv > My_comp_synteny_similarity.tsv.md5
    md5sum Synteny_matrix.tsv > Synteny_matrix.tsv.md5
    md5sum Summary_of_pairwise_comparisons.withDistance.tsv > Summary_of_pairwise_comparisons.withDistance.tsv.md5

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
