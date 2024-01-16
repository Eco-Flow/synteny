process SCORE_TREE {

    label 'process_low'
    tag "$sample_id"
    container = 'ecoflowucl/chopgo:r-4.3.2_python-3.10'
    publishDir "$params.outdir/Summary" , mode: "copy", pattern:"My_scores.tsv"
    publishDir "$params.outdir/Summary" , mode: "copy", pattern:"My_sim_cores.tsv"
    publishDir "$params.outdir/Summary" , mode: "copy", pattern:"My_comp_synteny_similarity.tsv"
    publishDir "$params.outdir/Summary" , mode: "copy", pattern:"My_pair_synteny_identity.pdf"
    publishDir "$params.outdir/Summary" , mode: "copy", pattern:"Synteny_matrix.tsv"
    publishDir "$params.outdir/Synt_gene_scores" , mode: "copy", pattern:"*geneScore.tsv"
    publishDir "$params.outdir/Summary" , mode: "copy", pattern:"Trans_location_version.out.txt"
    publishDir "$params.outdir/Summary" , mode: "copy", pattern:"*-all_treesort.pdf"

    input:
    path(anchors)
    path(simularity)
    path(gffs)
    path(tree)

    output:
    path("My_scores.tsv"), emit: score_combine
    path("My_sim_cores.tsv"), emit: simil_combine
    path("My_pair_synteny_identity.pdf"), emit: pairwiseplot
    path("My_comp_synteny_similarity.tsv"), emit: pairdata
    path("Synteny_matrix.tsv"), emit:synmat
    path("*geneScore.tsv"), emit: pairedgenescores
    path("*SpeciesScoreSummary.txt"), emit:speciesSummary
    path("Trans_location_version.out.txt"), emit:trans_inver_summary
    path("*-all_treesort.pdf"), emit:emeline_plots

    script:
    """
    #If gff files are compressed, decompress them (useful in testing)
    for gff in *.gz; do zcat \$gff > \${gff%.gz}; done

    #Run Score for each gene on how close it is to the edge of the syntenic block

    #Run score for genome X in terms of size of syntenic blacks to species Y.

    summarise_anchors.pl

    summarise_similarity.pl

    syntenicVSsimilarity.pl

    Synteny_gene_score.pl

    SyntenyScoreSummary.pl

    #Now we run newick tools to get the correct species order, based on the phylogenetic tree in nw format.

    #nw_prune $tree # This is the file called pruned_tree from Goatee

    nw_labels $tree | grep -v 'N[0-9]' > species_order

    Trans_location_Inversion_score_treeSort.pl

    Rscript "${projectDir}/bin/plotting-inversions-treeSort.R" > R_output.txt

    md5sum My_scores.tsv > My_scores.tsv.md5
    md5sum My_sim_cores.tsv > My_sim_cores.tsv.md5
    md5sum My_comp_synteny_similarity.tsv > My_comp_synteny_similarity.tsv.md5
    md5sum Synteny_matrix.tsv > Synteny_matrix.tsv.md5
    md5sum Trans_location_version.out.txt > Trans_location_version.out.txt.md5
    """
}
