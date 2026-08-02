process DATE_TREE {

    label 'process_single'
    tag "date_tree"
    container 'quay.io/ecoflowucl/cafe:r-4.3.1'
    publishDir "$params.outdir/algo/species_tree" , mode: "${params.publish_dir_mode}"

    input:
    path(tree_newick)
    path(calibrations)

    output:
    path("SpeciesTree_dated.nwk"), emit: dated_tree
    path("dating_calibrations.tsv"), emit: calibrations
    path("dating_qc.tsv"), emit: qc
    path "versions.yml", emit: versions

    script:
    def model     = params.chronos_model ?: 'discrete'
    def lambda    = params.chronos_lambda ?: 1
    def rate_cats = params.chronos_rate_categories ?: 10
    """
    Rscript ${projectDir}/bin/date_tree.R \\
        ${tree_newick} \\
        ${calibrations} \\
        ${model} \\
        ${lambda} \\
        ${rate_cats}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ape version: \$(Rscript -e 'cat(as.character(packageVersion("ape")))')
    END_VERSIONS
    """
}
