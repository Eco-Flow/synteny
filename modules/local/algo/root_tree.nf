process ROOT_TREE {

    label 'process_single'
    tag "root_tree"
    container 'python:3.12'
    publishDir "$params.outdir/algo/species_tree" , mode: "${params.publish_dir_mode}"

    input:
    path(tree_newick)

    output:
    path("SpeciesTree_rooted.nwk"), emit: tree
    path "versions.yml", emit: versions

    script:
    def outgroup = params.iqtree_outgroup ? "-g '${params.iqtree_outgroup}'" : ''
    """
    python3 ${projectDir}/bin/root_tree.py \\
        -i ${tree_newick} \\
        -o SpeciesTree_rooted.nwk \\
        ${outgroup}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Python version: \$(python3 --version | sed 's/[^0-9]*//')
    END_VERSIONS
    """
}
