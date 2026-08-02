process STRIP_TREE_SUPPORT {

    label 'process_single'
    tag "strip_tree_support"
    container 'python:3.12'

    input:
    path(tree)

    output:
    path("tree.no_support.nwk"), emit: tree
    path "versions.yml", emit: versions

    script:
    """
    python3 ${projectDir}/bin/strip_tree_support.py \\
        -i ${tree} \\
        -o tree.no_support.nwk

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Python version: \$(python3 --version | sed 's/[^0-9]*//')
    END_VERSIONS
    """
}
