process SYNGRAPH {

    label 'process_low'
    tag "syngraph"
    container "${params.syngraph_container}"
    publishDir "$params.outdir/algo/syngraph" , mode: "${params.publish_dir_mode}"

    input:
    path(filtered_tables)
    path(tree)

    output:
    path("algo.rearrangements.tsv"), emit: rearrangements
    path("algo.table.tsv"), emit: table
    path("algo*.pickle"), emit: pickles
    path "versions.yml", emit: versions

    script:
    """
    mkdir syngraph_input
    for f in *.filtered.tsv; do
        species=\$(basename \$f .filtered.tsv)
        cut -f1-4 \$f > syngraph_input/\${species}.tsv
    done

    syngraph build -d syngraph_input -o algo
    syngraph infer -g algo.pickle -t ${tree} -r ${params.syngraph_r} -m ${params.syngraph_m} -s ${params.syngraph_reference} -o algo
    syngraph tabulate -g algo.with_ancestors.pickle -o algo

    md5sum "algo.rearrangements.tsv" > "algo.rearrangements.tsv.md5"
    md5sum "algo.table.tsv" > "algo.table.tsv.md5"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Syngraph version: \$(syngraph --version 2>&1 | head -n1 || echo "unknown")
    END_VERSIONS
    """
}
