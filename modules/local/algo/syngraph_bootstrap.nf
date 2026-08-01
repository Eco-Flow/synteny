process SYNGRAPH_BOOTSTRAP {

    label 'process_low'
    tag "replicate_${replicate}"
    container 'quay.io/ecoflowucl/syngraph:v1.0'

    input:
    tuple val(replicate), path(resampled_dir)
    path(tree)

    output:
    tuple val(replicate), path("algo.replicate_${replicate}.table.tsv"), emit: table
    path "versions.yml", emit: versions

    script:
    """
    mkdir syngraph_input
    for f in ${resampled_dir}/*.filtered.tsv; do
        species=\$(basename \$f .filtered.tsv)
        cut -f1-4 \$f > syngraph_input/\${species}.tsv
    done

    syngraph build -d syngraph_input -o algo
    syngraph infer -g algo.pickle -t ${tree} -r ${params.syngraph_r} -m ${params.syngraph_m} -s ${params.syngraph_reference} -o algo
    syngraph tabulate -g algo.with_ancestors.pickle -o algo
    mv algo.table.tsv algo.replicate_${replicate}.table.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Syngraph version: \$(syngraph --version 2>&1 | head -n1 || echo "unknown")
    END_VERSIONS
    """
}
