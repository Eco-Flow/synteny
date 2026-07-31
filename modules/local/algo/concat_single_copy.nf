process CONCAT_SINGLE_COPY {

    label 'process_single'
    tag "concat_single_copy"
    container 'python:3.12'
    publishDir "$params.outdir/algo/species_tree/supermatrix" , mode: "${params.publish_dir_mode}"

    input:
    path(orthofinder_dir)
    path(alignment_dir)

    output:
    path("supermatrix.faa"), emit: alignment
    path("partitions.txt"), emit: partitions
    path "versions.yml", emit: versions

    script:
    def model = params.iqtree_partition_model ?: 'AA'
    """
    python3 ${projectDir}/bin/concat_single_copy.py \\
        -m ${alignment_dir} \\
        -g ${orthofinder_dir}/Orthogroups/Orthogroups.tsv \\
        -o supermatrix.faa \\
        -p partitions.txt \\
        --model '${model}'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Python version: \$(python3 --version | sed 's/[^0-9]*//')
    END_VERSIONS
    """
}
