process EXTRACT_SINGLE_COPY {

    label 'process_single'
    tag "extract_single_copy"
    container 'python:3.12'

    input:
    path(orthofinder_dir)
    path(proteomes)

    output:
    path("single_copy_orthogroups"), emit: orthogroups
    path "versions.yml", emit: versions

    script:
    """
    mkdir proteomes
    for f in ${proteomes}; do
        cp "\$f" proteomes/
    done

    python3 ${projectDir}/bin/extract_single_copy.py \\
        -g ${orthofinder_dir}/Orthogroups/Orthogroups.tsv \\
        -p proteomes \\
        -o single_copy_orthogroups

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Python version: \$(python3 --version | sed 's/[^0-9]*//')
    END_VERSIONS
    """
}
