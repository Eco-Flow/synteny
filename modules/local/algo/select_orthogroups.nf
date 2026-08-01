process SELECT_ORTHOGROUPS {

    label 'process_single'
    tag "select_orthogroups"
    container 'python:3.12'

    input:
    path(orthogroup_dir)

    output:
    path("selected_orthogroups"), emit: orthogroups
    path "versions.yml", emit: versions

    script:
    """
    python3 ${projectDir}/bin/select_orthogroups.py \\
        --max ${params.max_orthogroups} \\
        --indir ${orthogroup_dir} \\
        --outdir selected_orthogroups

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Python version: \$(python3 --version | sed 's/[^0-9]*//')
    END_VERSIONS
    """
}
