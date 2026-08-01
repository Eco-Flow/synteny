process RESAMPLE_MARKERS {

    label 'process_single'
    tag "replicate_${replicate}"
    container 'python:3.12'

    input:
    val(replicate)
    path(filtered_tables)

    output:
    tuple val(replicate), path("replicate_${replicate}"), emit: resampled
    path "versions.yml", emit: versions

    script:
    """
    mkdir replicate_${replicate}
    python3 ${projectDir}/bin/resample_markers.py \\
        --seed ${replicate} \\
        --outdir replicate_${replicate} \\
        ${filtered_tables}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Python version: \$(python3 --version | sed 's/[^0-9]*//')
    END_VERSIONS
    """
}
