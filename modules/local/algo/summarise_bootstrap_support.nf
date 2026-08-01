process SUMMARISE_BOOTSTRAP_SUPPORT {

    label 'process_single'
    tag "summarise_bootstrap_support"
    container 'python:3.12'
    publishDir "$params.outdir/algo/tables" , mode: "${params.publish_dir_mode}"

    input:
    path(reference_table)
    path(replicate_tables)

    output:
    path("bootstrap_support.tsv"), emit: support
    path "versions.yml", emit: versions

    script:
    """
    python3 ${projectDir}/bin/summarise_bootstrap_support.py \\
        --reference ${reference_table} \\
        --output bootstrap_support.tsv \\
        ${replicate_tables}

    md5sum "bootstrap_support.tsv" > "bootstrap_support.tsv.md5"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Python version: \$(python3 --version | sed 's/[^0-9]*//')
    END_VERSIONS
    """
}
