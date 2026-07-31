process SUMMARISE_ALG_TABLE {

    label 'process_single'
    tag "summarise_alg_table"
    container 'python:3.12'
    publishDir "$params.outdir/algo/tables" , mode: "${params.publish_dir_mode}"

    input:
    path(table)

    output:
    path("alg_summary.tsv"), emit: alg_summary
    path("alg_status.tsv"), emit: status_summary
    path "versions.yml", emit: versions

    script:
    """
    python3 ${projectDir}/bin/summarise_alg_table.py \\
        --table ${table} \\
        --alg-summary alg_summary.tsv \\
        --status-summary alg_status.tsv

    md5sum "alg_summary.tsv" > "alg_summary.tsv.md5"
    md5sum "alg_status.tsv" > "alg_status.tsv.md5"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Python version: \$(python3 --version | sed 's/[^0-9]*//')
    END_VERSIONS
    """
}
