process FRAGMENTATION_INDEX {

    label 'process_single'
    tag "fragmentation_index"
    container "${params.agora_container}"
    publishDir "$params.outdir/algo/tables" , mode: "${params.publish_dir_mode}"

    input:
    path(filtered_tables)
    path(ancestral_output)

    output:
    path("fragmentation_index.tsv"), emit: table
    path "versions.yml", emit: versions

    script:
    """
    tables_args=""
    for f in *.filtered.tsv; do
        species=\$(basename \$f .filtered.tsv)
        tables_args="\$tables_args \${species}=\$f"
    done

    python3 ${projectDir}/bin/fragmentation_index.py \\
        --tables \$tables_args \\
        --ancgenome-glob "${ancestral_output}/${params.agora_ancgenome_glob}" \\
        --output fragmentation_index.tsv

    md5sum "fragmentation_index.tsv" > "fragmentation_index.tsv.md5"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Python version: \$(python3 --version | sed 's/[^0-9]*//')
    END_VERSIONS
    """
}
