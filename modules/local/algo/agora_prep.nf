process AGORA_PREP {

    label 'process_single'
    tag "agora_prep"
    container "${params.agora_container}"

    input:
    path(filtered_tables)
    path(tree)

    output:
    path("agora_input"), emit: agora_input
    path "versions.yml", emit: versions

    script:
    """
    tables_args=""
    for f in *.filtered.tsv; do
        species=\$(basename \$f .filtered.tsv)
        tables_args="\$tables_args \${species}=\$f"
    done

    python3 ${projectDir}/bin/busco_to_agora.py --tables \$tables_args --tree ${tree} --outdir agora_input

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Python version: \$(python3 --version | sed 's/[^0-9]*//')
    END_VERSIONS
    """
}
