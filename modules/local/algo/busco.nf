process BUSCO {

    label 'process_medium'
    tag "$sample_id"
    container "${params.busco_container}"

    input:
    tuple val(sample_id), path(fasta)

    output:
    tuple val(sample_id), path("${sample_id}.full_table.tsv"), emit: full_table
    path "versions.yml", emit: versions

    script:
    """
    busco -i ${fasta} -o ${sample_id} -l ${params.busco_lineage} -m genome --metaeuk -c ${task.cpus} -f

    cp ${sample_id}/run_${params.busco_lineage}/full_table.tsv ${sample_id}.full_table.tsv

    md5sum "${sample_id}.full_table.tsv" > "${sample_id}.full_table.tsv.md5"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        BUSCO version: \$(busco --version | sed 's/BUSCO //')
    END_VERSIONS
    """
}
