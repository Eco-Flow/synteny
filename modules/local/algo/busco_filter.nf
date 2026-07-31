process BUSCO_FILTER {

    label 'process_single'
    tag "$sample_id"
    container "${params.busco_container}"
    publishDir "$params.outdir/algo/busco_filtered" , mode: "${params.publish_dir_mode}", pattern: "*.filtered.tsv"

    input:
    tuple val(sample_id), path(full_table)
    path(exclude_scaffolds)

    output:
    tuple val(sample_id), path("${sample_id}.filtered.tsv"), emit: filtered
    path "versions.yml", emit: versions

    script:
    def exclude_arg = exclude_scaffolds.name != 'NO_FILE' ? "--exclude-scaffolds ${exclude_scaffolds}" : ''
    """
    python3 ${projectDir}/bin/busco_filter.py ${sample_id} ${full_table} ${sample_id}.filtered.tsv ${exclude_arg}

    md5sum "${sample_id}.filtered.tsv" > "${sample_id}.filtered.tsv.md5"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Python version: \$(python3 --version | sed 's/[^0-9]*//')
    END_VERSIONS
    """
}
