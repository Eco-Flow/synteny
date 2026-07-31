process ALIGN_SINGLE_COPY {

    label 'process_medium'
    tag "align_single_copy"
    container 'quay.io/biocontainers/mulled-v2-12eba4a074f913c639117640936668f5a6a01da6:425707898cf4f85051b77848be253b88f1d2298a-0'

    input:
    path(orthogroup_dir)

    output:
    path("single_copy_alignments"), emit: alignments
    path "versions.yml", emit: versions

    script:
    def args = params.mafft_args ?: '--auto'
    """
    mkdir single_copy_alignments

    find -L ${orthogroup_dir} -name '*.fa' | sort > single_copy_files.txt

    if [ ! -s single_copy_files.txt ]; then
        echo "ERROR: no single-copy orthogroup sequences found in ${orthogroup_dir}" >&2
        exit 1
    fi

    xargs -a single_copy_files.txt -P ${task.cpus} -I {} \\
        sh -c 'mafft ${args} --quiet --anysymbol "\$1" > "single_copy_alignments/\$(basename "\$1")"' _ {}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mafft version: \$(mafft --version 2>&1 | sed 's/^v//;s/ .*//')
    END_VERSIONS
    """
}
