process SAMTOOLS_FAIDX {
    tag "$fasta"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'quay.io/biocontainers/samtools:1.21--h50ea8bc_0' }"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta2), path(fai)
    val get_sizes

    output:
    tuple val(meta), path("*.{fa,fasta}"), emit: fa, optional: true
    tuple val(meta), path("*.sizes"),      emit: sizes, optional: true
    tuple val(meta), path("*.fai"),        emit: fai, optional: true
    tuple val(meta), path("*.gzi"),        emit: gzi, optional: true
    path "versions.yml",                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def base = fasta.getBaseName()
    def uncompressed = base.endsWith('.gz') ? base.replaceFirst(/\.gz$/, '') : base
    def get_sizes_command = get_sizes ? "cut -f 1,2 ${uncompressed}.fai > ${uncompressed}.sizes" : ''

    """
    # Check if input is gzipped
    if [[ "$fasta" == *.gz ]]; then
        echo "Input is gzipped. Decompressing..."
        gunzip -c $fasta > ${uncompressed}
    else
        cp $fasta ${uncompressed}
    fi

    samtools faidx ${uncompressed}

    ${get_sizes_command}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1 | awk '{print \$2}')
    END_VERSIONS
    """
    
    stub:
    def uncompressed = fasta.getBaseName().replaceFirst(/\.gz$/, '')
    def get_sizes_command = get_sizes ? "touch ${uncompressed}.sizes" : ''
    """
    touch ${uncompressed}.fai
    if [[ "${fasta.extension}" == "gz" ]]; then
        touch ${uncompressed}.gzi
    fi

    ${get_sizes_command}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: stub
    END_VERSIONS
    """
}

