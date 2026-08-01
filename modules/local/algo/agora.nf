process AGORA {

    label 'process_medium'
    tag "agora"
    container 'quay.io/ecoflowucl/agora:v1.0'
    publishDir "$params.outdir/algo/agora" , mode: "${params.publish_dir_mode}"

    input:
    path(agora_input)

    output:
    path("ancestral_output"), emit: ancestral_output
    path "versions.yml", emit: versions

    script:
    // Confirmed against a real run of the vendored AGORA commit (containers/agora/Dockerfile):
    // orthologyGroups is a "%s"-templated path, one file per ancestor node
    // (doc/HowTo.md), not a single combined file -- AGORA_PREP writes
    // orthologyGroups/orthologyGroups.<ancestor>.list per bin/busco_to_agora.py.
    """
    mkdir -p ancestral_output
    agora-generic.py \\
        ${agora_input}/species_tree.nwk \\
        "${agora_input}/orthologyGroups/orthologyGroups.%s.list" \\
        "${agora_input}/genes/genes.%s.list" \\
        -workingDir=ancestral_output

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        AGORA container: ${task.container}
    END_VERSIONS
    """
}
