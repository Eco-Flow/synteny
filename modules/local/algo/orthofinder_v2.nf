process ORTHOFINDER_V2 {

    label 'process_high'
    tag "orthofinder_v2"
    // OrthoFinder 2.5.5, an opt-in alternative to the default vendored v3.x module
    // (--orthofinder_v2, see subworkflows/local/species_tree/main.nf) for machines
    // where the v3 biocontainers image doesn't run (confirmed on real hardware --
    // e.g. arm64). Matches the version already relied on elsewhere
    // (Eco-Flow/excon's ORTHOFINDER_V2 module). Only Orthogroups.tsv is used
    // downstream, and that output is unchanged between the two major versions.
    container 'quay.io/biocontainers/orthofinder:2.5.5--hdfd78af_2'

    input:
    path(proteomes, stageAs: 'input/')

    output:
    path("algo_species_tree"), emit: orthofinder
    path("algo_species_tree/Orthogroups/Orthogroups.tsv"), emit: orthologues
    path "versions.yml", emit: versions

    script:
    // Default (-M dendroblast, see --orthofinder_args) skips OrthoFinder's own
    // MSA-based gene-tree/species-tree inference: only Orthogroups.tsv is used below,
    // and dendroblast is both faster and avoids a famsa dependency issue seen with the
    // MSA-based default method.
    """
    mkdir temp_pickle

    orthofinder \\
        ${params.orthofinder_args} \\
        -t ${task.cpus} \\
        -a ${[task.cpus, 4].min()} \\
        -p temp_pickle \\
        -f input \\
        -n algo_species_tree

    mv input/OrthoFinder/Results_algo_species_tree algo_species_tree

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        OrthoFinder version: \$(orthofinder -h | grep -m1 -oE '[0-9]+\\.[0-9]+\\.[0-9]+')
    END_VERSIONS
    """
}
