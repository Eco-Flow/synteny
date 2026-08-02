//
// Reconstruct ancestral linkage groups (Syngraph) and ancestral gene order (AGORA)
// from BUSCO single-copy orthologues and a species tree, following the core method of
// Maulana et al. 2026 (bioRxiv 2026.07.17.739156)
//

include { BUSCO } from '../../../modules/local/algo/busco.nf'
include { BUSCO_FILTER } from '../../../modules/local/algo/busco_filter.nf'
include { STRIP_TREE_SUPPORT } from '../../../modules/local/algo/strip_tree_support.nf'
include { SYNGRAPH } from '../../../modules/local/algo/syngraph.nf'
include { SUMMARISE_ALG_TABLE } from '../../../modules/local/algo/summarise_alg_table.nf'
include { RESAMPLE_MARKERS } from '../../../modules/local/algo/resample_markers.nf'
include { SYNGRAPH_BOOTSTRAP } from '../../../modules/local/algo/syngraph_bootstrap.nf'
include { SUMMARISE_BOOTSTRAP_SUPPORT } from '../../../modules/local/algo/summarise_bootstrap_support.nf'
include { AGORA_PREP } from '../../../modules/local/algo/agora_prep.nf'
include { AGORA } from '../../../modules/local/algo/agora.nf'
include { FRAGMENTATION_INDEX } from '../../../modules/local/algo/fragmentation_index.nf'

workflow ANCESTRAL_RECONSTRUCTION {

    take:
    genome_fastas     // channel: [ id, fasta ]
    tree              // channel: path (rooted Newick species tree)
    exclude_scaffolds // channel: path (species<TAB>scaffold_name TSV, or the NO_FILE placeholder)

    main:
    ch_versions = Channel.empty()

    BUSCO ( genome_fastas )
    ch_versions = ch_versions.mix(BUSCO.out.versions.first())

    BUSCO_FILTER ( BUSCO.out.full_table, exclude_scaffolds.first() )
    ch_versions = ch_versions.mix(BUSCO_FILTER.out.versions.first())

    filtered_tables = BUSCO_FILTER.out.filtered.map { id, tsv -> tsv }.collect()

    // Syngraph's tree loader (ete3) parses an internal-node label as a single plain
    // support value and rejects IQ-TREE's combined SH-aLRT/UFBoot "100/100" format
    // outright. Syngraph doesn't use support values for anything (only topology and
    // branch lengths), so strip them for Syngraph's own consumption -- AGORA_PREP
    // still gets the original tree (unaffected: it already treats a support-only
    // label the same as a blank one when assigning ancestor names).
    STRIP_TREE_SUPPORT ( tree )
    ch_versions = ch_versions.mix(STRIP_TREE_SUPPORT.out.versions)
    syngraph_tree = STRIP_TREE_SUPPORT.out.tree.first()

    SYNGRAPH ( filtered_tables, syngraph_tree )
    ch_versions = ch_versions.mix(SYNGRAPH.out.versions)

    // Total ALG count per reconstructed ancestral node, plus per-species
    // intact/split/fused status per ALG -- derived from Syngraph's own
    // per-marker table, which it produces but doesn't summarise itself.
    SUMMARISE_ALG_TABLE ( SYNGRAPH.out.table )
    ch_versions = ch_versions.mix(SUMMARISE_ALG_TABLE.out.versions)

    // Bootstrap support for Syngraph's ALG calls (opt-in via --syngraph_bootstraps):
    // rerun build -> infer -> tabulate on N marker-resampled replicates, then, per
    // ancestral node, report what fraction of replicates place each marker back with
    // the same group of markers as the reference (unresampled) run above.
    if (params.syngraph_bootstraps > 0) {
        replicate_ids = Channel.of(1..params.syngraph_bootstraps)

        RESAMPLE_MARKERS ( replicate_ids, filtered_tables )
        ch_versions = ch_versions.mix(RESAMPLE_MARKERS.out.versions.first())

        SYNGRAPH_BOOTSTRAP ( RESAMPLE_MARKERS.out.resampled, syngraph_tree )
        ch_versions = ch_versions.mix(SYNGRAPH_BOOTSTRAP.out.versions.first())

        SUMMARISE_BOOTSTRAP_SUPPORT (
            SYNGRAPH.out.table,
            SYNGRAPH_BOOTSTRAP.out.table.map { replicate, table -> table }.collect()
        )
        ch_versions = ch_versions.mix(SUMMARISE_BOOTSTRAP_SUPPORT.out.versions)

        bootstrap_support = SUMMARISE_BOOTSTRAP_SUPPORT.out.support
    } else {
        bootstrap_support = Channel.empty()
    }

    AGORA_PREP ( filtered_tables, tree.first() )
    ch_versions = ch_versions.mix(AGORA_PREP.out.versions)

    AGORA ( AGORA_PREP.out.agora_input )
    ch_versions = ch_versions.mix(AGORA.out.versions)

    FRAGMENTATION_INDEX ( filtered_tables, AGORA.out.ancestral_output )
    ch_versions = ch_versions.mix(FRAGMENTATION_INDEX.out.versions)

    emit:
    rearrangements      = SYNGRAPH.out.rearrangements                  // channel: path (algo.rearrangements.tsv)
    alg_summary         = SUMMARISE_ALG_TABLE.out.alg_summary          // channel: path (total ALG count per ancestral node)
    alg_status          = SUMMARISE_ALG_TABLE.out.status_summary       // channel: path (per-species/per-ALG intact/split/fused)
    bootstrap_support   = bootstrap_support                            // channel: path (bootstrap_support.tsv), empty unless --syngraph_bootstraps is set
    ancestral_output    = AGORA.out.ancestral_output                   // channel: path (AGORA CARs directory)
    fragmentation_index = FRAGMENTATION_INDEX.out.table                // channel: path (fragmentation_index.tsv)
    versions            = ch_versions
}
