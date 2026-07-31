//
// Reconstruct ancestral linkage groups (Syngraph) and ancestral gene order (AGORA)
// from BUSCO single-copy orthologues and a species tree, following the core method of
// Maulana et al. 2026 (bioRxiv 2026.07.17.739156)
//

include { BUSCO } from '../../../modules/local/algo/busco.nf'
include { BUSCO_FILTER } from '../../../modules/local/algo/busco_filter.nf'
include { SYNGRAPH } from '../../../modules/local/algo/syngraph.nf'
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

    SYNGRAPH ( filtered_tables, tree.first() )
    ch_versions = ch_versions.mix(SYNGRAPH.out.versions)

    AGORA_PREP ( filtered_tables, tree.first() )
    ch_versions = ch_versions.mix(AGORA_PREP.out.versions)

    AGORA ( AGORA_PREP.out.agora_input )
    ch_versions = ch_versions.mix(AGORA.out.versions)

    FRAGMENTATION_INDEX ( filtered_tables, AGORA.out.ancestral_output )
    ch_versions = ch_versions.mix(FRAGMENTATION_INDEX.out.versions)

    emit:
    rearrangements      = SYNGRAPH.out.rearrangements       // channel: path (algo.rearrangements.tsv)
    ancestral_output    = AGORA.out.ancestral_output         // channel: path (AGORA CARs directory)
    fragmentation_index = FRAGMENTATION_INDEX.out.table      // channel: path (fragmentation_index.tsv)
    versions            = ch_versions
}
