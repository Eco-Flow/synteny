//
// Build a species tree with OrthoFinder + IQ-TREE from a set of annotated genomes,
// following the tree-building steps of
// https://github.com/Eco-Flow/excon/tree/tree_subsampling
//

include { LONGEST } from '../../../modules/local/longest.nf'
include { GFFREAD } from '../../../modules/local/gffread.nf'
include { EXTRACT_PROTEINS } from '../../../modules/local/algo/extract_proteins.nf'
include { ORTHOFINDER } from '../../../modules/nf-core/orthofinder/main'
include { ORTHOFINDER_V2 } from '../../../modules/local/algo/orthofinder_v2.nf'
include { EXTRACT_SINGLE_COPY } from '../../../modules/local/algo/extract_single_copy.nf'
include { SELECT_ORTHOGROUPS } from '../../../modules/local/algo/select_orthogroups.nf'
include { ALIGN_SINGLE_COPY } from '../../../modules/local/algo/align_single_copy.nf'
include { CONCAT_SINGLE_COPY } from '../../../modules/local/algo/concat_single_copy.nf'
include { IQTREE as IQTREE_SPECIES_TREE } from '../../../modules/nf-core/iqtree/main'
include { ROOT_TREE } from '../../../modules/local/algo/root_tree.nf'

workflow SPECIES_TREE {

    take:
    genome_annotations // channel: [ id, fasta, gff ]

    main:
    ch_versions = Channel.empty()

    genome_fastas = genome_annotations.map { id, fasta, gff -> [id, fasta] }

    LONGEST ( genome_annotations )
    ch_versions = ch_versions.mix(LONGEST.out.versions.first())

    GFFREAD ( LONGEST.out.longest_proteins )
    ch_versions = ch_versions.mix(GFFREAD.out.versions.first())

    genome_fastas.join( GFFREAD.out.proteins.map { id, nucl, gff -> [id, gff] } )
        .set { protein_extraction_input }

    EXTRACT_PROTEINS ( protein_extraction_input )
    ch_versions = ch_versions.mix(EXTRACT_PROTEINS.out.versions.first())

    proteomes = EXTRACT_PROTEINS.out.proteins.map { id, fa -> fa }.collect()

    // OrthoFinder is run with its default (fast) gene-tree method: only the
    // orthogroup calls (Orthogroups.tsv) are used below, not OrthoFinder's own
    // species tree or per-orthogroup gene trees.
    //
    // --orthofinder_v2 switches to the vendored v3.x module's OrthoFinder 2.5.5
    // alternative (modules/local/algo/orthofinder_v2.nf) for machines where v3's
    // biocontainers image doesn't run (confirmed on real hardware -- e.g. arm64).
    if (params.orthofinder_v2) {
        ORTHOFINDER_V2 ( proteomes )
        ch_versions = ch_versions.mix(ORTHOFINDER_V2.out.versions)

        orthofinder_dir = ORTHOFINDER_V2.out.orthofinder
    } else {
        ORTHOFINDER (
            proteomes.map { files -> [[id: 'algo_species_tree'], files] },
            [[:], []]
        )
        // Not mixed into ch_versions: ORTHOFINDER.out.versions_orthofinder uses the newer
        // nf-core "topic: versions" convention (a [process, tool, version] tuple), a
        // different shape from the plain path("versions.yml") every other module here
        // emits -- mixing the two shapes into one channel breaks consumers that expect a
        // single type. It's still captured automatically via Channel.topic('versions') if
        // ever needed.

        orthofinder_dir = ORTHOFINDER.out.orthofinder.map { meta, dir -> dir }
    }

    EXTRACT_SINGLE_COPY ( orthofinder_dir, proteomes )
    ch_versions = ch_versions.mix(EXTRACT_SINGLE_COPY.out.versions)

    // Strictly single-copy orthogroups get far more numerous the more closely related
    // the input species are (the "single copy in every species" filter is easier to
    // satisfy), so an uncapped run can land on wildly different IQ-TREE partition counts
    // -- and runtime, since -m MFP model selection runs per partition -- purely as a
    // side effect of how divergent the species happen to be. --max_orthogroups caps
    // that, keeping the longest orthogroups (more sites, more phylogenetic signal).
    // Off by default -- unlimited, matching prior behaviour.
    if (params.max_orthogroups) {
        SELECT_ORTHOGROUPS ( EXTRACT_SINGLE_COPY.out.orthogroups )
        ch_versions = ch_versions.mix(SELECT_ORTHOGROUPS.out.versions)

        orthogroups_for_alignment = SELECT_ORTHOGROUPS.out.orthogroups
    } else {
        orthogroups_for_alignment = EXTRACT_SINGLE_COPY.out.orthogroups
    }

    ALIGN_SINGLE_COPY ( orthogroups_for_alignment )
    ch_versions = ch_versions.mix(ALIGN_SINGLE_COPY.out.versions)

    CONCAT_SINGLE_COPY ( orthofinder_dir, ALIGN_SINGLE_COPY.out.alignments )
    ch_versions = ch_versions.mix(CONCAT_SINGLE_COPY.out.versions)

    // Edge-proportional partition model (-spp), per-partition ModelFinder (-m MFP),
    // 1000 ultrafast bootstrap + 1000 SH-aLRT replicates -- set via the
    // IQTREE_SPECIES_TREE ext.args closure in conf/modules.config.
    IQTREE_SPECIES_TREE (
        CONCAT_SINGLE_COPY.out.alignment.map { aln -> [ [id: 'algo_species_tree'], aln, [] ] },
        [], [], [], [],
        CONCAT_SINGLE_COPY.out.partitions,
        [], [], [], [], [], [], []
    )
    // Not mixed into ch_versions: IQTREE_SPECIES_TREE.out.versions_iqtree uses the
    // newer nf-core "topic: versions" convention (a [process, tool, version] tuple), a
    // different shape from the plain path("versions.yml") every other module here
    // emits -- mixing the two shapes into one channel breaks consumers that expect a
    // single type. It's still captured automatically via Channel.topic('versions') if
    // ever needed.

    // IQ-TREE returns an unrooted tree; ANCESTRAL_RECONSTRUCTION needs a rooted one.
    ROOT_TREE ( IQTREE_SPECIES_TREE.out.phylogeny.map { meta, phylo -> phylo } )
    ch_versions = ch_versions.mix(ROOT_TREE.out.versions)

    emit:
    tree     = ROOT_TREE.out.tree     // channel: path (rooted Newick tree)
    versions = ch_versions
}
