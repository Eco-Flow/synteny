#!/usr/bin/env nextflow

include { GENOME_ACQUISITION } from '../subworkflows/local/genome_acquisition/main.nf'
include { SPECIES_TREE } from '../subworkflows/local/species_tree/main.nf'
include { ANCESTRAL_RECONSTRUCTION } from '../subworkflows/local/ancestral_reconstruction/main.nf'

def algoErrorMessage() {
    log.info"""
    ==================
    algo entry error
    ==================
    You must supply --input, --busco_lineage and --syngraph_reference, plus exactly
    one of --species_tree (a pre-built Newick tree) or --iqtree_species_tree (build
    one with OrthoFinder + IQ-TREE from the same genomes).
    Example, supplying a tree:
      nextflow run main.nf -profile docker --mode algo \\
        --input data/Example-accession.csv \\
        --species_tree tree.nwk \\
        --busco_lineage coleoptera_odb12 \\
        --syngraph_reference SomeSpecies
    Example, building the tree instead:
      nextflow run main.nf -profile docker --mode algo \\
        --input data/Example-accession.csv \\
        --iqtree_species_tree \\
        --busco_lineage coleoptera_odb12 \\
        --syngraph_reference SomeSpecies
    The pipeline has exited with error status 1.
    """.stripIndent()
    exit 1
}

workflow ALGO {

    log.info """\
     =========================================

     nf-synteny :: ALGO (ancestral linkage groups)

     -----------------------------------------

     SPECIES_TREE (--species_tree, or OrthoFinder + IQ-TREE via --iqtree_species_tree)
           -> ANCESTRAL_RECONSTRUCTION (BUSCO -> Syngraph -> AGORA -> fragmentation index)

     =========================================""".stripIndent()

    def have_tree_input = params.species_tree as boolean
    def build_tree       = params.iqtree_species_tree as boolean

    // Exactly one of the two tree sources is required -- neither, or both, is an error.
    if (!params.input || !params.busco_lineage || !params.syngraph_reference || (have_tree_input == build_tree)) {
        algoErrorMessage()
    }

    ch_versions = Channel.empty()

    in_file = Channel.fromPath(params.input)
    exclude_ch = params.exclude_scaffolds ? Channel.fromPath(params.exclude_scaffolds) : Channel.fromPath("${projectDir}/assets/NO_FILE")

    // Same 2-col (NCBI accession) / 3-col (local path) samplesheet convention, download,
    // and validation as the default pairwise workflow (main.nf) -- shared via this subworkflow
    // rather than duplicated. BUSCO itself only needs the fasta (metaeuk calls genes ab
    // initio), but SPECIES_TREE needs the annotation too, to build proteomes.
    GENOME_ACQUISITION ( in_file )
    ch_versions = ch_versions.mix(GENOME_ACQUISITION.out.versions)

    genome_annotations = GENOME_ACQUISITION.out.genome_annotations
    genome_fastas      = GENOME_ACQUISITION.out.genome_fastas

    // Species tree: either supplied directly, or built from these same genomes' proteomes.
    if (build_tree) {
        SPECIES_TREE ( genome_annotations )
        tree = SPECIES_TREE.out.tree
        ch_versions = ch_versions.mix(SPECIES_TREE.out.versions)
    } else {
        tree = Channel.fromPath(params.species_tree)
    }

    ANCESTRAL_RECONSTRUCTION ( genome_fastas, tree, exclude_ch )
    ch_versions = ch_versions.mix(ANCESTRAL_RECONSTRUCTION.out.versions)

    def wf = workflow
    def outdir = params.outdir
    wf.onComplete {
        println(wf.success ? "\nALGO done! Check results in ${outdir}/algo/ \n" : "Hmmm .. something went wrong\n")
    }
}
