#!/usr/bin/env nextflow

log.info """\
 =========================================

 nf-synteny (v4.0.0)

 -----------------------------------------

 Authors:
   - Chris Wyatt <c.wyatt@ucl.ac.uk>
   - Simon Murray

 -----------------------------------------

 Copyright (c) 2022

 =========================================""".stripIndent()


def errorMessage() {
    log.info"""
    =============
    synteny error
    =============
    You failed to provide the input parameter
    Please provide this as follows:
      --input /full/path/to/sample/file
    or use the test profile:
      --profile test
    The pipeline has exited with error status 1.
    """.stripIndent()
    exit 1
}

include { DOWNLOAD_NCBI } from './modules/local/download_ncbi.nf'
include { GFFREAD } from './modules/local/gffread.nf'
include { LONGEST } from './modules/local/longest.nf'
include { JCVI } from './modules/local/jcvi.nf'
include { SYNTENY } from './modules/local/synteny.nf'
include { CHROMOPAINT } from './modules/local/chromo.nf'
include { SCORE } from './modules/local/score.nf'
include { SCORE_TREE } from './modules/local/score_tree.nf'
include { GO } from './modules/local/go.nf'
include { GO_JUNCTIONS } from './modules/local/go_junctions.nf'
include { GO_SUMMARISE } from './modules/local/go_summarise.nf'
include { FASTAVALIDATOR } from './modules/nf-core/fastavalidator/main'
include { SEQKIT_STATS } from './modules/nf-core/seqkit/stats/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from './modules/nf-core/custom/dumpsoftwareversions'
include { validateParameters; paramsHelp; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'
include { SCORE_PLOTS } from './modules/local/plot_score.nf'
include { SCORE_TREE_PLOTS } from './modules/local/plot_tree.nf'
include { SUMMARISE_PLOTS } from './modules/local/summarise_plots.nf'
include { RIBBON } from './modules/local/ribbon.nf'

// Set colours for figures.
Channel
    .fromPath(params.hex)
    .set { in_hex }

// Add input params cutoff for go
our_cutoff = Channel.from(params.cutoff)
                    .splitCsv(header: false, sep: ",")
                    .flatten()

def flattenCutoffGroups(groups) {
    return groups.collect { tuple ->
        def (key, nestedList) = tuple
        def flattenedList = nestedList.collectMany { it }  // collectMany flattens one level
        [key, flattenedList]
    }
}

// Caluclate buffer size, to split GO summarise results correctly.
import java.util.concurrent.atomic.AtomicInteger
import java.nio.file.Files
import java.nio.file.Paths
import groovy.io.FileType

workflow {

    // Print help message, supply typical command line usage for the pipeline
    if (params.help) {
       log.info paramsHelp("nextflow run main.nf -profile docker,local -resume --input /path/to/input/csv")
       exit 0
    }

    // Validate input parameters
    validateParameters()

    // Print summary of supplied parameters
    log.info paramsSummaryLog(workflow)
    
    // Make a channel for version outputs:
    ch_versions = Channel.empty()

    // Check if input is provided
    in_file = params.input != null ? Channel.fromPath(params.input) : errorMessage()

    in_file
        .splitCsv()
        .branch {
            ncbi: it.size() == 2
            path: it.size() == 3
        }
        .set { input_type }

    DOWNLOAD_NCBI ( input_type.ncbi )
    ch_versions = ch_versions.mix(DOWNLOAD_NCBI.out.versions.first())

    // Checks if paths are S3 objects if not ensures absolute paths are used for user inputted fasta and gff files
    input_type.path.map{ name, fasta , gff -> 
        def full_fasta = fasta =~ /^s3:\/\// ? fasta : new File(fasta).getAbsolutePath()
        def full_gff = gff =~ /^s3:\/\// ? gff : new File(gff).getAbsolutePath()
        [name, full_fasta, full_gff]
    }.set { local_full_tuple }

    // Split channel into 2, keep tuple the same for gffread and take just sample id and fasta for fastavalidator
    DOWNLOAD_NCBI.out.genome.mix(local_full_tuple)
        .multiMap { it ->
            gffread: it
            tuple: [[ id: it[0]], it[1]]
        }
        .set { fasta_inputs }
    FASTAVALIDATOR ( fasta_inputs.tuple )
    ch_versions = ch_versions.mix(FASTAVALIDATOR.out.versions.first())

    // Manipulate successful and error logs of fasta validator to be saved into output directory
    FASTAVALIDATOR.out.success_log.map { speciesname, logfile -> [speciesname.id, logfile] }
        .collectFile(name: { it[0] }, storeDir: "${params.outdir}/input_validation/fasta_validator/successful")
    FASTAVALIDATOR.out.error_log.map { speciesname, logfile -> [speciesname.id, logfile] }
        .collectFile(name: { it[0] }, storeDir: "${params.outdir}/input_validation/fasta_validator/error")

    SEQKIT_STATS( fasta_inputs.tuple )
    ch_versions = ch_versions.mix(SEQKIT_STATS.out.versions.first())

    // Manipulate seqkit_stats tsv to be saved into output directory
    SEQKIT_STATS.out.stats.map { speciesname, tsv -> [speciesname.id, tsv] }
        .collectFile(name: { it[0] }, storeDir: "${params.outdir}/input_validation/seqkit_stats")
 
    LONGEST ( fasta_inputs.gffread )
    ch_versions = ch_versions.mix(LONGEST.out.versions.first())

    GFFREAD ( LONGEST.out.longest_proteins )
    ch_versions = ch_versions.mix(GFFREAD.out.versions.first())

    JCVI ( GFFREAD.out.proteins )
    ch_versions = ch_versions.mix(JCVI.out.versions.first())

    // Do a pairwise combination of each species' JCVI output but filter out combinations of the same species
    SYNTENY ( JCVI.out.new_format.combine(JCVI.out.new_format).filter { it[0] != it[3] } )
    ch_versions = ch_versions.mix(SYNTENY.out.versions.first())

    // Use name of anchors file to identify the 2 species involved and create a tuple with these species as strings
    SYNTENY.out.anchors.map { it -> 
        def (sample1, sample2) = it.baseName.toString().split("\\.")
        [sample1, sample2, it]
    }.set { labelled_anchors }

    // Combine hex path with each tuple of species and anchor files for parallelisation of process
    in_hex.combine(labelled_anchors).set { hex_labelled_anchors }

    if (params.chromopaint) {
        CHROMOPAINT ( hex_labelled_anchors, JCVI.out.beds.collect() )
        ch_versions = ch_versions.mix(CHROMOPAINT.out.versions.first())
    }

    // To run the ribbon plot for each species given
    if (params.ribbon) {
        ribbonIn = Channel.from(params.ribbon)
        RIBBON ( SYNTENY.out.anchors_notlifted.collect(), JCVI.out.beds.collect(), ribbonIn )
    }

    // Code to measure a synteny scores for each species
    if (params.score) {
        if (params.lifted) {
            if (params.tree) {
                treeIn = Channel.fromPath(params.tree)
                SCORE_TREE ( SYNTENY.out.anchors.collect(), SYNTENY.out.percsim.collect(), GFFREAD.out.gff.collect(), JCVI.out.beds.collect(), SYNTENY.out.last.collect(), SYNTENY.out.unfilteredlast.collect(), treeIn )
                ch_versions = ch_versions.mix(SCORE_TREE.out.versions)
                SCORE_TREE_PLOTS(SCORE_TREE.out.filec, SCORE_TREE.out.species_order)
            } else {
                SCORE ( SYNTENY.out.anchors.collect(), SYNTENY.out.percsim.collect(), GFFREAD.out.gff.collect(), JCVI.out.beds.collect(), SYNTENY.out.last.collect(), SYNTENY.out.unfilteredlast.collect() )
                ch_versions = ch_versions.mix(SCORE.out.versions)
                SCORE_PLOTS(SCORE.out.filec)
            }
        } else {
            if (params.tree) {
                treeIn = Channel.fromPath(params.tree)
                SCORE_TREE ( SYNTENY.out.anchors_notlifted.collect(), SYNTENY.out.percsim.collect(), GFFREAD.out.gff.collect(), JCVI.out.beds.collect(), SYNTENY.out.last.collect(), SYNTENY.out.unfilteredlast.collect(), treeIn )
                ch_versions = ch_versions.mix(SCORE_TREE.out.versions)
                SCORE_TREE_PLOTS(SCORE_TREE.out.filec, SCORE_TREE.out.species_order)
            } else {
                SCORE ( SYNTENY.out.anchors_notlifted.collect(), SYNTENY.out.percsim.collect(), GFFREAD.out.gff.collect(), JCVI.out.beds.collect(), SYNTENY.out.last.collect(), SYNTENY.out.unfilteredlast.collect() )
                ch_versions = ch_versions.mix(SCORE.out.versions)
                SCORE_PLOTS(SCORE.out.filec)
            }
        }
    }

    // If you choose to run go
    if ( params.go && params.score ) {
        go_folder = Channel.fromPath(params.go)
        //Checks if SCORE_TREE output is not null and uses it, if it is null then SCORE was run instead and use that output

        if ( params.tree ){
            species_summary = SCORE_TREE.out.speciesSummary 
        }
        else{
            species_summary = SCORE.out.speciesSummary
        }

        //creating 3 instances of a channel with the GO hash files and species summary files 
        go_folder.combine(species_summary.flatten()).set{ go_and_summary }

        //Add input params cutoff for go
        our_cutoff = Channel.of(params.cutoff)

        // Split the params.cutoff string into a list of separate entries
        cutoffValues = Channel.from(params.cutoff.split(','))

        // Combine channels using cross to pair each element of go_and_summary with each element of cutoffValues
        mergedChannel = go_and_summary.combine(cutoffValues)

        //Get the inversion translocation scores in a new channel

        if ( params.tree ){
            species_inver = SCORE_TREE.out.geneinverdistancescores.flatten()
            species_trans = SCORE_TREE.out.genetransdistancescores.flatten()
        }
        else{
            species_inver = SCORE.out.geneinverdistancescores.flatten()
            species_trans = SCORE.out.genetransdistancescores.flatten()
        }

        //creating 3 instances of a channel with the GO hash files and species summary files 
        go_folder2 = Channel.fromPath(params.go)
        go_folder2.combine(species_inver.flatten()).set{ go_and_summary2 }

        // Combine the channels of species files with the cutoff values
        mergedChannel2 = go_and_summary2.combine(cutoffValues)

        GO ( mergedChannel , JCVI.out.beds.collect() )
        GO_JUNCTIONS ( mergedChannel2 , JCVI.out.beds.collect() )
        //GO_JUNCTIONS_TRANS ( species_trans , JCVI.out.beds.collect() )

        ch_versions = ch_versions.mix(GO.out.versions.first())


        GO_SUMMARISE ( GO.out.go_table.groupTuple() )
        ch_versions = ch_versions.mix(GO_SUMMARISE.out.versions)
        SUMMARISE_PLOTS(GO_SUMMARISE.out.go_summary_table)
    }


    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.collectFile(name: 'collated_versions.yml')
    )
}

workflow.onComplete {
    println(workflow.success ? "\nDone! Check results in $params.outdir/ \n" : "Hmmm .. something went wrong\n")
}
