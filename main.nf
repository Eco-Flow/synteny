#!/usr/bin/env nextflow

log.info """\
 =========================================

 nf-synteny (v4.0.0)

 -----------------------------------------

 Authors:
   - Chris Wyatt <c.wyatt@ucl.ac.uk>
   - Simon Murray <simon.murray@ucl.ac.uk>

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
                SCORE_TREE ( SYNTENY.out.anchors.collect(), SYNTENY.out.percsim.collect(), GFFREAD.out.gff.collect(), JCVI.out.beds.collect(), SYNTENY.out.last.collect(), treeIn )
                ch_versions = ch_versions.mix(SCORE_TREE.out.versions)
                SCORE_TREE_PLOTS(SCORE_TREE.out.filec, SCORE_TREE.out.species_order)
            } else {
                SCORE ( SYNTENY.out.anchors.collect(), SYNTENY.out.percsim.collect(), GFFREAD.out.gff.collect(), JCVI.out.beds.collect(), SYNTENY.out.last.collect())
                ch_versions = ch_versions.mix(SCORE.out.versions)
                SCORE_PLOTS(SCORE.out.filec)
            }
        } else {
            if (params.tree) {
                treeIn = Channel.fromPath(params.tree)
                SCORE_TREE ( SYNTENY.out.anchors_notlifted.collect(), SYNTENY.out.percsim.collect(), GFFREAD.out.gff.collect(), JCVI.out.beds.collect(), SYNTENY.out.last.collect(), treeIn )
                ch_versions = ch_versions.mix(SCORE_TREE.out.versions)
                SCORE_TREE_PLOTS(SCORE_TREE.out.filec, SCORE_TREE.out.species_order)
            } else {
                SCORE ( SYNTENY.out.anchors_notlifted.collect(), SYNTENY.out.percsim.collect(), GFFREAD.out.gff.collect(), JCVI.out.beds.collect(), SYNTENY.out.last.collect())
                ch_versions = ch_versions.mix(SCORE.out.versions)
                SCORE_PLOTS(SCORE.out.filec)
            }
        }
    }

    // If you choose to run go
    if (params.go && params.score) {
        goFolder = Channel.fromPath(params.go)
        // Checks if SCORE_TREE output is not null and uses it, if it is null then SCORE was run instead and use that output

        def speciesSummary = params.tree ? SCORE_TREE.out.speciesSummary : SCORE.out.speciesSummary

        // Creating 3 instances of a channel with the GO hash files and species summary files 
        goAndSummary = goFolder.combine(speciesSummary.flatten())
        goAndSummaryCut = goAndSummary.combine(our_cutoff)

        GO ( goAndSummaryCut, JCVI.out.beds.collect() )
        ch_versions = ch_versions.mix(GO.out.versions.first())

        //Calculate size of input and run next piece of code.
        buffer_val = new AtomicInteger()
        Channel
            .fromPath(params.input)
            .splitCsv()
            .count()
            .subscribe { count ->
                println "The pipeline will run with ${count} species"
                count_sp = count
                count_sp_6 = count_sp * 6
                buffer_val.set(count_sp_6 + 1)
                println "Buffer is ${buffer_val.get()}"
                cutoffGroups = GO.out.go_table.groupTuple()
                .flatten()
                .buffer(size: buffer_val.get())
                //cutoffGroups.view()
                //GO_SUMMARISE ( cutoffGroups )
            }

        def filePath = params.input
        // Read file length and save as integer
        def fileLength = Files.size(Paths.get(filePath)).intValue()

        // Print the file length for verification
        println "File length: $fileLength"

        //println "Length of the channel: $channelLength"

        //def buffer_val = ${count_sp_6} + 1

        //flattenedGroups = flattenCutoffGroups(cutoffGroups)



        //println list.size()


        
        //GO_SUMMARISE ( cutoffGroups )
    }

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.collectFile(name: 'collated_versions.yml')
    )
}

workflow.onComplete {
    println(workflow.success ? "\nDone! Check results in $params.outdir/ \n" : "Hmmm .. something went wrong\n")
}
