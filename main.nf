#!/usr/bin/env nextflow

log.info """\
 =========================================

 nf-synteny (v4.0.0)

 -----------------------------------------

 Authors:
   - Chris Wyatt <c.wyatt@ucl.ac.uk>
   - Simon Murray

 -----------------------------------------

 Copyright (c) 2024

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
include { SCORE3 } from './modules/local/score3.nf'
include { SCORE_TREE } from './modules/local/score_tree.nf'
include { GO } from './modules/local/go.nf'
include { GO_JUNCTIONS_INVER } from './modules/local/go_junctions_inver_score3.nf'
include { GO_JUNCTIONS_INTER } from './modules/local/go_junctions_inter_score3.nf'
include { GO_JUNCTIONS_INDEL_LARGE } from './modules/local/go_junctions_indel_large_score3.nf'
include { GO_JUNCTIONS_INDEL_TINY } from './modules/local/go_junctions_indel_tiny_score3.nf'
include { GO_JUNCTIONS_INDEL_SMALL } from './modules/local/go_junctions_indel_small_score3.nf'
include { GO_JUNCTIONS_TRANS } from './modules/local/go_junctions_trans.nf'
include { GO_JUNCTIONS_OTHER } from './modules/local/go_junctions_other.nf'
include { GO_SUMMARISE } from './modules/local/go_summarise.nf'
include { GO_SUMMARISE_INVER } from './modules/local/go_summarise_inver.nf'
include { GO_SUMMARISE_INTER } from './modules/local/go_summarise_inter.nf'
include { GO_SUMMARISE_INDEL_LARGE } from './modules/local/go_summarise_indel_large.nf'
include { GO_SUMMARISE_INDEL_TINY } from './modules/local/go_summarise_indel_tiny.nf'
include { GO_SUMMARISE_INDEL_SMALL } from './modules/local/go_summarise_indel_small.nf'
include { GO_SUMMARISE_TRANS } from './modules/local/go_summarise_trans.nf'
include { GO_SUMMARISE_OTHER } from './modules/local/go_summarise_other.nf'
include { SUMMARISE_PLOTS } from './modules/local/summarise_plots.nf'
include { SUMMARISE_PLOTS_INVER } from './modules/local/summarise_plots_inver.nf'
include { SUMMARISE_PLOTS_INTER } from './modules/local/summarise_plots_inter.nf'
include { SUMMARISE_PLOTS_INDEL_LARGE } from './modules/local/summarise_plots_indel_large.nf'
include { SUMMARISE_PLOTS_INDEL_TINY } from './modules/local/summarise_plots_indel_tiny.nf'
include { SUMMARISE_PLOTS_INDEL_SMALL } from './modules/local/summarise_plots_indel_small.nf'

include { SUMMARISE_PLOTS_TRANS } from './modules/local/summarise_plots_trans.nf'
include { SUMMARISE_PLOTS_OTHER } from './modules/local/summarise_plots_other.nf'
include { FASTAVALIDATOR } from './modules/nf-core/fastavalidator/main'
include { SEQKIT_STATS } from './modules/nf-core/seqkit/stats/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from './modules/nf-core/custom/dumpsoftwareversions'
include { validateParameters; paramsHelp; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'
include { SCORE_PLOTS } from './modules/local/score_plot.nf'
include { SCORE_PLOTS_2 } from './modules/local/score_plot2.nf'
include { SCORE_PLOTS_3 } from './modules/local/score_plot3.nf'
include { SCORE_PLOT_TREE } from './modules/local/score_plot_tree.nf'
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
                SCORE3 ( SYNTENY.out.anchors, SYNTENY.out.percsim.collect(), GFFREAD.out.gff.collect(), JCVI.out.beds.collect(), SYNTENY.out.last.collect(), SYNTENY.out.unfilteredlast.collect() )
                SCORE_PLOTS_3 ( SCORE_TREE.out.trans_inver_summary,  SCORE3.out.junction_summary.collect() )

            } else {
                SCORE ( SYNTENY.out.anchors.collect(), SYNTENY.out.percsim.collect(), GFFREAD.out.gff.collect(), JCVI.out.beds.collect(), SYNTENY.out.last.collect(), SYNTENY.out.unfilteredlast.collect() )
                ch_versions = ch_versions.mix(SCORE.out.versions)
                SCORE3 ( SYNTENY.out.anchors, SYNTENY.out.percsim.collect(), GFFREAD.out.gff.collect(), JCVI.out.beds.collect(), SYNTENY.out.last.collect(), SYNTENY.out.unfilteredlast.collect() )
                SCORE_PLOTS_3 ( SCORE.out.trans_inver_summary,  SCORE3.out.junction_summary.collect() )
            }
        } else {
            if (params.tree) {
                treeIn = Channel.fromPath(params.tree)
                SCORE_TREE ( SYNTENY.out.anchors_notlifted.collect(), SYNTENY.out.percsim.collect(), GFFREAD.out.gff.collect(), JCVI.out.beds.collect(), SYNTENY.out.last.collect(), SYNTENY.out.unfilteredlast.collect(), treeIn )
                ch_versions = ch_versions.mix(SCORE_TREE.out.versions)
                SCORE3 ( SYNTENY.out.anchors_notlifted, SYNTENY.out.percsim.collect(), GFFREAD.out.gff.collect(), JCVI.out.beds.collect(), SYNTENY.out.last.collect(), SYNTENY.out.unfilteredlast.collect() )
                SCORE_PLOTS_3 ( SCORE_TREE.out.trans_inver_summary,  SCORE3.out.junction_summary.collect() )

            } else {
                SCORE ( SYNTENY.out.anchors_notlifted.collect(), SYNTENY.out.percsim.collect(), GFFREAD.out.gff.collect(), JCVI.out.beds.collect(), SYNTENY.out.last.collect(), SYNTENY.out.unfilteredlast.collect() )
                ch_versions = ch_versions.mix(SCORE.out.versions)
                SCORE3 ( SYNTENY.out.anchors_notlifted, SYNTENY.out.percsim.collect(), GFFREAD.out.gff.collect(), JCVI.out.beds.collect(), SYNTENY.out.last.collect(), SYNTENY.out.unfilteredlast.collect() )
                SCORE_PLOTS_3 ( SCORE.out.trans_inver_summary, SCORE3.out.junction_summary.collect() )
            }
        }
    }

    // Add go data folder:
    go_folder = Channel.fromPath(params.go)
    // Split the params.cutoff string into a list of separate entries
    cutoffValues = Channel.from(params.cutoff.split(','))

    //Combine all the scores for each pairwise comparison:
    species_inver = SCORE3.out.inver_distancescores.groupTuple()
    species_inter = SCORE3.out.inter_distancescores.groupTuple()
    species_indel_large = SCORE3.out.large_indel_distancescores.groupTuple()
    species_indel_tiny  = SCORE3.out.tiny_indel_distancescores.groupTuple()
    species_indel_small = SCORE3.out.small_indel_distancescores.groupTuple()

    //For each junction type, combine the junctions with the go hashes and cutoff values (default 10,20), and run GO analysis.
    go_and_summary_inver = go_folder.combine(species_inver)
    mergedChannel_inver = go_and_summary_inver.combine(cutoffValues)
    GO_JUNCTIONS_INVER ( mergedChannel_inver , JCVI.out.beds.collect() )

    go_and_summary_inter = go_folder.combine(species_inter)
    mergedChannel_inter = go_and_summary_inter.combine(cutoffValues)
    GO_JUNCTIONS_INTER ( mergedChannel_inter , JCVI.out.beds.collect() )

    go_and_summary_indel_large = go_folder.combine(species_indel_large)
    mergedChannel_indel_large = go_and_summary_indel_large.combine(cutoffValues)
    GO_JUNCTIONS_INDEL_LARGE ( mergedChannel_indel_large , JCVI.out.beds.collect() )

    go_and_summary_indel_tiny = go_folder.combine(species_indel_tiny)
    mergedChannel_indel_tiny = go_and_summary_indel_tiny.combine(cutoffValues)
    GO_JUNCTIONS_INDEL_TINY ( mergedChannel_indel_tiny , JCVI.out.beds.collect() )

    go_and_summary_indel_small = go_folder.combine(species_indel_small)
    mergedChannel_indel_small = go_and_summary_indel_small.combine(cutoffValues)
    GO_JUNCTIONS_INDEL_SMALL ( mergedChannel_indel_small , JCVI.out.beds.collect() )

    ch_versions = ch_versions.mix(GO_JUNCTIONS_INVER.out.versions.first())

    //Summarise all the go tables in a single table
    GO_SUMMARISE_INVER ( GO_JUNCTIONS_INVER.out.go_table.groupTuple() )
    GO_SUMMARISE_INTER ( GO_JUNCTIONS_INTER.out.go_table.groupTuple() )
    GO_SUMMARISE_INDEL_LARGE ( GO_JUNCTIONS_INDEL_LARGE.out.go_table.groupTuple() )
    GO_SUMMARISE_INDEL_TINY ( GO_JUNCTIONS_INDEL_TINY.out.go_table.groupTuple() )
    GO_SUMMARISE_INDEL_SMALL ( GO_JUNCTIONS_INDEL_SMALL.out.go_table.groupTuple() )
    ch_versions = ch_versions.mix(GO_SUMMARISE_INVER.out.versions)

    //Plot all summary go plots
    SUMMARISE_PLOTS_INVER(GO_SUMMARISE_INVER.out.go_summary_table)
    SUMMARISE_PLOTS_INTER(GO_SUMMARISE_INTER.out.go_summary_table)
    SUMMARISE_PLOTS_INDEL_LARGE(GO_SUMMARISE_INDEL_LARGE.out.go_summary_table)
    SUMMARISE_PLOTS_INDEL_TINY(GO_SUMMARISE_INDEL_TINY.out.go_summary_table)
    SUMMARISE_PLOTS_INDEL_SMALL(GO_SUMMARISE_INDEL_SMALL.out.go_summary_table)
    ch_versions = ch_versions.mix(SUMMARISE_PLOTS_INVER.out.versions)

}

workflow.onComplete {
    println(workflow.success ? "\nDone! Check results in $params.outdir/ \n" : "Hmmm .. something went wrong\n")
}
