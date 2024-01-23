/*
 * Copyright (c) 2022
 */

/*
 * Authors:
 * - Chris Wyatt <c.wyatt@ucl.ac.uk>
 * - Simon Murray <simon.murray@ucl.ac.uk>
 */

/*
 * Default pipeline parameters (on test data). They can be overriden on the command line eg.
 * given `params.name` specify on the run command line `--name My_run_v1`.
 */

log.info """\
 ===================================

 nf-synteny (v3.0.0)

 ===================================
 input file                           : ${params.input}
 output directory                     : ${params.outdir}
 running go ?                         : ${params.go}
 """

//================================================================================
// Include modules
//================================================================================

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
include { JCVI } from './modules/local/jcvi.nf'
include { SYNTENY } from './modules/local/synteny.nf'
include { CHROMOPAINT } from './modules/local/chromo.nf'
include { SCORE } from './modules/local/score.nf'
include { SCORE_TREE } from './modules/local/score_tree.nf'
include { GO } from './modules/local/go.nf'
include { GO_SUMMARISE } from './modules/local/go_summarise.nf'
include { FASTAVALIDATOR } from './modules/nf-core/fastavalidator/main'

Channel
    .fromPath(params.hex)
    .set { in_hex }

workflow {

    //Check if input is provided
    in_file = params.input != null ? Channel.fromPath(params.input) : errorMessage()

    in_file
        .splitCsv()
        .branch {
            ncbi: it.size() == 2
            local: it.size() == 3
        }
        .set { input_type }

    DOWNLOAD_NCBI ( input_type.ncbi )

    //Ensures absolute paths are used if user inputs fasta and gff files
    input_type.local.map{ name, fasta , gff -> full_fasta = new File(fasta).getAbsolutePath(); full_gff = new File(gff).getAbsolutePath(); [name, full_fasta, full_gff] }.set{ local_full_tuple }

    //Split channel into 2, keep tuple the same for gffread and take just sample id and fasta for fastavalidator
    DOWNLOAD_NCBI.out.genome.mix(local_full_tuple)
         .multiMap { it ->
             gffread: it
             fastqc: [[ id: it[0]], it[1]]
          }
          .set { fasta_inputs }

    FASTAVALIDATOR(fasta_inputs.fastqc)

    GFFREAD ( fasta_inputs.gffread )

    JCVI ( GFFREAD.out.proteins )

    //Do a pairwise combination of each species' JCVI output but filter out combinations of the same species
    SYNTENY ( JCVI.out.new_format.combine(JCVI.out.new_format).filter{ it[0] != it[3] } )

    //Use name of anchors file to identify the 2 species involved and create a tuple with these species as strings
    SYNTENY.out.anchors.map{ it -> def(sample1, sample2) = it.baseName.toString().split("\\."); [sample1, sample2, it] }.set{ labelled_anchors }

    //Combine hex path with each tuple of species and anchor files for parallelisation of process
    in_hex.combine(labelled_anchors).set{ hex_labelled_anchors }

    CHROMOPAINT ( hex_labelled_anchors, JCVI.out.beds.collect() )

    if (params.tree) {
        tree_in = Channel.fromPath(params.tree)
        SCORE_TREE ( SYNTENY.out.anchors.collect(), SYNTENY.out.percsim.collect(), GFFREAD.out.gff.collect(), tree_in )
    }
    else {
        SCORE ( SYNTENY.out.anchors.collect(), SYNTENY.out.percsim.collect(), GFFREAD.out.gff.collect() )
    }
    if (params.go) {
        go_folder = Channel.fromPath(params.go)
        //Checks if SCORE_TREE output is not null and uses it, if it is null then SCORE was run instead and use that output
        species_summary = params.tree != null ? SCORE_TREE.out.speciesSummary : SCORE.out.speciesSummary
        //creating 3 instances of a channel with the GO hash files and species summary files 
        go_folder.combine(species_summary.flatten()).set{ go_and_summary }
        GO ( go_and_summary, JCVI.out.beds.collect() )
        GO_SUMMARISE ( GO.out.go_table.collect() )
    }
}

workflow.onComplete {
    println ( workflow.success ? "\nDone! Check results in $params.outdir/ \n" : "Hmmm .. something went wrong\n" )
}
