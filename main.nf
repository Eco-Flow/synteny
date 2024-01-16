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

include { GFFREAD } from './modules/gffread.nf'
include { JCVI } from './modules/jcvi.nf'
include { SYNTENY } from './modules/synteny.nf'
include { DOWNLOAD_NCBI } from './modules/download_ncbi.nf'
include { CHROMOPAINT } from './modules/chromo.nf'
include { SCORE } from './modules/score.nf'
include { GO } from './modules/go.nf'
include { SCORE_TREE } from './modules/score_tree.nf'
include { GO_SUMMARISE } from './modules/go_summarise.nf'

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

    GFFREAD ( DOWNLOAD_NCBI.out.genome.mix(local_full_tuple) )

    JCVI ( GFFREAD.out.proteins )

    //Do a pairwise combination of each species' JCVI output but filter out combinations of the same species
    SYNTENY ( JCVI.out.new_format.combine(JCVI.out.new_format).filter{ it[0] != it[3] } )

    CHROMOPAINT ( in_hex , SYNTENY.out.anchors , JCVI.out.beds.collect() )

    if (params.tree) {
        tree_in = Channel.fromPath(params.tree)
        SCORE_TREE ( SYNTENY.out.anchors.collect() , SYNTENY.out.percsim.collect() , GFFREAD.out.gff.collect() , tree_in )
    }
    else {
        SCORE ( SYNTENY.out.anchors.collect() , SYNTENY.out.percsim.collect() , GFFREAD.out.gff.collect() )
    }
    if (params.go) {
        //Checks if SCORE_TREE output is not null and uses it, if it is null then SCORE was run instead and use that output
        ch_go = params.tree != null ? SCORE_TREE.out.speciesSummary : SCORE.speciesSummary
        ch_go.view()
        GO ( go_datasets.collect() , ch_go.flatten(), JCVI.out.beds.collect() )
        GO_SUMMARISE ( GO.out.go_table.collect() )
    }
}

workflow.onComplete {
    println ( workflow.success ? "\nDone! Check results in $params.outdir/ \n" : "Hmmm .. something went wrong\n" )
}
