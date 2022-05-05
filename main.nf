/*
 * Copyright (c) 2022
 */

 /*
 * Authors:
 * - Chris Wyatt <chris.wyatt@seqera.io>
 */

/*
 * Default pipeline parameters (on test data). They can be overriden on the command line eg.
 * given `params.name` specify on the run command line `--name My_run_v1`.
 */
 
params.outdir = "Results"
params.input = "data/Example.csv"
params.seqids = "data/default1"
params.layout = "data/default2"
params.test=0


log.info """\
 ===================================
 input file                           : ${params.input}
 output directory                     : ${params.outdir}
 seqids file (optional)               : ${params.seqids}
 layout file (optional)               : ${params.layout}
 """

//================================================================================
// Include modules
//================================================================================

include { GFFREAD } from './modules/gffread.nf'
include { JCVI } from './modules/jcvi.nf'
include { SYNTENY } from './modules/synteny.nf'
include { MACRO } from './modules/macro.nf'
include { CONFIG } from './modules/default_config.nf'
include { DOWNLOAD_NCBI } from './modules/download_ncbi.nf'

Channel
	.fromPath(params.input)
    .splitCsv()
    .set { in_file }

Channel
    .fromPath(params.input)
    .splitCsv()
    .collect()
    .view()
    .set { in_file_config }

Channel
    .fromPath(params.seqids)
    .set { in_seqids }

Channel
    .fromPath(params.layout)
    .set { in_layout }   


workflow {

    GFFREAD ( in_file )

	JCVI ( GFFREAD.out.proteins )
    //JCVI.out.collect().view()

    CONFIG ( in_seqids , in_layout , in_file_config )

    SYNTENY ( JCVI.out.collect() )
    MACRO ( CONFIG.out.seqids_out , CONFIG.out.layout_out , SYNTENY.out.anchors, JCVI.out.collect() )

}

workflow.onComplete {
    println ( workflow.success ? "\nDone! Check results in $params.outdir/ \n" : "Hmmm .. something went wrong\n" )
}
