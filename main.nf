/*
 * Copyright (c) 2021
 */
 

 /*
 * Authors:
 * - Chris Wyatt <chris.wyatt@seqera.io>
 */

/* 
 * enable modules 
 */
nextflow.enable.dsl = 2

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
include { CREATE_SEQIDS } from './modules/create_seqids.nf'

Channel
	.fromPath(params.input)
    .splitCsv()
    .set { in_file }

Channel
    .fromPath(params.seqids)
    .set { in_seqids }

Channel
    .fromPath(params.layout)
    .set { in_layout }   


workflow {
    GFFREAD ( in_file )
	JCVI ( GFFREAD.out.proteins , in_seqids )
    JCVI.out.reformatted.collect().view()
    


    SYNTENY ( JCVI.out.reformatted.collect() )
    MACRO ( JCVI.out.seqids_out.collect() , in_layout , SYNTENY.out.anchors, JCVI.out.reformatted.collect() )
}

workflow.onComplete {
    println ( workflow.success ? "\nDone! Check results in $params.outdir/ \n" : "Hmmm .. something went wrong\n" )
}
