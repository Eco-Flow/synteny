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
params.seqids = "data/seqids"
params.layout = "data/layout"
params.test=0


log.info """\
 ===================================
 input file                           : ${params.input}
 out directory                        : ${params.outdir}
 """

//================================================================================
// Include modules
//================================================================================

include { GFFREAD } from './modules/gffread.nf'
include { JCVI } from './modules/jcvi.nf'
include { SYNTENY } from './modules/synteny.nf'
include { MACRO } from './modules/macro.nf'

Channel
	.fromPath(params.input)
    .splitCsv()
    .view () { row -> "${row[0]},${row[1]}" }
    .set { in_file }

Channel
    .fromPath(params.seqids)
    .set { in_seqids }

Channel
    .fromPath(params.layout)
    .set { in_layout }

workflow {
    GFFREAD ( in_file )
	JCVI ( GFFREAD.out.proteins )
    SYNTENY ( JCVI.out.collect() )
    MACRO ( in_seqids , in_layout , SYNTENY.out.anchors, JCVI.out.collect() )
}

workflow.onComplete {
	println ( workflow.success ? "\nDone! Open your report in your browser --> $params.outdir/report.html (if you added -with-report flag)\n" : "Hmmm .. something went wrong" )
}
