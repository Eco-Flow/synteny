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
include { CHROMOPAINT } from './modules/chromo.nf'

Channel
	.fromPath(params.input)
    .splitCsv()
    .set { in_file }

Channel
    .fromPath(params.input)
    .splitCsv()
    .collect()
    .set { in_file_config }

Channel
    .fromPath(params.seqids)
    .set { in_seqids }

Channel
    .fromPath(params.layout)
    .set { in_layout }   


Channel
    .fromPath(params.input)
    .splitCsv()
    .branch { 
        ncbi: it.size() == 2 
        local: it.size() == 3
    }
    .set { input_type }
    
// input_type.ncbi.view { "$it is small" }
// input_type.local.view { "$it is large" }



workflow {

    DOWNLOAD_NCBI ( input_type.ncbi )

    GFFREAD ( DOWNLOAD_NCBI.out.genome.mix(input_type.local) )

    JCVI ( GFFREAD.out.proteins )

    //JCVI.out.new_format.combine(JCVI.out.new_format).filter{ it[0] != it[3] }.view()

    //CONFIG ( in_seqids , in_layout , DOWNLOAD_NCBI.out.genome.mix(input_type.local).collect() )

    SYNTENY ( JCVI.out.new_format.combine(JCVI.out.new_format).filter{ it[0] != it[3] } )

    //MACRO ( CONFIG.out.seqids_out , CONFIG.out.layout_out , SYNTENY.out.anchors, JCVI.out.collect() )

    CHROMOPAINT ( SYNTENY.out.anchors , JCVI.out.collect() )
    
}

workflow.onComplete {
    println ( workflow.success ? "\nDone! Check results in $params.outdir/ \n" : "Hmmm .. something went wrong\n" )
}

