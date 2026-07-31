#!/usr/bin/env nextflow

include { PAIRWISE_SYNTENY } from './workflows/pairwise_synteny.nf'
include { ALGO } from './workflows/algo.nf'

workflow {

    // The strict parser used by this pipeline does not support `-entry`, so both
    // named workflows are dispatched here via --mode instead: the default pairwise
    // synteny workflow, or ALGO (ancestral linkage groups / ancestral gene order).
    if (params.mode == 'algo') {
        ALGO()
    } else {
        PAIRWISE_SYNTENY()
    }
}
