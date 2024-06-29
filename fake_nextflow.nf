#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Example of input params
params.cutoff = '10,20'

// Define go_and_summary with file paths
go_and_summary = Channel.of(
    ['/workspace/synteny/data/go_input/hash_files', '/workspace/synteny/work/22/2442c476ac7d5a6cded7e48e7178c4/Drosophila_santomea.SpeciesScoreSummary.txt'],
    ['/workspace/synteny/data/go_input/hash_files', '/workspace/synteny/work/22/2442c476ac7d5a6cded7e48e7178c4/Drosophila_simulans.SpeciesScoreSummary.txt'],
    ['/workspace/synteny/data/go_input/hash_files', '/workspace/synteny/work/22/2442c476ac7d5a6cded7e48e7178c4/Drosophila_yakuba.SpeciesScoreSummary.txt']
)

// Ensure that params.cutoff is a string
cutoffString = params.cutoff.toString()

// Split the params.cutoff string into a list of separate entries
cutoffValues = Channel.from(cutoffString.split(','))

// Combine channels using cross to pair each element of go_and_summary with each element of cutoffValues
mergedChannel = go_and_summary.cross(cutoffValues).map { pair -> 
    def (go_and_summary_entry, cutoffValue) = pair
    go_and_summary_entry + [cutoffValue]
}

// View the merged channel results
mergedChannel.view()
