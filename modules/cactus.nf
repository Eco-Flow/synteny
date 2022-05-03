process CACTUS {
    label 'cactus'
    publishDir "$params.outdir/Cactus_results/"
    //stageInMode 'copy'
    
    input:
	path genomes
               
    output:
        
	path("results_folder") , emit: cactus_result

    script:
    """
	cactus -genomes results
    """
}
