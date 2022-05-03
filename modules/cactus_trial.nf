process CACTUS_TRIAL {
    label 'cactus'
    publishDir "$params.outdir/Cactus_results/"
    container = 'quay.io/comparative-genomics-toolkit/cactus:v2.0.5'
    //stageInMode 'copy'
             
    input:

    path("*")

    output:
        
	path("*") , emit: cactus_result

    script:
    """
	cactus jobStore evolverMammals.txt evolverMammals.hal --root mr
    """
}

