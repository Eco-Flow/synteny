process SCORES {
    label 'score'
    tag "$sample_id"
    container = 'chriswyatt/jcvi'
             
    input:

        tuple val(sample_id), path(fasta), path(gff)

    output:
        
        tuple val(sample_id), path( "${sample_id}.cds" ), path( "${sample_id}.bed" ) , emit: new_format

    script:
    """
    #Run Score for each gene on how close it is to the edge of the syntenic block


    #Run score for genome X in terms of size of syntenic blacks to species Y.

    

    """
}

