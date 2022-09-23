process DOWNLOAD_NCBI {
    label 'download'
    tag "$sample_id via $accension_id"
             
    input:
        tuple val(sample_id), val(accension_id)

    output:    
        tuple val(sample_id), path("${sample_id}.genome.fna"), path("${sample_id}.genomic.gff"), emit: genome

    script:
    """

    #Get a genome and GFF assembly from NCBI using their datasets scripts

    datasets download genome accession ${accension_id}
    unzip ncbi_dataset.zip 

    download_collect.pl ${accension_id} ${sample_id}

    """
}

