process DOWNLOAD_NCBI {
    label 'download'
    tag "$sample_id via $accension_id"
             
    input:
        tuple val(sample_id), val(accension_id)

    output:    
        tuple val(sample_id), path("genome.fna"), path("genomic.gff"), emit: genome

    script:
    """

    #Get a genome and GFF assembly from NCBI using their datasets scripts

    datasets download genome accession ${accension_id}
    unzip ncbi_dataset.zip 

    cat ncbi_dataset/data/${accension_id}/chr*.fna > genome.fna
    cat ncbi_dataset/data/${accension_id}/unplaced.scaf.fna >> genome.fna 
    cat ncbi_dataset/data/${accension_id}/genomic.gff > genomic.gff

    """
}

