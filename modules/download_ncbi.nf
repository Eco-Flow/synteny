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

    
    if ls ncbi_dataset/data/${accension_id}/chr*.fna 1> /dev/null 2>&1; then
        cat ncbi_dataset/data/${accension_id}/chr*.fna > ${sample_id}.genome.fna
    fi
    if ls ncbi_dataset/data/${accension_id}/unplaced.scaf.fna 1> /dev/null 2>&1; then
        cat ncbi_dataset/data/${accension_id}/unplaced.scaf.fna >> ${sample_id}.genome.fna 
    fi
    
    cat ncbi_dataset/data/${accension_id}/genomic.gff > ${sample_id}.genomic.gff

    """
}

