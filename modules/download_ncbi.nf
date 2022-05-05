process DOWNLOAD_NCBI {
    label 'download'
    tag 'download'
             
    input:
    val(accension_id)

    output:    
    tuple val(accension_id), path("genome.fna"), path("genomic.gff")

    script:
    """
    #Get a genome and GFF assembly from NCBI using their datasets scripts

    conda create -n ncbi_datasets
    conda activate ncbi_datasets
    conda install -c conda-forge ncbi-datasets-cli

    datasets download genome accession GCF_003254395.2
    unzip ncbi_dataset.zip 
    cd ncbi_dataset/data/${accension_id}

    cat chr*.fna > genome.fna
    cat unplaced.scaf.fna >> genome.fna 

    """
}

