process DOWNLOAD_NCBI {

    label 'download'
    tag "$sample_id via $accension_id"
    container "${ params.architecture == 'arm' ? 'ecoflowucl/ncbi_download:v16.1.2-arm64' : 'ecoflowucl/ncbi_download:v16.1.2-amd64' }"

    input:
    tuple val(sample_id), val(accension_id)

    output:
    tuple val(sample_id), path("${sample_id}.genome.fna"), path("${sample_id}.genomic.gff"), emit: genome

    script:
    """
    #Get a genome and GFF assembly from NCBI using their datasets scripts
    datasets download genome accession ${accension_id} --include genome,gff3
    unzip ncbi_dataset.zip
    
    if [ -f ncbi_dataset/data/${accension_id}/chr*.fna ]; then
        cat ncbi_dataset/data/${accension_id}/chr*.fna > ${sample_id}.genome.fna
    elif [ -f ncbi_dataset/data/${accension_id}/unplaced.scaf.fna ]; then
        cat ncbi_dataset/data/${accension_id}/unplaced.scaf.fna >> ${sample_id}.genome.fna
    elif [ -f ncbi_dataset/data/${accension_id}/${accension_id}*_genomic.fna ]; then
        cat ncbi_dataset/data/${accension_id}/${accension_id}*_genomic.fna >> ${sample_id}.genome.fna
    fi
    
    cat ncbi_dataset/data/${accension_id}/genomic.gff > ${sample_id}.genomic.gff
    """
}
