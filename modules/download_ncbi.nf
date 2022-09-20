process DOWNLOAD_NCBI {
    label 'download'
    tag "$sample_id via $accension_id" 
    //errorStrategy { task.exitStatus == 2 ? 'ignore' }        
    errorStrategy 'ignore'

    input:
        tuple val(sample_id), val(accension_id)

    output:    
        tuple val(sample_id), path("${sample_id}.genome.fna"), path("${sample_id}.genomic.gff"), emit: genome

    script:
    """

    #Get a genome and GFF assembly from NCBI using their datasets scripts

        datasets download genome accession ${accension_id} --exclude-rna --exclude-genomic-cds --dehydrated --filename ncbi_dataset.zip
        unzip ncbi_dataset.zip 
        datasets rehydrate 

        ls ncbi_dataset/data/${accension_id}/chr*.fna > tmp.fna
        #Remove lines that are not in contiguous chromosomes.
        grep -v unlocalized tmp.fna > tmp.noloc.fna
        while read line; do cat \$line;  done < tmp.noloc.fna > ${sample_id}.init.genome.fna
        #Remove lines with ? for strand
        awk -F'\t' '!(\$7 == "\?")' ${sample_id}.init.genome.fna > ${sample_id}.genome.fna

        cat ncbi_dataset/data/${accension_id}/genomic.gff > ${sample_id}.genomic.gff

    """
}

