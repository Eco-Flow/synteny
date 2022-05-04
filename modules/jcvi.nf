process JCVI {
    label 'jcvi'
    tag 'jcvi'
    container = 'chriswyatt/jcvi'
             
    input:

    tuple val(sample_id), path(fasta), path(gff)
    path(seqids)

    output:
        
    tuple val(sample_id), path( "${sample_id}.cds" ), path( "${sample_id}.bed" ) , emit: reformatted
    path("seqids_default"), emit: seqids_out

    script:
    """
    #Run the basic transformation of gff to bed and fasta to cds conversions.

    python -m jcvi.formats.gff bed --type=mRNA --key=ID ${gff} -o ${sample_id}.bed
    python -m jcvi.formats.fasta format ${fasta} ${sample_id}.cds

    #Make a default seqids file

    if [ $params.seqids = 'data/default1' ]; then 
        grep '>' $fasta  | head -n 10 | cut -c 2- | tr '\n' ',' | sed 's/.\$//' > seqids_default       
    else
        cat $seqids > seqids_default     
    fi
    """
}

