process CREATE_SEQIDS {
    label 'createseqids'
    tag 'createseqids'
    container = 'chriswyatt/jcvi'
             
    input:
    tuple val(sample_id), path(fasta), path(gff), val(sample_id2), path(fasta2), path(gff2)

    output:
        
    path("seqids_default"), emit: scaffolds

    script:
    """
        grep '>' $fasta  | head -n 10 | cut -c 2- | tr '\n' ',' | sed 's/.\$//' >> seqids_default
        grep '>' $fasta2 | head -n 10 | cut -c 2- | tr '\n' ',' | sed 's/.\$//' >> seqids_default
    """
}

