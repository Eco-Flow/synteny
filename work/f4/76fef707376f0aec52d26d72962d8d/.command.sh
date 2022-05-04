#!/bin/bash -ue
#Run the basic transformation of gff to bed and fasta to cds conversions.

    python -m jcvi.formats.gff bed --type=mRNA --key=ID B_impatiens.gff3 -o B_impatiens.bed
    python -m jcvi.formats.fasta format B_impatiens.prot.fa B_impatiens.cds

    #Make a default seqids file

    if [ data/default1 = 'data/default1' ]; then 
        grep '>' B_impatiens.prot.fa  | head -n 10 | cut -c 2- | tr '
' ',' | sed 's/.$//' > seqids_default       
    else
        cat default1 > seqids_default     
    fi
