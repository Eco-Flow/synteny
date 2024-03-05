process LONGEST {
 
    label 'process_single'
    label 'process_long'
    tag "$sample_id"
    container = 'biocontainers/agat:0.8.0--pl5262hdfd78af_0'

    input:
    tuple val (sample_id), path(fasta), path(gff)

    output:
    tuple val (sample_id), path( "${fasta}" ), path( "${sample_id}.longest.gff3" ), emit: longest_proteins
    path "versions.yml", emit: versions

    script:
    """
    #if gzipped, un gzip and then run agat to find longest orf for each gene
    if [[ (${gff} =~ .gff.gz) || (${gff} =~ .gff3.gz) ]]; then 
        gunzip -c ${gff} > ${gff}.norm; agat_sp_keep_longest_isoform.pl -gff ${gff}.norm -o ${sample_id}.longest.gff3;
    else 
        agat_sp_keep_longest_isoform.pl -gff ${gff} -o ${sample_id}.longest.gff3;
    fi

    #un gzip the genome so it can be used in gffread
    if [[ ${fasta} =~ .fna.gz ]]; then 
        gunzip -c ${fasta};
    fi

    md5sum "${sample_id}.longest.gff3" > "${sample_id}.longest.gff3.md5"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Perl version: \$(perl --version | grep "version" | sed 's/.*(//g' | sed 's/[)].*//')
    END_VERSIONS
    """

}
