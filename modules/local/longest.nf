process LONGEST {

    label 'process_single'
    label 'process_med_memory'
    tag "$sample_id"
    container = 'biocontainers/agat:1.3.0--pl5321hdfd78af_0'

    input:
    tuple val (sample_id), path(fasta), path(gff)

    output:
    tuple val (sample_id), path( fasta ), path( "${sample_id}.longest.gff3" ), emit: longest_proteins
    path "versions.yml", emit: versions

    script:
    """
    # Run agat to find longest orf for each gene 
    agat_sp_keep_longest_isoform.pl -gff ${gff} -o ${sample_id}.longest.gff3;

    md5sum "${sample_id}.longest.gff3" > "${sample_id}.longest.gff3.md5"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Perl version: \$(perl --version | grep "version" | sed 's/.*(//g' | sed 's/[)].*//')
    END_VERSIONS
    """

}
