process LONGEST {

    label 'process_single'
    label 'process_med_memory'
    tag "$sample_id"
    container = 'biocontainers/agat:1.3.0--pl5321hdfd78af_0'
    publishDir "$params.outdir/output_data/longest" , mode: "${params.publish_dir_mode}", pattern:"*.txt"

    input:
    tuple val (sample_id), path(fasta), path(gff)

    output:
    tuple val (sample_id), path( fasta ), path( "${sample_id}.longest.gff3" ), emit: longest_proteins
    tuple val (sample_id), path( "${sample_id}.stat.original.txt" ), emit: agat_summary_original
    tuple val (sample_id), path( "${sample_id}.stat.long.txt" ), emit: agat_summary_longest
    path "versions.yml", emit: versions

    script:
    """
    # Run agat to find longest orf for each gene 
    agat_sp_keep_longest_isoform.pl -gff ${gff} -o ${sample_id}.longest.gff3
    
    # Run a few summarisation scripts to report the actual genes being considered.
    agat_sp_functional_statistics.pl --gff ${gff} -o ${sample_id}.stat.original.txt
    agat_sp_functional_statistics.pl --gff ${sample_id}.longest.gff3 -o ${sample_id}.stat.long.txt
    
    md5sum "${sample_id}.longest.gff3" > "${sample_id}.longest.gff3.md5"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Perl version: \$(perl --version | grep "version" | sed 's/.*(//g' | sed 's/[)].*//')
    END_VERSIONS
    """

}
