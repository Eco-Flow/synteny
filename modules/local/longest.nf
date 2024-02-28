process LONGEST {
 
    label 'process_single'
    tag "$sample_id"
    container = 'quay.io/ecoflowucl/bioseqio:perl-5.34.0'
    publishDir "$params.outdir/longest" , mode: "copy", pattern: "*${sample_id}*"

    input:
    tuple val(sample_id), path(fasta), path(gff)

    output:
    tuple val (sample_id), path( "${sample_id}.largestIsoform.fa" ), path( "${sample_id}.gff_for_jvci.gff3" ), emit: longest_proteins
    path "versions.yml", emit: versions

    script:
    """
    perl ${projectDir}/bin/ncbi_gffread_to_gene_universal.pl ${fasta} ${sample_id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Perl version: \$(perl --version | grep "version" | sed 's/.*(//g' | sed 's/[)].*//')
    END_VERSIONS
    """
}
