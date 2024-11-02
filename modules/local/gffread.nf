process GFFREAD {

    label 'process_single'
    tag "$sample_id"
    container = 'community.wave.seqera.io/library/gffread:0.12.7--33b95f1cfcc0e572'
    publishDir "$params.outdir/output_data/gffread" , mode: "${params.publish_dir_mode}"

    input:
    tuple val(sample_id), path(fasta), path(gff)

    output:
    tuple val(sample_id), path( "${sample_id}.nucl.fa" ), path( "${sample_id}.gff_for_jvci.gff3" ), emit: proteins
    path( "${sample_id}.gff_for_jvci.gff3" ), emit: gff
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    perl ${projectDir}/bin/gffread_unzip.pl ${sample_id} ${fasta} ${gff} ${args}
    md5sum "${sample_id}.nucl.fa" > "${sample_id}.nucl.fa.md5"
    md5sum "${sample_id}.gff_for_jvci.gff3" > "${sample_id}.gff_for_jvci.gff3.md5"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread version: \$(gffread --version)
        Perl version: \$(perl --version | grep "version" | sed 's/.*(//g' | sed 's/[)].*//')
    END_VERSIONS
    """
}
