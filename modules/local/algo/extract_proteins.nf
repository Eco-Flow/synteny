process EXTRACT_PROTEINS {

    label 'process_single'
    tag "$sample_id"
    container 'community.wave.seqera.io/library/gffread:0.12.7--33b95f1cfcc0e572'
    publishDir "$params.outdir/algo/proteomes" , mode: "${params.publish_dir_mode}"

    input:
    tuple val(sample_id), path(genome_input, stageAs: 'staged_genome_input'), path(gff)

    output:
    tuple val(sample_id), path("${sample_id}.fa"), emit: proteins
    path "versions.yml", emit: versions

    script:
    // A user's genome fasta can legitimately be named "<sample_id>.fa" -- the same
    // name this process's own output uses. Nextflow stages `path` inputs under their
    // original filename by default, so without the explicit stageAs above, that
    // name would collide with the output and gffread would end up writing its
    // result over the very file it's reading from (truncating the source genome).
    // Detected by content (`gzip -t`), not by filename, since stageAs above already
    // discards whatever suffix the original file had.
    """
    if gzip -t staged_genome_input 2>/dev/null; then
        zcat staged_genome_input > genome_for_extraction.fa
    else
        cp staged_genome_input genome_for_extraction.fa
    fi

    # -J: only emit mRNAs with a complete CDS (start + stop codon present, no
    # in-frame stop). Without it, a premature stop in a gene model -- routine
    # in real (e.g. BRAKER) annotations -- gets translated to a literal "."
    # in the protein, which diamond's parser inside OrthoFinder rejects with
    # "Error: Invalid character in sequence: '.'" and aborts the whole run.
    gffread -J -y ${sample_id}.fa -g genome_for_extraction.fa ${gff}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread version: \$(gffread --version)
    END_VERSIONS
    """
}
