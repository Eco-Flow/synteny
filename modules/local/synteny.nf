process SYNTENY {

    label 'process_single'
    tag "${sample_id}_VS_${sample_id2}"
    publishDir "$params.outdir/Jcvi_results" , mode: "copy", pattern: "*.anchors"
    publishDir "$params.outdir/Jcvi_results" , mode: "copy", pattern: "*.pdf"
    publishDir "$params.outdir/Jcvi_results" , mode: "copy", pattern: "*.percent.similarity"
    publishDir "$params.outdir/Karyotype_results" , mode: "copy", pattern: "*.karyotype.pdf"
    publishDir "$params.outdir/Last" , mode: "copy", pattern: "*last.filtered"
    container = 'ecoflowucl/jcvi:python-3.10_last-1522'

    input:
    tuple val(sample_id), path(fasta), path(gff), val(sample_id2), path(fasta2), path(gff2)

    output:
    path("${sample_id}.${sample_id2}.anchors"), emit: anchors
    tuple path("${sample_id}.${sample_id2}.pdf"), path("${sample_id}.${sample_id2}.depth.pdf"), path("${sample_id}.${sample_id2}.karyotype.pdf"), emit: pdf
    path("${sample_id}.${sample_id2}.percent.similarity"), emit: percsim
    path("${sample_id}.${sample_id2}.last.filtered"), emit: last
    path "versions.yml", emit: versions

    script:
    """
    python -m jcvi.compara.catalog ortholog ${sample_id} ${sample_id2} --no_strip_names
    python -m jcvi.compara.synteny depth --histogram ${sample_id}.${sample_id2}.anchors
    cut -f 3 ${sample_id}.${sample_id2}.last.filtered  | awk '{ sum += \$1; n++ } END { if (n > 0) print sum / n; }' > ${sample_id}.${sample_id2}.percent.similarity
    python -m jcvi.compara.synteny screen --minspan=30 --simple ${sample_id}.${sample_id2}.anchors ${sample_id}.${sample_id2}.anchors.new
    echo "# y, xstart, xend, rotation, color, label, va,  bed\n.6,     .1,    .8,       0,      , ${sample_id}, top, ${sample_id}.bed\n.4,     .1,    .8,       0,      , ${sample_id2}, top, ${sample_id2}.bed\n# edges\ne, 0, 1, ${sample_id}.${sample_id2}.anchors.simple" > layout
    syntenous_chromosomes.pl ${sample_id}.bed ${sample_id2}.bed ${sample_id}.${sample_id2}.anchors.new
    python -m jcvi.graphics.karyotype seqids_karyotype.txt layout
    mv karyotype.pdf ${sample_id}.${sample_id2}.karyotype.pdf

    md5sum "${sample_id}.${sample_id2}.anchors" > "${sample_id}.${sample_id2}.anchors.md5"
    md5sum "${sample_id}.${sample_id2}.percent.similarity" > "${sample_id}.${sample_id2}.percent.similarity.md5"
    md5sum "${sample_id}.${sample_id2}.last.filtered" > "${sample_id}.${sample_id2}.last.filtered.md5"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Python version: \$(python --version  | sed 's/[^0-9]*//')
        JCVI \$(pip show jcvi | grep "Version:")
    END_VERSIONS
    """
}
