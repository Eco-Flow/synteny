process SYNTENY {

    label 'process_single'
    tag "${sample_id}_VS_${sample_id2}"
    publishDir "$params.outdir/output_data/anchors" , mode: "copy", pattern: "*.anchors"
    publishDir "$params.outdir/figures/dotplots" , mode: "copy", pattern: "*dotplot.pdf"
    publishDir "$params.outdir/figures/depth_plots" , mode: "copy", pattern: "*depth.pdf"
    publishDir "$params.outdir/figures/karyotype_plots" , mode: "copy", pattern: "*.karyotype.pdf"
    publishDir "$params.outdir/output_data/last" , mode: "copy", pattern: "*last.filtered"
    container = 'ecoflowucl/jcvi:python-3.10_last-1522'

    input:
    tuple val(sample_id), path(fasta), path(gff), val(sample_id2), path(fasta2), path(gff2)

    output:
    path("${sample_id}.${sample_id2}.lifted.anchors"), emit: anchors
    tuple path("${sample_id}.${sample_id2}.dotplot.pdf"), path("${sample_id}.${sample_id2}.depth.pdf"), path("${sample_id}.${sample_id2}.karyotype.pdf"), emit: pdf
    path("${sample_id}.${sample_id2}.percent.similarity"), emit: percsim
    path("${sample_id}.${sample_id2}.last.filtered"), emit: last
    path "versions.yml", emit: versions

    script:
    """
    #main code block
    python -m jcvi.compara.catalog ortholog ${sample_id} ${sample_id2} --no_strip_names
    python -m jcvi.compara.synteny depth --histogram ${sample_id}.${sample_id2}.anchors
    cut -f 3 ${sample_id}.${sample_id2}.last.filtered  | awk '{ sum += \$1; n++ } END { if (n > 0) print sum / n; }' > ${sample_id}.${sample_id2}.percent.similarity
    python -m jcvi.compara.synteny screen --minspan=30 --simple ${sample_id}.${sample_id2}.anchors ${sample_id}.${sample_id2}.anchors.new
    

    #Echo a layout format based on the input file name
    echo "# y, xstart, xend, rotation, color, label, va,  bed" > layout
    echo ".6,     .1,    .8,       0,      , ${sample_id}, top, ${sample_id}.bed" >> layout
    echo ".4,     .1,    .8,       0,      , ${sample_id2}, top, ${sample_id2}.bed" >> layout
    echo "# edges" >> layout
    echo "e, 0, 1, ${sample_id}.${sample_id2}.anchors.simple" >> layout
    

    perl ${projectDir}/bin/syntenous_chromosomes.pl ${sample_id}.bed ${sample_id2}.bed ${sample_id}.${sample_id2}.anchors.new
    python -m jcvi.graphics.karyotype seqids_karyotype.txt layout
    mv karyotype.pdf ${sample_id}.${sample_id2}.karyotype.pdf
    mv ${sample_id}.${sample_id2}.pdf ${sample_id}.${sample_id2}.dotplot.pdf


    #md5 sums for unit testing
    md5sum "${sample_id}.${sample_id2}.anchors" > "${sample_id}.${sample_id2}.anchors.md5"
    md5sum "${sample_id}.${sample_id2}.percent.similarity" > "${sample_id}.${sample_id2}.percent.similarity.md5"
    md5sum "${sample_id}.${sample_id2}.last.filtered" > "${sample_id}.${sample_id2}.last.filtered.md5"

    #Versions print to file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Python version: \$(python --version  | sed 's/[^0-9]*//')
        Perl version: \$(perl --version | grep "version" | sed 's/.*(//g' | sed 's/[)].*//')
    END_VERSIONS
    """
}
