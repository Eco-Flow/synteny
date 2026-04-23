process SYNTENY {

    stageInMode = 'copy'
    label 'process_single'
    label 'process_long'
    label 'process_med_memory'
    tag "${sample_id}_VS_${sample_id2}"
    publishDir "$params.outdir/output_data/anchors" , mode: "${params.publish_dir_mode}", pattern: "*.anchors"
    publishDir "$params.outdir/figures/dotplots" , mode: "${params.publish_dir_mode}", pattern: "*dotplot.pdf"
    publishDir "$params.outdir/figures/depth_plots" , mode: "${params.publish_dir_mode}", pattern: "*depth.pdf"
    publishDir "$params.outdir/figures/karyotype_plots" , mode: "${params.publish_dir_mode}", pattern: "*.karyotype.pdf"
    publishDir "$params.outdir/figures/karyotype_plots" , mode: "${params.publish_dir_mode}", pattern: "*.karyotype.flipped.pdf"
    publishDir "$params.outdir/figures/karyotype_plots" , mode: "${params.publish_dir_mode}", pattern: "*.karyotype.svg"
    publishDir "$params.outdir/figures/karyotype_plots" , mode: "${params.publish_dir_mode}", pattern: "*.karyotype.flipped.svg"
    publishDir "$params.outdir/output_data/last" , mode: "${params.publish_dir_mode}", pattern: "*last.filtered"
    publishDir "$params.outdir/output_data/percentidentity" , mode: "${params.publish_dir_mode}", pattern: "*.percent.similarity"
    container = 'quay.io/ecoflowucl/jcvi:python-3.10_last-1522_StatisticsBasic'

    input:
    tuple val(sample_id), path(fasta), path(gff), val(sample_id2), path(fasta2), path(gff2)

    output:
    path("${sample_id}.${sample_id2}.lifted.anchors"), emit: anchors
    path("${sample_id}.${sample_id2}.anchors"), emit: anchors_notlifted
    tuple path("${sample_id}.${sample_id2}.dotplot.pdf"), path("${sample_id}.${sample_id2}.depth.pdf"), path("${sample_id}.${sample_id2}.karyotype.pdf"), emit: pdf
    path("${sample_id}.${sample_id2}.karyotype.flipped.pdf"), emit: pdf_flipped
    path("${sample_id}.${sample_id2}.karyotype.svg"), emit: svg
    path("${sample_id}.${sample_id2}.karyotype.flipped.svg"), emit: svg_flipped
    path("${sample_id}.${sample_id2}.percent.similarity"), emit: percsim
    path("${sample_id}.${sample_id2}.last.filtered"), emit: last
    path("${sample_id}.${sample_id2}.last"), emit: unfilteredlast
    path "versions.yml", emit: versions

    script:
    def prop_chrom = params.proportional_chromosomes
    """
    #main code block
    python -m jcvi.compara.catalog ortholog ${sample_id} ${sample_id2} ${params.jcvi_ortholog_arguments}
    python -m jcvi.compara.synteny depth --histogram ${sample_id}.${sample_id2}.anchors
    cut -f 3 ${sample_id}.${sample_id2}.last.filtered  | awk '{ sum += \$1; n++ } END { if (n > 0) print sum / n; }' > ${sample_id}.${sample_id2}.percent.similarity
    python -m jcvi.compara.synteny screen ${params.jcvi_screen_arguments} --simple ${sample_id}.${sample_id2}.anchors ${sample_id}.${sample_id2}.anchors.new
    mv ${sample_id}.${sample_id2}.pdf ${sample_id}.${sample_id2}.dotplot.pdf

    # Determine chromosome order and produce the flipped BED for species 2
    perl ${projectDir}/bin/syntenous_chromosomes.pl ${sample_id}.bed ${sample_id2}.bed ${sample_id}.${sample_id2}.anchors.new

    # Shorten chromosome names to last N chars (min 3) ensuring uniqueness across both species
    perl ${projectDir}/bin/shorten_chromnames.pl seqids_karyotype.txt ${sample_id}.bed ${sample_id2}.bed ${sample_id2}.bed.flipped.bed

    # Compute xend values — proportional to genome size (xstart=.2, xend_max=.8), or fixed at .8
    if [ "${prop_chrom}" = "true" ]; then
        size1=\$(awk '{if (\$3 > lens[\$1]) lens[\$1]=\$3} END {total=0; for (chr in lens) total+=lens[chr]; print total}' ${sample_id}.bed)
        size2=\$(awk '{if (\$3 > lens[\$1]) lens[\$1]=\$3} END {total=0; for (chr in lens) total+=lens[chr]; print total}' ${sample_id2}.bed)
        if [ "\$size1" -gt "\$size2" ]; then maxsize="\$size1"; else maxsize="\$size2"; fi
        xend1=\$(echo "\$size1 \$maxsize" | awk '{printf "%.3f", 0.2 + (\$1/\$2) * 0.6}')
        xend2=\$(echo "\$size2 \$maxsize" | awk '{printf "%.3f", 0.2 + (\$1/\$2) * 0.6}')
    else
        xend1=".8"
        xend2=".8"
    fi

    # Layout using shortened chromosome names (xstart=.2 leaves room for species labels)
    echo "# y, xstart, xend, rotation, color, label, va,  bed" > layout
    echo ".6,     .2,    \${xend1},       0,      , ${sample_id}, top, ${sample_id}.bed.short" >> layout
    echo ".4,     .2,    \${xend2},       0,      , ${sample_id2}, top, ${sample_id2}.bed.short" >> layout
    echo "# edges" >> layout
    echo "e, 0, 1, ${sample_id}.${sample_id2}.anchors.simple" >> layout

    python ${projectDir}/bin/karyotype_svg.py seqids_karyotype.txt.short layout --keep-chrlabels --outfile ${sample_id}.${sample_id2}.karyotype.pdf
    python ${projectDir}/bin/karyotype_svg.py seqids_karyotype.txt.short layout --keep-chrlabels --format svg --outfile ${sample_id}.${sample_id2}.karyotype.svg

    echo "# y, xstart, xend, rotation, color, label, va,  bed" > layout_flip
    echo ".6,     .2,    \${xend1},       0,      , ${sample_id}, top, ${sample_id}.bed.short" >> layout_flip
    echo ".4,     .2,    \${xend2},       0,      , ${sample_id2}, top, ${sample_id2}.bed.flipped.bed.short" >> layout_flip
    echo "# edges" >> layout_flip
    echo "e, 0, 1, ${sample_id}.${sample_id2}.anchors.simple" >> layout_flip

    python ${projectDir}/bin/karyotype_svg.py seqids_karyotype.txt.short layout_flip --keep-chrlabels --outfile ${sample_id}.${sample_id2}.karyotype.flipped.pdf
    python ${projectDir}/bin/karyotype_svg.py seqids_karyotype.txt.short layout_flip --keep-chrlabels --format svg --outfile ${sample_id}.${sample_id2}.karyotype.flipped.svg


    #md5 sums for unit testing
    md5sum "${sample_id}.${sample_id2}.anchors" > "${sample_id}.${sample_id2}.anchors.md5"
    md5sum "${sample_id}.${sample_id2}.percent.similarity" > "${sample_id}.${sample_id2}.percent.similarity.md5"
    md5sum "${sample_id}.${sample_id2}.last.filtered" > "${sample_id}.${sample_id2}.last.filtered.md5"

    #Versions print to file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Python version: \$(python --version  | sed 's/[^0-9]*//')
        Perl version: \$(perl --version | grep "version" | sed 's/.*(//g' | sed 's/[)].*//')
        Last version: \$(lastal --version | sed 's/[^0-9]*//')
    END_VERSIONS
    """
}
