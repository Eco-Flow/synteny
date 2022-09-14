process SYNTENY {
    label 'syn'
    tag "${sample_id}_VS_${sample_id2}"
    publishDir "$params.outdir/Jcvi_results" , mode: "copy"
    container = 'chriswyatt/jcvi'
             
    input:

    tuple val(sample_id), path(fasta), path(gff), val(sample_id2), path(fasta2), path(gff2)

    output:
        
    path("${sample_id}.${sample_id2}.anchors"), emit: anchors
    path("*.pdf"), emit: pdf
    path("${sample_id}.${sample_id2}.percent.similarity"), emit: percsim

    script:
    """
        python -m jcvi.compara.catalog ortholog ${sample_id} ${sample_id2} --no_strip_names
        python -m jcvi.compara.synteny depth --histogram ${sample_id}.${sample_id2}.anchors

        cut -f 3 work/c2/5a7eb6cb839c0a496008025bb4dfc5/Vespula_vulgaris.Dolichovespula_media.last.filtered  | awk '{ sum += \$1; n++ } END { if (n > 0) print sum / n; }' > ${sample_id}.${sample_id2}.percent.similarity
    """
}

