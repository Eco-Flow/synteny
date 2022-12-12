process SYNTENY {
    label 'syn'
    tag "${sample_id}_VS_${sample_id2}"
    publishDir "$params.outdir/Jcvi_results" , mode: "copy", pattern: "*.anchors"
    publishDir "$params.outdir/Jcvi_results" , mode: "copy", pattern: "*.pdf"
    publishDir "$params.outdir/Jcvi_results" , mode: "copy", pattern: "*.percent.similarity"
    publishDir "$params.outdir/Last" , mode: "copy", pattern: "*last.filtered"
    container = 'chriswyatt/jcvi'
             
    input:

    tuple val(sample_id), path(fasta), path(gff), val(sample_id2), path(fasta2), path(gff2)

    output:
        
    path("${sample_id}.${sample_id2}.anchors"), emit: anchors
    path("*.pdf"), emit: pdf
    path("${sample_id}.${sample_id2}.percent.similarity"), emit: percsim
    path("*.last.filtered"), emit: last

    script:
    """
        python -m jcvi.compara.catalog ortholog ${sample_id} ${sample_id2} --no_strip_names
        python -m jcvi.compara.synteny depth --histogram ${sample_id}.${sample_id2}.anchors

        cut -f 3 ${sample_id}.${sample_id2}.last.filtered  | awk '{ sum += \$1; n++ } END { if (n > 0) print sum / n; }' > ${sample_id}.${sample_id2}.percent.similarity
    """
}

