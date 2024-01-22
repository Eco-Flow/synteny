process CHROMOPAINT {

    label 'process_single'
    tag "${sample1} with ${sample2}"
    publishDir "$params.outdir/Chromosome_plots" , mode: "copy"
    container = 'ecoflowucl/jcvi:python-3.10_last-1522'

    input:
    tuple path (hex), val(sample1), val(sample2), path(anchors)
    path ('*')

    output:
    path("*.pdf"), emit: pdf
    path "versions.yml", emit: versions

    script:
    """
    echo '${anchors}' | rev | cut -d'/' -f 1 | rev > Name
    A="\$(cut -d'.' -f1 Name)"
    B="\$(cut -d'.' -f2 Name)"
    anchor.pl \$A.bed \$B.bed ${anchors} ${hex}
    python -m jcvi.graphics.chromosome Chromopaint.txt colour.idmap
    mv Chromopaint.pdf "\$A.\$B.chromo.pdf"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Python version: \$(python --version  | sed 's/[^0-9]*//')
        JCVI \$(pip show jcvi | grep "Version:")
        Perl version: \$(perl --version | grep "version" | sed 's/.*(//g' | sed 's/[)].*//')
    END_VERSIONS
    """
}
