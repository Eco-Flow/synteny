process RIBBON {

    label 'process_single'
    tag "$target"
    container = 'ecoflowucl/jcvi:python-3.10_last-1522_StatisticsBasic' 

    input:
    path(anchors)
    path(beds)
    val(target)

    output:
    path("Ribbon.pdf"), emit: ribbonplot
    path "versions.yml", emit: versions

    script:
    """
    echo ${target} > species.csv

    echo '${params.jcvi_screen_arguments}' > my_arguments

    #Perl script to make simple anchor files for each comparison.
    jcvisimple.pl my_arguments

    #Perl script to wrote input layout for the species chosen.
    ribbon.pl

    python -m jcvi.graphics.karyotype seqids_karyotype_all.txt layout_all

    mv karyotype.pdf Ribbon.pdf

    #Versions print to file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Python version: \$(python --version  | sed 's/[^0-9]*//')
        Perl version: \$(perl --version | grep "version" | sed 's/.*(//g' | sed 's/[)].*//')
    END_VERSIONS
    """
}
