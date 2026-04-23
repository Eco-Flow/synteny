process RIBBON {

    label 'process_single'
    tag "$target"
    container 'quay.io/ecoflowucl/jcvi:python-3.10_last-1522_StatisticsBasic'

    input:
    path(anchors)
    path(beds)
    val(target)

    output:
    path("Ribbon.pdf"), emit: ribbonplot
    path("Ribbon.svg"), emit: ribbonsvg
    path "versions.yml", emit: versions

    script:
    """
    echo ${target} > species.csv

    echo '${params.jcvi_screen_arguments}' > my_arguments
    echo "${params.proportional_chromosomes}" > proportional_flag

    #Perl script to make simple anchor files for each comparison.
    jcvisimple.pl my_arguments

    #Perl script to write input layout for the species chosen.
    ribbon.pl

    # Shorten chromosome names in seqids and all referenced BED files
    readarray -t _bed_files < <(awk -F',' '/^ /{gsub(/^[[:space:]]+/, "", \$NF); print \$NF}' layout_all | sort -u)
    shorten_chromnames.pl seqids_karyotype_all.txt "\${_bed_files[@]}"

    # Update layout to reference the .short BED files
    awk 'BEGIN{FS=OFS=","} /^ /{sub(/[[:space:]]*\$/, ".short", \$NF)} {print}' layout_all > layout_all.short

    python ${projectDir}/bin/karyotype_svg.py seqids_karyotype_all.txt.short layout_all.short --keep-chrlabels --outfile Ribbon.pdf
    python ${projectDir}/bin/karyotype_svg.py seqids_karyotype_all.txt.short layout_all.short --keep-chrlabels --format svg --outfile Ribbon.svg

    #Versions print to file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Python version: \$(python --version  | sed 's/[^0-9]*//')
        Perl version: \$(perl --version | grep "version" | sed 's/.*(//g' | sed 's/[)].*//')
    END_VERSIONS
    """
}
