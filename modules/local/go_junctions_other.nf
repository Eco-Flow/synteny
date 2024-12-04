process GO_JUNCTIONS_OTHER {

    label 'process_single'
    tag "Run $species with $cutoff percent cutoff"
    container = 'quay.io/ecoflowucl/chopgo:r-4.3.2_python-3.10_perl-5.38'
    publishDir "$params.outdir/output_data/go_results/individual/other" , mode: "${params.publish_dir_mode}", pattern:"*.tab"
    publishDir "$params.outdir/figures/go_results/individual/other" , mode: "${params.publish_dir_mode}", pattern:"*.pdf"
    publishDir "$params.outdir/output_data/go_results/individual/other/input_txt" , mode: "${params.publish_dir_mode}", pattern:"*.txt"

    input:
    tuple path(go, stageAs: 'Go'), val(species), path(junction_score), val (cutoff)
    path(beds)

    output:
    path( "*.pdf" ), emit: go_pdf, optional:true
    path( "*ALL.tab" ), emit: go_pvals, optional:true
    tuple val(cutoff), path( "*results_ALL.tab" ), emit: go_table, optional:true
    path( "*txt" ), emit: go_data, optional:true
    path "versions.yml", emit: versions

    script:
    """
    # ${task.process}

    perl "${projectDir}/bin/Junction_go.pl" ${cutoff} ${species}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R version: \$(R --version | grep "R version" | sed 's/[(].*//' | sed 's/ //g' | sed 's/[^0-9]*//')
        GO version: \$(Rscript -e "as.data.frame(installed.packages())[ ,c(1,3)]" | grep topGO | sed 's/[^0-9]*//')
        Perl version: \$(perl --version | grep "version" | sed 's/.*(//g' | sed 's/[)].*//')
    END_VERSIONS
    """
    
}
