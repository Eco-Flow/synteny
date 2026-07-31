//
// Download/resolve genomes+annotations from a samplesheet CSV (2-col NCBI accession /
// 3-col local path) and validate them. Shared by the default pairwise workflow
// (main.nf) and the ALGO subworkflow (workflows/algo.nf), which otherwise each
// duplicated this exact CSV-parsing / DOWNLOAD_NCBI / FASTAVALIDATOR / SEQKIT_STATS
// logic.
//

include { DOWNLOAD_NCBI } from '../../../modules/local/download_ncbi.nf'
include { FASTAVALIDATOR } from '../../../modules/nf-core/fastavalidator/main'
include { SEQKIT_STATS } from '../../../modules/nf-core/seqkit/stats/main'
include { SAMTOOLS_FAIDX } from '../../../modules/nf-core/samtools/faidx/main'

workflow GENOME_ACQUISITION {

    take:
    input_csv // channel: path to the samplesheet CSV (2-col NCBI accession / 3-col local path)

    main:
    ch_versions = Channel.empty()

    input_csv
        .splitCsv()
        .branch {
            ncbi: it.size() == 2
            path: it.size() == 3
        }
        .set { input_type }

    DOWNLOAD_NCBI ( input_type.ncbi )
    ch_versions = ch_versions.mix(DOWNLOAD_NCBI.out.versions.first())

    // Checks if paths are S3 objects if not ensures absolute paths are used for user inputted fasta and gff files
    input_type.path.map{ name, fasta, gff ->
        def full_fasta = fasta =~ /^s3:\/\// ? fasta : new File(fasta).getAbsolutePath()
        def full_gff   = gff   =~ /^s3:\/\// ? gff   : new File(gff).getAbsolutePath()
        [name, full_fasta, full_gff]
    }.set { local_annotations }

    genome_annotations = DOWNLOAD_NCBI.out.genome.mix(local_annotations)

    fasta_inputs = genome_annotations.map { name, fasta, gff -> [[id: name], fasta] }

    FASTAVALIDATOR ( fasta_inputs )
    ch_versions = ch_versions.mix(FASTAVALIDATOR.out.versions.first())

    // Manipulate successful and error logs of fasta validator to be saved into output directory
    FASTAVALIDATOR.out.success_log.map { meta, logfile -> [meta.id, logfile] }
        .collectFile(name: { it[0] }, storeDir: "${params.outdir}/input_validation/fasta_validator/successful")
    FASTAVALIDATOR.out.error_log.map { meta, logfile -> [meta.id, logfile] }
        .collectFile(name: { it[0] }, storeDir: "${params.outdir}/input_validation/fasta_validator/error")

    SEQKIT_STATS ( fasta_inputs )
    ch_versions = ch_versions.mix(SEQKIT_STATS.out.versions.first())

    // Manipulate seqkit_stats tsv to be saved into output directory
    SEQKIT_STATS.out.stats.map { meta, tsv -> [meta.id, tsv] }
        .collectFile(name: { it[0] }, storeDir: "${params.outdir}/input_validation/seqkit_stats")

    if (params.chromo_lengths) {
        SAMTOOLS_FAIDX ( fasta_inputs, [[],[]], true )
        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions.first())

        // Save chromosome lengths to output directory
        SAMTOOLS_FAIDX.out.sizes.map { meta, sizes -> [meta.id, sizes] }
            .collectFile(name: { it[0] }, storeDir: "${params.outdir}/output_data/chromosome_lengths")
    }

    emit:
    genome_annotations = genome_annotations                                     // channel: [ id, fasta, gff ]
    genome_fastas      = genome_annotations.map { id, fasta, gff -> [id, fasta] } // channel: [ id, fasta ]
    versions           = ch_versions
}
