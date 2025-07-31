/*
========================================================================================
    SUBWORKFLOW: TT-seq processing workflow
========================================================================================
*/

include { SAMPLESHEET_CHECK      } from '../modules/check_samplesheet'
include { FASTQC                 } from '../modules/qc'
include { MULTIQC                } from '../modules/qc'
include { STAR_ALIGN             } from '../modules/align'
include { STAR_ALIGN_SPIKE       } from '../modules/align'
include { SAMTOOLS_SORT          } from '../modules/samtools'
include { SAMTOOLS_INDEX         } from '../modules/samtools'
include { SAMTOOLS_STATS         } from '../modules/samtools'
include { PICARD_MARKDUPLICATES  } from '../modules/samtools'
include { BAMCOVERAGE_STRANDED   } from '../modules/coverage'
include { EXTRACT_SPIKE_COUNTS   } from '../modules/spike_in_norm'
include { CALCULATE_SIZE_FACTORS } from '../modules/spike_in_norm'

/*
========================================================================================
    SUBWORKFLOW: TT-seq analysis workflow
========================================================================================
*/

workflow TTSEQ {
    
    take:
    ch_samplesheet // channel: samplesheet read in from --sample_sheet

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // SUBWORKFLOW: Read and validate samplesheet
    //
    SAMPLESHEET_CHECK (
        ch_samplesheet
    )
    ch_versions = ch_versions.mix(SAMPLESHEET_CHECK.out.versions)

    //
    // Create reads channel from samplesheet
    //
    SAMPLESHEET_CHECK.out.csv
        .splitCsv(header:true, sep:',')
        .map { create_fastq_channel(it) }
        .set { ch_reads }

    //
    // MODULE: Run FastQC
    //
    if (!params.skip_fastqc) {
        FASTQC (
            ch_reads
        )
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
        ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    }

    //
    // MODULE: Align to primary genome with STAR
    //
    ch_genome_index = Channel.value([
        file(params.star_index, checkIfExists: true)
    ])
    ch_gtf = params.gtf ? Channel.value(file(params.gtf, checkIfExists: true)) : Channel.value([])
    
    STAR_ALIGN (
        ch_reads,
        ch_genome_index,
        ch_gtf
    )
    ch_versions = ch_versions.mix(STAR_ALIGN.out.versions.first()) 
    ch_multiqc_files = ch_multiqc_files.mix(STAR_ALIGN.out.log_final.collect{it[1]})

    //
    // MODULE: Conditional spike-in alignment
    //
    if (params.spike_in_genome) {
        ch_spike_index = Channel.value([
            file(params.spike_in_index, checkIfExists: true)
        ])
        ch_spike_gtf = params.spike_in_gtf ? Channel.value(file(params.spike_in_gtf, checkIfExists: true)) : Channel.value([])
        
        STAR_ALIGN_SPIKE (
            ch_reads,
            ch_spike_index, 
            ch_spike_gtf
        )
        ch_versions = ch_versions.mix(STAR_ALIGN_SPIKE.out.versions.first())
        
        //
        // MODULE: Extract spike-in counts
        //
        ch_bam_pairs = STAR_ALIGN.out.bam_sorted
            .join(STAR_ALIGN_SPIKE.out.bam, by: [0])
        
        EXTRACT_SPIKE_COUNTS (
            ch_bam_pairs
        )
        ch_versions = ch_versions.mix(EXTRACT_SPIKE_COUNTS.out.versions.first())
        
        //
        // MODULE: Calculate size factors from spike-in
        //
        CALCULATE_SIZE_FACTORS (
            EXTRACT_SPIKE_COUNTS.out.counts.collect{it[1]},
            EXTRACT_SPIKE_COUNTS.out.counts.collect{it[0].id + "_spike_stats.txt"}
        )
        ch_versions = ch_versions.mix(CALCULATE_SIZE_FACTORS.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(CALCULATE_SIZE_FACTORS.out.summary)

        //
        // Join BAM files with scaling factors
        //
        ch_scaling_factors = CALCULATE_SIZE_FACTORS.out.factors
            .splitCsv(header:true, sep:'\\t')
            .map { row -> [row.sample, row.scaling_factor] }

        ch_bam_with_factors = STAR_ALIGN.out.bam_sorted
            .map { meta, bam -> [meta.id, meta, bam] }
            .join(ch_scaling_factors, by: [0])
            .map { sample_id, meta, bam, scaling_factor -> [meta, bam, scaling_factor] }
    } else {
        //
        // No spike-in normalization - use BAM files without scaling factors
        //
        ch_bam_with_factors = STAR_ALIGN.out.bam_sorted
            .map { meta, bam -> [meta, bam, null] }
    }

    //
    // MODULE: Index BAM files
    //
    ch_bam_for_index = ch_bam_with_factors.map { meta, bam, scaling -> [meta, bam] }
    
    SAMTOOLS_INDEX (
        ch_bam_for_index
    )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    //
    // MODULE: Mark duplicates with Picard
    //
    PICARD_MARKDUPLICATES (
        ch_bam_for_index
    )
    ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(PICARD_MARKDUPLICATES.out.metrics.collect{it[1]})

    //
    // MODULE: Generate alignment statistics
    //
    ch_fasta = params.fasta ? Channel.value(file(params.fasta, checkIfExists: true)) : Channel.value([])
    
    SAMTOOLS_STATS (
        PICARD_MARKDUPLICATES.out.bam.join(PICARD_MARKDUPLICATES.out.bai, by: [0]),
        ch_fasta
    )
    ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS.out.stats.collect{it[1]})

    //
    // MODULE: Generate strand-specific BigWig files
    //
    ch_final_bams = PICARD_MARKDUPLICATES.out.bam
        .join(PICARD_MARKDUPLICATES.out.bai, by: [0])
        .join(
            ch_bam_with_factors.map { meta, bam, scaling -> [meta, scaling] }, 
            by: [0]
        )
        .map { meta, bam, bai, scaling -> [meta, bam, bai, scaling] }

    BAMCOVERAGE_STRANDED (
        ch_final_bams
    )
    ch_versions = ch_versions.mix(BAMCOVERAGE_STRANDED.out.versions.first())

    //
    // MODULE: MultiQC
    //
    if (!params.skip_multiqc) {
        MULTIQC (
            ch_multiqc_files.collect()
        )
        ch_versions = ch_versions.mix(MULTIQC.out.versions)
    }

    emit:
    multiqc_report    = params.skip_multiqc ? Channel.empty() : MULTIQC.out.report
    star_log         = STAR_ALIGN.out.log_final
    bam              = PICARD_MARKDUPLICATES.out.bam
    bam_index        = PICARD_MARKDUPLICATES.out.bai  
    bigwig_all       = BAMCOVERAGE_STRANDED.out.bigwig_all
    bigwig_forward   = BAMCOVERAGE_STRANDED.out.bigwig_forward
    bigwig_reverse   = BAMCOVERAGE_STRANDED.out.bigwig_reverse
    size_factors     = params.spike_in_genome ? CALCULATE_SIZE_FACTORS.out.factors : Channel.empty()
    norm_report      = params.spike_in_genome ? CALCULATE_SIZE_FACTORS.out.report : Channel.empty()
    versions         = ch_versions
}

/*
========================================================================================
    FUNCTIONS
========================================================================================
*/

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id           = row.sample
    meta.single_end   = row.containsKey("single_end") ? row.single_end.toBoolean() : false
    meta.strandedness = row.containsKey("strandedness") ? row.strandedness : 'unstranded'
    
    // add path(s) of the fastq file(s) to the meta map
    def fastq_meta = []
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\\n${row.fastq_1}"
    }
    if (meta.single_end) {
        fastq_meta = [ meta, [ file(row.fastq_1) ] ]
    } else {
        if (!file(row.fastq_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\\n${row.fastq_2}"
        }
        fastq_meta = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
    }
    return fastq_meta
}