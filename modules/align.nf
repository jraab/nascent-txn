/*
========================================================================================
    STAR alignment processes
========================================================================================
*/

process STAR_ALIGN {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::star=2.7.9a bioconda::samtools=1.15.1 conda-forge::gawk=5.1.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:59cdd445419f14abac76b31dd0d71217994cbcc9-0' :
        'quay.io/biocontainers/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:59cdd445419f14abac76b31dd0d71217994cbcc9-0' }"

    input:
    tuple val(meta), path(reads)
    path  index
    path  gtf

    output:
    tuple val(meta), path('*d.out.bam')       , emit: bam
    tuple val(meta), path('*Log.final.out')   , emit: log_final
    tuple val(meta), path('*Log.out')         , emit: log_out
    tuple val(meta), path('*Log.progress.out'), emit: log_progress
    path  "versions.yml"                      , emit: versions

    tuple val(meta), path('*sortedByCoord.out.bam')  , optional:true, emit: bam_sorted
    tuple val(meta), path('*toTranscriptome.out.bam'), optional:true, emit: bam_transcript
    tuple val(meta), path('*Aligned.unsort.out.bam') , optional:true, emit: bam_unsorted
    tuple val(meta), path('*fastq.gz')                , optional:true, emit: fastq
    tuple val(meta), path('*.tab')                    , optional:true, emit: tab
    tuple val(meta), path('*.SJ.out.tab')             , optional:true, emit: spl_junc_tab
    tuple val(meta), path('*.ReadsPerGene.out.tab')   , optional:true, emit: read_per_gene_tab

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ignore_gtf = params.star_ignore_sjdbgtf ? '' : "--sjdbGTFfile $gtf"
    def seq_platform = meta.seq_platform ? "--outSAMattrRGline ID:$prefix 'SM:$prefix' 'PL:$meta.seq_platform'" : "--outSAMattrRGline ID:$prefix 'SM:$prefix'"
    def seq_center = meta.seq_center ? "'CN:$meta.seq_center'" : ""
    def out_sam_attributes = seq_center ? "$seq_platform $seq_center" : seq_platform
    def mv_unsorted_bam = ""
    if (args.contains('--outSAMtype BAM Unsorted')) {
        mv_unsorted_bam = "mv ${prefix}.Aligned.out.bam ${prefix}.Aligned.unsort.out.bam"
    }
    
    if (meta.single_end) {
        """
        STAR \\
            --genomeDir $index \\
            --readFilesIn ${reads[0]} \\
            --runThreadN $task.cpus \\
            --outFileNamePrefix $prefix. \\
            $out_sam_attributes \\
            $ignore_gtf \\
            $args

        $mv_unsorted_bam

        if [ -f ${prefix}.Aligned.sortedByCoord.out.bam ]; then
            samtools index ${prefix}.Aligned.sortedByCoord.out.bam
        fi

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            star: \$(STAR --version | sed -e "s/STAR_//g")
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
            gawk: \$(echo \$(gawk --version 2>&1) | sed 's/^.*GNU Awk //; s/, .*\$//')
        END_VERSIONS
        """
    } else {
        """
        STAR \\
            --genomeDir $index \\
            --readFilesIn ${reads[0]} ${reads[1]} \\
            --runThreadN $task.cpus \\
            --outFileNamePrefix $prefix. \\
            $out_sam_attributes \\
            $ignore_gtf \\
            $args

        $mv_unsorted_bam

        if [ -f ${prefix}.Aligned.sortedByCoord.out.bam ]; then
            samtools index ${prefix}.Aligned.sortedByCoord.out.bam
        fi

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            star: \$(STAR --version | sed -e "s/STAR_//g")
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
            gawk: \$(echo \$(gawk --version 2>&1) | sed 's/^.*GNU Awk //; s/, .*\$//')
        END_VERSIONS
        """
    }
}

process STAR_ALIGN_SPIKE {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::star=2.7.9a bioconda::samtools=1.15.1 conda-forge::gawk=5.1.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:59cdd445419f14abac76b31dd0d71217994cbcc9-0' :
        'quay.io/biocontainers/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:59cdd445419f14abac76b31dd0d71217994cbcc9-0' }"

    input:
    tuple val(meta), path(reads)
    path  spike_index
    path  spike_gtf

    output:
    tuple val(meta), path('*spike.bam')         , emit: bam
    tuple val(meta), path('*spike.Log.final.out'), emit: log_final
    path  "versions.yml"                        , emit: versions

    when:
    params.spike_in_genome && (task.ext.when == null || task.ext.when)

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ignore_gtf = params.star_ignore_sjdbgtf ? '' : "--sjdbGTFfile $spike_gtf"
    def seq_platform = meta.seq_platform ? "--outSAMattrRGline ID:${prefix}_spike 'SM:${prefix}_spike' 'PL:$meta.seq_platform'" : "--outSAMattrRGline ID:${prefix}_spike 'SM:${prefix}_spike'"
    
    if (meta.single_end) {
        """
        STAR \\
            --genomeDir $spike_index \\
            --readFilesIn ${reads[0]} \\
            --runThreadN $task.cpus \\
            --outFileNamePrefix ${prefix}.spike. \\
            $seq_platform \\
            $ignore_gtf \\
            $args

        mv ${prefix}.spike.Aligned.sortedByCoord.out.bam ${prefix}.spike.bam
        samtools index ${prefix}.spike.bam

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            star: \$(STAR --version | sed -e "s/STAR_//g")
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
    } else {
        """
        STAR \\
            --genomeDir $spike_index \\
            --readFilesIn ${reads[0]} ${reads[1]} \\
            --runThreadN $task.cpus \\
            --outFileNamePrefix ${prefix}.spike. \\
            $seq_platform \\
            $ignore_gtf \\
            $args

        mv ${prefix}.spike.Aligned.sortedByCoord.out.bam ${prefix}.spike.bam
        samtools index ${prefix}.spike.bam

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            star: \$(STAR --version | sed -e "s/STAR_//g")
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
    }
}