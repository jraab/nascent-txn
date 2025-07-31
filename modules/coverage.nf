/*
========================================================================================
    BigWig coverage track generation
========================================================================================
*/

process BAMCOVERAGE {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::deeptools=3.5.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/deeptools:3.5.1--py_0' :
        'quay.io/biocontainers/deeptools:3.5.1--py_0' }"

    input:
    tuple val(meta), path(bam), path(bai), val(scaling_factor)

    output:
    tuple val(meta), path("*.bw"), emit: bigwig
    path  "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def scale_arg = scaling_factor ? "--scaleFactor ${scaling_factor}" : ""
    def effective_genome_size = params.effective_genome_size ?: 2913022398
    
    """
    bamCoverage \\
        --bam $bam \\
        --outFileName ${prefix}.bw \\
        --outFileFormat bigwig \\
        --binSize ${params.bigwig_binsize} \\
        --numberOfProcessors ${task.cpus} \\
        --effectiveGenomeSize ${effective_genome_size} \\
        $scale_arg \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(bamCoverage --version | sed -e "s/bamCoverage //g")
    END_VERSIONS
    """
}

process BAMCOVERAGE_STRANDED {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::deeptools=3.5.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/deeptools:3.5.1--py_0' :
        'quay.io/biocontainers/deeptools:3.5.1--py_0' }"

    input:
    tuple val(meta), path(bam), path(bai), val(scaling_factor)

    output:
    tuple val(meta), path("*_all.bw")     , emit: bigwig_all
    tuple val(meta), path("*_forward.bw") , emit: bigwig_forward  
    tuple val(meta), path("*_reverse.bw") , emit: bigwig_reverse
    path  "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def scale_arg = scaling_factor ? "--scaleFactor ${scaling_factor}" : ""
    def effective_genome_size = params.effective_genome_size ?: 2913022398
    
    """
    # All reads
    bamCoverage \\
        --bam $bam \\
        --outFileName ${prefix}_all.bw \\
        --outFileFormat bigwig \\
        --binSize ${params.bigwig_binsize} \\
        --numberOfProcessors ${task.cpus} \\
        --effectiveGenomeSize ${effective_genome_size} \\
        $scale_arg \\
        $args

    # Forward strand (first in pair reverse strand, second in pair forward strand)
    bamCoverage \\
        --bam $bam \\
        --outFileName ${prefix}_forward.bw \\
        --outFileFormat bigwig \\
        --binSize ${params.bigwig_binsize} \\
        --numberOfProcessors ${task.cpus} \\
        --effectiveGenomeSize ${effective_genome_size} \\
        --filterRNAstrand forward \\
        $scale_arg \\
        $args

    # Reverse strand (first in pair forward strand, second in pair reverse strand)  
    bamCoverage \\
        --bam $bam \\
        --outFileName ${prefix}_reverse.bw \\
        --outFileFormat bigwig \\
        --binSize ${params.bigwig_binsize} \\
        --numberOfProcessors ${task.cpus} \\
        --effectiveGenomeSize ${effective_genome_size} \\
        --filterRNAstrand reverse \\
        $scale_arg \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(bamCoverage --version | sed -e "s/bamCoverage //g")
    END_VERSIONS
    """
}