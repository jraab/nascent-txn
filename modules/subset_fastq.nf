process SUBSET_FASTQ {
    tag "$meta.id"
    label 'process_low'

    publishDir "${params.outdir}/debug/subset_fastq", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_subset_{1,2}.fastq.gz"), emit: reads
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def num_reads = params.debug_num_reads ?: 10000
    
    """
    # Extract subset of reads maintaining paired-end structure
    # Using seqtk to sample the same reads from both files
    
    seqtk sample -s100 ${reads[0]} ${num_reads} | gzip > ${prefix}_subset_1.fastq.gz
    seqtk sample -s100 ${reads[1]} ${num_reads} | gzip > ${prefix}_subset_2.fastq.gz
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqtk: \$(echo \$(seqtk 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_subset_1.fastq.gz
    touch ${prefix}_subset_2.fastq.gz
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqtk: \$(echo \$(seqtk 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}