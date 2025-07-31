/*
========================================================================================
    Spike-in normalization processes
========================================================================================
*/

process EXTRACT_SPIKE_COUNTS {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::samtools=1.15.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.15.1--h1170115_0' :
        'quay.io/biocontainers/samtools:1.15.1--h1170115_0' }"

    input:
    tuple val(meta), path(primary_bam), path(spike_bam)

    output:
    tuple val(meta), path("${meta.id}_spike_counts.txt"), emit: counts
    path  "versions.yml"                                , emit: versions

    when:
    params.spike_in_genome && (task.ext.when == null || task.ext.when)

    script:
    """
    # Count mapped reads in spike-in genome
    samtools view -c -F 4 $spike_bam > ${meta.id}_spike_counts.txt
    
    # Also get total reads for QC
    samtools view -c $spike_bam > ${meta.id}_spike_total.txt
    
    # Calculate mapping rate
    MAPPED=\$(cat ${meta.id}_spike_counts.txt)
    TOTAL=\$(cat ${meta.id}_spike_total.txt)
    if [ \$TOTAL -gt 0 ]; then
        RATE=\$(echo "scale=4; \$MAPPED / \$TOTAL * 100" | bc)
        echo "Sample: ${meta.id}, Spike-in mapped: \$MAPPED, Total: \$TOTAL, Rate: \$RATE%" > ${meta.id}_spike_stats.txt
    else
        echo "Sample: ${meta.id}, Spike-in mapped: 0, Total: 0, Rate: 0%" > ${meta.id}_spike_stats.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}

process CALCULATE_SIZE_FACTORS {
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::bioconductor-deseq2=1.32.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-deseq2:1.32.0--r41h399db7b_0' :
        'quay.io/biocontainers/bioconductor-deseq2:1.32.0--r41h399db7b_0' }"

    input:
    path spike_counts_files
    path spike_stats_files

    output:
    path "size_factors.tsv"           , emit: factors
    path "normalization_report.html" , emit: report
    path "spike_in_summary.tsv"      , emit: summary
    path "versions.yml"              , emit: versions

    when:
    params.spike_in_genome && (task.ext.when == null || task.ext.when)

    script:
    """
    #!/usr/bin/env Rscript
    
    # Load required libraries
    suppressPackageStartupMessages({
        library(DESeq2)
        library(knitr)
        library(rmarkdown)
    })
    
    # Read spike-in counts for all samples
    count_files <- list.files(".", pattern = "_spike_counts.txt", full.names = TRUE)
    sample_names <- gsub("_spike_counts.txt", "", basename(count_files))
    
    # Create count matrix
    spike_counts <- data.frame(row.names = "spike_in_total")
    for (i in seq_along(count_files)) {
        counts <- as.integer(readLines(count_files[i]))
        spike_counts[1, sample_names[i]] <- counts
    }
    
    # Ensure we have at least one count > 0
    if (all(spike_counts == 0)) {
        stop("ERROR: All spike-in counts are zero. Check spike-in alignment.")
    }
    
    # Create DESeq2 object
    coldata <- data.frame(
        sample = colnames(spike_counts),
        condition = rep("spike_in", ncol(spike_counts)),
        row.names = colnames(spike_counts)
    )
    
    dds <- DESeqDataSetFromMatrix(
        countData = spike_counts,
        colData = coldata,
        design = ~ 1
    )
    
    # Calculate size factors using DESeq2
    dds <- estimateSizeFactors(dds)
    size_factors <- sizeFactors(dds)
    
    # Calculate reciprocal scaling factors for BigWig generation
    scaling_factors <- 1 / size_factors
    
    # Create output table
    factor_table <- data.frame(
        sample = names(size_factors),
        spike_counts = as.numeric(spike_counts[1, ]),
        size_factor = as.numeric(size_factors),
        scaling_factor = as.numeric(scaling_factors),
        row.names = NULL
    )
    
    # Write size factors table
    write.table(
        factor_table,
        "size_factors.tsv",
        sep = "\\t",
        row.names = FALSE,
        quote = FALSE
    )
    
    # Read spike-in stats for summary
    stats_files <- list.files(".", pattern = "_spike_stats.txt", full.names = TRUE)
    stats_data <- data.frame()
    for (file in stats_files) {
        line <- readLines(file)
        # Parse the stats line 
        sample <- sub(".*Sample: ([^,]+),.*", "\\\\1", line)
        mapped <- as.numeric(sub(".*mapped: ([0-9]+),.*", "\\\\1", line))
        total <- as.numeric(sub(".*Total: ([0-9]+),.*", "\\\\1", line))
        rate <- as.numeric(sub(".*Rate: ([0-9.]+)%.*", "\\\\1", line))
        
        stats_data <- rbind(stats_data, data.frame(
            sample = sample,
            spike_mapped = mapped,
            spike_total = total,
            spike_rate = rate
        ))
    }
    
    # Merge with size factors
    summary_table <- merge(factor_table, stats_data, by = "sample")
    
    # Write summary
    write.table(
        summary_table,
        "spike_in_summary.tsv", 
        sep = "\\t",
        row.names = FALSE,
        quote = FALSE
    )
    
    # Generate HTML report
    cat('
    ---
    title: "Spike-in Normalization Report"
    output: html_document
    ---
    
    ```{r setup, include=FALSE}
    knitr::opts_chunk\$set(echo=FALSE, warning=FALSE, message=FALSE)
    library(ggplot2)
    library(knitr)
    
    # Read data
    summary_data <- read.table("spike_in_summary.tsv", header=TRUE, sep="\\t")
    ```
    
    ## Spike-in Alignment Summary
    
    ```{r summary_table}
    kable(summary_data, caption = "Spike-in alignment and normalization factors")
    ```
    
    ## Spike-in Count Distribution
    
    ```{r count_plot}
    ggplot(summary_data, aes(x=sample, y=spike_counts)) +
        geom_bar(stat="identity", fill="steelblue") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = "Spike-in Read Counts by Sample",
             x = "Sample", y = "Spike-in Counts")
    ```
    
    ## Size Factors
    
    ```{r size_factor_plot}
    ggplot(summary_data, aes(x=sample, y=size_factor)) +
        geom_bar(stat="identity", fill="darkgreen") +
        geom_hline(yintercept=1, linetype="dashed", color="red") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = "DESeq2 Size Factors",
             x = "Sample", y = "Size Factor",
             subtitle = "Red line indicates no normalization (size factor = 1)")
    ```
    
    ## Scaling Factors for BigWig
    
    ```{r scaling_factor_plot}
    ggplot(summary_data, aes(x=sample, y=scaling_factor)) +
        geom_bar(stat="identity", fill="orange") +
        geom_hline(yintercept=1, linetype="dashed", color="red") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = "Scaling Factors for BigWig Generation",
             x = "Sample", y = "Scaling Factor (1/Size Factor)",
             subtitle = "These factors are applied to bamCoverage")
    ```
    
    ', file = "report.Rmd")
    
    # Render the report
    rmarkdown::render("report.Rmd", output_file = "normalization_report.html")
    
    # Write versions
    cat(paste0(
        '"', '${task.process}', '":\\n',
        '    r-base: ', R.version.string, '\\n',
        '    bioconductor-deseq2: ', packageVersion("DESeq2"), '\\n'
    ), file = "versions.yml")
    """
}