# TT-seq Pipeline Output

## Introduction

This document describes the output produced by the TT-seq pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

* [FastQC](#fastqc) - Raw read QC
* [STAR](#star) - Alignment to reference genome
* [SAMtools](#samtools) - Alignment result processing
* [Picard](#picard) - Duplicate marking
* [deepTools](#deeptools) - BigWig coverage track generation
* [Spike-in normalization](#spike-in-normalization) - Optional spike-in based normalization
* [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
* [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

## FastQC

<details markdown="1">
<summary>Output files</summary>

* `fastqc/`
    * `*_fastqc.html`: FastQC report containing quality metrics.
    * `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

</details>

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your raw reads. It provides information about the quality score distribution across your reads, the per base sequence content (%A/T/G/C). You get information about adapter contamination and other overrepresented sequences.

For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

![MultiQC - FastQC sequence counts plot](images/mqc_fastqc_counts.png)

![MultiQC - FastQC mean quality scores plot](images/mqc_fastqc_quality.png)

![MultiQC - FastQC adapter content plot](images/mqc_fastqc_adapter.png)

> **NB:** The FastQC plots displayed in the MultiQC report shows _untrimmed_ reads. They may contain adapter sequence and potentially regions with low quality.

## STAR

<details markdown="1">
<summary>Output files</summary>

* `star/`
    * `*.Aligned.sortedByCoord.out.bam`: STAR aligned BAM file, coordinate sorted.
    * `*.Aligned.sortedByCoord.out.bam.bai`: BAM index file.
    * `*.Log.final.out`: STAR alignment log with alignment statistics.
    * `*.Log.out`: STAR log file containing detailed information about the run.
    * `*.Log.progress.out`: STAR log file containing job progress information.
    * `*.SJ.out.tab`: STAR splice junctions file.
    * `*.ReadsPerGene.out.tab`: STAR read counts per gene.

</details>

[STAR](https://github.com/alexdobin/STAR) is a read aligner designed for RNA sequencing. STAR stands for Spliced Transcripts Alignment to a Reference, it produces a coordinate sorted BAM file of alignments.

![MultiQC - STAR alignment scores plot](images/mqc_star_alignment_plot.png)

## SAMtools

<details markdown="1">
<summary>Output files</summary>

* `samtools/`
    * `*.stats`: SAMtools stats output with detailed alignment statistics.

</details>

[SAMtools](http://samtools.sourceforge.net/) is a suite of programs for interacting with high-throughput sequencing data. The stats subcommand produces comprehensive statistics from alignment files.

## Picard

<details markdown="1">
<summary>Output files</summary>

* `picard/`
    * `*.marked.bam`: BAM file with duplicates marked.
    * `*.marked.bam.bai`: BAM index file.
    * `*.metrics.txt`: Picard MarkDuplicates metrics file.

</details>

[Picard MarkDuplicates](https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates) is used to mark duplicate reads in the aligned BAM files. This is important for TT-seq analysis as it helps to distinguish between genuine nascent transcripts and PCR duplicates.

![MultiQC - Picard MarkDuplicates plot](images/mqc_picard_dups.png)

## deepTools

<details markdown="1">
<summary>Output files</summary>

* `bigwig/`
    * `*_all.bw`: BigWig file containing coverage from all reads.
    * `*_forward.bw`: BigWig file containing coverage from forward strand reads.
    * `*_reverse.bw`: BigWig file containing coverage from reverse strand reads.

</details>

[deepTools](https://deeptools.readthedocs.io/en/develop/) bamCoverage is used to generate BigWig coverage files from the aligned BAM files. For TT-seq data, strand-specific BigWig files are particularly important as they allow visualization of nascent transcription directionality.

Key features:
* **Strand-specific tracks**: Separate forward and reverse strand BigWig files allow visualization of transcription direction
* **Normalization**: Coverage tracks are normalized using RPGC (Reads Per Genomic Content) method
* **Spike-in scaling**: When spike-in genomes are used, scaling factors from DESeq2 are applied to account for differences in library preparation efficiency

## Spike-in normalization

<details markdown="1">
<summary>Output files</summary>

* `normalization/` (only generated when `--spike_in_genome` is specified)
    * `size_factors.tsv`: Table containing DESeq2 size factors and scaling factors for each sample.
    * `spike_in_summary.tsv`: Summary table with spike-in alignment statistics and normalization factors.
    * `normalization_report.html`: HTML report with plots showing spike-in count distribution and normalization factors.

</details>

When a spike-in genome is specified, the pipeline performs normalization based on spike-in read counts:

1. **Dual alignment**: Reads are aligned to both the main genome and spike-in genome
2. **Count extraction**: Mapped read counts are extracted from spike-in alignments
3. **Size factor calculation**: DESeq2's `estimateSizeFactors` function calculates normalization factors using only spike-in counts
4. **BigWig scaling**: The reciprocal of size factors is applied during BigWig generation

This approach ensures that differences in library preparation efficiency between samples are corrected, allowing for accurate comparison of nascent transcription levels.

![Spike-in count distribution](images/spike_in_counts.png)

![Size factors](images/size_factors.png)

## MultiQC

<details markdown="1">
<summary>Output files</summary>

* `multiqc/`
    * `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
    * `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
    * `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC, STAR, Picard, SAMtools. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

## Pipeline information

<details markdown="1">
<summary>Output files</summary>

* `pipeline_info/`
    * Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
    * Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used.
    * Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.

## Data interpretation

### BigWig files

The main output of the TT-seq pipeline are strand-specific BigWig coverage files:

* **`*_all.bw`**: Contains coverage from all reads, useful for overall transcription visualization
* **`*_forward.bw`**: Contains reads mapping to the forward strand
* **`*_reverse.bw`**: Contains reads mapping to the reverse strand

These files can be loaded into genome browsers (e.g., IGV, UCSC Genome Browser) to visualize nascent transcription patterns.

### Spike-in normalization

When spike-in normalization is used:

1. **Size factors < 1**: Sample has higher spike-in content, will be scaled down
2. **Size factors > 1**: Sample has lower spike-in content, will be scaled up  
3. **Scaling factors**: Applied to BigWig generation (reciprocal of size factors)

The normalization ensures that differences in library preparation efficiency are corrected across samples.

### Quality control considerations

Key QC metrics to examine:

* **Alignment rate**: Should be >80% for good quality TT-seq data
* **Duplicate rate**: Higher rates may indicate PCR amplification issues
* **Spike-in alignment**: Should be consistent across samples (typically 1-10%)
* **Strand specificity**: Forward/reverse strand balance should reflect experimental design