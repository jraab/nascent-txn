# raab-lab/nascent_txn

**A Nextflow pipeline for TT-seq (transient transcriptome sequencing) data analysis**

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.3-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction

**raab-lab/nascent_txn** is a bioinformatics pipeline for processing TT-seq (transient transcriptome sequencing) data. The pipeline takes FASTQ files as input and produces normalized, strand-specific BigWig coverage tracks suitable for visualization and downstream analysis of nascent transcription.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute environments in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies.

## Pipeline summary

1. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Alignment to reference genome ([`STAR`](https://github.com/alexdobin/STAR))
3. Optional alignment to spike-in genome for normalization
4. Sort and index alignments ([`SAMtools`](https://samtools.github.io/))
5. Mark duplicate reads ([`Picard`](https://broadinstitute.github.io/picard/))
6. Calculate spike-in normalization factors ([`DESeq2`](https://bioconductor.org/packages/release/bioc/html/DESeq2.html))
7. Generate strand-specific BigWig coverage tracks ([`deepTools`](https://deeptools.readthedocs.io/))
8. Present QC for raw reads and alignment results ([`MultiQC`](http://multiqc.info/))

## Key Features

### Spike-in Normalization
- Optional spike-in genome support for accurate normalization between samples
- DESeq2-based size factor calculation using spike-in read counts
- Reciprocal scaling factors applied to BigWig generation
- Normalization factors exported for downstream analysis

### Strand-specific Analysis
- Generates separate BigWig tracks for forward and reverse strands
- Essential for nascent RNA analysis and transcription directionality
- Compatible with standard genome browsers (IGV, UCSC)

### Quality Control
- Comprehensive QC reporting with MultiQC
- Alignment statistics and duplicate marking
- Spike-in alignment assessment
- Library complexity evaluation

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.10.3`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) (you can follow [this tutorial](https://singularity-tutorial.github.io/01-installation/)), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(you can use [`Conda`](https://conda.io/miniconda.html) both to install Nextflow itself and also to manage software within pipelines. Please only use it within pipelines as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_.

3. Download the pipeline and test it on a minimal dataset with a single command:

   ```bash
   nextflow run raab-lab/nascent_txn -profile test,docker
   ```

   Note that some form of configuration will be needed so that Nextflow knows how to fetch the required software. This is usually done in the form of a config profile (`YOURPROFILE` in the example command above). You can chain multiple config profiles in a comma-separated string.

   > * The pipeline comes with config profiles called `docker`, `singularity`, `podman`, `shifter`, `charliecloud` and `conda` which instruct the pipeline to use the named tool for software management. For example, `-profile test,docker`.
   > * Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines on your institution's cluster is already available!
   > * For more information about reproducible compute environments that can be used by Nextflow, please see [this link](https://www.nextflow.io/docs/latest/container.html).

4. Start running your own analysis!

   ```bash
   nextflow run raab-lab/nascent_txn --sample_sheet samplesheet.csv --genome hg38 -profile docker
   ```

## Usage

### Basic analysis
```bash
nextflow run raab-lab/nascent_txn \\
    --sample_sheet samplesheet.csv \\
    --genome hg38 \\
    --outdir results \\
    -profile docker
```

### With spike-in normalization
```bash
nextflow run raab-lab/nascent_txn \\
    --sample_sheet samplesheet.csv \\
    --genome hg38 \\
    --spike_in_genome yeast \\
    --spike_in_index /path/to/yeast/star_index \\
    --outdir results \\
    -profile docker
```

See [usage docs](docs/usage.md) for all of the available options when running the pipeline.

### Samplesheet

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 3 columns, and a header row as shown in the examples below.

```csv
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
CONTROL_REP2,AEG588A2_S2_L002_R1_001.fastq.gz,AEG588A2_S2_L002_R2_001.fastq.gz
TREATMENT_REP1,AEG588A4_S4_L003_R1_001.fastq.gz,AEG588A4_S4_L003_R2_001.fastq.gz
```

| Column    | Description                                                                                                                 |
| --------- | --------------------------------------------------------------------------------------------------------------------------- |
| `sample`  | Custom sample name. Spaces in sample names are automatically converted to underscores (`_`).                               |
| `fastq_1` | Full path to FastQ file for Illumina short reads 1. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz". |
| `fastq_2` | Full path to FastQ file for Illumina short reads 2. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz". |

An [example samplesheet](assets/samplesheet.csv) has been provided with the pipeline.

## Pipeline output

To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/rnaseq/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the [output documentation](docs/output.md).

## Credits

This pipeline was developed by the Jesse Raab Lab at UNC Chapel Hill, based on the TT-seq workflow from [crickbabs/DRB_TT-seq](https://github.com/crickbabs/DRB_TT-seq) and structured following the [raab-lab/cut-n-run](https://github.com/raab-lab/cut-n-run) template.

## Citations

If you use raab-lab/nascent_txn for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX)

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Holger Hotz, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).