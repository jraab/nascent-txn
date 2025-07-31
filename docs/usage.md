# TT-seq Pipeline Usage

## Introduction

This pipeline processes TT-seq (transient transcriptome sequencing) data from FASTQ files to generate normalized BigWig coverage tracks and perform basic quality control analysis.

## Running the pipeline

### Quick start

The typical command for running the pipeline is as follows:

```bash
nextflow run raab-lab/nascent_txn --sample_sheet samplesheet.csv --genome hg38 -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

### Core Nextflow arguments

> **NB:** These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

#### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Conda) - see below.

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from <https://github.com/nf-core/configs> when it runs, making multiple config profiles for various institutional clusters available at runtime. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended.

* `docker`
  * A generic configuration profile to be used with [Docker](https://docker.com/)
* `singularity`
  * A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
* `conda`
  * A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter or Charliecloud.

#### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

#### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Pipeline parameters

### Input/output options

#### `--sample_sheet`

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 3 columns, and a header row. See [usage docs](https://github.com/raab-lab/nascent_txn/blob/main/docs/usage.md#samplesheet-input).

#### `--outdir`

The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure.

### Reference genome options

#### `--genome`

Name of iGenomes reference to use. Currently supported: `hg38`, `mm10`.

Alternatively, you can specify the following parameters individually:

#### `--star_index`

Path to directory containing STAR genome index.

#### `--gtf`

Path to GTF annotation file.

#### `--fasta`

Path to genome FASTA file.

### Spike-in normalization options

#### `--spike_in_genome`

Name of spike-in genome to use for normalization (e.g., 'yeast'). When specified, reads will be aligned to both the main genome and spike-in genome, and DESeq2 size factors will be calculated from spike-in read counts.

#### `--spike_in_index`

Path to directory containing STAR index for spike-in genome.

#### `--spike_in_gtf`

Path to GTF annotation file for spike-in genome.

### Analysis options

#### `--skip_fastqc`

Skip FastQC quality control analysis.

#### `--skip_multiqc`

Skip MultiQC report generation.

#### `--bigwig_binsize`

Bin size for BigWig file generation (default: 10).

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 3 columns, and a header row as shown in the examples below.

```csv
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
CONTROL_REP2,AEG588A2_S2_L002_R1_001.fastq.gz,AEG588A2_S2_L002_R2_001.fastq.gz
CONTROL_REP3,AEG588A3_S3_L002_R1_001.fastq.gz,AEG588A3_S3_L002_R2_001.fastq.gz
TREATMENT_REP1,AEG588A4_S4_L003_R1_001.fastq.gz,AEG588A4_S4_L003_R2_001.fastq.gz
TREATMENT_REP2,AEG588A5_S5_L003_R1_001.fastq.gz,AEG588A5_S5_L003_R2_001.fastq.gz
TREATMENT_REP3,AEG588A6_S6_L003_R1_001.fastq.gz,AEG588A6_S6_L003_R2_001.fastq.gz
```

| Column    | Description                                                                                                                                                                            |
| --------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`  | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Spaces in sample names are automatically converted to underscores (`_`). |
| `fastq_1` | Full path to FastQ file for Illumina short reads 1. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |
| `fastq_2` | Full path to FastQ file for Illumina short reads 2. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

## Output

For a description of the output files and reports, please refer to the [output documentation](output.md).

## Example commands

### Basic TT-seq analysis (human genome)

```bash
nextflow run raab-lab/nascent_txn \\
    --sample_sheet samplesheet.csv \\
    --genome hg38 \\
    --outdir results \\
    -profile docker
```

### TT-seq analysis with yeast spike-in normalization

```bash
nextflow run raab-lab/nascent_txn \\
    --sample_sheet samplesheet.csv \\
    --genome hg38 \\
    --spike_in_genome yeast \\
    --spike_in_index /path/to/yeast/star_index \\
    --spike_in_gtf /path/to/yeast/annotation.gtf \\
    --outdir results \\
    -profile docker
```

### Mouse genome analysis

```bash
nextflow run raab-lab/nascent_txn \\
    --sample_sheet samplesheet.csv \\
    --genome mm10 \\
    --outdir results \\
    -profile docker
```