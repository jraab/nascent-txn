# TT-seq Nextflow Pipeline Development

## Project Overview
This is a Nextflow pipeline for processing TT-seq (Transient Transcriptome sequencing) data, designed to convert FASTQ files to normalized BigWig coverage tracks. The pipeline is based on the workflow from [crickbabs/DRB_TT-seq](https://github.com/crickbabs/DRB_TT-seq) but restructured using the raab-lab standards from [raab-lab/cut-n-run](https://github.com/raab-lab/cut-n-run).

## Key Features

### Core Workflow
1. **Quality Control**: FastQC + MultiQC reporting
2. **Alignment**: STAR two-pass alignment with optimized parameters for nascent RNA
3. **Post-processing**: BAM sorting, indexing, and duplicate marking
4. **Coverage Generation**: Strand-specific BigWig files with optional spike-in normalization
5. **TT-seq Analysis**: Transcription wave analysis and elongation rate calculations

### Spike-in Normalization Support
- Optional spike-in genome parameter (`--spike_in_genome`)
- DESeq2-based size factor calculation using spike-in reads only
- Reciprocal scaling factors applied to BigWig generation
- Normalization factors exported for downstream analysis

## Project Structure
```
nascent_txn/
├── bin/                    # Executable scripts
├── docs/                   # Documentation
├── helper/                 # Helper scripts
├── modules/                # Nextflow module definitions
│   ├── qc.nf              # FastQC + MultiQC
│   ├── align.nf           # STAR alignment
│   ├── samtools.nf        # BAM processing
│   ├── coverage.nf        # BigWig generation
│   ├── spike_in_norm.nf   # Spike-in normalization
│   └── ttseq_analysis.nf  # TT-seq specific analysis
├── subworkflows/          # Composite workflows
│   └── ttseq.nf          # Main TT-seq workflow
├── main.nf               # Primary pipeline entry
├── nextflow.config       # Configuration
└── README.md            # Documentation
```

## Key Parameters

### Basic Parameters
- `sample_sheet`: Input sample sheet path
- `outdir`: Output directory (default 'results')
- `genome`: Genome version (hg38, mm10)

### Spike-in Parameters
- `spike_in_genome`: Path to spike-in genome (optional)
- `spike_in_index`: Path to spike-in STAR index
- `spike_in_gtf`: GTF file for spike-in genome
- `spike_in_name`: Name for spike-in (e.g., 'yeast')

## Critical STAR Parameters (from DRB_TT-seq)
- `--runThreadN`: Multi-threading
- `--runMode alignReads`
- `--twopassMode Basic`
- `--quantMode TranscriptomeSAM GeneCounts`
- `--readFilesCommand zcat`
- `--outSAMunmapped None`
- `--outSAMtype BAM SortedByCoordinate`

## BigWig Generation Parameters
- Strand-specific tracks: `--filterRNAstrand forward/reverse`
- Scaling factors from spike-in normalization (reciprocal of DESeq2 size factors)
- Thread allocation for optimal performance

## Development Status
- **Completed**: 
  - Repository analysis and comprehensive planning
  - Complete pipeline structure implementation
  - All core modules (QC, alignment, BAM processing, BigWig generation)
  - Spike-in normalization workflow with DESeq2 integration
  - Configuration files for hg38, mm10, and yeast spike-in
  - Comprehensive documentation (usage, output, citations)
  - Sample validation scripts and example samplesheet
- **Ready for**: Testing and deployment

## Usage Examples

**With spike-in normalization:**
```bash
nextflow run main.nf --sample_sheet samples.csv --spike_in_genome yeast \
    --spike_in_index /path/to/yeast_index --profile hg38_yeast
```

**Without spike-in (standard mode):**
```bash
nextflow run main.nf --sample_sheet samples.csv --profile hg38
```

## Pipeline Structure Implemented

### Core Files Created
- `main.nf`: Pipeline entry point with sample sheet validation
- `nextflow.config`: Main configuration with profile support
- `subworkflows/ttseq.nf`: Complete TT-seq processing workflow

### Modules Implemented
- `modules/qc.nf`: FastQC and MultiQC for quality control
- `modules/align.nf`: STAR alignment for main and spike-in genomes
- `modules/samtools.nf`: BAM processing, sorting, indexing, duplicate marking
- `modules/coverage.nf`: Strand-specific BigWig generation with scaling
- `modules/spike_in_norm.nf`: DESeq2-based normalization calculation
- `modules/check_samplesheet.nf`: Sample sheet validation
- `modules/functions.nf`: Utility functions

### Configuration Files
- `conf/base.config`: Resource allocation and process settings
- `conf/modules.config`: Module-specific parameters and publishing
- `conf/genomes/hg38.config`: Human genome configuration
- `conf/genomes/mm10.config`: Mouse genome configuration  
- `conf/genomes/yeast_spike.config`: Yeast spike-in configuration

### Documentation
- `README.md`: Complete project documentation with usage examples
- `docs/usage.md`: Detailed usage instructions and parameter descriptions
- `docs/output.md`: Comprehensive output file documentation
- `CITATIONS.md`: Complete citation list for all tools used

### Supporting Files
- `bin/check_samplesheet.py`: Python script for sample sheet validation
- `lib/WorkflowMain.groovy`: Groovy functions for workflow management
- `assets/samplesheet.csv`: Example sample sheet template

## Reference Repositories
- Source workflow: [crickbabs/DRB_TT-seq](https://github.com/crickbabs/DRB_TT-seq)
- Structure template: [raab-lab/cut-n-run](https://github.com/raab-lab/cut-n-run)