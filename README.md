# MobiCT - ctDNA Analysis Pipeline

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A520.04.0-23aa62.svg)](https://www.nextflow.io/)
[![Run with conda](https://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![Run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)

## Introduction

**MobiCT** is an analysis pipeline designed for detecting SNVs (Single Nucleotide Variants) and small InDels in circulating tumor DNA (ctDNA) obtained through non-invasive liquid biopsy. The pipeline serves diagnostic, prognostic, and therapeutic purposes in precision oncology.

The pipeline is built using Nextflow, a workflow tool to run tasks across multiple compute environments in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible.

## Pipeline Summary

MobiCT performs the following key steps:

1. **Quality Control**: Raw read quality assessment using FastP
2. **Alignment**: Read mapping to reference genome using BWA-MEM
3. **Deduplication**: PCR duplicate removal using Picard/fgbio
4. **Variant Calling**: SNV and InDel detection using VarDict
5. **Annotation**: Variant annotation using Ensembl VEP
6. **Quality Metrics**: Comprehensive QC metrics generation
7. **Reporting**: MultiQC report generation

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/install.html) (`>=20.04.0`)

2. Create a conda environment with required tools:
   ```bash
   conda create -n mobict -c conda-forge -c bioconda \
     gatk4 fgbio bwa fastp samtools picard vardict ensembl-vep
   ```

3. Download the pipeline:
   ```bash
   git clone https://github.com/SimCab-CHU/MobiCT.git
   cd MobiCT
   ```

4. Test the pipeline with minimal dataset:
   ```bash
   conda activate mobict
   nextflow run MobiCT.nf -c nextflow.config --input test_data
   ```

## Usage

### Typical command

```bash
nextflow -log /path/to/output/my.log run MobiCT.nf -c nextflow.config
```

### Configuration

Before running the pipeline, edit the `nextflow.config` file to specify:
- Input FASTQ files paths
- Output directory
- Reference genome path
- Target intervals/BED files
- Resource allocation

### Input Requirements

The pipeline expects:
- **FASTQ files**: Paired-end sequencing data from ctDNA samples
- **Reference genome**: Human reference genome (e.g., GRCh38)
- **Target intervals**: BED file defining regions of interest
- **VEP database**: Pre-downloaded VEP cache and databases

### Reference Data Preparation

1. **Download reference genome** (GRCh38 recommended):
   ```bash
   wget http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
   ```

2. **Download VEP databases** (see [VEP documentation](https://github.com/Ensembl/ensembl-vep)):
   ```bash
   vep_install -a cf -s homo_sapiens -y GRCh38 -c /path/to/vep/cache
   ```

## Output

Results are organized in the specified output directory:

```
outdir/
├── sample1/
│   ├── sample1.deduplicated.bam
│   ├── sample1.deduplicated.bam.bai
│   ├── sample1.annotated.vcf
│   ├── sample1.HsMetrics.1.txt
│   ├── sample1.HsMetrics.3.txt
│   └── sample1.QC.bcftools_stats.stats
├── sample2/
│   └── ...
└── multiqc/
    └── multiqc_report.html
```

### Output Files Description

- **`.deduplicated.bam`**: Aligned, deduplicated BAM file
- **`.annotated.vcf`**: Variant calls with functional annotations
- **`.HsMetrics.*.txt`**: Hybrid selection metrics from Picard
- **`.QC.bcftools_stats.stats`**: Variant calling statistics
- **`multiqc_report.html`**: Comprehensive quality control report

## Parameters

### Core Options

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--input` | Path to input FASTQ files | Required |
| `--outdir` | Output directory | `./results` |
| `--genome` | Reference genome path | Required |
| `--intervals` | Target intervals BED file | Required |

### Resource Options

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--max_cpus` | Maximum number of CPUs | 16 |
| `--max_memory` | Maximum memory allocation | '128.GB' |
| `--max_time` | Maximum time per job | '240.h' |

## Profiles

The pipeline supports different execution profiles:

- `conda`: Use Conda for dependency management
- `docker`: Use Docker containers
- `singularity`: Use Singularity containers
- `test`: Run with test dataset

Example:
```bash
nextflow run MobiCT.nf -profile conda,test
```

## Test Data

Raw sequencing data (FASTQ files) of commercial controls used in the study are available at:
**NCBI SRA**: [PRJNA1209006](https://www.ncbi.nlm.nih.gov/sra/PRJNA1209006)

## Citations

If you use MobiCT for your analysis, please cite:

> **MobiCT: ctDNA Analysis Pipeline**
> 
> [Publication DOI will be added]

### Tools Citations

This pipeline uses several bioinformatics tools. Please also cite:

- **Nextflow**: Paolo Di Tommaso, et al. Nextflow enables reproducible computational workflows. Nature Biotechnology 35, 316–319 (2017)
- **BWA**: Li H. and Durbin R. Fast and accurate short read alignment with Burrows-Wheeler transform. Bioinformatics 25, 1754-1760 (2009)
- **VarDict**: Zhongwu Lai, et al. VarDict: a novel and versatile variant caller for next-generation sequencing in cancer research. Nucleic Acids Research 44, e108 (2016)
- **VEP**: McLaren W, et al. The Ensembl Variant Effect Predictor. Genome Biology 17, 122 (2016)
- **MultiQC**: Philip Ewels, et al. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics 32, 3047-3048 (2016)

## Credits

MobiCT was developed by the **SimCab team at CHU**.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Support

For questions and support:
- Create an issue on [GitHub](https://github.com/SimCab-CHU/MobiCT/issues)
- Contact the development team

## Changelog

### Version 1.0.0
- Initial release
- Support for SNV and small InDel detection in ctDNA
- Integrated quality control and reporting
- VEP-based variant annotation
