# MobiCT
ctDNA Analysis pipeline. Version 1.0.0

## Table of Contents

* [Introduction](#introduction)
* [Quick Start](#quick-start)
* [Installation](#installation)
    * [Nextflow](#nextflow)
    * [Execution environment](#execution-environment)
        * [Container platform](#container-platform)
            * [Docker](#docker)
            * [Singularity](#singularity)
        * [Conda](conda)
* [Usage](#usage)
* [Pipeline output](#pipeline-output)
* [Data availability](#data-availability)

## Introduction
Mobict is an analysis pipeline designed for detecting SNV and small InDels in circulating tumor DNA (ctDNA) obtained through non-invasive liquid biopsy, serving diagnostic, prognostic, and therapeutic purposes.
The pipeline is built using Nextflow.
The associated publication is available here: /doi/.../


<img src="https://github.com/user-attachments/assets/f44e5d99-4a85-423e-bf13-17ee13b5420c" width= 30% height=30%>



# Quick Start

1. Install nextflow (https://www.nextflow.io/docs/latest/install.html).
2. Create a conda environment for MobiCT:

    `conda create -n myenv -c conda-forge -c bioconda gatk4 fgbio bwa fastp samtools picard vardict ensembl-vep`
4. Download the reference genome
5. Download the datasets needed by VEP (see https://github.com/Ensembl/ensembl-vep)
6. Edit the *.config* file with input and output files/paths
7. Run MobiCT on your Dataset

    `Nextflow -log /output_directory/my.log run MobiCT.nf -c nextflow.config`

## Installation

The prerequisites to run the pipeline are:  

  * [Nextflow](https://www.nextflow.io/)  >= 22.04.0
  * [Docker](https://www.docker.com) or [Singularity](https://sylabs.io/singularity/) or [conda](https://juke34.github.io/fix-anaconda-licensing-issues/en/pages/conda-distrib/)

### Nextflow 

  * Via conda 

    <details>
      <summary>See here</summary>
      
      ```bash
      conda create -n nextflow
      conda activate nextflow
      conda install bioconda::nextflow
      ```  
    </details>

  * Manually
    <details>
      <summary>See here</summary>
      Nextflow runs on most POSIX systems (Linux, macOS, etc) and can typically be installed by running these commands:

      ```bash
      # Make sure 11 or later is installed on your computer by using the command:
      java -version
      
      # Install Nextflow by entering this command in your terminal(it creates a file nextflow in the current dir):
      curl -s https://get.nextflow.io | bash 
      
      # Add Nextflow binary to your user's PATH:
      mv nextflow ~/bin/
      # OR system-wide installation:
      # sudo mv nextflow /usr/local/bin
      ```
    </details>

### Execution environment

#### Container platform

To run the workflow you will need a container platform: docker or singularity.

##### Docker

Please follow the instructions at the [Docker website](https://docs.docker.com/desktop/)

##### Singularity

Please follow the instructions at the [Singularity website](https://docs.sylabs.io/guides/latest/admin-guide/installation.html)

#### Conda

It is recommended to use a container platform. However, if you prefer to use Conda, please follow the instructions available [here](https://juke34.github.io/fix-anaconda-licensing-issues/en/pages/conda-installation/)

## Usage

### Profile

To run the workflow you must select a profile according to the container platform you want to use:   
- `singularity`, a profile using Singularity containers to run the softwares
- `docker`, a profile using Docker containers to run the softwares
- `conda`, a profile using conda environments to run the softwares


### From Release

#### Default usage
The command will look like that: 

```bash
nextflow run SimCab-CHU/MobiCT -r v1.0.0 -profile docker <rest of paramaters>
```

#### Test

```bash
nextflow run SimCab-CHU/MobiCT -r v1.0.0 -profile docker,test <rest of paramaters>
```

### From repo

```bash
# clone the workflow repository
git clone https://github.com/SimCab-CHU/MobiCT.git

# Move in it
cd MobiCT
```
#### Default usage

```bash
# Run MobiCT
nextflow run MobiCT.nf -profile docker <rest of paramaters>
```

#### Test

```bash
# Run MobiCT
nextflow run MobiCT.nf -profile docker,test
```

## Pipeline output
The output files are stored in the directory you specified using the **outdir** parameter in the .config file. The **outdir** contains a folder per sample plus a **multiQC** folder. Each sample folder contains a deduplicated and aligned *.bam* file and its index, an annotated *.vcf* file, 3 metrics files (*.HsMetrics.1.txt*, *.HsMetrics.3.txt* and *QC.bcftools_stats.stats*).

## Data availability
Raw sequencing data (FASTQ files) of commercial controls used in the study are available @ https://www.ncbi.nlm.nih.gov/sra/PRJNA1209006.
