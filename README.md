# MobiCT
ctDNA Analysis pipeline. Version 1.0.0

# Introduction
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
7. Activate your conda environment
8. Run MobiCT on your Dataset

    `Nextflow -log /output_directory/my.log run MobiCT.nf -c nextflow.config`

# Pipeline output
The output files are stored in the directory you specified using the **outdir** parameter in the .config file. The **outdir** contains a folder per sample plus a **multiQC** folder. Each sample folder contains a deduplicated and aligned *.bam* file and its index, an annotated *.vcf* file, 3 metrics files (*.HsMetrics.1.txt*, *.HsMetrics.3.txt* and *QC.bcftools_stats.stats*).

# Data availability
Raw sequencing data (FASTQ files) of commercial controls used in the study are available @ https://www.ncbi.nlm.nih.gov/sra/PRJNA1209006.
