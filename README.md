# MobiCT
ctDNA Analysis pipeline 

# Introduction
Mobict is an analysis pipeline designed for detecting SNV and small InDels in circulating tumor DNA (ctDNA) obtained through non-invasive liquid biopsy, serving diagnostic, prognostic, and therapeutic purposes.

The pipeline is built using Nextflow, a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses conda containers making results highly reproducible.

<picture>
 <img alt="the workflow" src="https://github.com/Oussamadl/ctSOM/blob/main/pipeline.png">
</picture>

# Quick Start
1. Install nextflow.
2. creat this conda environment gatk4: conda create -n gatk4 -c conda-forge -c bioconda gatk4.
3. activate the conda environment you create, then istall the following tools:
  - fgbio  :  conda install bioconda::fgbio
  - bwa v0.7.17   :  conda install bioconda::bwa
  - fastp v0.23.2 :  conda install bioconda::fastp
  - samtools v1.17 :  conda install bioconda::samtools
  - picard : conda install bioconda::picard
  - vardict :  conda install bioconda::vardict
  - VEP :  conda install bioconda::ensembl-vep
5. Download the reference genome, and extract all files related to it through indexing: hg38.fa.amb, hg38.fa.ann, hg38.fa.bwt, hg38.fa.fai, hg38.fa.pac, hg38.fa.sa
6. Download the dataset needed for VEP use (Home_sapiens.....fa.gz , and the cache: homo_sapiens_vep_......tar.gz)
7. On nextflow.config file, define all input and output paths of your pipeline 
8. Download the pipeline and Run it on your Dataset:
   * Nextflow -log /output_directory/my.log run MobiCT.nf -with-report /output_directory/report.html -with-timeline /output_directory/timeline.html -with-dag /output_directory/flowchart.dot

# Pipeline output
The outcomes of your execution are stored within the directory you specified using the --outdir parameter. The log, HTML, and DOT files, which are provided as optional outputs, offer you valuable insights into the progress of your run:
- my.log : Logs are important for debugging, tracking the execution progress, and identifying any errors that might occur during the execution of the pipeline
- report.html typically includes information about the execution of each process in the workflow, such as input and output files, execution times, and any errors encountered.
- timeline.html provides a visual representation of the execution timeline of processes in the workflow, showing when each process started and finished.
- flowchart.dot: The **.dot** format is a standard format for describing graphs. it visualizes the workflow's structure and dependencies between processes.
