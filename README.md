# WARNING THIS IS A DEVELOPMENT VERSION NOT TO BE USED YET
Release soon!


# MobiCT
ctDNA Analysis pipeline 

# Introduction
Mobict is an analysis pipeline designed for detecting SNV and small InDels in circulating tumor DNA (ctDNA) obtained through non-invasive liquid biopsy, serving diagnostic, prognostic, and therapeutic purposes.

The pipeline is built using Nextflow, a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses conda containers making results highly reproducible.

# Quick Start

1. Install nextflow.
2. Create this conda environment MobiCT:

    `conda create -n MobiCT -c conda-forge -c bioconda`

4. Activate this conda environment, then istall the following tools:
   
    - fgbio
    - bwa
    - fastp
    - samtools
    - picard
    - vardict
    - VEP
  
    `conda install bioconda::fgbio bioconda::bwa bioconda::fastp bioconda::samtools bioconda::picard bioconda::vardict bioconda::ensembl-vep`
      
5. Download the reference genome, its indexed and dict
6. Download the dataset needed for VEP use
7. On nextflow.config file, edit input and output paths
8. Run it on your Dataset:
   * Nextflow -log /output_directory/my.log run MobiCT.nf -with-report /output_directory/report.html -with-timeline /output_directory/timeline.html -with-dag /output_directory/flowchart.dot

# Pipeline output
The outcomes of your execution are stored within the directory you specified using the --outdir parameter. The log, HTML, and DOT files, which are provided as optional outputs, offer you valuable insights into the progress of your run:
- my.log : Logs are important for debugging, tracking the execution progress, and identifying any errors that might occur during the execution of the pipeline
- report.html typically includes information about the execution of each process in the workflow, such as input and output files, execution times, and any errors encountered.
- timeline.html provides a visual representation of the execution timeline of processes in the workflow, showing when each process started and finished.
- flowchart.dot: The **.dot** format is a standard format for describing graphs. it visualizes the workflow's structure and dependencies between processes.

# Tips & tricks

- Each file outputed from MobiCT is named as follows: *sampleName_process* (*i.e.* *SampleTest_consensusMerge.bam* is the bam file output of the consensusMerge process). The sample name is returned from the **replaceExtension** function defined in the first lines of the **MobiCT.nf** file. It splits the fastq file names on a "_" separator, this splitting character can not be changed. However, if the user does'nt want the file name to be splitted, he/she must remove the "_" from the file name prior to the analysis.
- ongoing paragraph on UMI
