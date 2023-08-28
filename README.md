# ctSOM
analysis pipeline for somatic diagnosis on ctDNA
# Introduction
ctSOM is a variant analysis pipeline designed for circulating tumor DNA (ctDNA) obtained through non-invasive liquid biopsy, serving diagnostic, prognostic, and therapeutic purposes.

The pipeline is built using Nextflow, a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses conda containers making installation trivial and results highly reproducible.

The strength of this pipeline is attributed to the strategic integration of Unique Molecular Identifiers (UMIs). These molecular barcodes effectively mitigate the challenges posed by PCR duplicates and sequencing errors, contributing to the pipeline's enhanced accuracy and sensitivity in identifying genetic variants

# Pipeline summary
the pipeline can currently perform the following
- Extraction of the UMI from insert reads and generate consensus reads from them (fgbio )
- Adapter and quality trimming (fastp).
- Map Reads to Reference (bwa mem).
- Process BAM file (gatk SortSam, gatk MergeBamAlignment)
- convertion of file format ( gatk SamToFastq, gatk FastqToSam, vcf2tsv)
- quality control (samtools depth, picard collectHsMetrcis, multiQc)
- variant calling (vardict)
- Annotation (VEP)
# add image of pipeline

# Quick Start
1. Install nextflow.
2. creat this conda environment gatk4: conda create -n gatk4 -c conda-forge -c bioconda gatk4.
3. activate the conda environment you create, then istall the following tools:
  - fgbio v2.1.0  : fulcrumgenomics.github.io/fgbio
  - bwa v0.7.17   :bio-bwa.sourceforge.net
  - fastp v0.23.2 :github.com/OpenGene/fastp
  - samtools v1.17 :htslib.org
  - picard
  - vardict
  - VEP
  - vcf2tsvpy
5. Download the reference genome
6. Download the dataset needed for VEP use (Home_sapiens.....fa.gz , and the cache: homo_sapiens_vep_......tar.gz)

7. Download the pipeline and Run it on your Dataset:
   * Nextflow -log <OUTDIR>/my.log run ctSOM.nf --input /path_to_your_dataset/ *.fastq.gz --outdir <OUTDIR> --ref /path_to_your_refrence_genome/ .fa  -with-report <OUTDIR>/report.html -with-timeline <OUTDIR>/timeline.html -with-dag <OUTDIR>/flowchart.dot

# Pipeline output
the results of your run are all in the directory you defined <OUTDIR> . the log, html and dot files are optional output to give you insight about your run where:
- my.log : Logs are important for debugging, tracking the execution progress, and identifying any errors that might occur during the execution of the pipeline
- The HTML report typically includes information about the execution of each process in the workflow, such as input and output files, execution times, and any errors encountered.
- The timeline provides a visual representation of the execution timeline of processes in the workflow, showing when each process started and finished.
- The **.dot** format is a standard format for describing graphs. it visualizes the workflow's structure and dependencies between processes.
